functions {
  vector diagSPD_EQ(real alpha, real rho, real L, int M){
    return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
  }

  matrix PHI(int N, int M, real L, vector x){
    return sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), linspaced_vector(M, 1, M)))/sqrt(L);
  }

/** Kronecker multivariate product
  *
  * Enables efficient sampling from a 2D multivariate normal distribution.
  *
  * @param A
  * @param B
  * @param V A matrix of N(0,1) distributed random variables
  */
  matrix kron_mvprod(matrix A, matrix B, matrix V){
    return (B*V) * transpose(A);
  }

  matrix hsgp(int A, int B, real alpha, real rho1, real rho2, real L1, real L2, int M1, int M2,
            matrix PHI1, matrix PHI2, matrix z)
  {
    vector[M1] sqrt_spd_1 = diagSPD_EQ(alpha, rho1, L1, M1);
    vector[M2] sqrt_spd_2 = diagSPD_EQ(alpha, rho2, L2, M2);

    matrix[A,B] f = kron_mvprod(
      diag_post_multiply( PHI2, sqrt_spd_2 ),
      diag_post_multiply( PHI1, sqrt_spd_1 ),
      z
    );

    return(f);
  }
}

data {
  int N_M; // Number of observations
  int N_F;
  int A; // Number of age groups
  int B; // Number of time points

  array[N_M] int Y_M; // Mortality counts
  array[N_M] int Y_F;

  vector[A] AGE_IDX;
  vector[B] TIME_IDX;

  vector[A*B] P_M;
  vector[A*B] P_F;

  row_vector[A] W; // Weights calculated from the standard population

  // HSGP parameters
  int<lower=1> M1;
  int<lower=1> M2;
  int<lower=0> C1;
  int<lower=0> C2;
}

transformed data {
  int M = 1;
  int F = 2;
  int G = 2;
  int N = N_M + N_F;

  matrix[A,B] LOG_P_M = to_matrix(log(P_M), A, B);
  matrix[A,B] LOG_P_F = to_matrix(log(P_F), A, B);

  // Standardise Input
  vector[A] AGE_IDX_STD = (AGE_IDX - mean(AGE_IDX)) / sd(AGE_IDX);
  vector[B] TIME_IDX_STD = (TIME_IDX - mean(TIME_IDX)) / sd(TIME_IDX);

  real L1 = C1 * max(AGE_IDX_STD);
  real L2 = C2 * max(TIME_IDX_STD);

  matrix[A,M1] PHI1 = PHI(A, M1, L1, AGE_IDX_STD);
  matrix[B,M2] PHI2 = PHI(B, M2, L2, TIME_IDX_STD);

  array[N] int Y = append_array(Y_M, Y_F);
}

parameters {
  real beta0;  // Baseline rate
  vector[A-1] beta_a; // Age effects

  // GP parameters
  array[G] real<lower=0> sigma; // scale parameters
  array[G] real<lower=0> rho1;  // lengthscale parameters 1
  array[G] real<lower=0> rho2;  // lengthscale parameters 2
  matrix[G*M1,M2] z; // GP mean
}

transformed parameters {
  array[G] matrix[A,B] log_lambda;
  { // Local scope
    array[G] matrix[A,B] f2;
    f2[M] = hsgp(A, B, sigma[M], rho1[M], rho2[M], L1, L2, M1, M2, PHI1, PHI2, z[1:M1,:]);
    f2[F] = hsgp(A, B, sigma[F], rho1[F], rho2[F], L1, L2, M1, M2, PHI1, PHI2, z[(M1+1):2*M1,:]);

    for (a in 1:A){
      if (a == 1){
        log_lambda[M][a,:] = beta0 + f2[M][a,:] + 1e-13;
        log_lambda[F][a,:] = beta0 + f2[F][a,:] + 1e-13;
      } else {
        log_lambda[M][a,:] = beta0 + beta_a[a-1] + f2[M][a,:] + 1e-13;
        log_lambda[F][a,:] = beta0 + beta_a[a-1] + f2[F][a,:] + 1e-13;
      }
    }
  }
}

model {
  target += normal_lpdf(beta0 | 0, 10);
  target += normal_lpdf(beta_a | 0, 1);

  target += inv_gamma_lpdf(sigma | 5, 5);
  target += inv_gamma_lpdf(rho1 | 5, 5);
  target += inv_gamma_lpdf(rho2 | 5, 10);
  target += normal_lpdf(to_vector(z) | 0, 1);

  { // Local scope
    array[G] vector[A*B] mu;
    mu[M] = exp( to_vector(log_lambda[M] + LOG_P_M) );
    mu[F] = exp( to_vector(log_lambda[F] + LOG_P_F) );

    target += poisson_lpmf(Y | append_row(mu[M], mu[F]));
  }
}

generated quantities {
  array[N] real log_lik;    // For LOO
  array[G,A*B] int Yhat;  // For posterior predictive checks
  array[G] vector[B] AMR;        // Smoothed age-standardised mortality rates

  { // Local scope
    for (b in 1:B){
      AMR[M][b] = W * exp(log_lambda[M][:,b]);
      AMR[F][b] = W * exp(log_lambda[F][:,b]);
    }

    array[G] vector[A*B] mu;
    mu[M] = exp( to_vector(log_lambda[M] + LOG_P_M) );
    mu[F] = exp( to_vector(log_lambda[F] + LOG_P_F) );

    Yhat[M] = poisson_rng(mu[M]);
    Yhat[F] = poisson_rng(mu[F]);

    vector[N] mu_long = append_row(mu[M], mu[F]);

    for (i in 1:N){
      log_lik[i] = poisson_lpmf(Y[i] | mu_long[i]);
    }
  }
}




