# Load required libraries
library(optparse)

option_list <- list(
  make_option("--seed", type = "integer", default = 123,
              help = "Random number seed [default %default]",
              dest = "seed"),
  make_option("--iter_warmup", type = "integer", default = 1000,
              help = "HMC warmup iterations [default %default]",
              dest = 'iter.warmup'),
  make_option("--iter_sampling", type = "integer", default = 2500,
              help = "HMC of sampling iterations iterations [default %default]",
              dest = 'iter.sampling'),
  make_option("--n_chains", type = "integer", default = 4,
              help = "Number of MCMC chains",
              dest = 'n.chains'),
  make_option(c("-d", "--data"), type = "character", default = NA,
              help = "The data to run the models on",
              dest = "data"),
  make_option(c("-m", "--model"), type = "character", default = "hsgp-nb-eq.stan",
              help = "Name of Stan model file",
              dest = 'model'),
  make_option("--repo_path", type = "character", default = "~/Stanford/CARE 2022/KA-mortality-t10",
              help = "Absolute file path to repository directory",
              dest = 'repo.path'),
  make_option("--out_path", type = "character", default = NA,
              help = "Path to the output directory",
              dest = "out.path")
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

renv::activate(args$repo.path)

library(readr)
library(stringr)
library(data.table)
library(cmdstanr)
library(posterior)

cat(paste("========== Running Analysis for:", args$data, "==========\n"))

##### Configurations #####
out.path <- file.path(args$out.path, str_remove(args$data, ".rds"))
if (!file.exists(out.path)) {
  cat(paste(" Making export directory:", out.path, "\n"))
  dir.create(out.path, recursive = TRUE)
  dir.create(file.path(out.path, "stan_fits"))
  dir.create(file.path(out.path, "diagnostics"))
  dir.create(file.path(out.path, "results"))
}

# Helpers
source(file.path(args$repo.path, "R", "stan-utility.R"))
source(file.path(args$repo.path, "R", "postproc-utility.R"))

# Load data
dt <- as.data.table(readRDS(file.path(args$repo.path, "data", args$data)))
dt.pop <- unique(dt[,.(YEAR, SEX, AGEGRP, POPULATION, POP_SMOOTH, WEIGHT)])

##### PREPROCESS DATA #####
dt <- melt(dt[!is.na(COUNT)],
           id.vars = c("YEAR", "SEX", "AGEGRP"),
           measure.vars = patterns("COUNT_"),
           variable.name = "CAUSE", value.name = "COUNT")

dt[is.na(COUNT), COUNT := 0]
dt[, CAUSE := str_remove(CAUSE, "COUNT_")]
dt[, CAUSE := str_to_title(CAUSE)]

dt <- merge(dt, dt.pop, by = c("YEAR", "SEX", "AGEGRP"), all.x = TRUE)

age.groups <- c("0-4", "5-9", "10-14", "15-19",
                "20-24", "25-29", "30-34", "35-39",
                "40-44", "45-49", "50-54", "55-59", "60-64",
                "65-69", "70-74", "75-79", "80-84", "85 +")

dt[, AGEGRP := factor(AGEGRP, levels = age.groups)]
setkey(dt, YEAR, SEX, AGEGRP)

# Calculate crude age-standardised mortality rates
dt.camr <- dt[!is.na(COUNT), .(AMR = sum(COUNT/POPULATION*WEIGHT) * 100000),
              by = .(YEAR, SEX, CAUSE)]

##### Compile and run Stan models #####
model <- cmdstan_model(file.path(args$repo.path, "stan_models", args$model), force_recompile = TRUE)

model_inits <- function(){
  list(beta0 = -log(mean(dt$POPULATION)) + rnorm(1),
       z = rnorm(1, 0, 0.1),
       theta = 0.4 + rnorm(1, 0, 0.001))
}

print(model_inits())

causes_of_death <- unique(dt$CAUSE)
for (cause in causes_of_death){ # loop over all causes of deaths
  cat(paste("Running Stan model for:", cause, "\n"))

  stan_data <- make_stan_data(dt, cause)

  fit <- model$sample(data = stan_data,
                      seed = args$seed,
                      chains = args$n.chains,
                      parallel_chains = args$n.chains,
                      iter_warmup = args$iter.warmup,
                      iter_sampling = args$iter.sampling,
                      adapt_delta = 0.9,
                      max_treedepth = 13,
                      init = model_inits)
  # Save stan model
  fit$save_object(file = file.path(file.path(out.path, "stan_fits"), paste0(cause, "-", args$model.name, ".rds")))

  # Configure output paths
  path.diagnostics <- file.path(out.path, "diagnostics", cause)
  path.results <- file.path(out.path, "results", cause)
  if(!file.exists(path.diagnostics)){
    dir.create(path.diagnostics)
  }
  if(!file.exists(path.results)){
    dir.create(path.results)
  }

  # Diagnostic stats
  dt.dst <- diagnostic_stats(fit, path.diagnostics)

  # Compute LOO
  dt.loo <- summarise_loo(fit, cause, path.diagnostics)

  # Posterior predictive checks
  dt.ppc <- ppc(fit, dt[CAUSE == cause], path.diagnostics)

  # Posterior MR
  dt.mr <- posterior_mr(fit, dt[CAUSE == cause], path.results)

  # Posterior AMR
  dt.amr <- posterior_amr(fit, dt.camr[CAUSE == cause], path.results)
}

cat("========== DONE! ==========\n")


