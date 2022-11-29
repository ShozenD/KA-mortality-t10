#' Make convergence diagnostic statistics
#'
#' @param fit A fitted CmdStanModel
#' @param outdir Directory to save outputs
#'
#' @return A summary table of estimates and diagnostics
diagnostic_stats <- function(fit, outdir = NA) {

  # Rhat and effective sample size
  fit_summary <- fit$summary(variables = NULL, "rhat", "ess_bulk")

  cat("\n The min and max ESS are: ", range(fit_summary$ess_bulk, na.rm = TRUE))
  cat("\n The min and max Rhat are: ", range(fit_summary$rhat, na.rm = TRUE))
  if(min(fit_summary$ess_bulk, na.rm = TRUE) < 400) cat('\n Minimum effective sample size smaller than 400\n')

  # Diagnostics
  sampler_diagnostics <- fit$diagnostic_summary()

  # Time of execution
  time <- fit$time()

  # save
  if(!is.na(outdir)){
    saveRDS(fit_summary, file = file.path(outdir, "fit_summary.rds"))
    saveRDS(sampler_diagnostics, file = file.path(outdir, "sampler_diagnostics.rds"))
    saveRDS(time, file = file.path(outdir, "time_elapsed.rds"))
  } else {
    warning("\n outdir is not given. Results were not saved.")
  }

  return(fit_summary)
}

summarise_loo <- function(fit, cause = NA, outdir = NA){
  loo <- fit$loo()
  dt.loo <- data.table(loo$estimates, keep.rownames = TRUE)
  dt.loo$cause <- cause
  colnames(dt.loo) <- toupper(colnames(dt.loo))

  if(!is.na(outdir)){
    saveRDS(dt.loo, file = file.path(outdir, "loo.rds"))
  } else {
    warning("\n outdir is not given. Results were not saved.")
  }

  return(dt.loo)
}

ppc <- function(fit, data, outdir = NA){
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  age.groups <- c("0-4", "5-9", "10-14", "15-19",
                  "20-24", "25-29", "30-34", "35-39",
                  "40-44", "45-49", "50-54", "55-59", "60-64",
                  "65-69", "70-74", "75-79", "80-84", "85 +")

  po <- fit$draws(variables = "Yhat", format = "draws_matrix")
  dt.po <- as.data.table(reshape2::melt(po))

  .pattern <- "Yhat\\[([0-9]+),([0-9]+)\\]"
  dt.po[, SEX := as.numeric(gsub(.pattern, "\\1", variable))]
  dt.po[, IDX := as.numeric(gsub(.pattern, "\\2", variable))]

  dt.po <- dt.po[, list( q=quantile(value, prob=ps, na.rm=T), q_label = p_labs), by=.(IDX, SEX)]
  dt.po <- data.table::dcast(dt.po, SEX + IDX ~ q_label, value.var = "q")

  dt.po[, SEX := ifelse(SEX == 1, "M", "F")]

  dt.po[, AGEGRP := age.groups[(IDX - 1) %% 18 + 1]]
  dt.po[, YEAR := IDX %/% 18 + 2005]

  dt.ppc <- merge(data, dt.po, by=c("YEAR", "SEX", "AGEGRP"))
  dt.ppc[, AGEGRP := factor(AGEGRP, levels = age.groups)]
  dt.ppc[, inside.CI := COUNT >= CL & COUNT <= CU]
  cat(paste(" Pct of points in 95% posterior predictive interval:", mean(dt.ppc$inside.CI)))

  if(!is.na(outdir)){
    saveRDS(dt.ppc, file = file.path(outdir, "ppc.rds"))
  } else {
    warning("\n outdir is not given. Results were not saved.")
  }

  return(dt.ppc)
}

posterior_mr <- function(fit, data, outdir = NA){
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c("M", "CL", "CU")
  age.groups <- c("0-4", "5-9", "10-14", "15-19",
                  "20-24", "25-29", "30-34", "35-39",
                  "40-44", "45-49", "50-54", "55-59", "60-64",
                  "65-69", "70-74", "75-79", "80-84", "85 +")

  po <- fit$draws("MR", format = "draws_matrix")

  dt.po <- as.data.table(reshape2::melt(po))
  .pattern <- "MR\\[([0-9]+),([0-9]+),([0-9]+)\\]"
  dt.po[, SEX := as.numeric(gsub(.pattern, "\\1", variable))]
  dt.po[, AGEGRP := as.numeric(gsub(.pattern, "\\2", variable))]
  dt.po[, YEAR := as.numeric(gsub(.pattern, "\\3", variable))]

  # Summarise draws
  dt.po <- dt.po[, list( q=quantile(value*1e5, prob=ps, na.rm=T), q_label = p_labs), by=.(YEAR, SEX, AGEGRP)]
  dt.po <- data.table::dcast(dt.po, YEAR + SEX + AGEGRP ~ q_label, value.var = "q")

  # Recover year, sex, and age groups
  dt.po[, YEAR := (YEAR - 1) / 5 + 2005]
  dt.po[, SEX := ifelse(SEX == 1, "M", "F")]
  dt.po[, AGEGRP := age.groups[AGEGRP]]

  dt.mr <- merge(dt.po, data, by = c("YEAR", "SEX", "AGEGRP"), all.x = TRUE)
  dt.mr[is.na(CAUSE), CAUSE := unique(data[!is.na(CAUSE), CAUSE])]
  setcolorder(dt.mr, c("YEAR", "SEX", "AGEGRP", "CAUSE"))

  if(!is.na(outdir)){
    saveRDS(dt.mr, file = file.path(outdir, "mortality-rates.rds"))
  } else {
    warning("\n outdir is not given. Results were not saved.")
  }

  return(dt.mr)
}

posterior_amr <- function(fit, data, outdir = NA){
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c("M", "CL", "CU")

  po <- fit$draws("AMR", format = "draws_matrix")

  dt.po <- as.data.table(reshape2::melt(po))
  .pattern <- "AMR\\[([0-9]+),([0-9]+)\\]"
  dt.po[, SEX := as.numeric(gsub(.pattern, "\\1", variable))]
  dt.po[, YEAR := as.numeric(gsub(.pattern, "\\2", variable))]

  # Summarise draws
  dt.po <- dt.po[, list( q=quantile(value*1e5, prob=ps, na.rm=T), q_label = p_labs), by=.(YEAR, SEX)]
  dt.po <- data.table::dcast(dt.po, YEAR + SEX ~ q_label, value.var = "q")

  # Recover year and sex
  dt.po[, YEAR := (YEAR - 1) / 5 + 2005]
  dt.po[, SEX := ifelse(SEX == 1, "M", "F")]

  dt.amr <- merge(dt.po, data, by = c("YEAR", "SEX"), all.x = TRUE)
  dt.amr[is.na(CAUSE), CAUSE := unique(data[!is.na(CAUSE), CAUSE])]
  setcolorder(dt.amr, c("YEAR", "SEX", "CAUSE"))

  if(!is.na(outdir)){
    saveRDS(dt.amr, file = file.path(outdir, "adjusted-mortality-rates.rds"))
  } else {
    warning("\n outdir is not given. Results were not saved.")
  }

  return(dt.amr)
}
