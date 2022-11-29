make_stan_data <- function(data, cause){
  age.groups <- c("0-4", "5-9", "10-14", "15-19",
                  "20-24", "25-29", "30-34", "35-39",
                  "40-44", "45-49", "50-54", "55-59", "60-64",
                  "65-69", "70-74", "75-79", "80-84", "85 +")

  # Make Grid
  dt.grid <- as.data.table(
    expand.grid(
      YEAR = unique(data$YEAR),
      SEX = c("M", "F"),
      AGEGRP = age.groups
    )
  )
  dt.grid[, YEAR_IDX := 1:.N, by = .(AGEGRP, SEX)]
  setkey(dt.grid, YEAR, SEX, AGEGRP)

  # Make grid
  tmp <- merge(dt.grid, dt[CAUSE == cause], by = c("YEAR", "SEX", "AGEGRP"), all.x = TRUE)
  tmp[, POPULATION := NULL]
  tmp[, POP_SMOOTH := NULL]
  tmp[, WEIGHT := NULL]
  tmp <- merge(tmp, dt.pop, by=c("YEAR", "SEX", "AGEGRP"), all.x = TRUE)
  tmp[, AGEGRP := factor(AGEGRP, levels = age.groups)]
  setkey(tmp, YEAR, SEX, AGEGRP)

  dtm <- tmp[SEX == "M"]
  dtf <- tmp[SEX == "F"]

  dtm[, IDX := 1:.N]
  dtf[, IDX := 1:.N]

  stan_data <- list(
    N_M = nrow(dtm[!is.na(COUNT)]),
    N_F = nrow(dtf[!is.na(COUNT)]),
    A = length(unique(dtm$AGEGRP)),
    B = length(unique(dtm$YEAR)),
    B_PRED = length(seq(1, 16, 1/5)),

    Y_M = dtm[!is.na(COUNT), COUNT],
    Y_F = dtf[!is.na(COUNT), COUNT],
    OBS_IDX_M = dtm[!is.na(COUNT), IDX],
    OBS_IDX_F = dtm[!is.na(COUNT), IDX],

    AGE_IDX = as.numeric(unique(dtm[,AGEGRP])),
    TIME_IDX = unique(dtm[,YEAR_IDX]),
    TIME_PRED_IDX = seq(1,16,1/5),

    P_M = dtm[,POPULATION],
    P_F = dtf[,POPULATION],

    W = dtm[YEAR == 2005, WEIGHT],

    M1 = 20,
    M2 = 25,
    C1 = 2.0,
    C2 = 2.0
  )

  return(stan_data)
}
