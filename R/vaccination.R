##' Create a historical vaccine schedule from data
##'
##' @title Create historical vaccine schedule
##'
##' @param data A data.frame with columns `date`, `age_band_min`,
##'   and numbered doses columns, e.g. if there are three doses
##'   these should be `dose1`, `dose2` and `dose3`. Values of
##'   `age_band_min` should be either multiples of 5 or NA - the
##'   latter means those doses are not age-specific and will be
##'   distributed across all ages according to priority after
##'   all age-specific doses have already been dealt with
##'
##' @param region Region to use to get total population numbers
##'
##' @param uptake A matrix of 19 rows, and number of columns equal to
##'   number of doses. The (i,j)th entry gives the fractional uptake
##'   of dose j for group i. Should be non-increasing across rows
##'
##' @return A [vaccine_schedule] object
##' @export
vaccine_schedule_from_data <- function(data, region, uptake) {
  assert_is(data, "data.frame")
  dose_cols <- grep("dose[0-9]", names(data), value = TRUE)
  n_doses <- length(dose_cols)

  required <- c("age_band_min", "date")
  msg <- setdiff(required, names(data))
  if (length(msg) > 0) {
    stop("Required columns missing from 'data': ",
         paste(squote(msg), collapse = ", "))
  }

  dose_expected <- paste0("dose", seq_len(n_doses))
  dose_msg <- setdiff(dose_expected, names(data))
  if (length(dose_msg) > 0) {
    stop(sprintf("There are %s dose columns so expected dose column names: %s)",
                 n_doses, paste(squote(dose_expected), collapse = ", ")))
  }

  err <- data$age_band_min[!is.na(data$age_band_min)] %% 5 != 0
  if (any(err)) {
    stop("Invalid values for data$age_band_min: ",
         paste(unique(data$age_band_min[err]), collapse = ", "))
  }

  age_start <- sircovid_age_bins()$start
  if (nrow(uptake) != length(age_start) + 2) {
    stop(sprintf("Expected uptake to have %s rows as there are %s groups",
                 length(age_start) + 2, length(age_start) + 2))
  }
  if (ncol(uptake) != n_doses) {
    stop(sprintf("Data has %s dose columns so expected uptake to have %s
                 columns", n_doses, n_doses))
  }

  assert_proportion(uptake)
  if (any(apply(uptake, 1, diff) > 0)) {
    stop("Uptake should not increase with dose number for any group")
  }

  ## TODO: tidy up later:
  stopifnot(!is.na(data$date),
            all(data[, dose_cols] >= 0 | is.na(data[, dose_cols])),
            all(data[, dose_cols] >= 0 | is.na(data[, dose_cols])))

  ## First aggregate all the 80+ into one group
  data$date <- as_sircovid_date(data$date)
  data$age_band_min <- pmin(data$age_band_min, 80)
  data$age_band_min[is.na(data$age_band_min)] <- Inf
  data <- stats::aggregate(data[dose_cols],
                           data[c("age_band_min", "date")],
                           sum)

  dates <- seq(min(data$date), max(data$date), by = 1)

  doses <- lapply(dose_cols, function(i)
    stats::reshape(data[c("date", "age_band_min", i)],
                   direction = "wide", timevar = "date",
                   idvar = "age_band_min"))
  for (i in seq_len(n_doses)) {
    stopifnot(identical(dim(doses[[1]]), dim(doses[[i]])))
  }

  i_agg <- which(doses[[1]]$age_band_min == Inf)
  j <- match(dates, sub("^dose[12]\\.", "", names(doses[[1]])))

  agg_doses <- array(
    unlist(lapply(doses, function(d) unname(as.matrix(d)[i_agg, j]))),
    c(length(dates), n_doses))
  agg_doses <- aperm(agg_doses, c(2, 1))
  agg_doses[is.na(agg_doses)] <- 0

  ## TODO: add a test for missing days
  i <- match(age_start, doses[[1]]$age_band_min)

  doses <- array(
    unlist(lapply(doses, function(d) unname(as.matrix(d)[i, j]))),
    c(length(age_start), length(dates), n_doses))
  doses <- aperm(doses, c(1, 3, 2))
  doses[is.na(doses)] <- 0

  ## We have 16 groups, 12 priority groups
  priority_population <-
    vapply(seq_len(n_doses),
           function(j) vaccine_priority_population(region, uptake[, j]),
           array(0, c(16, 12)))

  ## Now distribute age-aggregated doses (if any) and add in
  if (any(agg_doses > 0)) {
    population_left <- priority_population
    doses_given <- apply(doses, c(1, 2), sum)

    for (j in seq_len(dim(priority_population)[2])) {
      vaccinated <- pmin(population_left[, j, ], doses_given)
      population_left[, j, ] <- population_left[, j, ] - vaccinated
      doses_given <- doses_given - vaccinated
    }

    n_groups <- dim(population_left)[1]
    n_priority_groups <- dim(population_left)[2]
    n_doses <- dim(population_left)[3]
    n_days <- dim(agg_doses)[2]

    population_to_vaccinate_mat <- array(0, c(n_groups, n_priority_groups,
                                              n_doses, n_days))
    daily_doses_tt <- array(0, dim(agg_doses))

    f <- function(dose) {
      vaccination_schedule_exec(daily_doses_tt, agg_doses[dose, ],
                                population_left,
                                population_to_vaccinate_mat,
                                Inf,
                                0L, dose)
    }

    doses <- doses +
      apply(Reduce("+", lapply(seq_len(n_doses), f)), c(1, 3, 4), sum)
  }


  vaccine_schedule(dates[[1]], doses, n_doses)
}


vaccine_priority_population <- function(region, uptake) {
  p <- vaccine_priority_proportion(uptake)
  pop <- ZamCovid_parameters(1)$N_tot

  pop_mat <- matrix(rep(pop, ncol(p)), nrow = nrow(p))
  round(p * pop_mat)
}


vaccine_priority_proportion <- function(uptake) {
  n_groups <- 16
  uptake <- ifelse(length(uptake) == n_groups,
                   uptake, rep_len(uptake, n_groups))

  n_priority_groups <- 12
  p <- matrix(0, n_groups, n_priority_groups)

  ## Aged base priority list
  jcvi_priority <- list(
    ## the age groups targeted in each priority group (see comments above)
    18:19, 17, 16, 15, 14, NULL, 13, 12, 11, 9:10, 7:8, 1:6)

  ## 1. Start with non aged based priority:
  ## helper function
  add_prop_to_vacc <- function(j, idx, prop_to_vaccinate, p) {
    p[idx, j] <- prop_to_vaccinate[idx]
    p
  }

  ## Group 2 includes frontline health and social care workers
  p <- add_prop_to_vacc(j = 2,
                        idx = seq_len(jcvi_priority[[2]] - 1),
                        prop_to_vaccinate = prop_hcw,
                        p)

  ## Group 4 includes clinically extremely vulnerable individuals
  p <- add_prop_to_vacc(j = 4,
                        idx = seq_len(jcvi_priority[[4]] - 1),
                        prop_to_vaccinate = (1 - prop_hcw) *
                          prop_very_vulnerable,
                        p)

  ## Group 6 includes all individuals aged 16 years to 64 years with
  ## underlying health conditions which put them at higher risk of
  ## serious disease and mortality
  p <- add_prop_to_vacc(j = 6,
                        idx = 4:13,
                        prop_to_vaccinate = (1 - prop_hcw) *
                          prop_underlying_condition,
                        p)

  ## 2. Add aged base priority
  for (j in seq_along(jcvi_priority)) {
    if (!is.null(jcvi_priority[[j]])) {
      ## discount those already vaccinated as part of non age based priority
      p[jcvi_priority[[j]], j] <- 1 - rowSums(p)[jcvi_priority[[j]]]
    }
  }

  ## 3. Account for uptake
  uptake_mat <- matrix(rep(uptake, n_priority_groups),
                       nrow = n_groups)

  p * uptake_mat
}
