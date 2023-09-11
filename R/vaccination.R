##' Create a hypothetical historic vaccine schedule.
##'
##' @title Create historic vaccine schedule
##'
##' @param data A data.frame with columns `date`, `age_band_min`,
##'   and numbered doses columns, e.g. if there are two doses
##'   these should be `dose1` and `dose2`. Values of
##'   `age_band_min` must be multiples of 5 from 0 to 75.
##'
##' @param uptake A matrix of 16 rows, and number of columns equal to
##'   number of doses. The (i,j)th entry gives the fractional uptake
##'   of dose j for group i. Should be non-increasing across rows.
##'
##' @return A [vaccine_schedule] object
##' @export
vaccine_schedule_historic <- function(data = NULL, uptake = NULL,
                                      age_priority = NULL,
                                      pop_to_vaccinate = NULL,
                                      daily_doses_value = NULL,
                                      days_between_doses = NULL,
                                      start = NULL) {

  if (is.null(data)) {
    # If no data provided, we need all these to create hypothetical data
    stopifnot(!is.null(daily_doses_value) && !is.null(pop_to_vaccinate) &&
                !is.null(days_between_doses) && !is.null(start))

    n_groups <- dim(pop_to_vaccinate)[1]
    n_priority_groups <- dim(pop_to_vaccinate)[2]
    n_doses <- 2L
    n_days <- length(daily_doses_value)

    population_to_vaccinate_mat <-
      array(0, c(n_groups, n_priority_groups, n_doses, n_days))


    population_left <- array(rep(c(pop_to_vaccinate), n_doses),
                             c(n_groups, n_priority_groups, n_doses))

    daily_doses_prev <- matrix(0, n_doses, 0)
    # n_prev <- 0L
    daily_doses_date <- start

    daily_doses_tt <- cbind(daily_doses_prev, matrix(0, n_doses, n_days))

    population_to_vaccinate_mat <- vaccination_schedule_exec(
      daily_doses_tt, daily_doses_value, population_left,
      population_to_vaccinate_mat, days_between_doses, 1:2)

    doses <- apply(population_to_vaccinate_mat, c(1, 3, 4), sum)

    vaccine_schedule(daily_doses_date, doses, n_doses)


  } else {
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

    age_start <- model_age_bins()$start
    if (nrow(uptake) != length(age_start)) {
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

    ## First aggregate all the 75+ into one group
    data$date <- numeric_date(data$date)
    data$age_band_min <- pmin(data$age_band_min, 75)
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

    ## We have 16 groups, 10 priority groups
    priority_population <-
      vapply(seq_len(n_doses),
             function(j) vaccine_priority_population(pop, uptake[, j],
                                                     age_priority),
             array(0, c(16, 10)))

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

}


vaccine_schedule <- function (date, doses, n_doses = 2L, n_groups = 16) {
  assert_scalar(date)
  n_groups <- n_groups
  if (length(dim(doses)) != 3L) {
    stop("Expected a 3d array for 'doses'")
  }
  if (nrow(doses) != n_groups) {
    stop(sprintf("'doses' must have %d rows", n_groups))
  }
  if (ncol(doses) != n_doses) {
    stop(sprintf("'doses' must have %d columns", n_doses))
  }
  if (dim(doses)[[3]] == 0) {
    stop("'doses' must have at least one element in the 3rd dimension")
  }
  if (any(is.na(doses))) {
    stop("'doses' must all be non-NA")
  }
  if (any(doses < 0)) {
    stop("'doses' must all be non-negative")
  }
  ret <- list(date = date, doses = doses, n_doses = n_doses)
  class(ret) <- "vaccine_schedule"
  ret
}

vaccine_priority_population <- function(pop, uptake, age_priority) {
  p <- vaccine_priority_proportion(uptake, age_priority)

  pop_mat <- matrix(rep(pop, ncol(p)), nrow = nrow(p))
  round(p * pop_mat)
}


vaccine_priority_proportion <- function(uptake, age_priority) {
  n_groups <- 16
  if(length(uptake) == n_groups) {
    uptake <- uptake
  } else {
    uptake <- rep_len(uptake, n_groups)
  }


  n_priority_groups <- 10
  p <- matrix(0, n_groups, n_priority_groups)

  ## 2. Add aged base priority
  for (j in seq_along(age_priority)) {
    if (!is.null(age_priority[[j]])) {
      ## discount those already vaccinated as part of non age based priority
      p[age_priority[[j]], j] <- 1 - rowSums(p)[age_priority[[j]]]
    }
  }

  ## 3. Account for uptake
  uptake_mat <- matrix(rep(uptake, n_priority_groups),
                       nrow = n_groups)

  p * uptake_mat
}


vaccination_schedule_exec <- function(daily_doses_tt, daily_doses_value,
                                      population_left,
                                      population_to_vaccinate_mat,
                                      mean_days_between_doses, dose_index){
  n_priority_groups <- ncol(population_left)
  i1 <- dose_index[[1]]
  i2 <- if (length(dose_index) == 2) dose_index[[2]] else dose_index[[1]]
  for (t in seq_along(daily_doses_value)) {
    tt <- t #+ n_prev
    tt_dose_1 <- tt - mean_days_between_doses
    if (tt_dose_1 >= 1) {
      ## If we have promised more 2nd doses than we can deliver, we
      ## move our debt forward in time by one day. If doses fluctuate
      ## this will eventually be paid off.
      if (daily_doses_tt[i1, tt_dose_1] > daily_doses_value[t]) {
        daily_doses_tt[i1, tt_dose_1 + 1] <-
          daily_doses_tt[i1, tt_dose_1 + 1] +
          (daily_doses_tt[i1, tt_dose_1] - daily_doses_value[t])
      }
      daily_doses_tt[i2, tt] <- min(daily_doses_value[t],
                                    daily_doses_tt[i1, tt_dose_1])
      daily_doses_tt[i1, tt] <- daily_doses_value[t] - daily_doses_tt[i2, tt]
    } else {
      ## Only distribute first doses
      daily_doses_tt[i2, tt] <- 0
      daily_doses_tt[i1, tt] <- daily_doses_value[t]
    }
    daily_doses_today <- daily_doses_tt[, tt]

    for (dose in dose_index) {
      eligible <- colSums(population_left[, , dose])
      ## Vaccinate the entire of the top priority groups
      n_full_vacc <- findInterval(daily_doses_today[dose], cumsum(eligible))
      if (n_full_vacc > 0) {
        i_full_vacc <- seq_len(n_full_vacc)
        population_to_vaccinate_mat[, i_full_vacc, dose, t] <-
          population_left[, i_full_vacc, dose]
      }

      ## Then partially vaccinate the next priority group, if possible
      if (n_full_vacc < n_priority_groups) {
        if (n_full_vacc == 0) {
          remaining_eligible <- daily_doses_today[dose]
        } else {
          remaining_eligible <- daily_doses_today[dose] -
            cumsum(eligible)[n_full_vacc]
        }
        i_vacc <- n_full_vacc + 1L

        ## Split remaining doses according to age
        population_to_vaccinate_mat[, i_vacc, dose, t] <-
          round(remaining_eligible * population_left[, i_vacc, dose] /
                  sum(population_left[, i_vacc, dose]))
      }

      population_left[, , dose] <- population_left[, , dose] -
        population_to_vaccinate_mat[, , dose, t]
    }
  }
  population_to_vaccinate_mat
}


create_index_dose_inverse <- function(n_vacc_classes, index_dose) {
  index_dose_inverse <- integer(n_vacc_classes)
  index_dose_inverse[index_dose] <- seq_along(index_dose)
  index_dose_inverse
}


build_vaccine_progression_rate <- function(vaccine_progression_rate,
                                           n_vacc_classes,
                                           index_dose) {

  n_groups <- 16L
  # if NULL, set vaccine_progression_rate to 0
  if (is.null(vaccine_progression_rate)) {
    mat_vaccine_progression_rate <- matrix(0, n_groups, n_vacc_classes)
  } else {
    if (is.matrix(vaccine_progression_rate)) {
      if (nrow(vaccine_progression_rate) != n_groups) {
        stop(
          "'vaccine_progression_rate' must have as many rows as age groups")
      }
      if (ncol(vaccine_progression_rate) != n_vacc_classes) {
        stop(
          "'vaccine_progression_rate' must have 'n_vacc_classes' columns")
      }
      if (any(vaccine_progression_rate < 0)) {
        stop("'vaccine_progression_rate' must have only non-negative values")
      }
      mat_vaccine_progression_rate <- vaccine_progression_rate
    } else { # vaccine_progression_rate vector of length n_vacc_classes
      if (!is.vector(vaccine_progression_rate) ||
          length(vaccine_progression_rate) != n_vacc_classes) {
        m1 <- "'vaccine_progression_rate' must be either:"
        m2 <- "a vector of length 'n_vacc_classes' or"
        m3 <- "a matrix with 'n_groups' rows and 'n_vacc_classes' columns"
        stop(paste(m1, m2, m3))
      }
      if (any(vaccine_progression_rate < 0)) {
        stop("'vaccine_progression_rate' must have only non-negative values")
      }
      # create matrix by repeating vaccine_progression_rate for each age group
      mat_vaccine_progression_rate <-
        matrix(rep(vaccine_progression_rate, each = n_groups), nrow = n_groups)
    }
  }
  for (i in seq_along(index_dose)) {
    j <- index_dose[[i]]
    if (!all(mat_vaccine_progression_rate[, j] == 0)) {
      stop(sprintf(
        "Column %d of 'vaccine_progression_rate' must be zero (dose %d)",
        j, i))
    }
  }
  mat_vaccine_progression_rate
}
