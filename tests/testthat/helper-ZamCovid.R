helper_data <- function(data, start_date, dt) {
  assert_numeric_date(start_date)
  data$date <- numeric_date(data$date)

  expected <- c(hosp_admissions = NA_real_, deaths_hosp = NA_real_,
                sero_pos_all = NA_real_, deaths_all = NA_real_,
                sero_pos_over15 = NA_real_, sero_tot_over15 = NA_real_,
                sero_pos_15_19 = NA_real_, sero_tot_15_19 = NA_real_,
                sero_pos_20_29 = NA_real_, sero_tot_20_29 = NA_real_,
                sero_pos_30_39 = NA_real_, sero_tot_30_39 = NA_real_,
                sero_pos_40_49 = NA_real_, sero_tot_40_49 = NA_real_,
                sero_pos_50_plus = NA_real_, sero_tot_50_plus = NA_real_)

  msg <- setdiff(names(expected), names(data))
  if (length(msg) > 0) {
    message(sprintf("Adding empty %s to data: %s",
                    ngettext(length(msg), "column", "columns"),
                    paste(squote(msg), collapse = ", ")))
    for (nm in msg) {
      data[[nm]] <- expected[[nm]]
    }
  }

  rate <- 1 / dt
  data <- mcstate::particle_filter_data(data, "date", rate, start_date)
  data
}


helper_particle_filter <- function(data, n_particles, n_threads = 1L,
                                   seed = NULL) {
  ZamCovid_check_data(data)
  mcstate::particle_filter$new(
    data,
    ZamCovid,
    n_particles,
    ZamCovid_compare,
    ZamCovid_index,
    ZamCovid_initial,
    n_threads = n_threads,
    seed = seed)
}


test_vaccine_schedule <- function(daily_doses = 250, n_days = 100,
                                  dose_waste = 0.1, start = 0,
                                  days_between_doses = 12 * 7,
                                  uptake = NULL, age_priority = NULL,
                                  population = NULL) {
  if (is.null(uptake)) {
    uptake <- c(rep(0, 3), # no vaccination in <15
                (2 / 5) * 0.75, # no vaccination in 15-17yo
                rep(0.75, 12))
  }
  if (is.null(age_priority)) {
    age_priority <- list(16, 15, 14, 13, 12, 11, 9:10, 7:8, 1:6)
  }
  if (is.null(population)) {
    population <- c(39460, 36388, 32525, 28413, 23396, 19197, 15983, 13281,
                    10824, 8111, 5807, 4388, 3181, 2227, 1510, 1512)
  }

  vaccine_schedule_historic(data = NULL,
                            pop_to_vaccinate = population, uptake = uptake,
                            age_priority = age_priority,
                            daily_doses_value = daily_doses,
                            days_between_doses = days_between_doses,
                            start = start, n_days = n_days)

}


expect_vector_equal <- function(x, y, digits = 100, tol = 0) {
  if (is.numeric(x) && is.numeric(y)) {
    expect_true(all(abs(round(x, digits) - round(y, digits)) <= tol),
                sprintf("\nNot all '%s' equal to '%s' (tol %s)",
                        deparse(substitute(x)), deparse(substitute(y)), tol))
  } else {
    expect_true(all(x == y))
  }
}


expect_approx_equal <- function(x1, x2, rel_tol = 0.05) {
  x1_zeros <- x1 == 0
  x2_zeros <- x2 == 0
  expect_true(all(abs(x1[!x1_zeros] - x2[!x1_zeros]) / x1[!x1_zeros] < rel_tol))
  expect_true(all(abs(x1[x1_zeros & !x2_zeros] - x2[x1_zeros & !x2_zeros]) /
                    x2[x1_zeros & !x2_zeros] < rel_tol))
}
