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
