context("ZamCovid")

test_that("ZamCovid severity pathways work as expected", {

  start_date <- numeric_date("2020-02-01")
  n_particles <- 10

  ## With default severity
  data <- read_csv(ZamCovid_file("extdata/severity_default.csv"))
  severity <- ZamCovid_parameters_severity(1, data)
  p <- ZamCovid_parameters(start_date, 1L, severity = severity)

  mod <- ZamCovid$new(p, 0, n_particles)
  initial <- ZamCovid_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)
  info <- mod$info()
  index <- ZamCovid_index(info)

  end_date <- numeric_date("2020-09-30") / p$dt

  # Check all infected and severity pathway compartments get populated
  t <- seq(4, end_date)
  res <- mod$simulate(t)

  comp <- c("E", "I_A", "I_P", "I_C_1", "I_C_2",
            "H_R_conf", "H_R_unconf", "H_D_conf", "H_D_unconf",
            "G_D", "R", "D", "D_hosp", "D_non_hosp")
  for (i in comp) {
    tmp <- res[info$index[[i]], , ]
    for (j in seq_along(nrow(tmp))) {
      chk <- apply(tmp[j, , ], 1, sum)
      expect_true(mean(chk) > 0,
                  label = paste0("Condition not TRUE for ", i, "_", j))
    }
  }

  ## Zero p_H_D and check people go into but don't die
  sev <- severity
  dim <- dim(sev$p_H_step)
  sev$p_H_D_step <- array(0, dim)
  p <- ZamCovid_parameters(start_date, 1L, severity = sev)
  mod <- ZamCovid$new(p, 0, n_particles)
  initial <- ZamCovid_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)
  res <- mod$simulate(t)

  hosp_D <- c("H_D_conf", "H_D_unconf", "D_hosp")
  for (i in comp) {
    tmp <- res[info$index[[i]], , ]
    for (j in seq_along(nrow(tmp))) {
      chk <- apply(tmp[j, , ], 1, sum)

      if (i %in% hosp_D) {
        expect_true(all(chk == 0),
                    label = paste0("hosp not TRUE for ", i, "_", j))
      } else {
        expect_true(mean(chk) > 0,
                    label = paste0("chk not TRUE for ", i, "_", j))
      }
    }
  }


  ## Zero p_H and check people never go into hospital
  sev <- severity
  dim <- dim(sev$p_H_step)
  sev$p_H_step <- array(0, dim)
  p <- ZamCovid_parameters(start_date, 1L, severity = sev)
  mod <- ZamCovid$new(p, 0, n_particles)
  initial <- ZamCovid_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)
  res <- mod$simulate(t)

  hosp <- c("H_D_conf", "H_D_unconf", "H_R_conf", "H_R_unconf")
  for (i in hosp) {
    tmp <- res[info$index[[i]], , ]
    for (j in seq_along(nrow(tmp))) {
      chk <- apply(tmp[j, , ], 1, sum)
      expect_true(all(chk == 0),
                  label = paste0("hosp not TRUE for ", i, "_", j))

    }
  }


  ## With p_H zeroed, IFR should be equal to p_C * p_G_D
  ifr_input <- as.vector(round(sev$p_C_step * sev$p_G_D_step, 3))
  ifr_output <- round(res[info$index$ifr_age, 6, dim(res)[3]], 3)
  expect_true(all(ifr_input == ifr_output))


  ## Also, all deaths should happen in the community (G_D)
  deaths <- c("G_D", "D_hosp", "D_non_hosp")
  for (i in deaths) {
    tmp <- res[info$index[[i]], , ]
    for (j in seq_along(nrow(tmp))) {
      chk <- apply(tmp[j, , ], 1, sum)

      if (i == "D_hosp") {
        expect_true(all(chk == 0),
                    label = paste0("D_hosp not TRUE for ", j))
      } else {
        expect_true(mean(chk) > 0,
                    label = paste0("chk not TRUE for ", i, "_", j))
      }
    }
  }
})


test_that("ZamCovid healthcare parallel flow works as expected", {

  start_date <- numeric_date("2020-02-01")
  n_particles <- 10

  ## With default severity
  data <- read_csv(ZamCovid_file("extdata/severity_default.csv"))
  severity <- ZamCovid_parameters_severity(1, data)
  p <- ZamCovid_parameters(start_date, 1L, severity = severity)

  mod <- ZamCovid$new(p, 0, n_particles)
  initial <- ZamCovid_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)
  info <- mod$info()
  index <- ZamCovid_index(info)

  end_date <- numeric_date("2020-09-30") / p$dt

  # Check all healthcare parallel flow compartments get populated
  t <- seq(4, end_date)
  res <- mod$simulate(t)

  comp <- c("P_H_sev", "P_H_rec", "P_mild")
  for (i in comp) {
    tmp <- res[info$index[[i]], , ]
    for (j in seq_along(nrow(tmp))) {
      chk <- apply(tmp[j, , ], 1, sum)
      expect_true(mean(chk) > 0,
                  label = paste0("Condition not TRUE for ", i, "_", j))
    }
  }

  ## Zero p_H_sev and check there are no svere cases
  sev <- severity
  dim <- dim(sev$p_sev_step)
  sev$p_sev_step <- array(0, dim)
  p <- ZamCovid_parameters(start_date, 1L, severity = sev)
  mod <- ZamCovid$new(p, 0, n_particles)
  initial <- ZamCovid_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)
  res <- mod$simulate(t)

  for (i in comp) {
    tmp <- res[info$index[[i]], , ]
    for (j in seq_along(nrow(tmp))) {
      chk <- apply(tmp[j, , ], 1, sum)

      if (i == "P_mild") {
        expect_true(mean(chk) > 0,
                    label = paste0("Not TRUE for ", i, "_", j))
      } else {
        expect_true(all(chk == 0),
                    label = paste0("Not TRUE for ", i, "_", j))
      }
    }
  }


  ## Zero gamma_H_sev and check severe cases never recover
  p <- ZamCovid_parameters(start_date, 1L)
  p$gamma_H_sev_step <- 0
  mod <- ZamCovid$new(p, 0, n_particles)
  initial <- ZamCovid_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)
  res <- mod$simulate(t)

  x <- res[info$index$P_H_rec, , ]
  y <- res[info$index$P_H_sev, , ]

  for (i in seq_along(nrow(x))) {
    chk <- apply(x[i, , ], 1, sum)
    expect_true(all(chk == 0),
                label = paste0("Not TRUE for ", i, "_", j))
  }

  for (i in seq_along(nrow(y))) {
    chk <- apply(y[i, , ], 1, sum)
    expect_true(all(chk > 0),
                label = paste0("Not TRUE for ", i, "_", j))
  }
})
