context("ZamCovid")


test_that("can run ZamCovid model", {
  start_date <- numeric_date("2020-02-01")
  p <- ZamCovid_parameters(start_date)
  mod <- ZamCovid$new(p, 0, 5, seed = 1L)
  end <- numeric_date("2020-09-30") / p$dt

  info <- mod$info()
  initial <- ZamCovid_initial(info, 10, p)
  mod$update_state(state = initial)

  res <- mod$run(end)

  expect_equal(c(672, 5), dim(res))

  ## TODO: add test with expected space observations given baseline parameters
})


test_that("can seed infections", {

  ## See one big lump of infections
  start_date <- numeric_date("2020-02-01")
  n_particles <- 10
  p <- ZamCovid_parameters(start_date)
  mod <- ZamCovid$new(p, 4, n_particles, seed = 1L)
  end <- numeric_date("2020-02-28") / p$dt

  initial <- ZamCovid_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)

  t <- seq(4, end)
  res <- mod$simulate(t)

  info <- mod$info()
  n_E <- res[info$index$E, , ]
  i <- apply(n_E[4, , ] > 0, 1, function(x) min(which(x)))
  expect_equal(t[i],
               rep(start_date * p$steps_per_day + 1, n_particles))

  # Total infections through seeding are plausible
  n <- mean(n_E[4, , i[[1]]])
  expect_gt(ppois(n, 10), 0.05)

  # No infections in any other age group before seed
  expect_true(all(n_E[-4, , seq_len(i[[1]])] == 0))


  ## Now seed infections over a vector of steps
  start_date <- numeric_date("2020-02-01") + 0.123
  n_particles <- 20
  pattern <- rep(1, 4) # over a 1 day window
  p <- ZamCovid_parameters(start_date,
                           initial_seed_size = 10,
                           initial_seed_pattern = pattern)

  expect_equal(p$seed_step_start, 128)
  expect_equal(p$seed_value, c(1.27, 2.5, 2.5, 2.5, 1.23))

  mod <- ZamCovid$new(p, 4, n_particles, seed = 1L)
  end <- numeric_date("2020-02-28") / p$dt

  initial <- ZamCovid_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)

  t <- seq(4, end)
  res <- mod$simulate(t)

  info <- mod$info()
  n_E <- res[info$index$E, , ]
  i <- apply(n_E[4, , ] > 0, 1, function(x) min(which(x)))
  expect_equal(min(t[i]), floor(start_date * p$steps_per_day) + 1)
  expect_gte(diff(range(t[i])), 1)
})


test_that("infections spread to other age groups after seed", {

  start_date <- numeric_date("2020-02-01") + 0.123
  n_particles <- 10
  pattern <- rep(1, 4) # over a 1 day window
  p <- ZamCovid_parameters(start_date,
                           initial_seed_size = 10,
                           initial_seed_pattern = pattern)

  mod <- ZamCovid$new(p, 4, n_particles, seed = 1L)
  end <- numeric_date("2020-09-30") / p$dt

  initial <- ZamCovid_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)

  t <- seq(4, end)
  res <- mod$simulate(t)

  info <- mod$info()
  n_E <- res[info$index$E, , ]

  i_donor <- apply(n_E[4, , ] > 0, 1, function(x) min(which(x)))
  i_recipient <- NULL
  for (i in c(seq(1, 3, 1), seq(5, 16))) {
    tmp <- suppressWarnings(
      apply(n_E[i, , ] > 0, 1, function(x) min(which(x)))
    )
    i_recipient <- rbind(i_recipient, tmp)
  }

  expect_true(all(i_recipient > mean(i_donor)))
})


test_that("people sero-convert", {
  start_date <- numeric_date("2020-02-01") + 0.123
  n_particles <- 5
  pattern <- rep(1, 4) # over a 1 day window
  p <- ZamCovid_parameters(start_date,
                           initial_seed_size = 10,
                           initial_seed_pattern = pattern)

  mod <- ZamCovid$new(p, 4, n_particles, seed = 1L)
  end <- numeric_date("2020-09-30") / p$dt

  initial <- ZamCovid_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)

  t <- seq(4, end)
  res <- mod$simulate(t)

  info <- mod$info()

  sero_pos <- sum(res[info$index$sero_pos_over15, , dim(res)[3]])

  expect_true(sero_pos > 0)
})


test_that("can run the particle filter on the model", {

  start_date <- numeric_date("2020-02-01")
  pars <- ZamCovid_parameters(start_date)

  data <- read_csv(ZamCovid_file("extdata/example.csv"))
  data <- helper_data(data, 0, pars$dt)

  pf <- helper_particle_filter(data, 1)
  expect_s3_class(pf, "particle_filter")

  pf$run(pars)
})


test_that("can run basic model", {

  p <- basic_parameters(model_end = 200)
  mod <- basic$new(p, 0, 5)
  info <- mod$info()

  y <- mod$run(200)

  expect_equal(c(104, 5), dim(y))
})
