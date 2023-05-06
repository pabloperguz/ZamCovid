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
  expect_equal(c(725, 5), dim(res))

  index <- ZamCovid_index(info)$run

  mod$set_index(index)
  res <- mod$run(end)

  expected <- rbind(
    time =             c(274,     274,  274,     274,     274),
    admitted_inc =     c(0,       1,    0,       1,       0),
    deaths_hosp_inc =  c(0,       1,    0,       2,       2),
    deaths_comm_inc =  c(39,      54,    0,      25,      81),
    base_death_inc =   c(0,       0,    0,       0,       0),
    sero_pos_all =     c(2165173, 2149797,    0, 1587775, 2526417),
    sero_pos_over15 =  c(1258374, 1250039,    0,  922442, 1469582),
    sero_pos_15_19 =   c(258108,  255407,    0,  188805,  300134),
    sero_pos_20_29 =   c(374755,  372427,    0,  274911,  437320),
    sero_pos_30_39 =   c(276178,  274552,    0,  202239,  322025),
    sero_pos_40_49 =   c(178880,  178160,    0,  131595,  210484),
    sero_pos_50_plus = c(170453,  169493,    0,  124892,  199619))

  expect_equal(res, expected)
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


test_that("Can model baseline deaths", {
  start_date <- numeric_date("2020-01-01")
  p <- ZamCovid_parameters(start_date)
  mod <- ZamCovid$new(p, 0, 5, seed = 1L)
  end <- numeric_date("2020-12-31") / p$dt

  info <- mod$info()
  initial <- ZamCovid_initial(info, 10, p)
  mod$update_state(state = initial)

  t <- seq(4, end)
  res <- mod$simulate(t)

  ## No baseline deaths with default parameters
  base_deaths <- res[info$index$base_death_inc, , ]
  expect_true(sum(base_deaths) == 0)

  ## Can introduce constant value of baseline deaths
  base_death_value <- 4
  base_death_date <- 0
  p <- ZamCovid_parameters(start_date, base_death_date = base_death_date,
                           base_death_value = base_death_value)
  mod <- ZamCovid$new(p, 0, 5, seed = 1L)
  mod$update_state(state = initial)
  res <- mod$simulate(t)
  base_deaths <- colMeans(res[info$index$base_death_inc, , ])
  expect_true(all(base_deaths == 1))


  ## Can introduce time-varying baseline deaths
  base_death_value <- c(4, 5, 3)
  base_death_date <- c(0, numeric_date(c("2020-06-01", "2020-10-01")))
  p <- ZamCovid_parameters(start_date, base_death_date = base_death_date,
                           base_death_value = base_death_value)
  mod <- ZamCovid$new(p, 0, 5, seed = 1L)
  mod$update_state(state = initial)
  res <- mod$simulate(t)
  base_deaths <- colMeans(res[info$index$base_death_inc, , ])
  expect_true(all(base_deaths %in% c(4 / 4, 5 / 4, 3 / 4)))

})


test_that("can run the particle filter on the model", {

  set.seed(2)
  start_date <- numeric_date("2020-02-01")
  pars <- ZamCovid_parameters(start_date)

  data <- read_csv(ZamCovid_file("extdata/example.csv"))
  data <- helper_data(data, 0, pars$dt)

  np <- 10
  pf <- helper_particle_filter(data, np, seed = 2)
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
