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
  expect_equal(c(932, 5), dim(res))

  index <- ZamCovid_index(info)$run

  mod$set_index(index)
  res <- mod$run(end)

  expected <- rbind(
    time =                  c(274,      274,      274,      274,      274),
    infections_inc =       c(8435,      251,     3298,      228,      252),
    reinfections_inc =      c(587,       20,      242,       11,       12),
    hosp_inc =              c(304,        9,       95,        6,        5),
    admitted_inc =           c(12,        0,        5,        0,        0),
    deaths_hosp_inc =        c(50,        0,       20,        2,        1),
    deaths_comm_inc =      c(1861,       51,      686,       47,       49),
    base_death_inc =          c(0,        0,        0,        0,        0),
    deaths_all_inc =       c(1438,       35,      517,       42,       38),
    deaths_cum_hosp =     c(13701,    14697,    14500,    14709,    14661),
    deaths_cum_comm =    c(970545,   987774,   986790,   990061,   992018),
    sero_pos_all =      c(5183117,  1998387,  4388864,  2039129,  2092334),
    sero_pos_over15 =   c(3031613,  1160838,  2560974,  1184134,  1215783),
    sero_pos_15_19 =     c(619864,   239112,   523476,   242727,   248705),
    sero_pos_20_29 =     c(893945,   345388,   758181,   352520,   361808),
    sero_pos_30_39 =     c(670298,   254398,   565229,   261177,   267346),
    sero_pos_40_49 =     c(436732,   165261,   367256,   167976,   173343),
    sero_pos_50_plus =   c(410774,   156679,   346832,   159734,   164581),
    inf_cum_all =      c(29275139, 29471071, 29439930, 29519482, 29513891),
    inf_cum_over15 =   c(17172239, 17284775, 17264676, 17314559, 17306170),
    inf_cum_15_19 =     c(3499062,  3525771,  3516340,  3531173,  3526723),
    inf_cum_20_29 =     c(5046915,  5080656,  5076318,  5088291,  5084757),
    inf_cum_30_39 =     c(3805080,  3828396,  3823880,  3835279,  3837186),
    inf_cum_40_49 =     c(2482226,  2496447,  2496716,  2500531,  2501333),
    inf_cum_50_plus =   c(2338956,  2353505,  2351422,  2359285,  2356171),
    immune_S_vacc =           c(0,        0,        0,        0,        0),
    immune_R_vacc =           c(0,        0,        0,        0,        0),
    immune_R_unvacc =   c(9652708,  9358293,  9631820,  9381036,  9386076))

  expect_equal(floor(res), floor(expected))
})


test_that("can seed infections", {

  ## See one big lump of infections
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
  expect_true(all(base_deaths == base_death_value))


  ## Can introduce time-varying baseline deaths
  base_death_value <- c(4, 5, 3)
  base_death_date <- c(0, numeric_date(c("2020-06-01", "2020-10-01")))
  p <- ZamCovid_parameters(start_date, base_death_date = base_death_date,
                           base_death_value = base_death_value)
  mod <- ZamCovid$new(p, 0, 5, seed = 1L)
  mod$update_state(state = initial)
  res <- mod$simulate(t)
  base_deaths <- colMeans(res[info$index$base_death_inc, , ])
  expect_true(all(base_deaths %in% base_death_value))

})


test_that("Can calculate YLL", {
  start_date <- numeric_date("2020-02-01")
  p <- ZamCovid_parameters(start_date)
  mod <- ZamCovid$new(p, 0, 5, seed = 1L)
  end <- numeric_date("2020-09-30") / p$dt

  info <- mod$info()
  initial <- ZamCovid_initial(info, 10, p)
  mod$update_state(state = initial)

  index <- ZamCovid_index(info)$state
  mod$set_index(index)

  res <- mod$run(end)

  expect_true(sum(res["yll_tot", ]) > 0)

  expect_true(
    all(round(res["D_65", ] * p$life_exp[14]) == round(res["yll_age_65", ])))

})


test_that("can re-seed other waves", {

  ## Only one wave with initial seed
  start_date <- numeric_date("2020-01-01")
  n_particles <- 10
  p <- ZamCovid_parameters(start_date)
  mod <- ZamCovid$new(p, 4, n_particles)
  end <- numeric_date("2021-12-31") / p$dt

  info <- mod$info()

  initial <- ZamCovid_initial(info, n_particles, p)
  mod$update_state(state = initial)

  t <- seq(4, end, by = 4)
  res <- mod$simulate(t)

  n_E <- res[info$index$E, , ]
  n_E <- apply(n_E, 3, mean)

  expect_true(sum(n_E[360:length(n_E)]) < 1e1)


  ## Can seed a second wave
  p <- ZamCovid_parameters(start_date,
                           re_seed_date = c(0, 360),
                           re_seed_value = c(0, 10))

  mod <- ZamCovid$new(p, 4, n_particles)
  end <- numeric_date("2021-12-31") / p$dt

  initial <- ZamCovid_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)

  t <- seq(4, end, by = 4)
  res <- mod$simulate(t)

  n_E <- res[info$index$E, , ]
  n_E <- apply(n_E, 3, mean)

  expect_true(sum(n_E[360:length(n_E)]) > 1e3)


  ## Can seed multiple waves
  p <- ZamCovid_parameters(start_date,
                           re_seed_date = c(0, 250, 500),
                           re_seed_value = c(0, 10, 10))

  mod <- ZamCovid$new(p, 4, n_particles)
  end <- numeric_date("2021-12-31") / p$dt

  initial <- ZamCovid_initial(mod$info(), n_particles, p)
  mod$update_state(state = initial)

  t <- seq(4, end, by = 4)
  res <- mod$simulate(t)

  n_E <- res[info$index$E, , ]
  n_E <- apply(n_E, 3, mean)

  expect_true(sum(n_E[500:length(n_E)]) > 1e3)
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
