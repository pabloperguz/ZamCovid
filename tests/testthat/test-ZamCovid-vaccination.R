context("ZamCovid (vaccination)")


test_that("Can create a hypothetical vaccination schedule", {

  population <- c(39460, 36388, 32525, 28413, 23396, 19197, 15983, 13281,
                  10824, 8111, 5807, 4388, 3181, 2227, 1510, 1512)
  uptake <- c(rep(0, 3), # no vaccination in <15
              (2 / 5) * 0.75, # no vaccination in 15-17yo
              rep(0.75, 12))
  age_priority <- list(16, 15, 14, 13, 12, 11, 9:10, 7:8, 1:6)
  daily_doses <- rep(250, 100) # a vector of number of doses to give each day
  mean_days_between_doses <- 12 * 7

  n <- vaccine_priority_population(population, uptake, age_priority)

  schedule <- vaccine_schedule_historic(pop_to_vaccinate = n,
                                        daily_doses_value = daily_doses,
                                        days_between_doses =
                                          mean_days_between_doses,
                                        start = 0)

  doses <- schedule$doses
  n_to_vaccinate1 <- schedule$doses[, 1, ]
  n_to_vaccinate2 <- schedule$doses[, 2, ]

  ## initially only first doses are delivered and they add up to daily target
  phase1 <- seq_len(mean_days_between_doses)
  expect_vector_equal(
    colSums(n_to_vaccinate1[, phase1]), daily_doses[phase1])
  expect_vector_equal(colSums(n_to_vaccinate2[, phase1]), 0)

  ## during phase 2 only second doses are delivered
  phase2 <- mean_days_between_doses + seq_len(mean_days_between_doses)
  phase2 <- phase2[which(phase2 <= dim(n_to_vaccinate2)[[2]])]
  expect_vector_equal(
    colSums(n_to_vaccinate2[, phase2]), daily_doses[phase2])
  expect_vector_equal(colSums(n_to_vaccinate1[, phase2]), 0)

  ## check that n_to_vaccinate1 + n_to_vaccinate2 = daily_doses
  ## not exactly equal because of some rounding
  expect_vector_equal(
    colSums(n_to_vaccinate1 + n_to_vaccinate2), daily_doses, tol = 5)

})


test_that("Can vaccinate susceptible individuals", {
  ## Tests that:
  ## Every susceptible moves to vaccinated and stays there if
  ## there is an infinite supply of vaccines and 100% uptake
  population <- data.frame(
    n = c(39460, 36388, 32525, 28413, 23396, 19197, 15983, 13281,
          10824, 8111, 5807, 4388, 3181, 2227, 1510, 1512))

  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            uptake = 1,
                                            days_between_doses = 1000)
  start_date <- numeric_date("2020-02-01")

  p <- ZamCovid_parameters(start_date,
                           population = population,
                           beta_value = c(0, 0, 1),
                           beta_date = c(0, 4, 5),
                           rel_susceptibility = c(1, 0),
                           rel_p_sympt = c(1, 1),
                           rel_p_hosp_if_sympt = c(1, 1),
                           rel_p_death = c(1, 1),
                           vaccine_schedule = vaccine_schedule,
                           vaccine_index_dose2 = 2L)

  mod <- ZamCovid$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  mod$update_state(state = ZamCovid_initial(info, 1, p))

  y <- mod$transform_variables(drop(mod$simulate(seq(0, 400, by = 4))))
  i <- 4:length(model_age_bins()$start)

  ## Predefined schedule means we may end up not vaccinating exactly everyone
  ## hence the approximate comparison
  expect_true(all(abs(y$S[i, 1, 1] - y$S[i, 2, 2]) / y$S[i, 1, 1] < 0.01))
  expect_equal(y$S[i, , 101], y$S[i, , 2])
})


test_that("Vaccination of other eligible compartments works", {
  ## Every exposed moves to vaccinated and stays there if everyone
  ## quickly gets vaccinated with a perfect vaccine.
  population <- data.frame(
    n = c(39460, 36388, 32525, 28413, 23396, 19197, 15983, 13281,
          10824, 8111, 5807, 4388, 3181, 2227, 1510, 1512))
  vaccine_schedule <- test_vaccine_schedule(daily_doses = Inf,
                                            days_between_doses = 1000,
                                            uptake = 1)
  p <- ZamCovid_parameters(0,
                           population = population,
                           rel_susceptibility = c(1, 0),
                           rel_p_sympt = c(1, 1),
                           rel_p_hosp_if_sympt = c(1, 1),
                           rel_p_death = c(1, 1),
                           vaccine_schedule = vaccine_schedule,
                           vaccine_index_dose2 = 2L)

  # Perfect vaccine: stop disease progression after E
  p$gamma_E_step <- 0
  mod <- ZamCovid$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  state <- ZamCovid_initial(info, 1, p)

  index_E <- array(info$index$E, info$dim$E)
  index_S <- array(info$index$S, info$dim$S)
  state[index_E[, , 1]] <- round(state[index_S] / 2)
  state[index_E[, , 2]] <- round(state[index_S] / 2)
  state[index_S] <- 0

  mod$update_state(state = state)
  mod$set_index(info$index$E)
  e <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of e
  expect_equal(length(e), prod(info$dim$E) * 101)
  e <- array(e, c(info$dim$E, 101))

  ## Every E moves from unvaccinated to vaccinated between steps 1 and 2
  E_compartment_idx <- 1
  unvacc_idx <- 1
  vacc_idx <- 2
  i <- 4:length(model_age_bins()$start)

  expect_equal(e[i, unvacc_idx, E_compartment_idx, 1],
               e[i, vacc_idx, E_compartment_idx, 2])
  ## And then they don't move anymore
  expect_equal(e[i, vacc_idx, E_compartment_idx, 2],
               e[i, vacc_idx, E_compartment_idx, 101])
  ## same for second E compartment
  E_compartment_idx <- 2
  expect_equal(e[i, unvacc_idx, E_compartment_idx, 1],
               e[i, vacc_idx, E_compartment_idx, 2])
  expect_equal(e[i, vacc_idx, E_compartment_idx, 2],
               e[i, vacc_idx, E_compartment_idx, 101])



  ## Now the same tests for infected asymptomatic (I_A) ----
  p <- ZamCovid_parameters(0,
                           population = population,
                           rel_susceptibility = c(1, 0),
                           rel_p_sympt = c(1, 1),
                           rel_p_hosp_if_sympt = c(1, 1),
                           rel_p_death = c(1, 1),
                           vaccine_schedule = vaccine_schedule,
                           vaccine_index_dose2 = 2L)

  #stop disease progression after I_A
  p$gamma_A_step <- 0
  mod <- ZamCovid$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  state <- ZamCovid_initial(info, 1, p)

  index_I_A <- array(info$index$I_A, info$dim$I_A)
  index_S <- array(info$index$S, info$dim$S)
  state[index_I_A] <- state[index_S]
  state[index_S] <- 0

  mod$update_state(state = state)
  mod$set_index(info$index$I_A)
  i_A <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of i_A
  expect_equal(length(i_A), prod(info$dim$I_A) * 101)
  i_A <- array(i_A, c(info$dim$I_A, 101))

  ## Every I_A moves from unvaccinated to vaccinated between steps 1 and 2
  expect_equal(i_A[i, unvacc_idx, 1],
               i_A[i, vacc_idx, 2])
  ## then they don't move anymore
  expect_equal(
    i_A[i, vacc_idx, 2],
    i_A[i, vacc_idx, 101])



  ## Now the same tests for infected pre-symptomatic (I_P) ----
  p <- ZamCovid_parameters(0,
                           population = population,
                           rel_susceptibility = c(1, 0),
                           rel_p_sympt = c(1, 1),
                           rel_p_hosp_if_sympt = c(1, 1),
                           rel_p_death = c(1, 1),
                           vaccine_schedule = vaccine_schedule,
                           vaccine_index_dose2 = 2L)

  #stop disease progression after I_P
  p$gamma_P_step <- 0
  mod <- ZamCovid$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  state <- ZamCovid_initial(info, 1, p)

  index_I_P <- array(info$index$I_P, info$dim$I_P)
  index_S <- array(info$index$S, info$dim$S)
  state[index_I_P] <- state[index_S]
  state[index_S] <- 0

  mod$update_state(state = state)
  mod$set_index(info$index$I_P)
  i_P <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of i_P
  expect_equal(length(i_P), prod(info$dim$I_P) * 101)
  i_P <- array(i_P, c(info$dim$I_P, 101))

  ## Every I_P moves from unvaccinated to vaccinated between steps 1 and 2
  expect_equal(i_P[i, unvacc_idx, 1],
               i_P[i, vacc_idx, 2])
  ## then they don't move anymore
  expect_equal(
    i_P[i, vacc_idx, 2],
    i_P[i, vacc_idx, 101])



  ## Now the same tests for infected pre-symptomatic (I_P) ----
  p <- ZamCovid_parameters(0,
                           population = population,
                           rel_susceptibility = c(1, 0),
                           rel_p_sympt = c(1, 1),
                           rel_p_hosp_if_sympt = c(1, 1),
                           rel_p_death = c(1, 1),
                           vaccine_schedule = vaccine_schedule,
                           vaccine_index_dose2 = 2L)

  p$gamma_R <- 0
  mod <- ZamCovid$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- ZamCovid_initial(info, 1, p)

  index_I_A <- array(info$index$I_A, info$dim$I_A)
  index_R <- array(info$index$R, info$dim$R)
  index_S <- array(info$index$S, info$dim$S)
  state[index_R] <- state[index_S]
  state[index_S] <- 0
  state[index_I_A] <- 0 # remove seeded infections

  mod$update_state(state = state)
  mod$set_index(info$index$R)
  r <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of r
  expect_equal(length(r), prod(info$dim$R) * 101)
  r <- array(r, c(info$dim$R, 101))

  ## every R moves from unvaccinated to vaccinated between
  ## time steps 1 and 2
  expect_approx_equal(r[i, unvacc_idx, 1], r[i, vacc_idx, 2])
  ## then they don't move anymore
  expect_equal(r[i, vacc_idx, 2], r[i, vacc_idx, 101])

})


test_that("No infections with perfect vaccine wrt rel_susceptibility", {

  ### No infections with perfect vaccine ----
  p <- ZamCovid_parameters(0,
                           rel_susceptibility = c(1, 0),
                           rel_p_sympt = c(1, 1),
                           rel_p_hosp_if_sympt = c(1, 1),
                           rel_p_death = c(1, 1))

  mod <- ZamCovid$new(p, 0, 1, seed = 1L)
  info <- mod$info()

  state <- ZamCovid_initial(info, 1, p)

  index_S <- array(info$index$S, info$dim$S)
  state[index_S[, 2]] <- state[index_S[, 1]]
  state[index_S[, 1]] <- 0

  mod$update_state(state = state)
  mod$set_index(info$index$S)
  s <- mod$simulate(seq(0, 400, by = 4))

  ## Reshape to show the full shape of s
  expect_equal(length(s), prod(info$dim$S) * 101)
  s <- array(s, c(info$dim$S, 101))

  ## None moves into unvaccinated
  expect_true(all(s[, 1, ] == 0))

  ## None changes compartment within the vaccinated individuals
  expect_true(all(s[, 2, ] == s[, 2, 1]))
})


test_that("Some infections with imperfect vaccine wrt rel_susceptibility", {

  ## Infections occur with a 50% effective vaccine
  p <- ZamCovid_parameters(0,
                           population = population,
                           rel_susceptibility = c(1, 0.5),
                           rel_p_sympt = c(1, 1),
                           rel_p_hosp_if_sympt = c(1, 1),
                           rel_p_death = c(1, 1),
                           rel_infectivity = c(1, 1))

  mod <- ZamCovid$new(p, 0, 1, seed = 1L)
  info <- mod$info()
  state <- ZamCovid_initial(info, 1, p)
  index_S <- array(info$index$S, info$dim$S)

  # Move everyone in S to vaccinated class
  state[index_S[, 2]] <- state[index_S[, 1]] / 2
  state[index_S[, 1]] <- state[index_S[, 1]] / 2

  mod$update_state(state = state)
  mod$set_index(info$index$S)
  y <- mod$simulate(seq(0, 400, 4))

  # Reshape simulated S compartment
  y <- array(y, c(info$dim$S, 101))

  # There have been infections across all age and vaccine classes
  expect_true(all(y[, , 1] > y[, , 101]))

  # More infections in the unvaccinated than in the vaccinated class
  expect_true(sum(y[, 1, 101]) < sum(y[, 2, 101]))



  ## Low protection
  p$rel_susceptibility <- array(rep(c(1, 0.1), each = 16), dim = c(16, 2))
  mod_1 <- ZamCovid$new(p, 0, 1, seed = 1L)

  state_1 <- ZamCovid_initial(info, 1, p)
  state_1[index_S[, 2]] <- state_1[index_S[, 1]] / 2
  state_1[index_S[, 1]] <- state_1[index_S[, 1]] / 2

  mod_1$update_state(state = state_1)
  mod_1$set_index(info$index$S)
  y_1 <- mod_1$simulate(seq(0, 400, 4))

  ## High protection
  p$rel_susceptibility <- array(rep(c(1, 0.9), each = 16), dim = c(16, 2))
  mod_2 <- ZamCovid$new(p, 0, 1, seed = 1L)

  state_2 <- ZamCovid_initial(info, 1, p)
  state_2[index_S[, 2]] <- state_2[index_S[, 1]] / 2
  state_2[index_S[, 1]] <- state_2[index_S[, 1]] / 2

  mod_2$update_state(state = state_2)
  mod_2$set_index(info$index$S)
  y_2 <- mod_2$simulate(seq(0, 400, 4))

  # Reshape to show the full shape of s
  y_1 <- array(y_1, c(info$dim$S, 101))
  y_2 <- array(y_2, c(info$dim$S, 101))


  expect_true(sum(y_1[, 2, 101]) > sum(y_2[, 2, 101]))
})
