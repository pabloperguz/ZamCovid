context("ZamCovid (parameters)")

test_that("ZamCovid_parameters returns a list of parameters", {
  date <- numeric_date("2020-02-01")
  p <- ZamCovid_parameters(date)

  expect_type(p, "list")
  expect_equal(p$beta_step, 0.1)

  ## Transmission matrix is of expected characteristics
  expect_equal(dim(p$m), c(16, 16))
  expect_identical(
    p$m,
    as.matrix(make_contact_matrix("extdata/matrix.csv", p$N_tot)))

  progression <- ZamCovid_parameters_progression(0.25)
  expect_identical(p[names(progression)], progression)

  vaccination <- ZamCovid_parameters_vaccination(
    p$dt, p$N_tot, p$n_groups, p$n_doses,
    p$rel_susceptibility, p$rel_p_sympt, p$rel_p_hosp_if_sympt,
    p$rel_p_death, p$rel_infectivity,
    p$vaccine_progression_rate, NULL, NULL)
  expect_identical(p[names(vaccination)], vaccination)

  severity <- ZamCovid_parameters_severity(0.25, NULL)
  expect_identical(p[names(severity)], severity)

  ## TODO: it'll be good to have a function for observation parameters
  # observation <- ZamCovid_parameters_observation(1e6)
  # expect_identical(p[names(observation)], observation)

  sens_and_spec <- ZamCovid_parameters_sens_and_spec()
  expect_identical(p[names(sens_and_spec)], sens_and_spec)

  expect_equal(p$N_tot_over15, sum(p$N_tot[4:16]))

  extra <- setdiff(names(p),
                   c("dt", "m", names(progression), names(severity),
                     names(vaccination), names(sens_and_spec)))
  expect_equal(
    extra,
    c("steps_per_day",       "beta_step",          "base_death_step",
      "cross_immunity_step", "n_groups",           "seed_size",
      "seed_age_band",       "seed_step_start",    "seed_value",
      "I_A_transmission",    "I_P_transmission",   "I_C_1_transmission",
      "I_C_2_transmission",  "hosp_transmission",  "G_D_transmission",
      "phi_death_all",       "kappa_death_all",    "phi_admitted",
      "kappa_admitted",      "phi_death_hosp",     "kappa_death_hosp",
      "exp_noise",           "N_tot",              "N_tot_all",
      "N_tot_over15",        "N_tot_15_19",        "N_tot_20_29",
      "N_tot_30_39",         "N_tot_40_49",        "N_tot_50_plus",
      "rel_p_H_D",           "rel_p_G_D",           "rel_p_R"))
})


test_that("Can input population data", {

  pop <- read_csv(ZamCovid_file("extdata/population.csv"))
  p <- ZamCovid_parameters(1, population = pop)
  expect_equal(p$N_tot, pop$n)

  m <- p$m
  expect_error(ZamCovid_parameters(1, contact_matrix = m))
  expect_error(ZamCovid_parameters(1, contact_matrix = as.list(m)))

  expect_true(
    is.list(ZamCovid_parameters(1, contact_matrix = m, population = pop)))

  pop[1L, 2] <- 0.5
  expect_error(
    ZamCovid_parameters(1, contact_matrix = m, population = pop))

  pop[1L, 2] <- -1
  expect_error(
    ZamCovid_parameters(1, contact_matrix = m, population = pop))

  pop <- read_csv(ZamCovid_file("extdata/population.csv"))[-16, ]
  expect_error(
    ZamCovid_parameters(1, contact_matrix = m, population = pop))
})


test_that("can compute severity for ZamCovid model", {

  # Get default values
  severity <- ZamCovid_parameters_severity(0.25, NULL)
  expect_equal(
    severity$p_star_step, array(0.2, c(1, 16)))

  # Can input severity data
  data <- read_csv(ZamCovid_file("extdata/severity_default.csv"))
  data[data$Name == "p_star", -1L] <- 0

  severity <- ZamCovid_parameters_severity(0.25, data)
  expect_equal(
    severity$p_star_step, array(0, c(1, 16)))

  # Compute time-varying severity
  dt <- 0.25
  p_star_date <- numeric_date(c("2020-02-01", "2020-05-01"))
  p_star_value <- c(0.2, 0.5)

  severity <-
    ZamCovid_parameters_severity(dt, NULL,
                                 p_star = list(date = p_star_date,
                                              value = p_star_value))
  p_star_step <- parameters_piecewise_linear(p_star_date, p_star_value, dt)

  expect_equal(severity$p_star_step[, 16], p_star_step)


  expect_error(
    ZamCovid_parameters_severity(dt,
                                 p_star = list(date = 1,
                                               value = 0.3)),
    "As 'p_star' has a single 'value', expected NULL or missing 'date'")

  expect_error(
    ZamCovid_parameters_severity(dt,
                                 p_star = list(date = c(1, 4, 5),
                                              value = c(0.2, 0.3))),
    "'date' and 'value' for 'p_star' must have the same length")

  expect_error(
    ZamCovid_parameters_severity(dt,
                                 p_star = list(date = c(1, 4),
                                              value = c(-1, 0.3))),
    "'p_star' must lie in [0, 1]", fixed = TRUE)

  expect_error(
    ZamCovid_parameters_severity(dt,
                                 p_star = list(date = c(1, 4),
                                              value = c(0.2, 3))),
    "'p_star' must lie in [0, 1]", fixed = TRUE)

})


test_that("All ZamCovid progression parameters are returned as default", {
  p <- ZamCovid_parameters_progression(dt = 1 / 4)
  expect_setequal(
    names(p),
    c("k_E",                    "k_H_R",                  "k_H_D",
      "k_G_D",                  "k_sero_pre",             "k_sero_pos",
      "k_PCR_pre",              "k_PCR_pos",              "gamma_R",
      "gamma_E_step",           "n_gamma_E_steps",        "gamma_A_step",
      "n_gamma_A_steps",        "gamma_P_step",           "n_gamma_P_steps",
      "gamma_C_1_step",         "n_gamma_C_1_steps",      "gamma_C_2_step",
      "n_gamma_C_2_steps",      "gamma_H_D_step",         "n_gamma_H_D_steps",
      "gamma_H_R_step",         "n_gamma_H_R_steps",      "gamma_G_D_step",
      "n_gamma_G_D_steps",      "gamma_U_step",           "n_gamma_U_steps",
      "gamma_PCR_pre_step",     "n_gamma_PCR_pre_steps",  "gamma_PCR_pos_step",
      "n_gamma_PCR_pos_steps",  "gamma_sero_pos_step",
      "n_gamma_sero_pos_steps", "gamma_sero_pre_step",
      "n_gamma_sero_pre_steps"))
})


# test_that("can compute time-varying progression parameters for ZamCovid
#           model", {
#             dt <- 0.25
#
#             gamma_H_R_value <- 0.3
#             gamma_H_D_date <- sircovid_date(c("2020-02-01", "2020-05-01"))
#             gamma_H_D_value <- c(0.2, 0.5)
#
#             progression <-
#               ZamCovid_parameters_progression(dt,
#                                               gamma_H_D =
#                                                 list(date = gamma_H_D_date,
#                                                      value = gamma_H_D_value),
#                                               gamma_H_R =
#                                                 list(value = gamma_H_R_value)
#               )
#
#             gamma_H_D_step <-
#               sircovid_parameters_piecewise_linear(gamma_H_D_date,
#                                                    gamma_H_D_value, dt)
#             expect_equal(progression$gamma_H_D_step, gamma_H_D_step)
#             expect_equal(progression$gamma_H_R_step, gamma_H_R_value)
#
#             expect_error(
#               ZamCovid_parameters_progression(dt,
#                                               gamma_E = list(date = 1,
#                                                              value = 3)),
#               "'gamma_E' has a single 'value', expected NULL or missing 'date'")
#
#             expect_error(
#               ZamCovid_parameters_progression(dt,
#                                               gamma_ICU_pre =
#                                                 list(date = c(1, 4, 5),
#                                                      value = c(2, 3))),
#               "'date' and 'value' for 'gamma_ICU_pre' must have the same length"
#             )
#
#             expect_error(
#               ZamCovid_parameters_progression(dt,
#                                               gamma_H_D =
#                                                 list(date = c(1, 4),
#                                                      value = c(-2, 3))),
#               "'gamma_H_D' must have only non-negative values")
#
#
#           })


test_that("ZamCovid_index is properly named", {
  p <- ZamCovid_parameters(1)
  mod <- ZamCovid$new(p, 0, 10)
  info <- mod$info()
  index <- ZamCovid_index(info)
  expect_false(any(is.na(names(index))))
})


test_that("Can compute initial conditions and control the seeding", {
  p <- ZamCovid_parameters(1)
  mod <- ZamCovid$new(p, 0, 10)
  info <- mod$info()

  initial <- ZamCovid_initial(info, 10, p)
  initial_y <- mod$transform_variables(initial)

  expect_equal(initial_y$N_tot_sero, sum(p$N_tot))
  expect_equal(initial_y$N_tot_PCR, sum(p$N_tot))
  expect_equal(initial_y$N_tot, p$N_tot)

  expect_equal(rowSums(initial_y$S), p$N_tot)

  ## We're at 674 compartments at the moment
  # 36 derived from:
  # S (16) + susceptible (1) + N_tot_PCR (1) + N_tot_sero (1) + N_tot(16)
  expect_equal(sum(initial != 0), 36)


  ## Change initial seed size
  ## TODO: this need properly checking!!
  p <- ZamCovid_parameters(1, initial_seed_size = 50)
  expect_equal(p$seed_size, 50)

  mod <- ZamCovid$new(p, 0, 10)
  info <- mod$info()

  initial <- ZamCovid_initial(info, 10, p)

  initial_y <- mod$transform_variables(initial)

  expect_equal(initial_y$N_tot_sero, sum(p$N_tot))
  expect_equal(initial_y$N_tot_PCR, sum(p$N_tot))
  expect_equal(initial_y$N_tot, p$N_tot)

  expect_equal(rowSums(initial_y$S) + drop(initial_y$I_A),
               p$N_tot)

  expect_equal(sum(initial != 0), 36)
})


test_that("ZamCovid check severity works as expected", {
  ## errors if params missing
  expect_error(
    check_severity(list(n_groups = 16)),
    "Parameter 'rel_p_sympt' is missing"
  )
  expect_error(
    check_severity(list(n_groups = 16, rel_p_sympt = 1)),
    "Parameter 'p_C_step' is missing"
  )

  ## errors if rel_ is not 1 or 4 cols
  # expect_error(check_severity(list(
  #   n_groups = 16,
  #   rel_p_sympt = array(1, c(16, 2, 3)),
  #   p_C_step = matrix(1, 1, 16)
  # )), "1 column")

  ## no error if 1 col
  expect_error(check_severity(list(
    n_groups = 16,
    rel_p_sympt = array(runif(16 * 3), c(16, 1, 3)),
    p_C_step = matrix(1, 1, 16)
  )), "Parameter 'rel_p_hosp_if_sympt' is missing", fixed = TRUE)

  ## check on required parameters
  steps <- c("p_C_step", "p_H_step", "p_H_D_step", "p_G_D_step", "p_R_step")
  rels <- c(
    "rel_p_sympt", "rel_p_hosp_if_sympt", "rel_p_H_D", "rel_p_G_D", "rel_p_R")

  p <- vector("list", 10)
  names(p) <- c(steps, rels)
  ## 16 groups, 3 vacc classes
  p[1:5] <- list(matrix(0.5, ncol = 16))
  p[6:10] <- list(array(1, c(16, 1, 3)))

  p$n_groups <- 16

  ## first check no errors
  expect_equal(check_severity(p), p)

  ## check errors when expected - we only check upper bound as it's safe to
  ##  assume negative probs won't be given in practice (but uses the same )

  ## rel_p_sympt Inf
  for (i in seq_along(steps)) {
    p[[rels[[i]]]][] <- -1
    expect_error(
      check_severity(p),
      sprintf("'%s' must have only non-negative values", rels[[i]]),
      fixed = TRUE
    )
    p[[rels[[i]]]][] <- 1

    p[[steps[[i]]]][] <- -1
    expect_error(
      check_severity(p),
      sprintf("'%s' must have only non-negative values", steps[[i]]),
      fixed = TRUE
    )
    p[[steps[[i]]]][] <- 1
  }
})


# test_that("ZamCovid_compare combines likelihood correctly", {
#   state <- rbind(
#     time = 45,
#     icu = 10:15,
#     general = 20:25,
#     deaths_carehomes_inc = 2:7,
#     deaths_comm_inc = 1:6,
#     deaths_comm_0_49_inc = 1:6,
#     deaths_comm_50_54_inc = 1:6,
#     deaths_comm_55_59_inc = 2:7,
#     deaths_comm_60_64_inc = 2:7,
#     deaths_comm_65_69_inc = 3:8,
#     deaths_comm_70_74_inc = 3:8,
#     deaths_comm_75_79_inc = 4:9,
#     deaths_comm_80_plus_inc = 4:9,
#     deaths_hosp_inc = 3:8,
#     deaths_hosp_0_49_inc = 1:6,
#     deaths_hosp_50_54_inc = 1:6,
#     deaths_hosp_55_59_inc = 2:7,
#     deaths_hosp_60_64_inc = 2:7,
#     deaths_hosp_65_69_inc = 3:8,
#     deaths_hosp_70_74_inc = 3:8,
#     deaths_hosp_75_79_inc = 4:9,
#     deaths_hosp_80_plus_inc = 4:9,
#     admitted_inc = 50:55,
#     diagnoses_inc = 60:65,
#     all_admission_0_9_inc = 1:6,
#     all_admission_10_19_inc = 1:6,
#     all_admission_20_29_inc = 1:6,
#     all_admission_30_39_inc = 1:6,
#     all_admission_40_49_inc = 9:14,
#     all_admission_50_59_inc = 12:17,
#     all_admission_60_69_inc = 17:22,
#     all_admission_70_79_inc = 20:25,
#     all_admission_80_plus_inc = 27:32,
#     sero_pos_1 = 4:9,
#     sero_pos_2 = 14:19,
#     sympt_cases_inc = 100:105,
#     sympt_cases_non_variant_inc = 70:75,
#     sympt_cases_over25_inc = 80:85,
#     sympt_cases_under15_inc = 5:10,
#     sympt_cases_15_24_inc = 5:10,
#     sympt_cases_25_49_inc = 19:24,
#     sympt_cases_50_64_inc = 19:24,
#     sympt_cases_65_79_inc = 19:24,
#     sympt_cases_80_plus_inc = 19:24,
#     sympt_cases_non_variant_over25_inc = 60:65,
#     ons_pos = 2:7,
#     react_pos = 2:7,
#     react_5_24_pos = 1:6,
#     react_25_34_pos = 1:6,
#     react_35_44_pos = 1:6,
#     react_45_54_pos = 1:6,
#     react_55_64_pos = 1:6,
#     react_65_plus_pos = 1:6)
#   observed <- list(
#     icu = 13,
#     general = 23,
#     hosp = 36,
#     deaths_carehomes = 4,
#     deaths_hosp = 5,
#     deaths_hosp_0_49 = 1,
#     deaths_hosp_50_54 = 2,
#     deaths_hosp_55_59 = 3,
#     deaths_hosp_60_64 = 4,
#     deaths_hosp_65_69 = 5,
#     deaths_hosp_70_74 = 6,
#     deaths_hosp_75_79 = 7,
#     deaths_hosp_80_plus = 8,
#     deaths_comm = 3,
#     deaths_comm_0_49 = 1,
#     deaths_comm_50_54 = 2,
#     deaths_comm_55_59 = 3,
#     deaths_comm_60_64 = 4,
#     deaths_comm_65_69 = 5,
#     deaths_comm_70_74 = 6,
#     deaths_comm_75_79 = 7,
#     deaths_comm_80_plus = 8,
#     deaths = 8,
#     deaths_non_hosp = 6,
#     admitted = 53,
#     diagnoses = 63,
#     all_admission = 116,
#     all_admission_0_9 = 2,
#     all_admission_10_19 = 4,
#     all_admission_20_29 = 5,
#     all_admission_30_39 = 5,
#     all_admission_40_49 = 10,
#     all_admission_50_59 = 15,
#     all_admission_60_69 = 20,
#     all_admission_70_79 = 25,
#     all_admission_80_plus = 30,
#     sero_pos_15_64_1 = 43,
#     sero_tot_15_64_1 = 83,
#     sero_pos_15_64_2 = 58,
#     sero_tot_15_64_2 = 98,
#     pillar2_pos = 35,
#     pillar2_tot = 600,
#     pillar2_cases = 35,
#     pillar2_over25_pos = 25,
#     pillar2_over25_tot = 500,
#     pillar2_over25_cases = 25,
#     pillar2_under15_cases = 8,
#     pillar2_15_24_cases = 8,
#     pillar2_25_49_cases = 20,
#     pillar2_50_64_cases = 20,
#     pillar2_65_79_cases = 20,
#     pillar2_80_plus_cases = 20,
#     pillar2_under15_pos = 8,
#     pillar2_15_24_pos = 8,
#     pillar2_25_49_pos = 20,
#     pillar2_50_64_pos = 20,
#     pillar2_65_79_pos = 20,
#     pillar2_80_plus_pos = 20,
#     pillar2_under15_tot = 160,
#     pillar2_15_24_tot = 160,
#     pillar2_25_49_tot = 400,
#     pillar2_50_64_tot = 400,
#     pillar2_65_79_tot = 400,
#     pillar2_80_plus_tot = 400,
#     ons_pos = 4,
#     ons_tot = 600,
#     react_pos = 3,
#     react_tot = 500,
#     react_5_24_pos = 1,
#     react_5_24_tot = 50,
#     react_25_34_pos = 1,
#     react_25_34_tot = 50,
#     react_35_44_pos = 1,
#     react_35_44_tot = 100,
#     react_45_54_pos = 1,
#     react_45_54_tot = 100,
#     react_55_64_pos = 1,
#     react_55_64_tot = 100,
#     react_65_plus_pos = 2,
#     react_65_plus_tot = 100,
#     strain_non_variant = 40,
#     strain_tot = 50,
#     strain_over25_non_variant = 20,
#     strain_over25_tot = 25)
#   date <- sircovid_date("2020-01-01")
#   pars <- ZamCovid_parameters(date, "uk", exp_noise = Inf)
#
#   observed_keep <- function(nms) {
#     observed[setdiff(names(observed), nms)] <- NA_real_
#     observed
#   }
#   observed_drop <- function(nms) {
#     observed[nms] <- NA_real_
#     observed
#   }
#
#   ## This function is more complicated to test than the basic model
#   ## because it's not a simple sum
#   nms_sero_1 <- c("sero_pos_15_64_1", "sero_tot_15_64_1")
#   nms_sero_2 <- c("sero_pos_15_64_2", "sero_tot_15_64_2")
#   nms_pillar2 <- c("pillar2_pos", "pillar2_tot")
#   nms_pillar2_over25 <- c("pillar2_over25_pos", "pillar2_over25_tot")
#   nms_pillar2_under15 <- c("pillar2_under15_pos", "pillar2_under15_tot")
#   nms_pillar2_15_24 <- c("pillar2_15_24_pos", "pillar2_15_24_tot")
#   nms_pillar2_25_49 <- c("pillar2_25_49_pos", "pillar2_25_49_tot")
#   nms_pillar2_50_64 <- c("pillar2_50_64_pos", "pillar2_50_64_tot")
#   nms_pillar2_65_79 <- c("pillar2_65_79_pos", "pillar2_65_79_tot")
#   nms_pillar2_80_plus <- c("pillar2_80_plus_pos", "pillar2_80_plus_tot")
#   nms_ons <- c("ons_pos", "ons_tot")
#   nms_react <- c("react_pos", "react_tot")
#   nms_react_5_24 <- c("react_5_24_pos", "react_5_24_tot")
#   nms_react_25_34 <- c("react_25_34_pos", "react_25_34_tot")
#   nms_react_35_44 <- c("react_35_44_pos", "react_35_44_tot")
#   nms_react_45_54 <- c("react_45_54_pos", "react_45_54_tot")
#   nms_react_55_64 <- c("react_55_64_pos", "react_55_64_tot")
#   nms_react_65_plus <- c("react_65_plus_pos", "react_65_plus_tot")
#   nms_strain <- c("strain_non_variant", "strain_tot")
#   nms_strain_over25 <- c("strain_over25_non_variant", "strain_over25_tot")
#   parts <- c(as.list(setdiff(names(observed),
#                              c(nms_sero_1, nms_sero_2, nms_pillar2,
#                                nms_pillar2_over25, nms_pillar2_under15,
#                                nms_pillar2_15_24, nms_pillar2_25_49,
#                                nms_pillar2_50_64, nms_pillar2_65_79,
#                                nms_pillar2_80_plus, nms_ons,
#                                nms_react, nms_react_5_24,
#                                nms_react_25_34, nms_react_35_44,
#                                nms_react_45_54, nms_react_55_64,
#                                nms_react_65_plus, nms_strain,
#                                nms_strain_over25))),
#              list(nms_sero_1), list(nms_sero_2), list(nms_pillar2),
#              list(nms_pillar2_over25), list(nms_pillar2_under15),
#              list(nms_pillar2_15_24), list(nms_pillar2_25_49),
#              list(nms_pillar2_50_64), list(nms_pillar2_65_79),
#              list(nms_pillar2_80_plus), list(nms_ons), list(nms_react),
#              list(nms_react_5_24), list(nms_react_25_34),
#              list(nms_react_35_44), list(nms_react_45_54),
#              list(nms_react_55_64), list(nms_react_65_plus),
#              list(nms_strain), list(nms_strain_over25))
#
#   ll_parts <- lapply(parts, function(x)
#     ZamCovid_compare(state, observed_keep(x), pars))
#
#   ## Extremely light testing, though this has already flushed out some
#   ## issues
#   expect_vector_equal(lengths(ll_parts), 6)
#   expect_equal(
#     ZamCovid_compare(state, observed, pars),
#     rowSums(do.call(cbind, ll_parts)))
#
#   ## Test that there are non-zero values for each log-likelihood part.
#   ## This helps make sure all parts contribute to the log-likelihood.
#   expect_true(all(sapply(ll_parts, function(x) any(x != 0))))
#
#
#   ## Check that weekend effect parameters work as expected. For each day
#   ## use same state values except time
#   state <- state[, 1, drop = FALSE]
#   time <- seq(20, 26, 1)
#   age_nms <- c("_under15", "_15_24", "_25_49",
#                "_50_64", "_65_79", "_80_plus")
#
#   for (i in age_nms) {
#     pars[[paste0("p_NC_weekend", i)]] <- pars[[paste0("p_NC", i)]]
#     pars[[paste0("phi_pillar2_cases_weekend", i)]] <-
#       pars[[paste0("phi_pillar2_cases", i)]]
#   }
#
#   helper <- function(parameter, age_nm) {
#
#     pars[[paste0(parameter, "_weekend", age_nm)]] <-
#       0.5 * pars[[paste0(parameter, age_nm)]]
#
#     ll_time <- vnapply(time, function(x) {
#       state["time", ] <- x
#       ZamCovid_compare(state, observed, pars)
#     })
#
#     ## First 5 days are weekdays, last 2 are weekend
#     expect_equal(grepl("^S", weekdays(sircovid_date_as_date(time))),
#                  c(rep(FALSE, 5), rep(TRUE, 2)))
#     ## Weekdays should yield the same values. Weekends should yield the same
#     ## values, but different to weekdays
#     expect_equal(ll_time[2:5], rep(ll_time[1], 4))
#     expect_equal(ll_time[6], ll_time[7])
#     expect_true(ll_time[1] != ll_time[6])
#   }
#
#   lapply(age_nms, function(x) helper("p_NC", x))
#   lapply(age_nms, function(x) helper("phi_pillar2_cases", x))
#
# })


# test_that("ZamCovid vaccination parameters", {
#   n_groups <- length(ZamCovid:::model_age_bins()$start)
#   # test default values
#   ntot <- rep(1000, n_groups)
#   p <- ZamCovid_parameters_vaccination(0.24, ntot, n_groups)
#   expect_setequal(
#     names(p),
#     c("rel_susceptibility", "rel_p_sympt", "rel_p_hosp_if_sympt", "rel_p_death",
#       "rel_infectivity",
#       "n_vacc_classes",
#       "vaccine_progression_rate_base", "vaccine_dose_step",
#       "vaccine_catchup_fraction",
#       "index_dose", "index_dose_inverse", "n_doses"))
#   expect_equal(nrow(p$rel_susceptibility), n_groups)
#   expect_equal(ncol(p$rel_susceptibility), 1)
#   expect_equal(nlayer(p$rel_susceptibility), 1)
#   expect_equal(nrow(p$vaccine_progression_rate_base), n_groups)
#   expect_equal(ncol(p$vaccine_progression_rate_base), 1)
#
#   # test when more vaccinated categories than default
#   rel_susceptibility <- c(1, 0.2, 0.1, 0.4)
#   rel_p_sympt <- c(1, 0.75, 0.5, 0.75)
#   rel_p_hosp_if_sympt <- c(1, 0.8, 0.6, 0.9)
#   vaccine_progression_rate <- c(0, 1, 0, 1)
#
#   region <- "london"
#   vaccine_daily_doses <- c(5000, 10000)
#   vaccine_daily_doses_date <- sircovid_date(c("2020-12-01", "2021-02-01"))
#   daily_doses <- c(rep(vaccine_daily_doses[1], diff(vaccine_daily_doses_date)),
#                    rep(vaccine_daily_doses[2], 200))
#   uptake <- test_example_uptake()
#   n <- vaccine_priority_population(region, uptake)
#
#   vaccine_schedule <- vaccine_schedule_future(
#     vaccine_daily_doses_date[1], daily_doses, 1e6, n)
#
#   p <- ZamCovid_parameters_vaccination(ntot,
#                                        dt = 0.25,
#                                        rel_susceptibility = rel_susceptibility,
#                                        rel_p_sympt = rel_p_sympt,
#                                        rel_p_hosp_if_sympt =
#                                          rel_p_hosp_if_sympt,
#                                        vaccine_progression_rate =
#                                          vaccine_progression_rate,
#                                        vaccine_schedule = vaccine_schedule,
#                                        vaccine_index_dose2 = 3L)
#   expect_setequal(
#     names(p),
#     c("rel_susceptibility", "rel_p_sympt", "rel_p_hosp_if_sympt", "rel_p_death",
#       "rel_infectivity",
#       "n_vacc_classes",
#       "vaccine_progression_rate_base", "vaccine_dose_step",
#       "vaccine_catchup_fraction",
#       "index_dose", "index_dose_inverse", "n_doses"))
#
#   expect_equal(nrow(p$rel_susceptibility), n_groups)
#   expect_equal(ncol(p$rel_susceptibility), 1)
#   expect_equal(nlayer(p$rel_susceptibility), length(rel_susceptibility))
#
#   expect_equal(nrow(p$rel_p_sympt), n_groups)
#   expect_equal(ncol(p$rel_p_sympt), 1)
#   expect_equal(dim(p$rel_p_sympt)[3], length(rel_p_sympt))
#
#   expect_equal(nrow(p$rel_p_hosp_if_sympt), n_groups)
#   expect_equal(ncol(p$rel_p_hosp_if_sympt), 1)
#   expect_equal(dim(p$rel_p_hosp_if_sympt)[3], length(rel_p_hosp_if_sympt))
#
#   expect_equal(nrow(p$rel_infectivity), n_groups)
#   expect_equal(ncol(p$rel_infectivity), 1)
#   expect_equal(dim(p$rel_infectivity)[3], length(rel_susceptibility))
#
#   expect_equal(nrow(p$vaccine_progression_rate_base), n_groups)
#   expect_equal(ncol(p$vaccine_progression_rate_base),
#                length(vaccine_progression_rate))
#
#   expect_equal(dim(p$vaccine_dose_step),
#                c(19, 2,
#                  (length(daily_doses) + vaccine_daily_doses_date[1]) * 4))
#   ## daily doses are as expected
#   expect_vector_equal(colSums(p$vaccine_dose_step[, 1, ]),
#                       c(rep(0, vaccine_daily_doses_date[1] * 4),
#                         rep(daily_doses / 4, each = 4)),
#                       digits = 0, tol = 2)
#   msg1 <- "rel_susceptibility, rel_p_sympt, rel_p_hosp_if_sympt, rel_p_death,"
#   msg2 <- "rel_infectivity should have the same dimension"
#   expect_error(
#     ZamCovid_parameters_vaccination(ntot,
#                                     rel_susceptibility = 1,
#                                     rel_p_sympt = c(1, 0.5, 0.25),
#                                     rel_p_hosp_if_sympt = c(1, 0.1),
#                                     rel_p_death = c(1, 0.2),
#                                     rel_infectivity = 1),
#     paste(msg1, msg2))
#   expect_error(ZamCovid_parameters_vaccination(ntot,
#                                                rel_susceptibility = c(1, 1),
#                                                rel_p_sympt = c(1, 0.5, 0.25),
#                                                rel_p_hosp_if_sympt = 1,
#                                                rel_p_death = 1,
#                                                rel_infectivity = 1),
#                paste(msg1, msg2))
#   expect_error(ZamCovid_parameters_vaccination(ntot,
#                                                rel_susceptibility = c(1, 1),
#                                                rel_p_sympt = c(1, 0.5),
#                                                rel_p_hosp_if_sympt = c(1, 1),
#                                                rel_p_death = c(1, 0.2),
#                                                rel_infectivity = c(1, 1, 0.5)),
#                paste(msg1, msg2))
#   expect_error(ZamCovid_parameters_vaccination(ntot,
#                                                vaccine_catchup_fraction = -1),
#                "'vaccine_catchup_fraction' must lie in [0, 1]",
#                fixed = TRUE)
#   expect_error(ZamCovid_parameters_vaccination(ntot,
#                                                vaccine_catchup_fraction = 1.5),
#                "'vaccine_catchup_fraction' must lie in [0, 1]",
#                fixed = TRUE)
#   expect_error(
#     ZamCovid_parameters_vaccination(ntot,
#                                     vaccine_catchup_fraction = c(1, 1)),
#     "'vaccine_catchup_fraction' must be a scalar")
# })

# test_that("can tune the noise parameter", {
#   p1 <- ZamCovid_parameters_observation(1e6)
#   p2 <- ZamCovid_parameters_observation(1e4)
#   expect_setequal(names(p1), names(p2))
#   v <- setdiff(names(p1), "exp_noise")
#   expect_mapequal(p1[v], p2[v])
#   expect_equal(p1$exp_noise, 1e6)
#   expect_equal(p2$exp_noise, 1e4)
# })
