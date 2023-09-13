##' ##' "Basic" SEIR model. This is a dust model.
##' @name basic
##' @title The basic SEIR model
##' @export basic
NULL

##' ##' Main ZamCovid SEIR model. This is a dust model.
##' @name ZamCovid
##' @title The main ZamCovid model
##' @export ZamCovid
NULL


#' Define ZamCovid model parameters
#'
#' @param start_date A positive integer, greated than `model_start` for the model
#'    end time.
#'
#' @param steps_per_day A positive integer for the number of steps per day to
#'    run the model.
#'
#' @param beta_date A vector of `numeric_date`s for constructing piece-wise
#'    linear interpolate function for beta (transmission scaling) parameter.
#'
#' @param beta_value A vector of `double` values for constructing piece-wise
#'    linear interpolate function for beta (transmission scaling) parameter.
#'
#' @param base_death_date A vector of `numeric_deaths`s for constructing
#'    piece-wise linear interpolate function of daily baseline deaths.
#'
#' @param base_death_value A vector of `double` values for constructing
#'    piece-wise linear interpolate function of daily baseline deaths.
#'
#' @param cross_immunity_date A vector of `numeric_deaths`s for constructing
#'    piece-wise linear interpolate function of cross-immunity vs infection.
#'
#' @param cross_immunity_value A vector of `double` values for constructing
#'    piece-wise linear interpolate function of cross-immunity vs infection.
#'
#' @param population Either `NULL` or a `data.frame` with columns `age_group`
#'    and `n` of same length (age groups) as `contact_matrix` supplied. If
#'    `NULL` defaults to nationally representative values for Zambia.
#'
#' @param contact_matrix Either `NULL` or a `data.frame` of a symmetric
#'    individual-to-individual contact rate matrix. If `NULL` defaults to
#'    nationally representative values for Zambia.
#'
#' @param N_tot Either `NULL` or a `vector` integer numbers of age-specific
#'    population of same length as age-groups in `contact_matri`. Must be
#'    supplied if `population` and `contact_matrix` are not `NULL`.
#'
#' @param n_groups Either `NULL` or an `integer` number specifying the number of
#'    age groups in the populsation to model. Must be supplied if `population`
#'    and `contact_matrix` are not `NULL`.
#'
#' @param severity Either `NULL` or a named `list` containing all elements in
#'    [ZamCovid::ZamCovid_parameters_severity].
#'
#' @param progression Either `NULL` or a named `list` containing all elements in
#'    [ZamCovid::ZamCovid_parameters_progression].
#'
#' @param observation Either `NULL` or a named `list` containing all elements in
#'    [ZamCovid::ZamCovid_parameters_observation].
#'
#' @param sens_and_spec Either `NULL` or a named `list` containing all elements
#'    in [ZamCovid::ZamCovid_parameters_sens_and_spec].
#'
#' @param vaccination Either `NULL` or a named `list` containing all elements in
#'    [ZamCovid::ZamCovid_parameters_vaccination]. Nota that, if `NULL` and
#'    modelling vaccination, a `vaccine_progression_rate`, `vaccine_schedule`,
#'    `vaccine_index_dose2`, `vaccine_catchup_fraction`
#'    and `n_doses` must be supplied.
#'
#' @return Return a `list` of parameters to run the basic ZamCovid model.
#'
#' @export
ZamCovid_parameters <- function(start_date,
                                steps_per_day = 4L,
                                beta_date = NULL,
                                beta_value = NULL,
                                base_death_date = NULL,
                                base_death_value = NULL,
                                cross_immunity_date = NULL,
                                cross_immunity_value = NULL,
                                population = NULL,
                                contact_matrix = NULL,
                                N_tot = NULL,
                                n_groups = NULL,
                                severity = NULL,
                                progression = NULL,
                                observation = NULL,
                                sens_and_spec = NULL,
                                vaccination = NULL,
                                vaccine_progression_rate = NULL,
                                vaccine_schedule = NULL,
                                vaccine_index_dose2 = NULL,
                                vaccine_catchup_fraction = NULL,
                                n_doses = 2L,
                                initial_seed_pattern = 1,
                                initial_seed_size = 10,
                                I_A_transmission = 0.223,
                                I_P_transmission = 1,
                                I_C_1_transmission = 1,
                                I_C_2_transmission = 0,
                                hosp_transmission = 0,
                                G_D_transmission = 0,
                                rel_susceptibility = 1,
                                rel_p_sympt = 1,
                                rel_p_hosp_if_sympt = 1,
                                rel_p_death = 1,
                                rel_infectivity = 1) {

  if (!is.numeric(start_date)) {
    stop("start_date must be numeric! Did you forget to do `numeric_date` on
         string data of format yyyy-mm-dd?")
  }

  if (!is.integer(steps_per_day)) {
    stop("Steps per day must be an integer!")
  }

  if (length(beta_date) != length(beta_value)) {
    stop("beta_date and beta_value must be of same length!")
  }

  if (!is.null(beta_date)) {
    if (!is.numeric(beta_date)) {
      stop("beta_date must be numeric! Did you forget to do `numeric_date`?")
    }
  }

  if (!is.null(contact_matrix)) {
    if (is.null(population)) {
      stop("A `population` must be provided with contact matrix!")
    }

    stopifnot(colnames(population) == c("age_group", "n"))
    N_tot <- population$n
    n_groups <- length(N_tot)

    if (!all(is.integer(population$n))) {
      stop("Population must all be integer")
    }

    if (!all(population$n >= 0)) {
      stop("Population must all be non-negative values")
    }

    stopifnot(is.matrix(contact_matrix))
    stopifnot(nrow(contact_matrix) == n_groups)

  } else {

    if (is.null(population)) {
      population <- read_csv(ZamCovid_file("extdata/population.csv"))
    } else {
      population <- population
    }

    N_tot <- as.numeric(population$n)
    contact_matrix <-
      make_contact_matrix("extdata/matrix.csv", N_tot)
    n_groups <- nrow(contact_matrix)
    contact_matrix <- as.matrix(contact_matrix)
  }

  dt <- 1 / steps_per_day

  ## First epidemic wave seed
  start_step <- start_date * steps_per_day
  initial_seed_value <-
    seed_over_steps(start_step, initial_seed_pattern) * initial_seed_size


  ## Make piece-wise linear beta_step
  beta_step <- parameters_piecewise_linear(beta_date,
                                           beta_value %||% 0.1, dt)

  ## Make piece-wise constant step-function of baseline deaths
  base_death_step <- parameters_piecewise_constant(base_death_date,
                                                   base_death_value %||% 0,
                                                   dt)

  ## Make piece-wise linear function of cross_immunity (proxy for VOC emergence)
  cross_immunity_step <-
    parameters_piecewise_linear(cross_immunity_date,
                                cross_immunity_value %||% 0.95,
                                dt)

  ## Default parameters
  ret <- list(
    steps_per_day = steps_per_day,
    dt = dt,
    beta_step = beta_step,
    base_death_step = base_death_step,
    cross_immunity_step = cross_immunity_step,

    m = contact_matrix,
    n_groups = n_groups, #  16 n_groups (5 year age bands and 75+)

    seed_size = initial_seed_size,
    seed_age_band = 4L, # seed first wave in 15-19yo
    seed_step_start = floor(start_step),
    seed_value = initial_seed_value,
    I_A_transmission = I_A_transmission,
    I_P_transmission = I_P_transmission,
    I_C_1_transmission = I_C_1_transmission,
    I_C_2_transmission = I_C_2_transmission,
    hosp_transmission = hosp_transmission,
    G_D_transmission = G_D_transmission
  )

  severity <- severity %||% ZamCovid_parameters_severity(dt, severity)

  progression <- progression %||% ZamCovid_parameters_progression(dt)

  sens_and_spec <- sens_and_spec %||% ZamCovid_parameters_sens_and_spec()

  observation <- observation %||% ZamCovid_parameters_observation()

  vaccination <-
    vaccination %||% ZamCovid_parameters_vaccination(dt,
                                                     N_tot,
                                                     n_groups,
                                                     n_doses,
                                                     rel_susceptibility,
                                                     rel_p_sympt,
                                                     rel_p_hosp_if_sympt,
                                                     rel_p_death,
                                                     rel_infectivity,
                                                     vaccine_progression_rate,
                                                     vaccine_schedule,
                                                     vaccine_index_dose2)

  ret <- c(ret, severity, progression, vaccination, sens_and_spec, observation)

  ## Parameters to support compare function
  ret$N_tot <- N_tot
  ret$N_tot_all <- sum(N_tot[1:length(N_tot)])
  ret$N_tot_over15 <- sum(N_tot[4:length(N_tot)])
  ret$N_tot_15_19 <- N_tot[4]
  ret$N_tot_20_29 <- sum(N_tot[5:6])
  ret$N_tot_30_39 <- sum(N_tot[7:8])
  ret$N_tot_40_49 <- sum(N_tot[9:10])
  ret$N_tot_50_plus <- sum(N_tot[11:length(N_tot)])

  ## Add some bespoke rel_p parameters
  rel_p_death <- build_rel_param(n_groups, rel_p_death,
                                 vaccination$n_vacc_classes, "rel_p_death")
  ret$rel_p_H_D <- rel_p_death
  ret$rel_p_G_D <- rel_p_death
  ret$rel_p_R <- array(1, c(ret$n_groups, vaccination$n_vacc_classes))

  check_severity(ret)
}


##' ZamCovid progression parameters.  The `k_` parameters are the
##' scaling parameters for Erlang distributed progressions, while
##' the `gamma_` parameters are either the gamma parameters of that
##' distribution or the mean of exponentially distributed progressions.
##'
##' @title ZamCovid progression parameters
##'
##' @param dt The step size
##'
##' @section Time-varying parameters:
##' Every time varying parameter has the same format, which can be `NULL` (in
##'   which case a default single value is used) or a list with `date` and
##'   `value` for changes in the parameter. If `value` is scalar then
##'   `date` can be `NULL` or missing. If `value` is a vector then `date`
##'   must be a vector of numeric dates of the same length as `value`.
##'
##' @param gamma_E Time-varying mean duration in the E (exposed) compartment.
##'
##' @param gamma_A Time-varying mean duration in the I_A (asymptomatic)
##'   compartment.
##'
##' @param gamma_P Time-varying mean duration in the I_P (pre-symptomatic)
##'   compartment.
##'
##' @param gamma_C_1 Time-varying mean duration in the I_C_1 first symptomatic
##'   compartment.
##'
##' @param gamma_C_2 Time-varying mean duration in the I_C_2 second symptomatic
##'   compartment.
##'
##' @param gamma_H_D Time-varying parameters for the Erlang rate parameter of
##'   the duration in the H_D (hospitalisation leading to death) compartment.
##'
##' @param gamma_H_R Time-varying parameters for the Erlang rate parameter of
##'   the duration in the H_R (hospitalisation leading to recovery) compartment.
##'
##' @param gamma_G_D Time-varying parameters for the Erlang rate parameter of
##'   the duration in the G_D (death in the community) compartment.
##'
##' @param gamma_R Time-varying mean duration in the R (recovered) compartment;
##'   effectively duration of infection-induced immunity.
##'
##' @param gamma_U Time-varying mean delay from hospital admission to diagnosis
##'   for those not confirmed on admission.
##'
##' @param gamma_sero_pre Time-varying parameters for the Erlang rate parameter
##'   of the time from exposed to becoming seropositive.
##'
##' @param gamma_sero_pos Time-varying parameters for the Erlang rate parameter
##'   of the duration of seropositivity.
##'
##' @param gamma_PCR_pre Time-varying parameters for the Erlang rate parameter
##'   of the time from infection exposure to becoming PCR positive.
##'
##' @param gamma_PCR_pos Time-varying parameters for the Erlang rate parameter
##'   of the duration of PCR positivity.
##'
##' @return A list of parameter values
##'
##' @export
ZamCovid_parameters_progression <- function(dt,
                                            gamma_E = NULL,
                                            gamma_A = NULL,
                                            gamma_P = NULL,
                                            gamma_C_1 = NULL,
                                            gamma_C_2 = NULL,
                                            gamma_H_D = NULL,
                                            gamma_H_R = NULL,
                                            gamma_G_D = NULL,
                                            gamma_R = NULL,
                                            gamma_U = NULL,
                                            gamma_sero_pre = NULL,
                                            gamma_sero_pos = NULL,
                                            gamma_PCR_pre = NULL,
                                            gamma_PCR_pos = NULL) {

  ## The k_ parameters are the shape parameters for the Erlang
  ## distribution, while the gamma parameters are the rate
  ## parameters of that distribution.
  ret <- list(k_E = 2,
              k_H_R = 1,
              k_H_D = 2,
              k_G_D = 2,
              k_sero_pre = 2,
              k_sero_pos = 2,
              k_PCR_pre = 2,
              k_PCR_pos = 2,

              gamma_E = 1 / 1.156069,
              gamma_A = 1 / 2.88,
              gamma_P = 1 / 1.68,
              gamma_C_1 = 1 / 2.14,
              gamma_C_2 = 1 / 1.86,
              gamma_H_D = 1 / 5.2,
              gamma_H_R = 1 / 10.7,
              gamma_G_D = 1 / (3 / 2),
              gamma_R = 1 / (3 * 365),
              gamma_U = 3 / 10,
              gamma_sero_pre = 1 / 10,
              gamma_sero_pos = 1 / 25,
              gamma_PCR_pre = 2 / 3,
              gamma_PCR_pos = 1 / 5
  )

  time_varying_gammas <- list(E = gamma_E,
                              A = gamma_A,
                              P = gamma_P,
                              C_1 = gamma_C_1,
                              C_2 = gamma_C_2,
                              H_D = gamma_H_D,
                              H_R = gamma_H_R,
                              G_D = gamma_G_D,
                              U = gamma_U,
                              PCR_pre = gamma_PCR_pre,
                              PCR_pos = gamma_PCR_pos,
                              sero_pos = gamma_sero_pos,
                              sero_pre = gamma_sero_pre)

  get_gamma_step <- function(x, name) {

    gamma_name <- paste0("gamma_", name)
    gamma <- x[[gamma_name]]
    time_vary <- time_varying_gammas[[name]]

    if (is.null(time_vary)) {
      gamma_value <- NULL
      gamma_date <- NULL
    } else {
      gamma_value <- time_vary$value
      if ("date" %in% names(time_vary)) {
        gamma_date <- time_vary$date
      } else {
        gamma_date <- NULL
      }
    }

    if (!is.null(gamma_value)) {
      assert_non_negative(gamma_value, gamma_name)
      if (length(gamma_value) == 1L) {
        if (length(gamma_date) != 0) {
          stop(sprintf(
            "As '%s' has a single 'value', expected NULL or missing 'date'",
            gamma_name))
        }
      } else if (length(gamma_date) != length(gamma_value)) {
        stop(sprintf("'date' and 'value' for '%s' must have the same length",
                     gamma_name))
      }
    }

    if (is.null(gamma_value)) {
      gamma_step <- gamma
    } else {
      gamma_step <- parameters_piecewise_linear(gamma_date, gamma_value, dt)
    }

    x[[paste0("gamma_", name, "_step")]] <- gamma_step
    x[[paste0("gamma_", name)]] <- NULL
    x[[paste0("n_gamma_", name, "_steps")]] <- length(gamma_step)

    x
  }

  for (name in names(time_varying_gammas)) {
    ret <- get_gamma_step(ret, name)
  }
  ret
}


##' ZamCovid diagnostic test performance parameters
##'
##' @title ZamCovid sensitivity and specificity parameters
##'
##' @return A list of parameter values
##'
##' @param sero_specificity Specificity of the serology test assay
##'
##' @param sero_sensitivity Sensitivity of the serology test assay
##'
##' @param pcr_specificity Specificity of the PCR test
##'
##' @param pcr_sensitivity Sensitivity of the PCR test
##'
##' @export
ZamCovid_parameters_sens_and_spec <- function(sero_specificity = 0.999,
                                              sero_sensitivity = 0.927,
                                              pcr_specificity = 0.99,
                                              pcr_sensitivity = 0.99) {

  ## Sero_sens of Abbott Architect IgG from PHE 2020:
  # Evaluation of sensitivity and specificity of four commercially available
  # SARS-CoV-2 antibody immunoassays

  ret <- list(
    sero_specificity = sero_specificity,
    sero_sensitivity = sero_sensitivity,
    pcr_specificity = pcr_specificity,
    pcr_sensitivity = pcr_sensitivity)

  lapply(seq_len(length(ret)),
         function(i) assert_proportion(ret[i], names(ret)[i]))

  ret
}


##' ZamCovid diagnostic test performance parameters
##'
##' @title ZamCovid sensitivity and specificity parameters
##'
##' @return A list of parameter values
##'
##' @param sero_specificity Specificity of the serology test assay
##'
##' @param sero_sensitivity Sensitivity of the serology test assay
##'
##' @param pcr_specificity Specificity of the PCR test
##'
##' @param pcr_sensitivity Sensitivity of the PCR test
##'
##' @export
ZamCovid_parameters_observation <- function(exp_noise = 1e6) {

  ## Note alphas should be fitted observation probabilities; so 1 / alpha in
  ## dnbinom will be `size` (i.e. we assume overdispersion, hence not using
  ## Poisson for observation counts).
  ## Also note we assume perfect reporting here (phis = 1), but "inflate" counts
  ## were relevant in ZamCovid-outputs/src/ZamCovid_parameters for now
  list(
    phi_death_all = 1,
    kappa_death_all = 2,
    # kappa_pcr_cases = 1 / x
    # phi_pcr_cases = 1
    phi_admitted = 1,
    kappa_admitted = 2,
    phi_death_hosp = 1,
    kappa_death_hosp = 2,
    exp_noise = exp_noise
  )
}

##' ZamCovid severity parameters
##'
##' @title ZamCovid severity parameters
##'
##' @param dt The step size
##'
##' @param severity [data.frame] object with severity data used to determine
##'   default severity values by age. First column is a string (names of
##'   severity parameter in ZamCovid) and rest of columns must be `real` numeric
##'   values for each age group (16). If `NULL`, default package parameters will
##'   be used.
##'
##' @section Time-varying parameters:
##' Every time varying parameter has the same format, which can be `NULL` (in
##'   which case the value from `severity` is used) or a list with `date` and
##'   `value` for changes in the parameter. If `value` is scalar then
##'   `date` can be `NULL` or missing. If `value` is a vector then `date`
##'   must be a vector of numeric dates of the same length as `value`.
##'
##' @param p_C Time-varying parameters for p_C (the probability of an infected
##'   individual becoming symptomatic). See Details.
##'
##' @param p_H Time-varying parameters for p_H (the probability of a symptomatic
##'   individual requiring hospitalisation). See Details.
##'
##' @param p_H_D Time-varying parameters for p_H_D (the probability of death in
##'   general beds).See Details.
##'
##' @param p_G_D Time-varying parameters for p_G_D (the probability of
##'   individuals requiring hospitalisation dying in the community or a care
##'   home). See Details.
##'
##' @param p_R Time-varying parameters for p_R (the probability of an non-
##'   fatally infected individual having immunity post-infection). See Details.
##'
##' @param p_star Time-varying parameters for p_star (the probability of
##'   patients being confirmed as covid on admission to hospital). See Details.
##'
##' @return A list of severity parameters
##'
##' @export
ZamCovid_parameters_severity <- function(dt,
                                         severity = NULL,
                                         p_C = NULL,
                                         p_H = NULL,
                                         p_H_D = NULL,
                                         p_G_D = NULL,
                                         p_R = NULL,
                                         p_star = NULL,
                                         p_sero_pos = NULL) {

  severity <- process_parameters_severity(severity)

  time_varying_severity <- list(C = p_C,
                                H = p_H,
                                H_D = p_H_D,
                                G_D = p_G_D,
                                R = p_R,
                                star = p_star)

  get_p_step <- function(x, name) {

    p_name <- paste0("p_", name)
    p <- x[[p_name]]
    time_vary <- time_varying_severity[[name]]

    if (is.null(time_vary)) {
      p_value <- NULL
      p_date <- NULL
    } else {
      p_value <- time_vary$value
      if ("date" %in% names(time_vary)) {
        p_date <- time_vary$date
      } else {
        p_date <- NULL
      }
    }

    if (!is.null(p_value)) {
      assert_proportion(p_value, p_name)
      if (length(p_value) == 1L) {
        if (length(p_date) != 0) {
          stop(sprintf(
            "As '%s' has a single 'value', expected NULL or missing 'date'",
            p_name))
        }
      } else if (length(p_date) != length(p_value)) {
        stop(sprintf("'date' and 'value' for '%s' must have the same length",
                     p_name))
      }
    }

    if (all(p == 0)) {
      psi <- p
    } else {
      psi <- p / max(p)
    }

    if (is.null(p_value)) {
      p_step <- max(p)
    } else {
      p_step <- parameters_piecewise_linear(p_date, p_value, dt)
    }

    p_step <- outer(p_step, psi)

    x[[paste0(p_name, "_step")]] <- p_step
    x[[paste0(p_name)]] <- NULL
    x[[paste0("n_", p_name, "_steps")]] <- dim(p_step)[1]

    x
  }

  for (name in names(time_varying_severity)) {
    severity <- get_p_step(severity, name)
  }

  severity
}


##' ZamCovid vaccination parameters
##'
##' @title ZamCovid vaccination parameters
##'
##' @param dt The step size
##'
##' @return A list of vaccination parameters
##'
##' @export
ZamCovid_parameters_vaccination <- function(dt,
                                            N_tot,
                                            n_groups,
                                            n_doses,
                                            rel_susceptibility,
                                            rel_p_sympt,
                                            rel_p_hosp_if_sympt,
                                            rel_p_death,
                                            rel_infectivity,
                                            vaccine_progression_rate,
                                            vaccine_schedule,
                                            vaccine_index_dose2,
                                            vaccine_catchup_fraction = 1) {

  stopifnot(length(N_tot) == n_groups)
  calc_n_vacc_classes <- function(x) {
    if (is.array(x)) nlayer(x) else length(x)
  }

  assert_proportion(rel_susceptibility)

  rel_params <- list(rel_susceptibility = rel_susceptibility,
                     rel_p_sympt = rel_p_sympt,
                     rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
                     rel_p_death = rel_p_death,
                     rel_infectivity = rel_infectivity)

  n <- vnapply(rel_params, calc_n_vacc_classes)

  if (any(n > 1) && length(unique(n[n > 1])) != 1) {
    msg1 <- paste(names(rel_params), collapse = ", ")
    msg2 <- "should have the same dimension"
    stop(paste(msg1, msg2))
  }
  n_vacc_classes <- max(n)

  ret <- Map(function(value, name)
    build_rel_param(n_groups, value, n_vacc_classes, name),
    rel_params, names(rel_params))

  if (is.null(vaccine_schedule)) {
    if (!is.null(vaccine_index_dose2) && vaccine_index_dose2 != 1L) {
      stop("'vaccine_index_dose2' set without schedule")
    }
    ret$vaccine_dose_step <- array(0, c(n_groups, n_doses, 1))
    ret$index_dose <- rep(1L, n_doses)

  } else {

    assert_is(vaccine_schedule, "vaccine_schedule")
    vaccine_index_dose2 <- vaccine_index_dose2 %||% 1L
    if (vaccine_index_dose2 > n_vacc_classes) {
      stop(sprintf(
        "Invalid value for 'vaccine_index_dose2', must be in [1, %d]",
        n_vacc_classes))
    }
    stopifnot(vaccine_schedule$n_doses == n_doses)

    n_days <- dim(vaccine_schedule$doses)[[3]]
    i <- rep(seq_len(n_days), each = 1 / dt)
    len <- vaccine_schedule$date / dt
    ret$index_dose <- c(1L, vaccine_index_dose2)

    ret$vaccine_dose_step <- mcstate::array_bind(
      array(0, c(n_groups, n_doses, len)),
      (vaccine_schedule$doses * dt)[, , i])

  }

  ret$index_dose_inverse <- create_index_dose_inverse(n_vacc_classes,
                                                      ret$index_dose)

  ret$n_vacc_classes <- n_vacc_classes
  ret$vaccine_progression_rate_base <- build_vaccine_progression_rate(
    vaccine_progression_rate, n_vacc_classes, ret$index_dose)


  assert_scalar(vaccine_catchup_fraction)
  assert_proportion(vaccine_catchup_fraction)
  ret$vaccine_catchup_fraction <- vaccine_catchup_fraction

  ret$n_doses <- n_doses

  ret
}



##' Create initial conditions for the ZamCovid model. This matches the
##' interface required for mcstate
##'
##' @title Initial conditions for the ZamCovid model
##'
##' @return A numeric vector of initial conditions
##' @export
##' @examples
##' p <- ZamCovid_parameters(numeric_date("2020-02-07"))
##' mod <- ZamCovid$new(p, 0, 10)
##' ZamCovid_initial(mod$info(), 10, p, 4L)
ZamCovid_initial <- function(info, n_particles, pars) {

  index <- info$index
  state <- numeric(info$len)

  index_S <- index[["S"]]
  index_N_tot <- index[["N_tot"]]
  index_susceptible <- index[["susceptible"]]
  index_S_no_vacc <- index_S[seq_len(length(pars$N_tot))]
  index_N_tot_sero <- index[["N_tot_sero"]][[1L]]
  index_N_tot_PCR <- index[["N_tot_PCR"]][[1L]]
  index_I_weighted <- index[["I_weighted"]][[1L]] + pars$seed_age_band - 1L

  ## S0 is the population totals
  initial_S <- pars$N_tot

  state[index_N_tot] <- pars$N_tot
  state[index_susceptible] <- sum(pars$N_tot)
  state[index_S_no_vacc] <- initial_S
  state[index_N_tot_sero] <- sum(pars$N_tot)
  state[index_N_tot_PCR] <- sum(pars$N_tot)
  state[index_I_weighted] <- 1

  state
}


##' Check that data for the particle filter has required columns and apply
##' constraints for fitting to it.
##'
##' @title Check data for particle filter
##' @param data A data.frame of data
##' @return Invisibly, a data frame, identical to `data`
##' @export
ZamCovid_check_data <- function(data) {
  ## TODO: develop constraints here, e.g. not allowing to fit to aggregated and
  ## age-disaggregated data, etc.
  invisible(data)
}


#' Define basic model parameters
#'
#' @param model_end A positive integer, greated than `model_start` for the model
#'    end time.
#'
#' @param steps_per_day A positive integer for the number of steps per day to
#'    run the model.
#'
#' @param model_start Either 0 or a positive integer for the initial start time
#'    for the model.
#'
#' @param contact_matrix Either `NULL` or a `data.frame` of a symmetric
#'    individual-to-individual contact rate matrix. If `NULL` defaults to
#'    nationally representative values.
#'
#' @param population Either `NULL` or a `data.frame` with columns `age_group`
#'    and `n` of same length (age groups) as `contact_matrix` supplied. If
#'    `NULL` defaults to nationally representative values.
#'
#' @param seed_pars Either `NULL` or a `list` with elements `seed_group`,
#'    `seed_number` and `seed_time` for the `age_group` in which initial
#'    infections will be seeded, the number of infection to seed and the model
#'    time in which they will be seeded, respectively. If `NULL` will
#'    default to seeding 5 infection in `age_group` 5 at `t = 0`.
#'
#' @param pars_list Either `NULL` or a `list` of all other model parameters to
#'    fix. If `NULL` default values will be supplied.
#'
#' @return Return a `list` of parameters to run the basic ZamCovid model.
#'
#' @export
#'
#' @examples basic_parameters(model_end = 200)
basic_parameters <- function(model_end,
                                steps_per_day = 1,
                                model_start = 0,
                                contact_matrix = NULL,
                                population = NULL,
                                seed_pars = NULL,
                                par_list = NULL) {

  if(is.null(steps_per_day)) {
    stop("A number of steps per day must be provided!")
  }

  if(is.null(model_end)) {
    stop("A model_end value must be provided!")
  }

  if (!is.null(contact_matrix)) {
    if(is.null(population)) {
      stop("A `population` must be provided with contact matrix!")
    }

  } else {
    population <- read_csv(ZamCovid_file("extdata/population.csv"))
    contact_matrix <- read_csv(ZamCovid_file("extdata/matrix.csv"))
  }

  pop <- as.numeric(population$n)
  contact_matrix <- contact_matrix / rep(pop, each = ncol(contact_matrix))

  N_age <- nrow(contact_matrix)

  # Initialise model population
  S_ini <- as.numeric(pop)
  E_ini <- rep(0, N_age)
  E2_ini <- rep(0, N_age)
  R_ini <- rep(0, N_age)

  # Seed infections
  if (is.null(seed_pars)) {
    I_ini <- rep(0, N_age)
    I_ini[5] <- 5
    S_ini[5] <- S_ini[5] - 5
  } else {
    I_ini <- rep(0, N_age)
    I_ini[seed_pars$seed_group] <- seed_pars$seed_number
    S_ini[seed_pars$seed_group] <- S_ini[seed_pars$seed_group] -
      seed_pars$seed_number
  }



  # If no parameters provided, use the following as default values
  if (is.null(par_list)){
    # Transition parameters
    gamma_E <- 1 / 3
    gamma_I <- 1 / 8
    eta_R <- 1 / (365 * 3)

    # Observation parameters
    rho <- 0.3
    kappa <- 0.2

    # Time time-varying parameters
    beta_date <- c(model_start, 20, 50, 100, 150, model_end)
    beta_value <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.8)


  } else {
    # Specified parameters will be fixed in the model

    gamma_E <- ifelse(!is.null(par_list$gamma_E),
                      par_list$gamma_E,
                      1 / 3)
    gamma_I <- ifelse(!is.null(par_list$gamma_I),
                      par_list$gamma_I,
                      1 / 8)
    eta_R <- ifelse(!is.null(par_list$eta_R),
                    par_list$eta_R,
                    1 / (365 * 3))
    rho <- ifelse(!is.null(par_list$rho),
                  par_list$rho,
                  0.3)
    kappa <- ifelse(!is.null(par_list$kappa),
                    par_list$kappa,
                    0.2)

    beta_date <- ifelse(!is.null(par_list$beta_date),
                        par_list$beta_date,
                        c(model_start, 20, 50, 100, 150, model_end))
    beta_value <- ifelse(!is.null(par_list$beta_value),
                         par_list$beta_value,
                         c(0.1, 0.2, 0.3, 0.4, 0.5, 0.8))
  }

  # Create piece-wise linear function for time time-varying parameters
  beta_step <- parameters_piecewise_linear(
    beta_date, beta_value, 1 / steps_per_day)

  # Final list of all user-defined parameters, vectors and arrays
  list(steps_per_day = steps_per_day,
       N_age = N_age,
       m = as.matrix(contact_matrix),
       beta_step = beta_step,
       S_ini = S_ini,
       E_ini = E_ini,
       E2_ini = E2_ini,
       I_ini = I_ini,
       R_ini = R_ini,
       gamma_E = gamma_E,
       gamma_I = gamma_I,
       eta_R = eta_R,
       rho = rho,
       kappa = kappa)
}


parameters_piecewise_linear <- function (date, value, dt) {
  if (!inherits(value, "matrix")) {
    value <- matrix(value, ncol = 1)
  }
  if (is.null(date)) {
    if (nrow(value) != 1L) {
      stop("As 'date' is NULL, expected single value")
    }
    if (ncol(value) == 1) {
      value <- as.numeric(value)
    }
    return(value)
  }
  if (length(date) != nrow(value)) {
    stop("'date' and 'value' must have the same length")
  }
  if (length(date) < 2) {
    stop("Need at least two dates and values for a varying piecewise linear")
  }
  if (!is.numeric(date)) {
    stop("'date' must be numeric - did you forget numeric_date()?")
  }
  if (any(date < 0)) {
    stop("Negative dates, numeric_date likely applied twice")
  }
  if (any(diff(date) <= 0)) {
    stop("Dates must be strictly increasing")
  }
  if (date[[1]] != 0) {
    date <- c(0, date)
    value <- rbind(value[1, ], value)
  }
  value <- apply(value, 2, function(x) {
    stats::approx(date, x, seq(0, date[[length(date)]], by = dt))$y
  })
  if (ncol(value) == 1) {
    value <- as.numeric(value)
  }
  value
}


parameters_piecewise_constant <- function (date, value, dt) {

  if (is.null(date)) {
    if (length(value) != 1L) {
      stop("As 'date' is NULL, expected single value")
    }
    return(value)
  }

  if (length(date) != length(value)) {
    stop("'date' and 'value' must have the same length")
  }

  if (!is.null(date)) {
    if (date[1L] != 0) {
      stop("As 'date' is not NULL, first date should be 0")
    }
  }

  times <- seq(0, date[[length(date)]], by = dt)
  stats::approx(date, value, method = "constant",
                xout = times)$y
}



ZamCovid_parameters_expand_step <- function(step, value_step) {
  # Convert between the values passed to parameters_piecewise_linear()
  # and the actual values for a given set of steps.
  value_step[pmin(step, length(value_step) - 1L) + 1L]
}


#' Recover fitted parameters model samples
#'
#' @param samples A `samples` object from a pMCMC model run.
#' @param priors `list` of priors used for model run.
#'
#' @return A `list` of new model priors and vcv matrix.
#' @export
#'
#' @examples basic_new_pars(samples, priors)
basic_new_pars <- function(samples, priors) {

  i <- which.max(samples$probabilities[, "log_posterior"])
  initial <- samples$pars[i, ]

  vcv <- cov(samples$pars)
  dimnames(vcv) <- NULL

  for (i in names(priors)) {
    nm <- priors[[i]]$name
    priors[[i]]$initial <- initial[[nm]]
  }

  list(vcv = vcv,
       priors = priors)
}


numeric_date <- function(date) {
  days_into_2020 <- as.numeric(as.Date(date) - as.Date("2019-12-31"))
  if (any(days_into_2020 < 0)) {
    stop("Negative dates, numeric_date likely applied twice")
  }
  days_into_2020
}


numeric_date_as_date <- function(date) {
  assert_numeric_date(date)
  as.Date("2019-12-31") + date
}


assert_numeric_date <- function(date) {
  if (!is.numeric(date)) {
    stop("'date' must be numeric - did you forget sircovid_date()?")
  }
  if (any(date < 0)) {
    stop("Negative dates, sircovid_date likely applied twice")
  }
  date
}


seed_over_steps <- function(start_step, weights) {
  ## The weights vector must over steps, not dates
  weights <- weights / sum(weights)
  p <- start_step %% 1
  if (p == 0) {
    ret <- weights
  } else{
    ret <- p * c(0, weights) + (1 - p) * c(weights, 0)
  }
  ret
}


process_parameters_severity <- function(params) {

  if (is.null(params)) {
    params <- read_csv(ZamCovid_file("extdata/severity_default.csv"))
  }

  ## Transpose so columns are parameters, rownames are age groups
  data <- t(as.matrix(params[-1L]))
  colnames(data) <- params[[1]]
  data <- cbind(age = rownames(data),
                data.frame(data, check.names = FALSE),
                stringsAsFactors = FALSE)
  rownames(data) <- NULL

  required <- c(
    p_C = "p_C",
    p_H = "p_H",
    p_H_D = "p_H_D",
    p_sero_pos = "p_sero_pos",
    p_G_D = "p_G_D",
    p_R = "p_R",
    p_star = "p_star")
  data <- rename(data, required, names(required))

  list(
    p_star = data[["p_star"]],
    p_C = data[["p_C"]],
    p_G_D = data[["p_G_D"]],
    p_H_D = data[["p_H_D"]],
    p_sero_pos = data[["p_sero_pos"]],
    p_H = data[["p_H"]],
    p_R = data[["p_R"]])
}


check_severity <- function(pars) {
  check_parameters <- function(p_step, rel_p) {
    if (is.null(pars[[rel_p]])) {
      stop(sprintf("Parameter '%s' is missing", rel_p))
    } else if (is.null(pars[[p_step]])) {
      stop(sprintf("Parameter '%s' is missing", p_step))
    }

    if (!(ncol(pars[[rel_p]]) %in% c(1, 2, 5))) {
      stop(sprintf("%s should have 1 column", rel_p))
    }

    assert_non_negative(pars[[rel_p]], rel_p)
    assert_non_negative(pars[[p_step]], p_step)
  }

  step_pars <- c("p_C_step", "p_H_step", "p_H_D_step", "p_G_D_step", "p_R_step")
  rel_pars <- c("rel_p_sympt", "rel_p_hosp_if_sympt", "rel_p_H_D",
                "rel_p_G_D", "rel_p_R")

  Map(check_parameters,
      p_step = step_pars,
      rel_p = rel_pars
  )
  invisible(pars)
}


build_rel_param <- function(n_groups, rel_param, n_vacc_classes, name_param) {

  if (length(rel_param) == 1) {
    mat_rel_param <- array(rel_param, c(n_groups, n_vacc_classes))
  } else if (is.array(rel_param)) {
    if (length(dim(rel_param)) != 2) {
      stop(paste(name_param, "should be a two dimensional array with",
                 "dimensions: age groups, vaccine classes"))
    }
    if (nrow(rel_param) != n_groups) {
      stop(paste(name_param, "should have as many rows as age groups"))
    }
    if (dim(rel_param)[2] != n_vacc_classes) {
      stop(paste(name_param,
                 "should have number of vaccine classes as 3rd dimension"))
    }
    mat_rel_param <- rel_param
  } else { # create array by repeating rel_param for each age group and strain
    mat_rel_param <-
      array(rep(rel_param, each = n_groups),
            dim = c(n_groups, n_vacc_classes))
  }
  mat_rel_param
}



make_contact_matrix <- function(file, N_tot) {

  m <- read_csv(ZamCovid_file(file))
  m <- m / rep(N_tot, each = ncol(m))
  age_groups <- paste0("[", seq(0, 75, 5), ",", c(seq(5, 75, 5), 120), ")")
  colnames(m) <- age_groups
  rownames(m) <- age_groups
  m
}
