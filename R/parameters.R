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
#' @return Return a `list` of parameters to run the basic ZamCovid model.
#'
#' @export
ZamCovid_parameters <- function(start_date,
                                steps_per_day = 4L,
                                beta_date = NULL,
                                beta_value = NULL,
                                population = NULL,
                                contact_matrix = NULL,
                                N_tot = NULL,
                                n_groups = NULL,
                                severity = NULL,
                                progression = NULL,
                                vaccination = NULL,
                                initial_seed_pattern = 1.0,
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
                                rel_infectivity = 1,
                                vaccine_progression_rate = NULL,
                                vaccine_schedule = NULL,
                                vaccine_index_dose2 = NULL,
                                vaccine_index_booster = NULL,
                                vaccine_catchup_fraction = 1,
                                n_doses = 2L) {

  if (!is.character(start_date)) {
    stop("start_date must be character string in format yyyy-mm-dd!")
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

    if (is.null(N_tot)) {stop("Expected N_tot vector!")}
    if (is.null(n_groups)) {stop("Expected n_groups!")}

  } else {
    population <- read_csv(ZamCovid_file("extdata/population.csv"))
    N_tot <- as.numeric(population$n)
    contact_matrix <-
      make_contact_matrix("extdata/matrix.csv", N_tot)
    n_groups <- nrow(contact_matrix)
  }

  dt <- 1 / steps_per_day

  ## First epidemic wave seed
  start_step <- numeric_date(start_date) * steps_per_day
  initial_seed_value <-
    seed_over_steps(start_step, initial_seed_pattern) * initial_seed_size


  ## Make piece-wise linear beta_step
  # beta_date <- numeric_date(beta_date)
  beta_step <- parameters_piecewise_linear(beta_date,
                                           beta_value %||% 0.1, dt)


  ## Default parameters
  ret <- list(
    steps_per_day = steps_per_day,
    dt = dt,
    beta_step = beta_step,

    m = as.matrix(contact_matrix),
    n_groups = n_groups, #  16 n_groups (5 year age bands and 75+)

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
                                                     vaccine_index_dose2,
                                                     vaccine_index_booster)

  ret <- c(ret, severity, progression, vaccination)

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
              k_H_R = 2,
              k_H_D = 2,
              k_G_D = 2,
              k_sero_pre = 2,
              k_sero_pos = 2,
              k_PCR_pre = 2,
              k_PCR_pos = 2,

              gamma_E = 1 / (3.42 / 2),
              gamma_A = 1 / 2.88,
              gamma_P = 1 / 1.68,
              gamma_C_1 = 1 / 2.14,
              gamma_C_2 = 1 / 1.86,
              gamma_H_D = 2 / 5,
              gamma_H_R = 2 / 10,
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
                                         p_star = NULL) {

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
                                            vaccine_index_booster,
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
    stop("Let's cross that bridge when we get to it!")
  }

  ret$index_dose_inverse <- 2L

  ret$n_vacc_classes <- n_vacc_classes
  ret$vaccine_progression_rate_base <- matrix(0, n_groups, n_vacc_classes)

  assert_scalar(vaccine_catchup_fraction)
  assert_proportion(vaccine_catchup_fraction)
  ret$vaccine_catchup_fraction <- vaccine_catchup_fraction

  ret$n_doses <- n_doses

  ret
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


#' Index function for basic model
#'
#' @param info A `list` object containing `info` from generated ZamCovid model.
#'
#' @return A `list` with elements `run` and `state`
#'
#' @export
#'
#' @examples basic_index(info)
basic_index <- function(info) {

  index <- info$index

  # An index of all model states required for the particle filter
  index_core <- c(cases_under_15 = index[["cases_under_15"]],
                  cases_15_19 = index[["cases_15_19"]],
                  cases_20_29 = index[["cases_20_29"]],
                  cases_30_39 = index[["cases_30_39"]],
                  cases_40_49 = index[["cases_40_49"]],
                  cases_50_plus = index[["cases_50_plus"]])

  # An index of only incidence versions for the likelihood function
  # For now, this is basically the same as index_core, but we will need to
  # modify this when we introduce new data to fit the model to, like serology
  index_run <- c(time = index[["time"]], index_core)

  # Similarly, we might introduce an additional index_save if we have other
  # variables we sometimes, but not always, want to save for post-processing,
  # such as model states by vaccination or strain classes
  # index_save <- NULL
  index_state <- index_core

  list(run = index_run, state = index_state)
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

    if (!(ncol(pars[[rel_p]]) %in% c(1, 4))) {
      stop(sprintf("%s should have 1 or 4 columns", rel_p))
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
      stop(paste(name_param, "should be a three dimensional array with",
                 "dimensions: age groups, strains, vaccine classes"))
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
