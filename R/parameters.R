#' Title
#'
#' @param date
#' @param value
#' @param dt
#'
#' @return
#' @export
#'
#' @examples
zamcovid_parameters <- function(steps_per_day,
                                date_start,
                                date_end,
                                contact_matrix = NULL,
                                population = NULL,
                                seed_group = 5,
                                seed_number = 5,
                                par_list = NULL) {

  if(is.null(steps_per_day)) {
    stop("A number of steps per day must be provided!")
  }

  if (!is.null(contact_matrix)) {
    if(is.null(population)) {
      stop("A `population` must be provided with contact matrix!")
    }

  } else {
    population <- read.csv("data/population.csv", stringsAsFactors = FALSE)
    contact_matrix <- read.csv("data/matrix.csv", stringsAsFactors = FALSE)
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
  I_ini <- rep(0, N_age)
  I_ini[seed_group] <- seed_number
  S_ini[seed_group] <- S_ini[seed_group] - seed_number


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
    beta_date <- c(date_start, 20, 50, 100, 150, date_end)
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
                        c(date_start, 20, 50, 100, 150, date_end))
    beta_value <- ifelse(!is.null(par_list$beta_value),
                         par_list$beta_value,
                         c(0.1, 0.2, 0.3, 0.4, 0.5, 0.8))
  }

  # Create piece-wise linear function for time time-varying parameters
  beta_step <- sircovid:::sircovid_parameters_piecewise_linear(
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


#' Title
#'
#' @param info
#'
#' @return
#' @export
#'
#' @examples
zamcovid_index <- function(info) {

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


piecewise_linear <- function (date, value, dt) {
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
    stop("'date' must be numeric - did you forget sircovid_date()?")
  }
  if (any(date < 0)) {
    stop("Negative dates, sircovid_date likely applied twice")
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


get_new_pars <- function(samples, priors) {

  i <- which.max(samples$probabilities[, "log_posterior"])
  initial <- samples$pars[i, ]

  vcv <- cov(samples$pars)
  dimnames(vcv) <- NULL

  for (i in 1:length(priors)) {
    nm <- priors[[i]]$name
    priors[[i]]$initial <- initial[[nm]]
  }

  list(vcv = vcv,
       priors = priors)
}
