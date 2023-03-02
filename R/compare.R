#' Compare function for basic ZamCovid model.
#'
#' @param state State variable outputs from particle filter of generated
#'    ZamCovid model.
#'
#' @param observed Space variables from data used for particle filter.
#'
#' @param pars A `list` of `mcstate::pmcmc_parameters` to fit. Note these cannot
#'    be directly supplied to the function, but rather as an argument to
#'    `mcstate::pmcmc()`.
#'
#' @return A numerical value for the negative log=likelihood of the
#'    state-space model.
#'
#' @export
#'
#' @examples
zamcovid_compare <- function(state, observed, pars = NULL) {

  ll_nbinom <- function(data, model, kappa, exp_noise = 1e6) {

    if (is.na(data)) {
      return(numeric(length(model)))
    }

    mu <- model + rexp(length(model), rate = exp_noise)

    dnbinom(data, kappa, mu = mu, log = TRUE)
  }

  model_cases_under_15 <- state["cases_under_15", , drop = TRUE]
  model_cases_15_19 <- state["cases_15_19", , drop = TRUE]
  model_cases_20_29 <- state["cases_20_29", , drop = TRUE]
  model_cases_30_39 <- state["cases_30_39", , drop = TRUE]
  model_cases_40_49 <- state["cases_40_49", , drop = TRUE]
  model_cases_50_plus <- state["cases_50_plus", , drop = TRUE]

  ll_cases_under_15 <- ll_nbinom(observed$cases_under_15,
                                 model_cases_under_15, pars$kappa)
  ll_cases_15_19 <- ll_nbinom(observed$cases_15_19,
                              model_cases_15_19, pars$kappa)
  ll_cases_20_29 <- ll_nbinom(observed$cases_20_29,
                              model_cases_20_29, pars$kappa)
  ll_cases_30_39 <- ll_nbinom(observed$cases_30_39,
                              model_cases_30_39, pars$kappa)
  ll_cases_40_49 <- ll_nbinom(observed$cases_40_49,
                              model_cases_40_49, pars$kappa)
  ll_cases_50_plus <- ll_nbinom(observed$cases_50_plus,
                                model_cases_50_plus, pars$kappa)

  ll_cases_under_15 + ll_cases_15_19 + ll_cases_20_29 + ll_cases_30_39 +
    ll_cases_40_49 + ll_cases_50_plus
}


#' Make transform function for ZamCovid basic model parameters
#'
#' @param baseline A `list` containing all baseline parameters.
#'
#' @return A `list` of `odin` parameters for the basic model.
#' @export
#'
#' @examples
zamcovid_transform <- function(baseline) {

  function(theta) {

    steps_per_day <- baseline$steps_per_day
    N_age <- baseline$N_age
    m <- baseline$m
    beta_step <- baseline$beta_step
    S_ini <- baseline$S_ini
    E_ini <- baseline$E_ini
    E2_ini <- baseline$E2_ini
    I_ini <- baseline$I_ini
    R_ini <- baseline$R_ini
    gamma_I <- baseline$gamma_I
    eta_R <- baseline$eta_R
    kappa <- baseline$kappa

    list(
      # Fixed parameters
      steps_per_day = steps_per_day,
      N_age = N_age,
      m = m,
      beta_step = beta_step,
      S_ini = S_ini,
      E_ini = E_ini,
      E2_ini = E2_ini,
      I_ini = I_ini,
      R_ini = R_ini,
      gamma_I = gamma_I,
      eta_R = eta_R,
      kappa = kappa,

      # Parameters to fit
      rho = theta[["rho"]],
      gamma_E = theta[["gamma_E"]])
  }
}
