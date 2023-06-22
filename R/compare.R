#' Compare function for ZamCovid model.
#'
#' @param state State variable outputs from particle filter of generated
#'    basic model.
#'
#' @param observed Space variables from data used for particle filter.
#'
#' @param pars A `list` of `mcstate::pmcmc_parameters` to fit. Note these cannot
#'    be directly supplied to the function, but rather as an argument to
#'    `mcstate::pmcmc()`.
#'
#' @return A numerical value for the negative log-likelihood of the
#'    state-space model.
#'
#' @export
#'
#' @examples basic_compare(state, observed)
ZamCovid_compare <- function(state, observed, pars) {

  model_admissions_conf <- state["admitted_inc", ]
  model_sero_pos_all <- state["sero_pos_all", ]
  model_sero_pos_over15 <- state["sero_pos_over15", ]
  model_sero_pos_15_19 <- state["sero_pos_15_19", ]
  model_sero_pos_20_29 <- state["sero_pos_20_29", ]
  model_sero_pos_30_39 <- state["sero_pos_30_39", ]
  model_sero_pos_40_49 <- state["sero_pos_40_49", ]
  model_sero_pos_50_plus <- state["sero_pos_50_plus", ]

  model_base_deaths <- state["base_death_inc", ]
  model_deaths_hosp <- state["deaths_hosp_inc", ]
  model_deaths_comm <- state["deaths_comm_inc", ]
  model_deaths_all <- state["deaths_all_inc", ]


  ## Serology assay
  # It is possible that model_sero_pos > pars$N_tot; this is capped here to
  # avoid probabilities > 1
  model_sero_capped_all <- pmin(model_sero_pos_all, pars$N_tot)
  model_sero_capped_over15 <- pmin(model_sero_pos_over15, pars$N_tot_over15)
  model_sero_capped_15_19 <- pmin(model_sero_pos_15_19, pars$N_tot_15_19)
  model_sero_capped_20_29 <- pmin(model_sero_pos_20_29, pars$N_tot_20_29)
  model_sero_capped_30_39 <- pmin(model_sero_pos_30_39, pars$N_tot_30_39)
  model_sero_capped_40_49 <- pmin(model_sero_pos_40_49, pars$N_tot_40_49)
  model_sero_capped_50_plus <- pmin(model_sero_pos_50_plus, pars$N_tot_50_plus)

  model_sero_prob_pos_all <-
    test_prob_pos(model_sero_capped_all,
                  pars$N_tot - model_sero_capped_all,
                  pars$sero_sensitivity,
                  pars$sero_specificity,
                  pars$exp_noise)

  model_sero_prob_pos_over15 <-
    test_prob_pos(model_sero_capped_over15,
                  pars$N_tot_over15 - model_sero_capped_over15,
                  pars$sero_sensitivity,
                  pars$sero_specificity,
                  pars$exp_noise)

  model_sero_prob_pos_15_19 <-
    test_prob_pos(model_sero_capped_15_19,
                  pars$N_tot_15_19 - model_sero_capped_15_19,
                  pars$sero_sensitivity,
                  pars$sero_specificity,
                  pars$exp_noise)

  model_sero_prob_pos_20_29 <-
    test_prob_pos(model_sero_capped_20_29,
                  pars$N_tot_20_29 - model_sero_capped_20_29,
                  pars$sero_sensitivity,
                  pars$sero_specificity,
                  pars$exp_noise)

  model_sero_prob_pos_30_39 <-
    test_prob_pos(model_sero_capped_30_39,
                  pars$N_tot_30_39 - model_sero_capped_30_39,
                  pars$sero_sensitivity,
                  pars$sero_specificity,
                  pars$exp_noise)

  model_sero_prob_pos_40_49 <-
    test_prob_pos(model_sero_capped_40_49,
                  pars$N_tot_40_49 - model_sero_capped_40_49,
                  pars$sero_sensitivity,
                  pars$sero_specificity,
                  pars$exp_noise)

  model_sero_prob_pos_50_plus <-
    test_prob_pos(model_sero_capped_50_plus,
                  pars$N_tot_50_plus - model_sero_capped_50_plus,
                  pars$sero_sensitivity,
                  pars$sero_specificity,
                  pars$exp_noise)


  ## Log-likelihood
  ll_admitted <- ll_nbinom(observed$hosp_admissions,
                           pars$phi_admitted * model_admissions_conf,
                           pars$kappa_admitted, pars$exp_noise)

  ll_deaths_hosp <- ll_nbinom(observed$deaths_hosp,
                              pars$phi_death_hosp * model_deaths_hosp,
                              pars$kappa_death_hosp, pars$exp_noise)

  ll_deaths_all <- ll_nbinom(observed$deaths_all,
                             pars$phi_death_all * model_deaths_all,
                             pars$kappa_death_all, pars$exp_noise)

  ll_serology_all <- ll_binom(observed$sero_pos_all,
                              observed$sero_tot_all,
                              model_sero_prob_pos_all)

  ll_serology_over15 <- ll_binom(observed$sero_pos_over15,
                                 observed$sero_tot_over15,
                                 model_sero_prob_pos_over15)

  ll_serology_15_19 <- ll_binom(observed$sero_pos_15_19,
                                 observed$sero_tot_15_19,
                                 model_sero_prob_pos_15_19)

  ll_serology_20_29 <- ll_binom(observed$sero_pos_20_29,
                                 observed$sero_tot_20_29,
                                 model_sero_prob_pos_20_29)

  ll_serology_30_39 <- ll_binom(observed$sero_pos_30_39,
                                 observed$sero_tot_30_39,
                                 model_sero_prob_pos_30_39)

  ll_serology_40_49 <- ll_binom(observed$sero_pos_40_49,
                                 observed$sero_tot_40_49,
                                 model_sero_prob_pos_40_49)

  ll_serology_50_plus <- ll_binom(observed$sero_pos_50_plus,
                                 observed$sero_tot_50_plus,
                                 model_sero_prob_pos_50_plus)

  ll_admitted + ll_deaths_all + ll_serology_all + ll_serology_over15 +
    ll_serology_15_19 + ll_serology_20_29 + ll_serology_30_39 +
    ll_serology_40_49 + ll_serology_50_plus
}


#' Compare function for basic model.
#'
#' @param state State variable outputs from particle filter of generated
#'    basic model.
#'
#' @param observed Space variables from data used for particle filter.
#'
#' @param pars A `list` of `mcstate::pmcmc_parameters` to fit. Note these cannot
#'    be directly supplied to the function, but rather as an argument to
#'    `mcstate::pmcmc()`.
#'
#' @return A numerical value for the negative log-likelihood of the
#'    state-space model.
#'
#' @export
#'
#' @examples basic_compare(state, observed)
basic_compare <- function(state, observed, pars = NULL) {

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


#' Make transform function for basic model parameters
#'
#' @param baseline A `list` containing all baseline parameters.
#'
#' @return A `list` of `odin` parameters for the basic model.
#' @export
#'
#' @examples basic_transform(baseline)
basic_transform <- function(baseline) {

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
