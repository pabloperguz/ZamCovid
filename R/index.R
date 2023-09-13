#' Index function for ZamCovid model
#'
#' @param info A `list` object containing `info` from generated ZamCovid model.
#'
#' @return A `list` with elements `run` and `state`
#'
#' @export
#'
#' @examples basic_index(info)
ZamCovid_index <- function(info) {

  index <- info$index

  # Model states required for the particle filter to run:
  index_core <- c(infections_inc = index[["infections_inc"]],
                  reinfections_inc = index[["reinfections_inc"]],
                  admitted_inc = index[["admit_conf_inc"]],
                  deaths_hosp_inc = index[["hosp_deaths_inc"]],
                  deaths_comm_inc = index[["comm_deaths_inc"]],
                  base_death_inc = index[["base_death_inc"]],
                  deaths_all_inc = index[["all_deaths_inc"]],
                  sero_pos_all = index[["sero_pos_all"]],
                  sero_pos_over15 = index[["sero_pos_over15"]],
                  sero_pos_15_19 = index[["sero_pos_15_19"]],
                  sero_pos_20_29 = index[["sero_pos_20_29"]],
                  sero_pos_30_39 = index[["sero_pos_30_39"]],
                  sero_pos_40_49 = index[["sero_pos_40_49"]],
                  sero_pos_50_plus = index[["sero_pos_50_plus"]],
                  inf_cum_all = index[["inf_cum_all"]],
                  inf_cum_over15 = index[["inf_cum_over15"]],
                  inf_cum_15_19 = index[["inf_cum_15_19"]],
                  inf_cum_20_29 = index[["inf_cum_20_29"]],
                  inf_cum_30_39 = index[["inf_cum_30_39"]],
                  inf_cum_40_49 = index[["inf_cum_40_49"]],
                  inf_cum_50_plus = index[["inf_cum_50_plus"]])

  # An index of only incidence versions for the likelihood function
  # For now, this is basically the same as index_core, but we will need to
  # modify this when we introduce new data to fit the model to, like serology
  index_run <- c(time = index[["time"]], index_core)

  # Other states we want to save for post-processing
  index_save <- c(hosp = index[["hosp_tot"]],
                  D_hosp_tot = index[["D_hosp_tot"]])


  # States disaggregated by age and/or vaccine class
  # Be careful of what disaggregated states you extract, as too many can easily
  # lead to memory issues when post-processing!
  age_suffix <- paste0("_", model_age_bins()$start)
  n_vacc_classes <- info$dim$S[[2]]

  index_S <- calculate_index(index, "S", list(n_vacc_classes), age_suffix)

  index_R <- calculate_index(index, "R", list(n_vacc_classes), age_suffix)

  index_D <- calculate_index(index, "D", list(n_vacc_classes), age_suffix)

  index_inf_inc_age <- calculate_index(index, "infections_inc_age",
                                       list(), age_suffix)

  index_reinf_inc_age <- calculate_index(index, "reinfections_inc_age",
                                       list(), age_suffix)

  index_severity <- c(ifr = index[["ifr"]],
                      calculate_index(index, "ifr_age",
                                      list(), age_suffix))

  index_state <- c(index_core, index_save, index_S, index_R, index_D,
                   index_severity, index_inf_inc_age, index_reinf_inc_age)

  list(run = index_run, state = index_state)
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


## Calculates the index of a given state and adds a suffix corresponding to
##  how the state is disaggregated (e.g. _V1, _V2 for two vaccine classes
calculate_index <- function(index, state, suffix_list, suffix0 = NULL,
                            state_name = state) {
  if (is.null(suffix0)) {
    suffixes <- list()
  } else {
    suffixes <- list(suffix0)
  }
  for (i in seq_along(suffix_list)) {
    nm <- names(suffix_list)[[i]]
    if (length(nm) == 0) {
      nm <- ""
    }
    suffixes <- c(suffixes,
                  list(c("", sprintf("_%s%s", nm,
                                     seq_len(suffix_list[[i]] - 1L)))))
  }
  suffixes <- expand.grid(suffixes)
  nms <- apply(suffixes, 1,
               function(x) sprintf("%s%s",
                                   state_name, paste0(x, collapse = "")))

  set_names <- function(x, nms) {
    names(x) <- nms
    x
  }

  set_names(index[[state]], nms)
}


## We will always use these age bands, so rather than detect them, we will
## check that things conform to them.
model_age_bins <- function() {
  end <- c(seq(4, 75, by = 5), 120)
  start <- c(0, end[-length(end)] + 1L)
  list(start = start, end = end)
}
