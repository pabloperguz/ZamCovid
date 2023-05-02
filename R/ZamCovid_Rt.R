##' Compute "Rt" for a single simulated trajectory and parameter set.
##'
##' @title Compute "Rt"
##'
##' @param time A vector of time steps that the model was run over
##'
##' @param S A (n groups x n vaccine classes) x steps matrix of "S"
##'   compartment counts
##'
##' @param p A [ZamCovid_parameters()] object
##'
##' @param type A character vector of possible Rt types to
##'   compute. Can be `eff_Rt_general` or `Rt_general`, for Rt accounting or
##'   not accounting for immunity in the population, respectively.
##'
##' @param eigen_method The eigenvalue method to use (passed to
##'   [eigen1::eigen1] as `method`)
##'
##' @param R A (n groups x n strains x n vaccine classes) x time steps matrix of
##'   "R" compartment counts, required for multi-strain models.
##'
##' @param weight_Rt If `TRUE` then computes the weighted average
##'  of the Rt for all strains, otherwise all calculations are returned with
##'  an additional dimension to index each strain.
##'
##' @return A list with elements `time`, `beta`, and any of the `type`
##'   values specified above.
##'
##' @export
ZamCovid_Rt <- function(time, S, p,
                        type = NULL,
                        eigen_method = "power_iteration", R = NULL,
                        weight_Rt = FALSE) {

  if (sum(p$hosp_transmission, p$ICU_transmission, p$G_D_transmission) > 0) {
    stop("Cannot currently compute Rt if any of 'hosp_transmission',
    'ICU_transmission' or 'G_D_transmission' are non-zero")
  }

  all_types <- c("eff_Rt_all", "eff_Rt_general", "Rt_all", "Rt_general")
  if (is.null(type)) {
    type <- all_types
  } else {
    err <- setdiff(type, all_types)
    if (length(err) > 0) {
      stop(sprintf("Unknown R type %s, must match %s",
                   paste(squote(err), collapse = ", "),
                   paste(squote(all_types), collapse = ", ")))
    }
  }

  n_strains <- 1
  n_real_strains <- 1
  R <- NULL

  if (nrow(S) != nlayer(p$rel_susceptibility) * nrow(p$m)) {
    stop(sprintf(
      "Expected 'S' to have %d rows = %d groups x %d vaccine classes",
      p$n_groups * nlayer(p$rel_susceptibility),
      p$n_groups,
      nlayer(p$rel_susceptibility)))
  }
  if (ncol(S) != length(time)) {
    stop(sprintf("Expected 'S' to have %d columns, following 'time'",
                 length(time)))
  }
  if (!is.null(R)) {
    if (nrow(R) != nlayer(p$rel_susceptibility) * nrow(p$m) * p$n_strains_R) {
      stop(sprintf(
        "Expected 'R' to have %d rows = %d groups x %d strains_R x %d vaccine
          classes",
        p$n_groups * nlayer(p$rel_susceptibility) * p$n_strains_R,
        p$n_groups, p$n_strains_R, nlayer(p$rel_susceptibility)))
    }
    if (ncol(R) != length(time)) {
      stop(sprintf("Expected 'R' to have %d columns, following 'time'",
                   length(time)))
    }
  }

  n_vacc_classes <- nlayer(p$rel_susceptibility)

  ### here mean_duration accounts for relative infectivity of
  ### different infection / vaccination stages
  beta <- ZamCovid_parameters_expand_step(time, p$beta_step)

  ages <- seq_len(p$n_groups)
  ## We only need to do this section for each different beta value,
  ## and then it's still a pretty straightforward scaling; we've
  ## applied beta to all *but* a fraction of the matrix (i.e., in the matrix
  ##
  ##   A B
  ##   C D
  ##
  ## Everything but D is scaled by beta
  ##
  ## When we have more than one vaccination group we replicate this to
  ## get block structure
  ##
  ##   A B A B
  ##   C D C D
  ##   A B A B
  ##   C D C D
  ##
  ## And scale all of A, B, C (not D) by beta
  m <- block_expand(unname(p$m), n_vacc_classes)
  mt <- m %o% beta

  n_time <- length(time)

  mean_duration <-
    ZamCovid_Rt_mean_duration_weighted_by_infectivity(time, p)

  compute_ngm <- function(x, S, rel_sus, R = 0, rel_sus_strain = 1) {
    n_groups <- p$n_groups
    len <- n_groups * n_vacc_classes
    Sw <- (S + R * rel_sus_strain) * c(rel_sus)
    mt * vapply(seq_len(n_time), function(t)
      tcrossprod(c(x[, , t]), Sw[, t]),
      matrix(0, len, len))
  }

  ## NOTE the signs on the exponents here is different! This gives
  ## good performance and reasonable accuracy to the point where
  ## this calculation is small in the profile.
  eigen <- function(m) {
    eigen1::eigen1(m, max_iterations = 1e5, tolerance = 1e-6,
                   method = eigen_method)
  }

  ret <- list(time = time,
              date = time * p$dt,
              beta = beta)

  n_groups <- nrow(p$m)

  ngm_computer <- function(x, effective) {
    compute_ngm(mean_duration, x, p$rel_susceptibility)
  }

  if ("eff_Rt_general" %in% type) {
    ngm <- ngm_computer(S, effective = TRUE)
    ret$eff_Rt_general <- array(eigen(ngm), dim = c(n_time, n_vacc_classes))
  }

  if ("Rt_general" %in% type) {
    N_tot_non_vacc <- array(p$N_tot, dim = c(p$n_groups, ncol(S)))
    N_tot_all_vacc_groups <- N_tot_non_vacc
    if (n_vacc_classes > 1) {
      for (i in 2:n_vacc_classes) {
        N_tot_all_vacc_groups <- rbind(N_tot_all_vacc_groups,
                                       0 * N_tot_non_vacc)
      }
    }

    ngm <- ngm_computer(N_tot_all_vacc_groups, effective = FALSE)
    ret$Rt_general <- array(eigen(ngm), dim = c(n_time, n_vacc_classes))
  }

  ## ensure backwards compatibility by dropping columns for single_strain and
  ## separating classes
  class(ret) <- c("Rt")
  is_single <- (inherits(ret$Rt_general, c("matrix", "array")) &&
                  is.null(ncol(ret$Rt_general))) ||
    (!inherits(ret$Rt_general, c("matrix", "array")) &&
       length(ret$Rt_general) == 1)


  if (is_single) {
    ## adding 'nocov' as this is a safety check that should never be hit
    class(ret) <- c("single_strain", class(ret)) # nocov
  } else if (isTRUE(ncol(ret$Rt_general) == 1)) {
    ret[intersect(all_types, names(ret))] <-
      lapply(ret[intersect(all_types, names(ret))], drop)
    class(ret) <- c("single_strain", class(ret))
  }

  ret
}


## Here we expect 'S' in order:
##
##   state x sample x time
##
## We expect 'pars' to be a list along sample (or a shared parameter set)
## We expect 'time' to be a vector along time

##' Compute "Rt" for a set of simulated trajectories (e.g., the result
##' of the `$iterate()` method of [ZamCovid], [mcstate::pmcmc()] or
##' [mcstate::pmcmc_predict()]. The trajectories may or may not share
##' parameters.
##'
##' @title Compute Rt for a set of trajectories
##'
##' @param time A vector of time steps
##'
##' @param S A 3d ((n groups x n vaccine classes) x n trajectories x n time
##'   steps) array of "S" compartment counts
##'
##' @param pars Either a single [ZamCovid_parameters()] object
##'   (shared parameters) or an unnamed list of
##'   [ZamCovid_parameters()] objects, the same length as `ncol(S)`.
##'
##' @param initial_time_from_parameters If `TRUE`, then `time[[1]]` is
##'   replaced by the value of `initial_time` from the parameters.
##'   This is usually what you want.
##'
##' @param R A 2d ((n groups x n strains x n vaccine classes) x
##'   n trajectories x n time steps) array of "R" compartment counts, required
##'   for multi-strain models.
##'
##' @inheritParams ZamCovid_Rt
##'
##' @return As for [ZamCovid_Rt()], except that every element is a
##'   matrix, not a vector.
##'
##' @export
ZamCovid_Rt_trajectories <- function(time, S, pars,
                                     initial_time_from_parameters = TRUE,
                                     type = c("eff_Rt_general", "Rt_general"),
                                     eigen_method = "power_iteration",
                                     R = NULL, weight_Rt = FALSE) {
  calculate_Rt_trajectories(
    calculate_Rt = ZamCovid_Rt, time = time,
    S = S, pars = pars,
    initial_time_from_parameters = initial_time_from_parameters,
    type = type,
    eigen_method = eigen_method,
    R = R,
    weight_Rt = weight_Rt)
}


ZamCovid_Rt_mean_duration_weighted_by_infectivity <- function(time, pars) {

  dt <- pars$dt
  n_time_steps <-
    length(ZamCovid_parameters_expand_step(time, pars$p_H_step))

  p_C <- combine_steps_groups(
    time, pars$n_groups, n_time_steps,
    n_strains = 1, n_vacc_classes = pars$n_vacc_classes,
    p_step = pars$p_C_step, rel_p = pars$rel_p_sympt,
    strain_rel_p = pars$strain_rel_p_sympt
  )

  ## Compute mean duration (in time steps) of each stage of infection,
  ## weighed by probability of going through that stage
  ## and by relative infectivity of that stage

  ## Note the mean duration (in time steps) of a compartment for
  ## a discretised Erlang(k, gamma) is k / (1 - exp(dt * gamma))
  calculate_mean <- function(transmission, prob, name) {
    gamma_step <-
      ZamCovid_parameters_expand_step(time,
                                      pars[[paste0("gamma_", name, "_step")]])
    rel_gamma <- 1
    k <- 1
    gamma <- aperm(outer(outer(gamma_step, rel_gamma),
                         array(1, pars$n_groups)),
                   c(3, 2, 1))
    transmission * k * prob / stats::pexp(gamma, dt)
  }

  mean_duration_I_A <- calculate_mean(pars$I_A_transmission, (1 - p_C), "A")
  mean_duration_I_P <- calculate_mean(pars$I_P_transmission, p_C, "P")
  mean_duration_I_C_1 <- calculate_mean(pars$I_C_1_transmission, p_C, "C_1")
  mean_duration_I_C_2 <- calculate_mean(pars$I_C_2_transmission, p_C, "C_2")

  mean_duration <- mean_duration_I_A + mean_duration_I_P +
    mean_duration_I_C_1 + mean_duration_I_C_2

  ## Account for different infectivity levels depending on vaccination stage
  mean_duration <- mean_duration * outer(pars$rel_infectivity,
                                         rep(1, n_time_steps))

  ## Multiply by dt to convert from time steps to days
  mean_duration <- dt * mean_duration

  mean_duration
}


calculate_Rt_trajectories <- function(calculate_Rt, time, S, pars,
                                      initial_time_from_parameters,
                                      type, R = NULL, ...) {
  if (length(dim(S)) != 3) {
    stop("Expected a 3d array of 'S'")
  }

  if (!is.null(R) && length(dim(R)) != 3) {
    stop("Expected a 3d array of 'R'")
  }


  if (!is.null(names(pars))) {
    stop("If not using shared parameters, expected a unnamed list for 'pars'")
  }
  if (length(pars) != ncol(S)) {
    stop(sprintf(
      "Expected 2nd dimension of 'S' to have length %d, following 'pars'",
      length(pars)))
  }


  if (dim(S)[[3]] != length(time)) {
    stop(sprintf(
      "Expected 3rd dimension of 'S' to have length %d, following 'time'",
      length(time)))
  }

  if (!is.null(R) && any(dim(R)[2:3] != dim(S)[2:3])) {
    stop(sprintf(
      "Expected 2nd and 3rd dimension of 'R' to be the same as 'S' (%d x %d)'",
      ncol(S), dim(S)[[3]]))
  }

  calculate_rt_one_trajectory <- function(i) {
    if (initial_time_from_parameters) {
      time[[1L]] <- pars[[i]]$initial_time %||% 0
    }

    rt_1 <- calculate_Rt(time, S[, i, ], pars[[i]], type = type,
                         R = R[, i, ], ...)
    rt_1
  }

  res <- lapply(seq_along(pars), calculate_rt_one_trajectory)

  ## These are stored in a list-of-lists and we convert to a
  ## list-of-matrices here
  collect <- function(nm) {
    if (nm %in% c("time", "date", "beta")) {
      matrix(unlist(lapply(res, "[[", nm)), length(time), length(res))
    } else {
      array(unlist(lapply(res, "[[", nm)),
            c(length(time), ncol(res[[1]][[nm]]), length(res)))
    }
  }
  nms <- names(res[[1]])
  ret <- set_names(lapply(nms, collect), nms)

  # Ensure backwards compatibility by dropping columns for single_strain and
  # separating classes
  all_types <- c("eff_Rt_general", "Rt_general")
  ret[intersect(all_types, names(ret))] <-
    lapply(ret[intersect(all_types, names(ret))], drop)
  class(ret) <- c("single_strain", "Rt_trajectories", "Rt")

  ret
}


nlayer <- function(x) {
  dim(x)[2L]
}

block_expand <- function(m, n) {
  if (n == 1L) {
    return(m)
  }
  len <- nrow(m) * n
  matrix(t(matrix(m, nrow(m), len)), len, len, byrow = TRUE)
}

combine_steps_groups <- function(step, n_groups, n_time_steps,
                                 n_strains = 1, n_vacc_classes,
                                 p_step, rel_p, strain_rel_p = 1) {

  ret <- vapply(
    seq_len(n_groups),
    function(i) {
      outer(
        ZamCovid_parameters_expand_step(step, p_step[, i]),
        rel_p[i, ]
      )
    },
    array(0, c(n_time_steps, n_vacc_classes))
  )

  ## If these dimensions are all 1 then the above reduces to a vector
  ## rather than an array, so here we just reshape into an array
  if (all(c(n_time_steps, n_strains, n_vacc_classes) == 1)) {
    ret <- array(ret, c(1, 1, 1, length(ret)))
  }

  ret <- pmin(ret, 1)
  ret <- aperm(ret, c(3, 2, 1))

  ret
}

set_names <- function(x, nms) {
  names(x) <- nms
  x
}
