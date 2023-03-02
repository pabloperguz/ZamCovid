#' Plot basci ZamCovid pMCMC traces.
#'
#' @param samples A `list` object containing the sample outputs of a basic
#'    model run.
#'
#' @return A base R plot of the pMCMC chains for fitted parameters and
#'    log-likelihood.
#'
#' @export
#'
#' @examples
zamcovid_plot_traces <- function(samples) {

  if (is.null(samples$chain)) {
    n_chains <- 1L
  } else {
    n_chains <- length(unique(samples$chain))
  }

  cols <- rev(viridisLite::viridis(n_chains))

  pars <- samples$pars

  nms <- colnames(pars)
  probs <- samples$probabilities

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  new_grid <- function(n, title) {
    par(mfrow = rep(ceiling(sqrt(n + 1)), 2),
        mar = c(3, 3, 2, 1),
        mgp = c(2, 0.5, 0),
        oma = c(1, 1, 1 + as.integer(title), 1))
  }

  plot_traces1 <- function(p, name) {
    traces <- matrix(p, ncol = n_chains)
    ess <- coda::effectiveSize(coda::as.mcmc(traces))

    if (name == "log_likelihood") {
      main <- ""
    } else {
      main <- sprintf("ess = %s", round(sum(ess)))
    }
    matplot(traces, type = "l", lty = 1,
            xlab = "Iteration", bty = "n",
            ylab = name, col = cols,
            main = main,
            font.main = 1)
    rug(samples$iteration[samples$chain == 1], ticksize = 0.1)
  }

  new_grid(length(nms), FALSE)
  for (nm in nms) {
    plot_traces1(samples$pars[, nm], nm)
  }
  plot_traces1(probs[, "log_likelihood"], "log_likelihood")
}


#' Plot basic ZamCovid model state outputs.
#'
#' @param samples A `list` of model state outputs.
#'
#' @param data A `data.frame` of the data used for model fitting.
#'
#' @return A plot of fitted model trajectories.
#'
#' @export
#'
#' @examples
zamcovid_plot_trajectories <- function(samples, data) {

  states <- samples$trajectories$state
  nms <- rownames(states)
  time <- data$date

  parse_traj <- function(n, tmp) {
    data.frame(
      time = time,
      state = n,
      mean = colMeans(tmp),
      lb = matrixStats::colQuantiles(tmp, probs = 0.025),
      ub = matrixStats::colQuantiles(tmp, probs = 0.975),
      data = data[, n]
    )
  }

  df <- NULL
  for (n in nms) {
    ret <- parse_traj(n, states[n, , -1])
    df <- rbind(df, ret)
  }

  df$state <- factor(df$state, levels = unique(df$state))

  ggplot(df, aes(x = time)) +
    geom_line(aes(y = mean, col = state)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = state), alpha = 0.4) +
    geom_point(aes(y = data, col = state), size = 0.7, alpha = 0.9) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~state) +
    theme_bw() +
    theme(legend.position = "none")
}
