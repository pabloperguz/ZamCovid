##' @importFrom stats rexp
test_prob_pos <- function(pos, neg, sensitivity, specificity, exp_noise) {

  ## We add some exponential noise to the number of positives and negatives
  ## to help ensure prob_pos is not 0 or 1. If e.g. prob_pos were 0 and there
  ## were individuals who tested positive, this would result in a weight of 0
  ## for a particle. The exponential noise produces small non-zero weights in
  ## these circumstances to prevent the particle filter from breaking.

  pos <- pos + rexp(length(pos), exp_noise)
  neg <- neg + rexp(length(neg), exp_noise)

  prob_pos <- (sensitivity * pos + (1 - specificity) * neg) / (pos + neg)
  prob_pos
}


##' @importFrom stats dnbinom rexp
ll_nbinom <- function(data, model, kappa, exp_noise) {
  if (is.na(data)) {
    return(numeric(length(model)))
  }
  mu <- model + rexp(length(model), rate = exp_noise)
  dnbinom(data, kappa, mu = mu, log = TRUE)
}


##' @importFrom stats dbinom
ll_binom <- function(data_x, data_size, model_prob) {
  if (is.na(data_x) || is.na(data_size)) {
    return(numeric(length(model_prob)))
  }
  dbinom(data_x, data_size, model_prob, log = TRUE)
}


##' @importFrom stats dpois
ll_dpois <- function(data, model_lambda) {
  if (is.na(data)) {
    return(numeric(length(model_lambda)))
  }
  dpois(data, model_lambda, log = TRUE)
}
