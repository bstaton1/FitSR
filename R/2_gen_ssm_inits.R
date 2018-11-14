#' Generate initial values for a state-space SRA fit
#'
#' Takes the output of the \code{SimFit} functions to prepare
#' for fitting a state-space spawner recruit model using JAGS
#'
#' @param params a list object created by \code{SimFit::init_sim()} or \code{kusko_data_prep()$params}.
#'   Contains the driving parameters and dimensional variables.
#' @param obs a list object created by the functions that generate observations in the
#'   \code{SimFit} package or the \code{kusko_data_prep()$obs} function.
#' @param n_chains a numeric vector of length 1: the number of MCMC chains
#'   to generate initial values for
#'
#' @export

gen_ssm_inits = function(params, obs, n_chains) {

  output = with(append(params, obs), {
    # fit a basic regression approach to get Umsy and Smsy for each substock
    lm_data = lm_data_prep(params = params, obs = obs)
    lm_alpha = NULL
    lm_beta = NULL
    for (s in 1:ns) {
      tmp_S = lm_data$S_obs[lm_data$stock == s]
      tmp_log_RPS = lm_data$obs_log_RPS_lm[lm_data$stock == s]

      tmp_fit = tryCatch({
        lm(tmp_log_RPS ~ tmp_S)
      }, error = function(e) NULL)

      if (is.null(tmp_fit)) {
        lm_alpha = c(lm_alpha, NA)
        lm_beta = c(lm_beta, NA)
      } else {
        lm_alpha = c(lm_alpha, unname(exp(coef(tmp_fit)[1])))
        lm_beta = c(lm_beta, unname(abs(coef(tmp_fit)[2])))
      }
    }

    # fix very small alpha values
    lm_alpha[lm_alpha <= 1] = 1.5

    # replace NAs with the mean of the others
    lm_alpha[lm_alpha > 20] = 20 # cap it at 20
    lm_alpha[is.na(lm_alpha)] = mean(lm_alpha, na.rm = T)
    lm_beta[is.na(lm_beta)] = mean(lm_beta, na.rm = T)

    # get the reference points of U_msy and S_msy
    lm_mgmt = SimSR::gen_lm_mgmt(lm_alpha, lm_beta)

    # randomly perturb these n_chains times and store
    inits = list()
    for (i in 1:n_chains) {
      inits[[i]] = list(
        U_msy = {
          y = rbeta(ns, lm_mgmt$U_msy * 100,(1 - lm_mgmt$U_msy) * 100)
          ifelse (y < 0.1, 0.15, y)
        },

        log_S_msy = log(rlnorm(ns, log(lm_mgmt$S_msy), 0.1)),
        log_R = apply(R_ys_obs, 2, function(x) { # loop over stocks
          mu = mean(x, na.rm = T)   # calculate mean when available
          x[is.na(x)] = mu    # fill in NA with the mean
          log(rlnorm(length(x), log(x), 0.2))  # perturb it
        }),
        D_scale = runif(1, 0.08, 0.12)
      )
    }
    inits
  })

  return(output)
}

