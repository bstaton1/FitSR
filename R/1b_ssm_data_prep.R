#' Prepare data for a state-space SRA fit
#'
#' Takes the output of the \code{SimFit} functions to prepare
#' for fitting a state-space spawner recruit model
#'
#' @param params a list object created by \code{SimFit::init_sim()} or \code{kusko_data_prep()$params}.
#'   Contains the driving parameters and dimensional variables.
#' @param obs a list object created by the functions that generate observations in the
#'   \code{SimFit} package or the \code{kusko_data_prep()$obs} function.
#' @param covariance a character vector of length 1: \code{"simple"} or \code{"complex"}
#'
#' @export

ssm_data_prep = function(params, obs, covariance) {

  output = with(append(params, obs), {

    # vectorize escapement observations
    # don't want to waste time looping over NAs
    S_ts_obs_m = S_ts_obs
    S_ts_obs_v = as.numeric(S_ts_obs_m)

    S_obs_s = rep(1:ns, each = nt)
    S_obs_t = rep(1:nt, ns)

    no_na_yrs = which(!is.na(S_ts_obs_v))

    S_obs = S_ts_obs_v[no_na_yrs]
    S_obs_s = S_obs_s[no_na_yrs]
    S_obs_t = S_obs_t[no_na_yrs]
    S_obs_n = length(S_obs)

    sig_S_obs_ts_v = as.numeric(sig_S_ts_obs)
    sig_S_obs = sig_S_obs_ts_v[no_na_yrs]

    # remove NAs from age comps: turn them to zeros
    x_tas_obs[is.na(x_tas_obs)] = 0

    # do the covariance matrix info
    if (covariance == "simple") {  # this is if common variance, common rho
      m = matrix(1:(ns^2), ns, ns, byrow = T)
      d = diag(m)

      vcov_list = list(
        vcov_row = rep(1:ns, each = ns),
        vcov_col = rep(1:ns, ns),
        vcov_ind = ifelse(1:ns^2 %in% d, 1, 2),
        vcov_N = ns * ns
      )
    }

    if (covariance == "complex") { # this is if full Sigma_R mat estimated as inverse wishart
      vcov_list = list(
        R_wish = diag(rep(1,ns)),
        df_wish = ns + 1
      )
    }

    # the base data needed regardless of covariance structure
    base = list(
      # dimension variables
      ns = ns,
      nt = nt,
      ny = ny,
      na = na,
      a_max = a_max,

      # observed harvest states
      C_tot_t_obs = C_tot_t_obs,
      tau_C_obs = 1/sig_C_t_obs^2,
      v = v,

      # vectorized observe escapement counts
      S_obs = S_obs, # the count
      S_obs_t = S_obs_t, # the year of the ith count
      S_obs_s = S_obs_s, # the stock of the ith count
      S_obs_n = S_obs_n, # the number of escapement observations
      tau_S_obs = 1/sig_S_obs^2,

      # observed age comp states
      x_tas_obs = x_tas_obs,
      ESS_ts = apply(x_tas_obs, 3, rowSums),
      age_stocks = age_comp_stocks,
      n_age_stocks = length(age_comp_stocks)
    )

    append(base, vcov_list)
  })

  return(output)
}


