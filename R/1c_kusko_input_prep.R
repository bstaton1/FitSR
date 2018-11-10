#' Prepare the Kuskokwim Substock Data for this analysis
#'
#' Prepare the output data files from the GitHub repository: bstaton1/kusko-multi-sra
#' into the format needed by the functions in the \code{FitSR} package.
#'
#' @param S_dat the output file for substock spawner abundance
#' @param H_dat the output file for drainage-wide
#'   aggregate harvest summed across all substocks
#' @param age_dat the output file for substock-specific age composition data
#'   for substocks that have these data
#'
#' @export

# arguments are the raw output files in the data-prep directories
kusko_prep = function(S_dat, H_dat, age_dat) {

  ### HANDLE THE ESCAPEMENT DATA ###
  # get total escapement each year from stocks in this analysis: used in getting U
  S_tot_t_obs = tapply(S_dat$mean, S_dat$year, sum)

  # drop out years not used in this analysis
  S_dat$mean[S_dat$obs == 0] = NA

  # calendar year dims
  years = unique(S_dat$year)
  nt = length(years)

  # stock dims
  stocks = unique(S_dat$stock)
  ns = length(stocks)

  # prepare escapement data
  S_ts_obs = round(reshape2::dcast(S_dat, year ~ stock, value.var = "mean")); rownames(S_ts_obs) = years; S_ts_obs = S_ts_obs[,-1]
  cv_S_ts_obs = reshape2::dcast(S_dat, year ~ stock, value.var = "cv"); rownames(cv_S_ts_obs) = years; cv_S_ts_obs = cv_S_ts_obs[,-1]
  sig_S_ts_obs = sqrt(log(cv_S_ts_obs^2+1))

  ### HANDLE THE HARVEST DATA ###
  C_tot_t_obs = H_dat$mean
  cv_C_tot_t_obs = H_dat$cv
  sig_C_tot_t_obs = sqrt(log(cv_C_tot_t_obs^2+1))

  U_t_obs = C_tot_t_obs/(S_tot_t_obs + C_tot_t_obs)

  ### SET VULNERABILITY ###
  region = c("mid", "mid", "mid", "mid", "lwr", "mid", "lwr", "mid", "upr", "mid", "upr", "upr", "lwr")
  v = ifelse(region == "lwr" | region == "mid", 0.95, 1)

  # age dimensions
  a_min = 4
  a_max = 7
  na = a_max - a_min + 1
  ages = a_min:a_max
  ny = nt + na - 1

  ### HANDLE THE AGE COMP DATA ###

  age_stocks = which(stocks %in% age_dat$stock)
  n_age_stocks = length(age_stocks)

  # extract the counts for each aged stock and place in the right spot
  x_tas_obs = array(NA, c(nt, na, n_age_stocks))
  for (j in 1:n_age_stocks) {
    x_tas_obs[,,j] = as.matrix(age_dat[age_dat$stock == stocks[age_stocks[j]],paste("a", a_min:a_max, sep = "")])
  }
  dimnames(x_tas_obs) = list(years, paste("a", a_min:a_max, sep = ""), stocks[age_stocks])
  ESS_ts = apply(x_tas_obs, 3, rowSums)

  obs = list(
    C_tot_t_obs = C_tot_t_obs,
    S_ts_obs = as.matrix(S_ts_obs),
    x_tas_obs = x_tas_obs,
    sig_S_ts_obs = as.matrix(sig_S_ts_obs),
    sig_C_t_obs = sig_C_tot_t_obs,
    U_t_obs = U_t_obs,
    age_comp_stocks = age_stocks
  )

  params = list(
    ns = ns,
    nt = nt,
    a_min = a_min,
    a_max = a_max,
    na = na,
    ages = a_min:a_max,
    ny = ny,
    stocks = stocks,
    v = v
  )

  list(
    obs = obs,
    params = params
  )
}
