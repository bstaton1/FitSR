#' Write the lm and lme model code
#'
#' @param dir_name the directory name in your working directory to place the file.
#'   If the directory does not exist, it will be created with a warning
#' @param file_name the file name to save the file as.
#' @param silent logical. Do you wish to suppress the printing of the location of the file?
#'
#' @export

write_lme_file = function(dir_name = "Model Files", file_name = "lme_model.txt", silent = F) {

  # error handle for directory name
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
    warning(paste("No directory named '", dir_name, "' found. One was created.", sep = ""))
  }

  # specify the JAGS model
  mod = function() {

    ### FIT THE MIXED EFFECTS VERSION ###
    sig_fit_lme ~ dunif(0, 5)
    tau_fit_lme <- 1/sig_fit_lme^2

    log_alpha_bar_lme ~ dunif(0,10)
    sig_log_alpha_lme ~ dunif(0,10)
    tau_log_alpha_lme <- 1/sig_log_alpha_lme^2

    # stock level parameters
    for (s in 1:ns) {
      log_alpha_lme[s] ~ dnorm(log_alpha_bar_lme, tau_log_alpha_lme)
      beta_lme[s] ~ dunif(0, 1)
      alpha_lme[s] <- exp(log_alpha_lme[s])
    }

    # likelihood
    for (i in 1:nobs) {
      obs_log_RPS_lme[i] ~ dnorm(pred_log_RPS_lme[i], tau_fit_lme)
      pred_log_RPS_lme[i] <- log_alpha_lme[stock[i]] - beta_lme[stock[i]] * S_obs[i]
    }

    ### FIT THE INDEPENDENT REGRESSIONS VERSION ###
    # stock level parameters
    for (s in 1:ns) {
      log_alpha_lm[s] ~ dunif(0, 5)
      alpha_lm[s] <- exp(log_alpha_lm[s])
      beta_lm[s] ~ dunif(0, 1)
      sig_fit_lm[s] ~ dunif(0, 5)
      tau_fit_lm[s] <- 1/sig_fit_lm[s]^2
    }

    # likelihood
    for (i in 1:nobs) {
      obs_log_RPS_lm[i] ~ dnorm(pred_log_RPS_lm[i], tau_fit_lm[stock[i]])
      pred_log_RPS_lm[i] <- log_alpha_lm[stock[i]] - beta_lm[stock[i]] * S_obs[i]
    }
  }

  R2OpenBUGS::write.model(model = mod, file.path(dir_name, file_name))
  if(!silent) return(file.path(dir_name, file_name))
}

#' Write the ssm 1 model code
#'
#' @param dir_name the directory name in your working directory to place the file.
#'   If the directory does not exist, it will be created with a warning
#' @param file_name the file name to save the file as.
#' @param silent logical. Do you wish to suppress the printing of the location of the file?
#'
#' @export

write_ssm_1_file = function(dir_name = "Model Files", file_name = "ssm_1_model.txt", silent = F) {

  # error handle for directory name
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
    warning(paste("No directory named '", dir_name, "' found. One was created.", sep = ""))
  }

  # specify JAGS model code
  mod = function() {
    # Priors on primary parameters
    # mean_sigma_R <- 0.71
    # mean_tau_R <- 1/mean_sigma_R^2
    # tau_R_red <- mean_tau_R * (1 - phi^2)
    # log_resid_0a ~ dnorm(0, tau_R_red)
    log_resid_0a <- 0

    phi ~ dunif(-0.99, 0.99)
    for (s in 1:ns) {
      U_msy[s] ~ dunif(0.01, 0.99)
      log_S_msy[s] ~ dnorm(0, 0.001) %_% I(1, 11.5)
      S_msy[s] <- exp(log_S_msy[s])

      alpha[s] <- exp(U_msy[s])/(1 - U_msy[s])
      log_alpha[s] <- log(alpha[s])
      beta[s] <- U_msy[s]/S_msy[s]
      log_resid_0[s] <- log_resid_0a
    }

    # build the Sigma_R[,] matrix
    sig.common ~ dunif(0,2)
    var.common <- sig.common^2
    rho.common ~ dunif(-0.05,1)
    rho.vec[1] <- 1
    rho.vec[2] <- rho.common

    for (i in 1:vcov_N) {
      rho_mat[vcov_row[i],vcov_col[i]] <- rho.vec[vcov_ind[i]]
      Sigma_R[vcov_row[i],vcov_col[i]] <- var.common * rho.vec[vcov_ind[i]]
    }

    Tau_R[1:ns,1:ns] <- inverse(Sigma_R[1:ns,1:ns])

    # white noise process sd for each substock
    for (s in 1:ns) {
      sigma_R[s] <- sqrt(Sigma_R[s,s])
    }

    # produce Ricker predictions
    for (s in 1:ns) {
      # for years without SR link: use unfished equilibrium recruitment
      R_eq[s] <- log_alpha[s]/beta[s]
      R0[s] <- R_eq[s]
      log_R0[s] <- log(R0[s])

      # log_R_mean1 = deterministic ricker; log_R_mean2 = time-corrected expectation
      log_R_mean1[1,s] <- log_R0[s]
      R_mean1[1,s] <- R0[s]
      log_R_mean2[1,s] <- log_R_mean1[1,s] + phi * log_resid_0[s]
      for (y in 2:a_max) {
        R_mean1[y,s] <- R0[s]
        log_R_mean1[y,s] <- log_R0[s]
        log_R_mean2[y,s] <- log_R_mean1[y,s] + phi * log_resid[y-1,s]
      }

      # for years with SR link
      for (y in (a_max+1):ny) {
        R_mean1[y,s] <- S[y-a_max,s] * exp(log_alpha[s] - beta[s] * S[y-a_max,s])
        log_R_mean1[y,s] <- log(R_mean1[y,s])
        log_R_mean2[y,s] <- log_R_mean1[y,s] + phi * log_resid[y-1,s]
      }
    }

    # draw true recruitment states
    for (y in 1:ny) {
      log_R[y,1:ns] ~ dmnorm(log_R_mean2[y,1:ns], Tau_R[1:ns,1:ns])
    }

    # calculate residuals
    for (y in 1:ny) {
      for (s in 1:ns) {
        R[y,s] <- exp(log_R[y,s])
        log_resid[y,s] <- log_R[y,s] - log_R_mean1[y,s]
      }
      # R_tot[y] <- sum(R[y,1:ns])
    }

    # maturity schedule
    prob[1] ~ dbeta(1,1)
    prob[2] ~ dbeta(1,1)
    prob[3] ~ dbeta(1,1)
    pi[1] <- prob[1]
    pi[2] <- prob[2] * (1 - pi[1])
    pi[3] <- prob[3] * (1 - pi[1] - pi[2])
    pi[4] <- 1 - pi[1] - pi[2] - pi[3]

    for (y in 1:ny) {
      p[y,1:na] <- pi[1:na]
    }

    # allocate R[y,s] to N[t,a,s] and create predicted calendar year states for each stock
    for (s in 1:ns) {
      for (t in 1:nt) {
        for (a in 1:na) {
          N_tas[t,a,s] <- R[t+na-a,s] * p[t+na-a,a]
        }
        N[t,s] <- sum(N_tas[t,1:na,s])
        S[t,s] <- N[t,s] * (1 - U[t] * v[s])
        C[t,s] <- N[t,s] * (U[t] * v[s])
      }
    }

    # create calendar year totals across stocks
    for (t in 1:nt) {
      U[t] ~ dbeta(1,1)

      # N_tot[t] <- sum(N[t,1:ns])
      # S_tot[t] <- sum(S[t,1:ns])
      C_tot[t] <- sum(C[t,1:ns])
    }

    # obtain calendar year age composition for each stock
    # use i because not looping over all stocks - only those that have have comps
    for (i in 1:n_age_stocks) {
      for (t in 1:nt) {
        for (a in 1:na) {
          q[t,a,i] <- N_tas[t,a,age_stocks[i]]/N[t,age_stocks[i]]
        }
      }
    }

    # observe calendar year total harvest
    for (t in 1:nt) {
      log_C_tot[t] <- log(C_tot[t])
      C_tot_t_obs[t] ~ dlnorm(log_C_tot[t], tau_C_obs[t])
    }

    # observe calendar year substock specific harvests
    # vectorized to avoid looping over many NAs (hence the i not s)
    for (i in 1:S_obs_n) {
      log_S[i] <- log(S[S_obs_t[i],S_obs_s[i]])
      S_obs[i] ~ dlnorm(log_S[i], tau_S_obs[i])
    }

    # observe age composition
    # only loop over stocks that have age comp data (hence the i not s)
    for (i in 1:n_age_stocks) {
      for (t in 1:nt) {
        x_tas_obs[t,1:na,i] ~ dmulti(q[t,1:na,i], ESS_ts[t,i])
      }
    }
  }

  R2OpenBUGS::write.model(model = mod, file.path(dir_name, file_name))

  if(!silent) return(file.path(dir_name, file_name))
}

#' Write the ssm 2 model code
#'
#' @param dir_name the directory name in your working directory to place the file.
#'   If the directory does not exist, it will be created with a warning
#' @param file_name the file name to save the file as.
#' @param silent logical. Do you wish to suppress the printing of the location of the file?
#'
#' @export

write_ssm_2_file = function(dir_name = "Model Files", file_name = "ssm_2_model.txt", silent = F) {

  # error handle for directory name
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
    warning(paste("No directory named '", dir_name, "' found. One was created.", sep = ""))
  }

  # specify JAGS model code
  mod = function() {
    # Priors on primary parameters
    # mean_sigma_R <- 0.71
    # mean_tau_R <- 1/mean_sigma_R^2
    # tau_R_red <- mean_tau_R * (1 - phi^2)
    # log_resid_0a ~ dnorm(0, tau_R_red)
    log_resid_0a <- 0

    phi ~ dunif(-0.99, 0.99)
    for (s in 1:ns) {
      U_msy[s] ~ dunif(0.01, 0.99)
      log_S_msy[s] ~ dnorm(0, 0.001) %_% I(1, 11.5)
      S_msy[s] <- exp(log_S_msy[s])

      alpha[s] <- exp(U_msy[s])/(1 - U_msy[s])
      log_alpha[s] <- log(alpha[s])
      beta[s] <- U_msy[s]/S_msy[s]
      log_resid_0[s] <- log_resid_0a
    }

    # estimate the covariance matrix on log(recruitment states
    Tau_R[1:ns,1:ns] ~ dwish(R_wish[1:ns,1:ns], df_wish)
    Sigma_R[1:ns,1:ns] <- inverse(Tau_R)

    # white noise process sd for each substock
    for (s in 1:ns) {
      sigma_R[s] <- sqrt(Sigma_R[s,s])
    }

    # get the pairwise correlation matrix
    for (i in 1:ns) {
      for (j in 1:ns) {
        rho_mat[i,j] <- Sigma_R[i,j]/(sigma_R[i] * sigma_R[j])
      }
    }

    # produce Ricker predictions
    for (s in 1:ns) {
      # for years without SR link: use unfished equilibrium recruitment
      R_eq[s] <- log_alpha[s]/beta[s]
      R0[s] <- R_eq[s]
      log_R0[s] <- log(R0[s])

      # log_R_mean1 = deterministic ricker; log_R_mean2 = time-corrected expectation
      log_R_mean1[1,s] <- log_R0[s]
      R_mean1[1,s] <- R0[s]
      log_R_mean2[1,s] <- log_R_mean1[1,s] + phi * log_resid_0[s]
      for (y in 2:a_max) {
        R_mean1[y,s] <- R0[s]
        log_R_mean1[y,s] <- log_R0[s]
        log_R_mean2[y,s] <- log_R_mean1[y,s] + phi * log_resid[y-1,s]
      }

      # for years with SR link
      for (y in (a_max+1):ny) {
        R_mean1[y,s] <- S[y-a_max,s] * exp(log_alpha[s] - beta[s] * S[y-a_max,s])
        log_R_mean1[y,s] <- log(R_mean1[y,s])
        log_R_mean2[y,s] <- log_R_mean1[y,s] + phi * log_resid[y-1,s]
      }
    }

    # draw true recruitment states
    for (y in 1:ny) {
      log_R[y,1:ns] ~ dmnorm(log_R_mean2[y,1:ns], Tau_R[1:ns,1:ns])
    }

    # calculate residuals
    for (y in 1:ny) {
      for (s in 1:ns) {
        R[y,s] <- exp(log_R[y,s])
        log_resid[y,s] <- log_R[y,s] - log_R_mean1[y,s]
      }
      # R_tot[y] <- sum(R[y,1:ns])
    }

    # maturity schedule
    prob[1] ~ dbeta(1,1)
    prob[2] ~ dbeta(1,1)
    prob[3] ~ dbeta(1,1)
    pi[1] <- prob[1]
    pi[2] <- prob[2] * (1 - pi[1])
    pi[3] <- prob[3] * (1 - pi[1] - pi[2])
    pi[4] <- 1 - pi[1] - pi[2] - pi[3]

    for (y in 1:ny) {
      p[y,1:na] <- pi[1:na]
    }

    # allocate R[y,s] to N[t,a,s] and create predicted calendar year states for each stock
    for (s in 1:ns) {
      for (t in 1:nt) {
        for (a in 1:na) {
          N_tas[t,a,s] <- R[t+na-a,s] * p[t+na-a,a]
        }
        N[t,s] <- sum(N_tas[t,1:na,s])
        S[t,s] <- N[t,s] * (1 - U[t] * v[s])
        C[t,s] <- N[t,s] * (U[t] * v[s])
      }
    }

    # create calendar year totals across stocks
    for (t in 1:nt) {
      U[t] ~ dbeta(1,1)

      # N_tot[t] <- sum(N[t,1:ns])
      # S_tot[t] <- sum(S[t,1:ns])
      C_tot[t] <- sum(C[t,1:ns])
    }

    # obtain calendar year age composition for each stock
    # use i because not looping over all stocks - only those that have have comps
    for (i in 1:n_age_stocks) {
      for (t in 1:nt) {
        for (a in 1:na) {
          q[t,a,i] <- N_tas[t,a,age_stocks[i]]/N[t,age_stocks[i]]
        }
      }
    }

    # observe calendar year total harvest
    for (t in 1:nt) {
      log_C_tot[t] <- log(C_tot[t])
      C_tot_t_obs[t] ~ dlnorm(log_C_tot[t], tau_C_obs[t])
    }

    # observe calendar year substock specific harvests
    # vectorized to avoid looping over many NAs (hence the i not s)
    for (i in 1:S_obs_n) {
      log_S[i] <- log(S[S_obs_t[i],S_obs_s[i]])
      S_obs[i] ~ dlnorm(log_S[i], tau_S_obs[i])
    }

    # observe age composition
    # only loop over stocks that have age comp data (hence the i not s)
    for (i in 1:n_age_stocks) {
      for (t in 1:nt) {
        x_tas_obs[t,1:na,i] ~ dmulti(q[t,1:na,i], ESS_ts[t,i])
      }
    }
  }

  R2OpenBUGS::write.model(model = mod, file.path(dir_name, file_name))
  if(!silent) return(file.path(dir_name, file_name))
}

#' Write the ssm 3 model code
#'
#' @param dir_name the directory name in your working directory to place the file.
#'   If the directory does not exist, it will be created with a warning
#' @param file_name the file name to save the file as.
#' @param silent logical. Do you wish to suppress the printing of the location of the file?
#'
#' @export

write_ssm_3_file = function(dir_name = "Model Files", file_name = "ssm_3_model.txt", silent = F) {

  # error handle for directory name
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
    warning(paste("No directory named '", dir_name, "' found. One was created.", sep = ""))
  }

  # specify JAGS model code
  mod = function() {
    # Priors on primary parameters
    # mean_sigma_R <- 0.71
    # mean_tau_R <- 1/mean_sigma_R^2
    # tau_R_red <- mean_tau_R * (1 - phi^2)
    # log_resid_0a ~ dnorm(0, tau_R_red)
    log_resid_0a <- 0

    phi ~ dunif(-0.99, 0.99)
    for (s in 1:ns) {
      U_msy[s] ~ dunif(0.01, 0.99)
      log_S_msy[s] ~ dnorm(0, 0.001) %_% I(1, 11.5)
      S_msy[s] <- exp(log_S_msy[s])

      alpha[s] <- exp(U_msy[s])/(1 - U_msy[s])
      log_alpha[s] <- log(alpha[s])
      beta[s] <- U_msy[s]/S_msy[s]
      log_resid_0[s] <- log_resid_0a
    }

    # build the Sigma_R[,] matrix
    sig.common ~ dunif(0,2)
    var.common <- sig.common^2
    rho.common ~ dunif(-0.05,1)
    rho.vec[1] <- 1
    rho.vec[2] <- rho.common

    for (i in 1:vcov_N) {
      rho_mat[vcov_row[i],vcov_col[i]] <- rho.vec[vcov_ind[i]]
      Sigma_R[vcov_row[i],vcov_col[i]] <- var.common * rho.vec[vcov_ind[i]]
    }

    Tau_R[1:ns,1:ns] <- inverse(Sigma_R[1:ns,1:ns])

    # white noise process sd for each substock
    for (s in 1:ns) {
      sigma_R[s] <- sqrt(Sigma_R[s,s])
    }

    # produce Ricker predictions
    for (s in 1:ns) {
      # for years without SR link: use unfished equilibrium recruitment
      R_eq[s] <- log_alpha[s]/beta[s]
      R0[s] <- R_eq[s]
      log_R0[s] <- log(R0[s])

      # log_R_mean1 = deterministic ricker; log_R_mean2 = time-corrected expectation
      log_R_mean1[1,s] <- log_R0[s]
      R_mean1[1,s] <- R0[s]
      log_R_mean2[1,s] <- log_R_mean1[1,s] + phi * log_resid_0[s]
      for (y in 2:a_max) {
        R_mean1[y,s] <- R0[s]
        log_R_mean1[y,s] <- log_R0[s]
        log_R_mean2[y,s] <- log_R_mean1[y,s] + phi * log_resid[y-1,s]
      }

      # for years with SR link
      for (y in (a_max+1):ny) {
        R_mean1[y,s] <- S[y-a_max,s] * exp(log_alpha[s] - beta[s] * S[y-a_max,s])
        log_R_mean1[y,s] <- log(R_mean1[y,s])
        log_R_mean2[y,s] <- log_R_mean1[y,s] + phi * log_resid[y-1,s]
      }
    }

    # draw true recruitment states
    for (y in 1:ny) {
      log_R[y,1:ns] ~ dmnorm(log_R_mean2[y,1:ns], Tau_R[1:ns,1:ns])
    }

    # calculate residuals
    for (y in 1:ny) {
      for (s in 1:ns) {
        R[y,s] <- exp(log_R[y,s])
        log_resid[y,s] <- log_R[y,s] - log_R_mean1[y,s]
      }
      # R_tot[y] <- sum(R[y,1:ns])
    }

    # maturity schedule
    prob[1] ~ dbeta(1,1)
    prob[2] ~ dbeta(1,1)
    prob[3] ~ dbeta(1,1)
    pi[1] <- prob[1]
    pi[2] <- prob[2] * (1 - pi[1])
    pi[3] <- prob[3] * (1 - pi[1] - pi[2])
    pi[4] <- 1 - pi[1] - pi[2] - pi[3]

    D_scale ~ dunif(0.03, 1)
    D_sum <- 1/D_scale^2
    for (a in 1:na) {
      dir_alpha[a] <- D_sum * pi[a]
      for (y in 1:ny) {
        g[y,a] ~ dgamma(dir_alpha[a],1)
        p[y,a] <- g[y,a]/sum(g[y,1:na])
      }
    }

    # allocate R[y,s] to N[t,a,s] and create predicted calendar year states for each stock
    for (s in 1:ns) {
      for (t in 1:nt) {
        for (a in 1:na) {
          N_tas[t,a,s] <- R[t+na-a,s] * p[t+na-a,a]
        }
        N[t,s] <- sum(N_tas[t,1:na,s])
        S[t,s] <- N[t,s] * (1 - U[t] * v[s])
        C[t,s] <- N[t,s] * (U[t] * v[s])
      }
    }

    # create calendar year totals across stocks
    for (t in 1:nt) {
      U[t] ~ dbeta(1,1)

      # N_tot[t] <- sum(N[t,1:ns])
      # S_tot[t] <- sum(S[t,1:ns])
      C_tot[t] <- sum(C[t,1:ns])
    }

    # obtain calendar year age composition for each stock
    # use i because not looping over all stocks - only those that have have comps
    for (i in 1:n_age_stocks) {
      for (t in 1:nt) {
        for (a in 1:na) {
          q[t,a,i] <- N_tas[t,a,age_stocks[i]]/N[t,age_stocks[i]]
        }
      }
    }

    # observe calendar year total harvest
    for (t in 1:nt) {
      log_C_tot[t] <- log(C_tot[t])
      C_tot_t_obs[t] ~ dlnorm(log_C_tot[t], tau_C_obs[t])
    }

    # observe calendar year substock specific harvests
    # vectorized to avoid looping over many NAs (hence the i not s)
    for (i in 1:S_obs_n) {
      log_S[i] <- log(S[S_obs_t[i],S_obs_s[i]])
      S_obs[i] ~ dlnorm(log_S[i], tau_S_obs[i])
    }

    # observe age composition
    # only loop over stocks that have age comp data (hence the i not s)
    for (i in 1:n_age_stocks) {
      for (t in 1:nt) {
        x_tas_obs[t,1:na,i] ~ dmulti(q[t,1:na,i], ESS_ts[t,i])
      }
    }
  }

  R2OpenBUGS::write.model(model = mod, file.path(dir_name, file_name))
  if(!silent) return(file.path(dir_name, file_name))
}

#' Write the ssm 4 model code
#'
#' @param dir_name the directory name in your working directory to place the file.
#'   If the directory does not exist, it will be created with a warning
#' @param file_name the file name to save the file as.
#' @param silent logical. Do you wish to suppress the printing of the location of the file?
#'
#' @export

write_ssm_4_file = function(dir_name = "Model Files", file_name = "ssm_4_model.txt", silent = F) {

  # error handle for directory name
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
    warning(paste("No directory named '", dir_name, "' found. One was created.", sep = ""))
  }

  # specify JAGS model code
  mod = function() {
    # Priors on primary parameters
    # mean_sigma_R <- 0.71
    # mean_tau_R <- 1/mean_sigma_R^2
    # tau_R_red <- mean_tau_R * (1 - phi^2)
    # log_resid_0a ~ dnorm(0, tau_R_red)
    log_resid_0a <- 0

    phi ~ dunif(-0.99, 0.99)
    for (s in 1:ns) {
      U_msy[s] ~ dunif(0.01, 0.99)
      log_S_msy[s] ~ dnorm(0, 0.001) %_% I(1, 11.5)
      S_msy[s] <- exp(log_S_msy[s])

      alpha[s] <- exp(U_msy[s])/(1 - U_msy[s])
      log_alpha[s] <- log(alpha[s])
      beta[s] <- U_msy[s]/S_msy[s]
      log_resid_0[s] <- log_resid_0a
    }

    # estimate the covariance matrix on log(recruitment states
    Tau_R[1:ns,1:ns] ~ dwish(R_wish[1:ns,1:ns], df_wish)
    Sigma_R[1:ns,1:ns] <- inverse(Tau_R)

    # white noise process sd for each substock
    for (s in 1:ns) {
      sigma_R[s] <- sqrt(Sigma_R[s,s])
    }

    # get the pairwise correlation matrix
    for (i in 1:ns) {
      for (j in 1:ns) {
        rho_mat[i,j] <- Sigma_R[i,j]/(sigma_R[i] * sigma_R[j])
      }
    }

    # produce Ricker predictions
    for (s in 1:ns) {
      # for years without SR link: use unfished equilibrium recruitment
      R_eq[s] <- log_alpha[s]/beta[s]
      R0[s] <- R_eq[s]
      log_R0[s] <- log(R0[s])

      # log_R_mean1 = deterministic ricker; log_R_mean2 = time-corrected expectation
      log_R_mean1[1,s] <- log_R0[s]
      R_mean1[1,s] <- R0[s]
      log_R_mean2[1,s] <- log_R_mean1[1,s] + phi * log_resid_0[s]
      for (y in 2:a_max) {
        R_mean1[y,s] <- R0[s]
        log_R_mean1[y,s] <- log_R0[s]
        log_R_mean2[y,s] <- log_R_mean1[y,s] + phi * log_resid[y-1,s]
      }

      # for years with SR link
      for (y in (a_max+1):ny) {
        R_mean1[y,s] <- S[y-a_max,s] * exp(log_alpha[s] - beta[s] * S[y-a_max,s])
        log_R_mean1[y,s] <- log(R_mean1[y,s])
        log_R_mean2[y,s] <- log_R_mean1[y,s] + phi * log_resid[y-1,s]
      }
    }

    # draw true recruitment states
    for (y in 1:ny) {
      log_R[y,1:ns] ~ dmnorm(log_R_mean2[y,1:ns], Tau_R[1:ns,1:ns])
    }

    # calculate residuals
    for (y in 1:ny) {
      for (s in 1:ns) {
        R[y,s] <- exp(log_R[y,s])
        log_resid[y,s] <- log_R[y,s] - log_R_mean1[y,s]
      }
      # R_tot[y] <- sum(R[y,1:ns])
    }

    # maturity schedule
    prob[1] ~ dbeta(1,1)
    prob[2] ~ dbeta(1,1)
    prob[3] ~ dbeta(1,1)
    pi[1] <- prob[1]
    pi[2] <- prob[2] * (1 - pi[1])
    pi[3] <- prob[3] * (1 - pi[1] - pi[2])
    pi[4] <- 1 - pi[1] - pi[2] - pi[3]

    D_scale ~ dunif(0.03, 1)
    D_sum <- 1/D_scale^2
    for (a in 1:na) {
      dir_alpha[a] <- D_sum * pi[a]
      for (y in 1:ny) {
        g[y,a] ~ dgamma(dir_alpha[a],1)
        p[y,a] <- g[y,a]/sum(g[y,1:na])
      }
    }

    # allocate R[y,s] to N[t,a,s] and create predicted calendar year states for each stock
    for (s in 1:ns) {
      for (t in 1:nt) {
        for (a in 1:na) {
          N_tas[t,a,s] <- R[t+na-a,s] * p[t+na-a,a]
        }
        N[t,s] <- sum(N_tas[t,1:na,s])
        S[t,s] <- N[t,s] * (1 - U[t] * v[s])
        C[t,s] <- N[t,s] * (U[t] * v[s])
      }
    }

    # create calendar year totals across stocks
    for (t in 1:nt) {
      U[t] ~ dbeta(1,1)

      # N_tot[t] <- sum(N[t,1:ns])
      # S_tot[t] <- sum(S[t,1:ns])
      C_tot[t] <- sum(C[t,1:ns])
    }

    # obtain calendar year age composition for each stock
    # use i because not looping over all stocks - only those that have have comps
    for (i in 1:n_age_stocks) {
      for (t in 1:nt) {
        for (a in 1:na) {
          q[t,a,i] <- N_tas[t,a,age_stocks[i]]/N[t,age_stocks[i]]
        }
      }
    }

    # observe calendar year total harvest
    for (t in 1:nt) {
      log_C_tot[t] <- log(C_tot[t])
      C_tot_t_obs[t] ~ dlnorm(log_C_tot[t], tau_C_obs[t])
    }

    # observe calendar year substock specific harvests
    # vectorized to avoid looping over many NAs (hence the i not s)
    for (i in 1:S_obs_n) {
      log_S[i] <- log(S[S_obs_t[i],S_obs_s[i]])
      S_obs[i] ~ dlnorm(log_S[i], tau_S_obs[i])
    }

    # observe age composition
    # only loop over stocks that have age comp data (hence the i not s)
    for (i in 1:n_age_stocks) {
      for (t in 1:nt) {
        x_tas_obs[t,1:na,i] ~ dmulti(q[t,1:na,i], ESS_ts[t,i])
      }
    }
  }

  R2OpenBUGS::write.model(model = mod, file.path(dir_name, file_name))
  if(!silent) return(file.path(dir_name, file_name))
}
