#' Quickly summarize the linear model posteriors
#'
#' @param post an object of class \code{mcmc.list}
#' @param params a list object created by \code{SimFit::init_sim()} or \code{kusko_data_prep()$params}.
#'   Contains the driving parameters and dimensional variables.
#' @param seed a numeric vector of length 1: represents a secondary identifier. Typically would be used to identify the seed.
#'   To exclude, set to the default: \code{NA}.
#' @param diag_plots logical. Do you wish to save traceplots and density plots for convergence diagnostics?
#'   Defaults to \code{FALSE}.
#' @param plot_dir character vector representing the directory or full path to the location
#'   where you wish to save the plot files, if \code{diag_plots = TRUE}.
#'   Defaults to \code{NULL}, which will place the files in the working directory.
#' @param return_post logical. Do you wish to return a list with elements equal to
#'   \code{$ests} (the standard estimate summary), \code{$post_lm} (the updated mcmc.list object for the
#'   LM method), and \code{$post_lme} (the updated mcmc.list object for the LME method)?
#'   Defaults to \code{FALSE}.
#'
#' @export

lme_summary = function(post, params, seed = NA, diag_plots = F, plot_dir = NULL, return_post = F) {

  # print message
  cat("  Summarizing LME Model Output", "\n", sep = "")

  # check if post is NULL. if TRUE, that means JAGS crashed.
  if (!is.null(post)) {

    # base parameters to extract
    p = c("alpha", "beta")

    # extract base posterior samples
    post_samps_lm = codaTools::filter_post(
      post = post, p = paste(p, "lm[", sep = "_"), format = "matrix", chains = T, iters = T
    )

    post_samps_lme = codaTools::filter_post(
      post = post, p = paste(p, "lme[", sep = "_"), format = "matrix", chains = T, iters = T
    )

    # number of stocks
    ns = (ncol(post_samps_lme) - 2)/2

    # number of posterior samples
    n_samp = nrow(post_samps_lm)

    # obtain posterior samples of U_msy and S_msy for both lm and lme
    U_msy_post_lm = matrix(NA, n_samp, ns); colnames(U_msy_post_lm) = paste("U_msy_lm[", 1:ns, "]", sep = "")
    S_msy_post_lm = matrix(NA, n_samp, ns); colnames(S_msy_post_lm) = paste("S_msy_lm[", 1:ns, "]", sep = "")
    U_msy_post_lme = matrix(NA, n_samp, ns); colnames(U_msy_post_lme) = paste("U_msy_lme[", 1:ns, "]", sep = "")
    S_msy_post_lme = matrix(NA, n_samp, ns); colnames(S_msy_post_lme) = paste("S_msy_lme[", 1:ns, "]", sep = "")

    a_lm_names = paste("alpha_lm[", 1:ns, "]", sep = "")
    b_lm_names = paste("beta_lm[", 1:ns, "]", sep = "")
    a_lme_names = paste("alpha_lme[", 1:ns, "]", sep = "")
    b_lme_names = paste("beta_lme[", 1:ns, "]", sep = "")
    for (s in 1:ns) {
      temp_lm_mgmt = SimSR::gen_lm_mgmt(
        alpha = post_samps_lm[,a_lm_names[s]],
        beta = post_samps_lm[,b_lm_names[s]])

      U_msy_post_lm[,s] = temp_lm_mgmt$U_msy
      S_msy_post_lm[,s] = temp_lm_mgmt$S_msy

      temp_lme_mgmt = SimSR::gen_lm_mgmt(
        alpha = post_samps_lme[,a_lme_names[s]],
        beta = post_samps_lme[,b_lme_names[s]])
      U_msy_post_lme[,s] = temp_lme_mgmt$U_msy
      S_msy_post_lme[,s] = temp_lme_mgmt$S_msy
    }

    # calculate drainage-wide reference points
    mgmt_post_lm = mgmt_post_lme = matrix(NA, n_samp, 4)
    colnames(mgmt_post_lm) = colnames(mgmt_post_lme) = c("S_obj", "U_obj", "S_MSY", "U_MSY")
    for (i in 1:n_samp) {
      mgmt_post_lm[i,] = SimSR::gen_mgmt(
        params = list(
          alpha = post_samps_lm[i,a_lm_names],
          beta = post_samps_lm[i,b_lm_names],
          U_msy = U_msy_post_lm[i,],
          S_msy = S_msy_post_lm[i,],
          U_range = seq(0,1,0.01),
          max_p_overfished = params$max_p_overfished,
          ns = ns)
      )$mgmt

      mgmt_post_lme[i,] = SimSR::gen_mgmt(
        params = list(
          alpha = post_samps_lme[i,a_lme_names],
          beta = post_samps_lme[i,b_lme_names],
          U_msy = U_msy_post_lme[i,],
          S_msy = S_msy_post_lme[i,],
          U_range = seq(0,1,0.01),
          max_p_overfished = params$max_p_overfished,
          ns = ns)
      )$mgmt
    }

    # combine all posterior samples into one big matrix
    post_samps_lm = cbind(
      post_samps_lm,
      U_msy_post_lm,
      S_msy_post_lm,
      mgmt_post_lm
    )

    post_samps_lme = cbind(
      post_samps_lme,
      U_msy_post_lme,
      S_msy_post_lme,
      mgmt_post_lme
    )

    # coerce to mcmc.list
    post_samps_lm = codaTools::matrix2mcmclist(post_samps_lm)
    post_samps_lme = codaTools::matrix2mcmclist(post_samps_lme)

    # extract posterior summaries
    post_summs_lm = codaTools::summ_post(
      post = post_samps_lm,
      p = paste("^", codaTools::get_nodes(post_samps_lm), sep = ""),
      ess = T, bgr = T
    )

    post_summs_lme = codaTools::summ_post(
      post = post_samps_lme,
      p = paste("^", codaTools::get_nodes(post_samps_lme), sep = ""),
      ess = T, bgr = T
    )

    # create the "id variables" for each output model
    id_lm = data.frame(
      seed = seed,
      param = stringr::str_remove(colnames(post_summs_lm), "\\_lm\\[.+\\]"),
      stock = c(rep(1:ns, 4),
                rep(NA, 4)),
      method = "lm",
      stringsAsFactors = F
    )

    id_lme = data.frame(
      seed = seed,
      param = stringr::str_remove(colnames(post_summs_lme), "\\_lme\\[.+\\]"),
      stock = c(rep(1:ns, 4),
                rep(NA, 4)),
      method = "lme",
      stringsAsFactors = F
    )

    # combine
    ests = rbind(
      cbind(id_lm, t(post_summs_lm)),
      cbind(id_lme, t(post_summs_lme))
    )

    # do the diagnostic plot if desired
    if (diag_plots) {


      file_lm = paste("lm", "_",
                      ifelse(is.na(seed), "", paste(seed, "_", sep = "")),
                      "diag_plots.pdf", sep = "")
      file_lme = paste("lme", "_",
                      ifelse(is.na(seed), "", paste(seed, "_", sep = "")),
                      "diag_plots.pdf", sep = "")

      if (!is.null(plot_dir)) {
        file_lm = file.path(plot_dir, file_lm)
        file_lme = file.path(plot_dir, file_lme)
      }

      codaTools::diag_plots(
        post = post_samps_lm,
        p = c("alpha_lm", "beta_lm", "U_msy_lm" ,"S_msy_lm",
              "U_obj", "S_obj", "S_MSY", "U_MSY"),
        save = T,
        file = file_lm
      )

      codaTools::diag_plots(
        post = post_samps_lme,
        p = c("alpha_lme", "beta_lme", "U_msy_lme" ,"S_msy_lme",
              "U_obj", "S_obj", "S_MSY", "U_MSY"),
        save = T,
        file = file_lme
      )
    }

  } else { # if is.null(post)

    ests = data.frame(seed = seed, param = NA, stock = NA,
                        method = c("lm", "lme"), mean = NA, sd = NA,
                        x1 = NA, x2 = NA, x3 = NA)

    colnames(ests)[(ncol(ests) - 2):ncol(ests)] = c("50%", "2.5%", "97.5%")
    ests$bgr = NA
    ests$ess = NA
  }

  # remove rownames from summary
  rownames(ests) = NULL

  if (return_post) {
    output = list(
      post_lm = post_samps_lm,
      post_lme = post_samps_lme,
      ests = ests
    )
  } else {
    output = ests
  }

  # return output
  return(output)
}
