#' Quickly summarize the state-space model posteriors
#'
#' @param post an object of class \code{mcmc.list}
#' @param params a list object created by \code{SimFit::init_sim()} or \code{kusko_data_prep()$params}.
#'   Contains the driving parameters and dimensional variables.
#' @param model the model identifier: e.g., a number or letter
#' @param maturity character vector of length 1: does the model have \code{"simple"} or \code{"complex"}
#'   maturity schedules? The extracted parameters depends on this
#' @param diag_plots logical. Do you wish to save traceplots and density plots for convergence diagnostics?
#'   Defaults to \code{FALSE}.
#' @param file character vector of length 1: optional file name of the saved diagnostic plots.
#' @param seed a numeric vector of length 1: represents a secondary identifier. Typically would be used to identify the seed.
#'   To exclude, set to the default: \code{NA}.
#' @param plot_dir character vector representing the directory or full path to the location
#'   where you wish to save the plot files, if \code{diag_plots = TRUE}.
#'   Defaults to \code{NULL}, which will place the files in the working directory.
#' @param return_post logical. Do you wish to return a list with elements equal to
#'   \code{$ests} (the standard estimate summary) and \code{$post} (the updated mcmc.list object)?
#'   Defaults to \code{FALSE}.
#' @export

ssm_summary = function(post, params, model, maturity, diag_plots = F, file = NULL, seed = NA, plot_dir = NULL, return_post = F) {

  # print message
  cat("  Summarizing SSM Model #", model, " Output", "\n", sep = "")

  # check if post is NULL. if TRUE, that means JAGS crashed.
  if (!is.null(post)) {

    # base parameters to extract
    p = c("alpha", "beta", "U_msy", "S_msy", "sigma_R", "pi", "phi", "R", "^S\\[", "C_tot", "^U\\[")

    # add on "D_sum" doing complex maturity
    if (maturity == "complex") p = c(p, "D_sum")

    # extract base posterior samples for parameters in all models
    post_samps = codaTools::filter_post(
      post = post, p = p, format = "matrix", chains = T, iters = T
    )

    # get mean rho posterior
    rho_mat_post = codaTools::filter_post(post, "rho_mat", format = "matrix")
    diag_names = paste("rho_mat[", 1:params$ns, ",", 1:params$ns, "]", sep = "")
    mean_rho_post = apply(rho_mat_post[,-which(colnames(rho_mat_post) %in% diag_names)], 1, mean)

    # get mean sigma R posterior
    sigma_R_post = codaTools::filter_post(post, p = "sigma_R", format = "matrix")
    mean_sigma_R_post = apply(sigma_R_post, 1, mean)

    # calculate drainage-wide reference points for each mcmc iteration
    a_names = paste("alpha[", 1:params$ns, "]", sep = "")
    b_names = paste("beta[", 1:params$ns, "]", sep = "")
    U_names = paste("U_msy[", 1:params$ns, "]", sep = "")
    S_names = paste("S_msy[", 1:params$ns, "]", sep = "")

    n_samp = length(mean_sigma_R_post)
    mgmt_post = matrix(NA, n_samp, 8)
    colnames(mgmt_post) = c("Sstar_0.1", "Sstar_0.3", "Sstar_0.5", "Ustar_0.1", "Ustar_0.3", "Ustar_0.5", "S_MSY", "U_MSY")
    for (i in 1:n_samp) {
      mgmt_post[i,] = SimSR::gen_mgmt(
        params = list(
          alpha = post_samps[i,a_names],
          beta = post_samps[i,b_names],
          U_msy = post_samps[i,U_names],
          S_msy = post_samps[i,S_names],
          U_range = seq(0,1,0.01),
          ns = params$ns)
      )$mgmt
    }

    # combine all posterior samples into one big matrix
    post_samps = cbind(
      post_samps,
      mean_sigma_R = mean_sigma_R_post,
      mean_rho = mean_rho_post,
      mgmt_post)

    # coerce to mcmc.list
    post_samps = codaTools::matrix2mcmclist(post_samps)

    # extract posterior summaries
    post_summs = codaTools::summ_post(
      post = post_samps,
      p = paste("^", codaTools::get_nodes(post_samps), sep = ""),
      ess = T, bgr = T
    )

    # combine output
    id = data.frame(
      seed = seed,
      param = stringr::str_remove(colnames(post_summs), "\\[.+\\]"),
      stringsAsFactors = F
    )
    id$stock = NA
    id$stock[id$param == "alpha"] = 1:params$ns
    id$stock[id$param == "beta"] = 1:params$ns
    id$stock[id$param == "U_msy"] = 1:params$ns
    id$stock[id$param == "S_msy"] = 1:params$ns
    id$stock[id$param == "sigma_R"] = 1:params$ns
    id$stock[id$param == "R"] = rep(1:params$ns, each = params$ny)
    id$stock[id$param == "S"] = rep(1:params$ns, each = params$nt)

    id$year = NA
    id$year[id$param == "R"] = rep(1:params$ny, params$ns)
    id$year[id$param == "S"] = rep(1:params$nt, params$ns)
    id$year[id$param == "C_tot"] = 1:params$nt
    id$year[id$param == "U"] = 1:params$nt

    id$method = paste("ssm", model, sep = "")

    ests = cbind(id, t(post_summs))

    # do the diagnostic plot if desired
    if (diag_plots) {

      # include a nice file name if not supplied
      if (is.null(file)) {
        file = paste("ssm_", model, "_",
              ifelse(is.na(seed), "", seed),
              "_diag_plots.pdf", sep = "")
      }

      if (!is.null(plot_dir)) {
        file = file.path(plot_dir, file)
      }

      if ("D_sum" %in% codaTools::get_nodes(post_samps)) {
        inc_D_sum = "D_sum"
      } else {
        inc_D_sum = NULL
      }

      codaTools::diag_plots(
        post = post_samps,
        p = c("U_MSY", "S_MSY", "mean_sigma_R",
              "mean_rho", "phi", inc_D_sum,
              "alpha", "beta", "U_msy", "S_msy"),
        save = T,
        file = file
      )
    }

  } else {  # if is.null(post)
    # just return a blank df
    ests = data.frame(seed = seed, param = NA, stock = NA, year = NA,
                        method = paste("ssm", model, sep = ""), mean = NA, sd = NA,
                        x1 = NA, x2 = NA, x3 = NA)

    colnames(ests)[(ncol(ests) - 2):ncol(ests)] = c("50%", "2.5%", "97.5%")
    ests$bgr = NA
    ests$ess = NA
  }

  # remove rownames from summary
  rownames(ests) = NULL

  if (return_post) {
    output = list(
      post = post_samps,
      ests = ests
    )
  } else {
    output = ests
  }

  # return outuput
  return(output)

}
