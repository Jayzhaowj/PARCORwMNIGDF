########################################################################
#### Multivariate parcor model with shrinkage prior of inverse gamma
########################################################################
build_model <- function(F1_fwd, F1_bwd, P, m, K, n_Total){
  ### the dimension of F1_fwd and F1_bwd: n_Total * K
  ###

  ## the forward PARCOR
  n_1 <- m + 1
  n_T <- n_Total
  data <- rep(list(NA), 2*K)
#  if(m == 2){
#    tmp_x <- cbind(F1_bwd[(n_1-1):(n_T-1), , drop = FALSE], F1_bwd[n_1:n_T-2, , drop = FALSE])
#  }else{
    tmp_x <- F1_bwd[(n_1-m):(n_T-m), , drop = FALSE]
#  }


  y <- F1_fwd[(n_1):n_T, , drop = FALSE]
  index <- 0
  #browser()
  for(k in 1:K){
    index <- index + 1
    if(k==1){
      data[[index]] <- data.frame(y = y[, k, drop = FALSE], tmp_x)
    }else{
      data[[index]] <- data.frame(y = y[, k, drop = FALSE],
                                  cbind(tmp_x,
                                        y[, 1:(k-1)]))
    }
  }
  #browser()

  n_1 <- 1
  n_T <- n_Total - m
#  if(m == 2){
#    tmp_x <- cbind(F1_fwd[(n_1+1):(n_T+1), , drop = FALSE], F1_fwd[(n_1+m):(n_T+m), , drop = FALSE])
#  }else{
    tmp_x <- F1_fwd[(n_1+m):(n_T+m), , drop = FALSE]
#  }

  y <- F1_bwd[n_1:n_T, , drop = FALSE]
  for(k in 1:K){
    index <- index + 1
    if(k == 1){
      data[[index]] <- data.frame(y = y[, k, drop = FALSE],
                                  tmp_x)
    }else{
      data[[index]] <- data.frame(y = y[, k, drop = FALSE],
                                  cbind(tmp_x,
                                        y[, 1:(k-1)]))
    }

  }
  #browser()
  return(data)
}



PARCOR_shrinkage <- function(Y, P,
                             cpus = 1,
                             niter = 10000,
                             nburn = round(niter / 2),
                             nthin = 1,
                             hyperprior_param,
                             display_progress = TRUE,
                             ret_beta_nc = FALSE,
                             SAVS = TRUE,
                             sv = FALSE,
                             sv_param,
                             simple = TRUE,
                             delta){
  nsave <- (niter - nburn)/nthin
  default_hyper <- list(c0 = 2.5,
                        g0 = 5,
                        G0 = 5 / (2.5 - 1),
                        eta2_d1 = 5,
                        eta2_d2 = 4,
                        p_d1 = 1,
                        p_d2 = 1,
                        par_c = 2.5/(10^5))

  if (missing(hyperprior_param)){
    hyperprior_param <- default_hyper
  }
  default_hyper_sv <- list(Bsigma_sv = 1,
                           a0_sv = 5,
                           b0_sv = 1.5,
                           bmu = 0,
                           Bmu = 1)
  if (missing(sv_param) | sv == FALSE){
    sv_param <- default_hyper_sv
  }
  ## number of time series
  K <- ncol(Y)
  ## number of time points
  n_Total <- nrow(Y)

  ## storage
  result_all <- list("median" = NA, "mean" = NA, "qtl" = NA, "qtu" = NA)
  PHI_fwd <- array(NA, dim = c(K^2, n_Total, P))
  PHI_bwd <- array(NA, dim = c(K^2, n_Total, P))
  PHI_star_fwd <- array(NA, dim = c(K^2, n_Total, P))
  PHI_star_bwd <- array(NA, dim = c(K^2, n_Total, P))
  u_inv_fwd <- array(0, dim = c((K^2-K)/2, n_Total, P))
  u_inv_bwd <- array(0, dim = c((K^2-K)/2, n_Total, P))
  theta_sr <- rep(list(NA), P)
  result <- list(PHI_fwd = PHI_fwd, PHI_bwd = PHI_bwd,
                 PHI_star_fwd = PHI_star_fwd, PHI_star_bwd = PHI_star_bwd,
                 u_inv_fwd = u_inv_fwd, u_inv_bwd = u_inv_bwd,
                 SIGMA = NA)

  result_all$median <- result
  result_all$qtl <- result
  result_all$qtu <- result

  data <- build_model(F1_fwd = Y, F1_bwd = Y, P = P,
                      m = 1, K = K, n_Total = n_Total)
  #browser()
  #res <- rep(list(NA), P)


  ##

  PHI_fwd_samp <- array(NA, dim = c(K^2, n_Total, P, nsave))
  PHI_bwd_samp <- array(NA, dim = c(K^2, n_Total, P, nsave))

  #browser()
  ### storage results

  #browser()
  sfInit(parallel = TRUE, cpus = cpus, type = "SOCK")
  sfLibrary(GIGrvg)
  sfLibrary(coda)
  sfLibrary(Rcpp)
  sfLibrary(RcppArmadillo)
  sfLibrary(stochvol)
  sfLibrary(PARCORwMNIGDF)
  #sfClusterEval(sourceCpp("shrinkTVP_mNIG.cpp"))
  #sfSource("shrinkTVP_mNIG.R")
  sfExportAll()
  for(i in 1:P){
    #sfExport("data")
    #browser()
    res_tmp <- sfLapply(1:length(data), function(x) shrinkTVP_mNIG(formula = y ~ .-1, data = data[[x]], K = K,
    #browser()

    #res_tmp <- shrinkTVP_mNIG(formula = y ~ .-1, data = data[[2]], K = K,
                                                         niter = niter,
                                                         nburn = nburn,
                                                         nthin = nthin,
                                                         hyperprior_param = hyperprior_param[[i]],
                                                         display_progress = FALSE,
                                                         ret_beta_nc = ret_beta_nc,
                                                         SAVS = SAVS,
                                                         sv = sv,
                                                         sv_param = sv_param,
                                                         simple = simple,
    #                                                     delta = delta[[2]])
                                                         delta = delta[[i]][[x]])
    )

    #print("hello \n")
    ### obtain median
    if(simple){
      res <- unpack_res(res = res_tmp, m = i, n_Total = n_Total, K = K, type = "beta_median")
      #browser()
      data <- build_model(F1_fwd = res$resid_fwd,
                          F1_bwd = res$resid_bwd, P = P, m = i+1,
                          K = K, n_Total = n_Total)
      #browser()
      result_all$median$PHI_fwd[, , i] <- res$phi_fwd
      result_all$median$PHI_bwd[, , i] <- res$phi_bwd
      result_all$median$PHI_star_fwd[, , i] <- res$phi_star_fwd
      result_all$median$PHI_star_bwd[, , i] <- res$phi_star_bwd
      result_all$median$u_inv_fwd[, , i] <- res$u_inv_fwd
      result_all$median$u_inv_bwd[, , i] <- res$u_inv_bwd
      result_all$median$theta_sr[[i]] <- res_tmp[[1]]$theta_sr
      if(i == P)
        result_all$median$SIGMA <- res$SIGMA


      ### obtain lower bound of 95% credible interval
      res <- unpack_res(res = res_tmp, m = i, n_Total = n_Total, K = K, type = "beta_qtl")

      result_all$qtl$PHI_fwd[, , i] <- res$phi_fwd
      result_all$qtl$PHI_bwd[, , i] <- res$phi_bwd
      result_all$qtl$PHI_star_fwd[, , i] <- res$phi_star_fwd
      result_all$qtl$PHI_star_bwd[, , i] <- res$phi_star_bwd
      result_all$qtl$u_inv_fwd[, , i] <- res$u_inv_fwd
      result_all$qtl$u_inv_bwd[, , i] <- res$u_inv_bwd
      result_all$qtl$theta_sr[[i]] <- res_tmp[[1]]$theta_sr
      if(i == P)
        result_all$qtl$SIGMA <- res$SIGMA


      ### obtain upper bound of 95% credible interval
      res <- unpack_res(res = res_tmp, m = i, n_Total = n_Total, K = K, type = "beta_qtu")

      result_all$qtu$PHI_fwd[, , i] <- res$phi_fwd
      result_all$qtu$PHI_bwd[, , i] <- res$phi_bwd
      result_all$qtu$PHI_star_fwd[, , i] <- res$phi_star_fwd
      result_all$qtu$PHI_star_bwd[, , i] <- res$phi_star_bwd
      result_all$qtu$u_inv_fwd[, , i] <- res$u_inv_fwd
      result_all$qtu$u_inv_bwd[, , i] <- res$u_inv_bwd
      result_all$qtu$theta_sr[[i]] <- res_tmp[[1]]$theta_sr
      if(i == P)
        result_all$qtu$SIGMA <- res$SIGMA

    }else{
      res <- unpack_res_uni(res = res_tmp, m = i, n_Total = n_Total, K = K, nsave = nsave)
      PHI_fwd_samp[, , i, ] <- res$phi_fwd
      PHI_bwd_samp[, , i, ] <- res$phi_bwd
      data <- build_model(F1_fwd = res$resid_fwd,
                          F1_bwd = res$resid_bwd, P = P, m = i+1,
                          K = K, n_Total = n_Total)
      if(i == P)
        SIGMA <- res$SIGMA
    }
    cat("Stage: ", i, "/", P, "\n")
  }
  if(!simple){
    result_all <- list(phi_fwd = PHI_fwd_samp, phi_bwd = PHI_bwd_samp, SIGMA = SIGMA)
  }
  sfStop()
  return(result_all)
}

unpack_res_uni <- function(res, m, n_Total, K, nsave){
  ## unpack PARCOR coefficients
  ### forward time index
  n_1_fwd <- m + 1
  n_T_fwd <- n_Total

  ### backward time index
  n_1_bwd <- 1
  n_T_bwd <- n_Total - m

  ## unpack residuals
  resid_tmp <- t(matrix(unlist(lapply(res, function(x) x$residuals)),
                        nrow = length(res), byrow = TRUE))

  ## unpack observational innovation
  tmp_sigma2 <- simplify2array(lapply(res, function(x) x$sigma2))
  sigma2 <- tmp_sigma2[, , 1:K]
  if(is.matrix(sigma2)){
    sigma2_mean <- apply(sigma2, 2, mean)
  }else if(is.array(sigma2)){
    sigma2_mean <- apply(sigma2, 2:3, mean)
  }else{
    sigma2_mean <- mean(sigma2)
  }
  resid_fwd <- matrix(NA, nrow = n_Total, ncol = K)
  resid_bwd <- matrix(NA, nrow = n_Total, ncol = K)
  phi_fwd <- array(NA, dim = c(n_Total, K^2, nsave))
  phi_bwd <- array(NA, dim = c(n_Total, K^2, nsave))

  #browser()
  resid_fwd[n_1_fwd:n_T_fwd, ] <- resid_tmp[, 1]
  resid_bwd[n_1_bwd:n_T_bwd, ] <- resid_tmp[, 2]
  phi_fwd[n_1_fwd:n_T_fwd, , ] <- res[[1]][["beta"]][-1, , ]
  phi_bwd[n_1_bwd:n_T_bwd, , ] <- res[[2]][["beta"]][-1, , ]

  phi_fwd <- aperm(phi_fwd, c(2, 1, 3))
  phi_bwd <- aperm(phi_bwd, c(2, 1, 3))
  return(list(phi_fwd = phi_fwd,
              phi_bwd = phi_bwd,
              resid_fwd = resid_fwd,
              resid_bwd = resid_bwd,
              SIGMA = sigma2_mean))
}


unpack_res <- function(res, m, n_Total, K, type){
  ## unpack residuals
  resid_tmp <- t(matrix(unlist(lapply(res, function(x) x$residuals)),
                        nrow = length(res), byrow = TRUE))
  #browser()
  ## unpack observational innovation
  tmp_sigma2 <- simplify2array(lapply(res, function(x) x$sigma2))
  sigma2 <- tmp_sigma2[, , 1:K]
  if(is.matrix(sigma2)){
    sigma2_mean <- apply(sigma2, 2, mean)
  }else if(is.array(sigma2)){
    sigma2_mean <- apply(sigma2, 2:3, mean)
  }else{
    sigma2_mean <- mean(sigma2)
  }

  ## storage
  resid_fwd <- matrix(NA, nrow = n_Total, ncol = K)
  resid_bwd <- matrix(NA, nrow = n_Total, ncol = K)
  phi_fwd <- matrix(NA, nrow = K^2, ncol = n_Total)
  phi_bwd <- matrix(NA, nrow = K^2, ncol = n_Total)
  phi_star_fwd <- matrix(NA, nrow = K^2, ncol = n_Total)
  phi_star_bwd <- matrix(NA, nrow = K^2, ncol = n_Total)
  u_inv_fwd <- matrix(0, nrow = (K^2-K)/2, ncol = n_Total)
  u_inv_bwd <- matrix(0, nrow = (K^2-K)/2, ncol = n_Total)
  index_fwd <- 1
  index_bwd <- 1

  ## unpack PARCOR coefficients
  ### forward time index
  n_1_fwd <- m + 1
  n_T_fwd <- n_Total

  ### backward time index
  n_1_bwd <- 1
  n_T_bwd <- n_Total - m

  ###
  for(k in 1:K){
    ### forward
    phi_star_fwd[seq(k, K^2, by = K), n_1_fwd:n_T_fwd] <- t(res[[k]][[type]][-1, 1:K])
    if(k > 1){
      n <- nrow(t(res[[k]][[type]][, 1:(k-1)+K]))
      u_inv_fwd[index_fwd:(index_fwd + n - 1), n_1_fwd:n_T_fwd] <- t(res[[k]][[type]][-1, 1:(k-1)+K])
      index_fwd <- index_fwd + n
    }
    ### backward
    phi_star_bwd[seq(k, K^2, by = K), n_1_bwd:n_T_bwd] <- t(res[[k+K]][[type]][-1, 1:K])

    if(k > 1){
      n <- nrow(t(res[[k+K]][[type]][, 1:(k-1)+K]))
      u_inv_bwd[index_bwd:(index_bwd + n - 1), n_1_bwd:n_T_bwd] <- t(res[[k+K]][[type]][-1, 1:(k-1)+K])
      index_bwd <- index_bwd + n
    }
  }

  SIGMA <- rep(list(NA), n_T_fwd - n_1_fwd + 1)
  ### transformation for forward residuals
  index <- 1
  for(i in n_1_fwd:n_T_fwd){
    u_inv <- diag(K)
    u_inv[lower.tri(u_inv)] <- -u_inv_fwd[, i]
    u <- solve(u_inv, diag(K))
    phi_fwd[, i] <- as.vector(u%*%matrix(phi_star_fwd[, i], nrow = K, ncol = K))
    #browser()
    resid_fwd[i, ] <- as.vector(u%*% matrix(resid_tmp[index, 1:K], ncol = 1))
    #browser()
    ### transformation for innovations
    if(length(sigma2_mean) > 1){
      SIGMA[[i]] <- u %*% diag(sigma2_mean)%*%t(u)
    }else{
      SIGMA[[i]] <- u %*% as.matrix(sigma2_mean)%*%t(u)
    }
    #browser()
    index <- index + 1
  }
  index <- 1
  ### transformation for backward residuals
  for(i in n_1_bwd:n_T_bwd){
    u_inv <- diag(K)
    u_inv[lower.tri(u_inv)] <- -u_inv_bwd[, i]
    u <- solve(u_inv, diag(K))
    phi_bwd[, i] <- as.vector(u%*%matrix(phi_star_bwd[, i], nrow = K, ncol = K))
    resid_bwd[i, ] <- as.vector(u%*% matrix(resid_tmp[index, (K+1):(2*K)], ncol = 1))
    index <- index + 1
  }

  return(list(phi_fwd = phi_fwd,
              phi_bwd = phi_bwd,
              phi_star_fwd = phi_star_fwd,
              phi_star_bwd = phi_star_bwd,
              u_inv_fwd = u_inv_fwd,
              u_inv_bwd = u_inv_bwd,
              resid_fwd = resid_fwd,
              resid_bwd = resid_bwd,
              SIGMA = SIGMA))
}





