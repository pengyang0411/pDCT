##---------------------------------
## Functions for running DCTs
## Author: Peng Yang
##---------------------------------

##-----------------------------
## Mean function for outcomets
##-----------------------------

fi <- function(t, phi){
  
  return(phi[1] + phi[2] * t^2)

  
}



##-----------------------------
## Function to generate data
##-----------------------------

sim_data <- function(n     = 50,  # number of subjects
                     n_T   = 5,   # number of observations per subject
                     n_p   = 5,   # number of baseline covariates
                     d_s   = 1,
                     sig_1 = 0.1, # variance for onsite visits
                     sig_2 = 5,   # variance for offsite visits
                     if_null = FALSE,
                     if_plot = FALSE) {
  
  # Time range
  t_all <- rep(seq(0, 1, length.out = n_T), n)
  
  # Treatment assignment
  X <- rbinom(n, size = 1, prob = 0.5)
  
  # Site assignment schedule
  S <- rep(1, n_T)
  S[seq(2, n_T - 1, by = 2)] <- 0
  
  # Expanded treatment/site vectors
  X_all <- rep(X, each = n_T)
  S_all <- rep(S, n)
  
  # Subject ID
  ID_all <- rep(1:n, each = n_T)
  
  # Baseline covariates
  Z <- matrix(runif(n * n_p, min = -0.5, max = 0.5), nrow = n)
  Z_scaled <- scale(Z, center = TRUE, scale = FALSE)
  
  # Visit-specific SD
  sig <- S * sig_1 + (1 - S) * sig_2
  
  # Subject-specific random effects
  phi <- cbind(
    rnorm(n = n, mean = 0, sd = 0.1),
    rnorm(n = n, mean = 0, sd = 0.1)
  )
  
  # Simulate potential outcomes
  fi_all <- c()
  Y0 <- Y1 <- c()
  
  for (i in 1:n) {
    tt <- t_all[((i - 1) * n_T + 1):(i * n_T)]
    Si <- S_all[((i - 1) * n_T + 1):(i * n_T)]
    
    ff   <- fx(tt)
    ffi  <- fi(tt, phi[i, ])
    bias <- Gamma_bias(tt)
    
    mu_Y0 <- ff + ffi +
      d_s * (1 - Si) +
      fx_a0(tt, n_p) %*% Z_scaled[i, ]^2 +
      fx_b0(tt, n_p) %*% Z_scaled[i, ]
    
    mu_Y1 <- ff + ffi +
      bias * (1 - Si) +
      d_s * (1 - Si) +
      fx_a1(tt, n_p) %*% Z_scaled[i, ]^2 +
      fx_b1(tt, n_p) %*% Z_scaled[i, ] +
      Delta(tt) * (1 - as.numeric(if_null))
    
    Y0 <- c(Y0, rnorm(n = n_T, mean = mu_Y0, sd = sig))
    Y1 <- c(Y1, rnorm(n = n_T, mean = mu_Y1, sd = sig))
    fi_all <- c(fi_all, ffi)
  }
  
  # Observed outcomes
  Y <- X_all * Y1 + (1 - X_all) * Y0
  
  if (if_plot) {
    plot(t_all, Y)
    lines(t_all[X_all == 0], Y[X_all == 0], type = "p", col = "orange")
    lines(t_all[X_all == 1], Y[X_all == 1], type = "p", col = "light blue")
    
    for (t in unique(t_all)) {
      lines(t, mean(Y[t_all == t & X_all == 1]),
            type = "p", pch = 3, cex = 3, col = "dark blue")
      lines(t, mean(Y[t_all == t & X_all == 0]),
            type = "p", pch = 3, cex = 3, col = "dark orange")
    }
  }
  
  out <- list(
    Z  = Z,
    Y  = Y,
    X  = X,
    S  = S,
    Y0 = Y0,
    Y1 = Y1,
    fi_all = fi_all,
    phi = phi,
    
    n = n,
    n_T = n_T,
    n_p = n_p,
    
    S_all = S_all,
    X_all = X_all,
    t_all = t_all,
    ID_all = ID_all
  )
  
  return(out)
}


##-----------------------------
## Functions to train
## 4 separate tsBART model
##-----------------------------


## Four seperate models conditoning on A
train_tsBART <- function(Y, t, id, X_all, S_all, Z_all,
                         nsim, ntree, if_all_onsite) {

  ## Test data set
  data_test <- data.frame(Y = Y,
                          t = t,
                          id = id)

  out <- list()

  count <- 1
  
  for(s in c(0,1)){
    
    for(a in c(0,1)){
      
      ## Training data
      data_train <- data.frame(
        Y = Y[which(X_all == x & S_all == s)],
        t = t[which(X_all == x & S_all == s)],
        id = id[which(X_all == x & S_all == s)]
      )
      
      zz_train <- Z_all[which(X_all == x & S_all == s), ]
      
      ##-------------------------------
      ## Hyperparameters setup
      ##-------------------------------
      
      # Calibrate BART's error variance a la CGM 2010 (method 2)
      df = data.frame(data_train, 't'=data_train$t, 'y'=data_train$Y)
      lmf = lm(y~.,df)
      sighat = sigma(lmf)
      
      # Hyperparameters
      nu = 3
      sigq = .9
      qchi = qchisq(1.0-sigq,nu)
      lambda = (sighat*sighat*qchi)/nu
      
      
      ## Tune length-scale via the expected number of crossings
      
      # Evaluate optimal number of crossings.
      ecross_candidates = seq(.5,3,by=.5)
      ecrossTune = tuneEcross(ecross_candidates,
                              y=data_train$Y, tgt=data_train$t,  x=zz_train,
                              nburn=100, nsim=100, ntree=200,
                              lambda=lambda, sigq=sigq, sighat=sighat, nu=nu,
                              base_tree=.95, power_tree=2,
                              probit=FALSE, yobs=NULL)
      
      # Set expected number of crossings.
      exp_cross = ecrossTune[["ecross_opt"]]
      
      # Print plot.
      # print(ecrossTune[["waic_plot"]])
      
      ## Fit Targeted Smooth BART Model
      fit <- tsbart(
        y = data_train$Y, tgt = data_train$t, tpred = data_test$t,
        x = zz_train, xpred = Z_all,
        nburn = 0, nsim = nsim, ntree = ntree,
        lambda = lambda, sigq = .9, nu = 3,
        ecross = exp_cross, base_tree = .95, power_tree = 2,
        probit = FALSE, yobs = NULL, verbose = FALSE
      )
      
      out[[count]] <- fit
      count <- count + 1
      
      
    }
    
  }
    
    

  

  return(out)

}

##-----------------------------
## Functions to train to 
## estimate the ATE
##-----------------------------

run_LDT <- function(data, 
                    K = 2, ## Order of time
                    L = 2, ## Number of knots
                    d = 2, ## Degree of knots
                    G = K + 1 + L,
                    nsamp = 1e3, ## Number of posterior draws
                    nburn = round(nsamp / 2),
                    ntree = 200, ## Number of trees for BART model
                    if_plot = FALSE,
                    if_site = TRUE,
                    if_digital_twins = TRUE,
                    if_covariates = TRUE,
                    if_model_diagnostic = FALSE, 
                    plot_path = NULL
){
  
  ## Output from data
  Z <- data$Z
  Y <- data$Y
  n <- data$n
  n_T <- data$n_T
  n_p <- data$n_p
  
  S_all <- data$S_all
  X_all <- data$X_all
  t_all <- data$t_all
  ID_all <- data$ID_all
  
  X <- data$X
  
  n_S0   <- sum(S_all == 0) ## Offsite
  n_S1   <- sum(S_all == 1) ## Onsite
  Indx_S <- which(S_all == 1) ## Indices of onsite visits
  
  if_all_onsite <- ifelse(n_S0 == 0, TRUE, FALSE)
  
  ##----------------------------
  ## Loading the design matrix
  ##----------------------------
  Q_f <- array(0, dim = c(n * n_T, G))
  
  kappa <- seq(0, 1, length.out = L + 2)
  kappa <- kappa[-c(1, L + 2)]
  
  ## Design matrix for f
  for(k in 0:K){
    Q_f[, k + 1] <- t_all^k
  }
  
  for(l in 1:L){
    Q_f[, K + 1 + l] <- ifelse(t_all <= kappa[l], (t_all - kappa[l])^d, 0)
  }
  
  ## Design matrix for tau
  Q_tau <- Q_f[, -1] * X_all
  
  ## Design matrix for gamma
  Q_gamma <- Q_f[, -1] * X_all * (1 - S_all)
  
  ## Design matrix to expand subject-level covariates to visit level
  Q_z <- array(0, dim = c(n * n_T, n))
  for(i in 1:n){
    Q_z[(i * n_T - (n_T - 1)):(i * n_T), i] <- 1
  }
  
  ## Design matrix for subject-specific random effects
  Q_f_all <- array(0, c(n * n_T, n * G))
  for(i in 1:n){
    Q_f_all[((i - 1) * n_T + 1):(i * n_T), ((i - 1) * G + 1):(i * G)] <- 
      Q_f[((i - 1) * n_T + 1):(i * n_T), ]
  }
  
  ##----------------------------
  ## Running Gibbs sampler
  ##----------------------------
  
  ##----------------------------
  ## Train tsBART model
  ##----------------------------
  Z_scaled <- Q_z %*% scale(Z, center = TRUE, scale = FALSE)
  
  if(if_digital_twins){
    
    obj_tsBART <- train_tsBART(
      Y = Y, t = t_all, id = ID_all, 
      X_all = X_all, S_all = S_all, Z_all = Q_z %*% Z,
      nsim = nsamp, ntree = ntree, if_all_onsite = if_all_onsite
    )
    
    ## Center the predicted DT by time
    for(t in unique(t_all)){
      obj_tsBART[[1]]$mcmcdraws_oos[, t_all == t] <- 
        t(scale(t(obj_tsBART[[1]]$mcmcdraws_oos[, t_all == t]), center = TRUE, scale = FALSE))
      obj_tsBART[[2]]$mcmcdraws_oos[, t_all == t] <- 
        t(scale(t(obj_tsBART[[2]]$mcmcdraws_oos[, t_all == t]), center = TRUE, scale = FALSE))
      if(!if_all_onsite){
        obj_tsBART[[3]]$mcmcdraws_oos[, t_all == t] <- 
          t(scale(t(obj_tsBART[[3]]$mcmcdraws_oos[, t_all == t]), center = TRUE, scale = FALSE))
        obj_tsBART[[4]]$mcmcdraws_oos[, t_all == t] <- 
          t(scale(t(obj_tsBART[[4]]$mcmcdraws_oos[, t_all == t]), center = TRUE, scale = FALSE))
      }
    }
    
    if(if_all_onsite){
      m_0 <- cbind(
        obj_tsBART[[1]]$mcmcdraws_oos[1, ],
        obj_tsBART[[2]]$mcmcdraws_oos[1, ]
      )
    }else{
      m_0 <- cbind(
        obj_tsBART[[1]]$mcmcdraws_oos[1, ],
        obj_tsBART[[2]]$mcmcdraws_oos[1, ],
        obj_tsBART[[3]]$mcmcdraws_oos[1, ],
        obj_tsBART[[4]]$mcmcdraws_oos[1, ]
      )
    }
    
    if(if_covariates){
      Z_aug <- cbind(Z_scaled, m_0)
    }else{
      Z_aug <- cbind(m_0)
    }
    
    if(if_site) Z_aug <- cbind(1 - S_all, Z_aug)
    
  }else{
    
    Z_aug <- c()
    
    if(if_covariates){
      Z_aug <- Z_scaled
    }
    
    if(if_site & !if_covariates) Z_aug <- matrix(1 - S_all, ncol = 1)
    if(if_site & if_covariates)  Z_aug <- cbind(1 - S_all, Z_aug)
  }
  
  n_p_all <- ncol(Z_aug)
  if(is.null(n_p_all)) n_p_all <- 0
  
  if(if_site){
    Q_fixed <- cbind(Q_f, Q_tau, Q_gamma, Z_aug)
  }else{
    Q_fixed <- cbind(Q_f, Q_tau, Z_aug)
  }
  
  ## Parameters
  post_fixed <- array(0, c(nsamp, ncol(Q_fixed)))
  post_f_i <- array(0, c(nsamp, n * G))
  post_gamma <- array(0, c(nsamp, G))
  
  ## Variance parameters
  post_sigma_e <- post_sigma_delta <- post_sigma_delta_tau <- post_sigma_phi_tau <- 
    post_sigma_delta_gamma <- post_sigma_phi_gamma <- rep(0, nsamp)
  post_sigma_delta_f_all <- post_sigma_phi_f_all <- array(0, c(nsamp, n))
  post_sigma_e <- array(0, c(nsamp, 2))
  
  ## Initialization
  f_i_s <- rnorm(n * G)
  
  sigma_e_s <- rep(1, 2)
  sigma_delta_s <- sigma_delta_tau_s <- sigma_phi_gamma_s <- sigma_delta_gamma_s <- 1
  sigma_delta_f_all_s <- sigma_phi_f_all_s <- rep(1, n)
  
  sigma_phi <- sigma_beta <- sigma_tau <- sigma_gamma <- 1e3
  sigma_gamma <- 1
  sigma_phi_tau_s <- 1
  a <- b <- 1e-2
  
  ## Construct covariance matrix inverse for subject-specific effects
  Sigma_f_all_inv <- diag(NA, n * G)
  
  Indx_K <- c()
  for(k in 1:(K + 1)) Indx_K <- c(Indx_K, seq(k, n * G, by = G))
  Indx_K <- sort(Indx_K)
  
  Indx_L <- c()
  for(l in 1:L) Indx_L <- c(Indx_L, seq(K + 1 + l, n * G, by = G))
  Indx_L <- sort(Indx_L)
  
  diag(Sigma_f_all_inv)[Indx_K] <- 1 / rep(sigma_phi_f_all_s, each = (K + 1))
  diag(Sigma_f_all_inv)[Indx_L] <- 1 / rep(sigma_delta_f_all_s, each = L)
  
  ##----------------------------
  ## Gibbs sampler
  ##----------------------------
  for(s in 1:nsamp){
    
    ##------------------------------
    ## Update digital twins
    ##------------------------------
    if(if_digital_twins){
      
      if(if_all_onsite){
        m_s <- cbind(
          obj_tsBART[[1]]$mcmcdraws_oos[s, ],
          obj_tsBART[[2]]$mcmcdraws_oos[s, ]
        )
      }else{
        m_s <- cbind(
          obj_tsBART[[1]]$mcmcdraws_oos[s, ],
          obj_tsBART[[2]]$mcmcdraws_oos[s, ],
          obj_tsBART[[3]]$mcmcdraws_oos[s, ],
          obj_tsBART[[4]]$mcmcdraws_oos[s, ]
        )
      }
      
      if(if_covariates){
        if(if_site){
          Z_aug <- cbind(1 - S_all, Z_scaled, m_s)
          Q_fixed <- cbind(Q_f, Q_tau, Q_gamma, Z_aug)
        }else{
          Z_aug <- cbind(Z_scaled, m_s)
          Q_fixed <- cbind(Q_f, Q_tau, Z_aug)
        }
      }else{
        if(if_site){
          Z_aug <- cbind(1 - S_all, m_s)
          Q_fixed <- cbind(Q_f, Q_tau, Q_gamma, Z_aug)
        }else{
          Z_aug <- cbind(m_s)
          Q_fixed <- cbind(Q_f, Q_tau, Z_aug)
        }
      }
    }
    
    ##---------------------------------
    ## Block 1: sample population-level
    ## mean trajectory, ATE trajectory,
    ## bias trajectory, and fixed effects
    ##---------------------------------
    if(if_site){
      Sigma_fixed_inv <- diag(c(
        rep(1 / sigma_phi, K + 1),       ## f(t): polynomial coefficients
        rep(1 / sigma_delta_s, L),       ## f(t): spline coefficients
        rep(1 / sigma_phi_tau_s, K),     ## tau(t): polynomial coefficients
        rep(1 / sigma_delta_tau_s, L),   ## tau(t): spline coefficients
        rep(1 / sigma_phi_gamma_s, K),   ## gamma(t): polynomial coefficients
        rep(1 / sigma_delta_gamma_s, L), ## gamma(t): spline coefficients
        rep(1 / sigma_beta, n_p_all)     ## fixed-effect coefficients
      ))
    }else{
      Sigma_fixed_inv <- diag(c(
        rep(1 / sigma_phi, K + 1),       ## f(t): polynomial coefficients
        rep(1 / sigma_delta_s, L),       ## f(t): spline coefficients
        rep(1 / sigma_phi_tau_s, K),     ## tau(t): polynomial coefficients
        rep(1 / sigma_delta_tau_s, L),   ## tau(t): spline coefficients
        rep(1 / sigma_beta, n_p_all)     ## fixed-effect coefficients
      ))
    }
    
    Sigma_e_inv <- diag(rep(1 / sigma_e_s[1], n * n_T) * (1 - S_all)) +
      diag(rep(1 / sigma_e_s[2], n * n_T) * S_all)
    
    V_fixed <- Sigma_fixed_inv + t(Q_fixed) %*% Sigma_e_inv %*% Q_fixed
    mu_fixed <- t(Q_fixed) %*% Sigma_e_inv %*% (Y - Q_f_all %*% f_i_s)
    
    ch_Q <- chol(V_fixed)
    fixed_s <- backsolve(
      ch_Q,
      forwardsolve(t(ch_Q), mu_fixed) + rnorm(ncol(Q_fixed))
    )
    
    f_s   <- fixed_s[1:G]
    tau_s <- fixed_s[(G + 1):(2 * G - 1)]
    if(if_site) gamma_s <- fixed_s[(2 * G):(3 * G - 2)]
    
    ##------------------------------
    ## Block 2: sample subject-level
    ## random effects
    ##------------------------------
    diag(Sigma_f_all_inv)[Indx_K] <- 1 / rep(sigma_phi_f_all_s, each = (K + 1))
    diag(Sigma_f_all_inv)[Indx_L] <- 1 / rep(sigma_delta_f_all_s, each = L)
    
    V_f_all <- Sigma_f_all_inv + t(Q_f_all) %*% Sigma_e_inv %*% Q_f_all
    mu_f_all <- t(Q_f_all) %*% Sigma_e_inv %*% (Y - Q_fixed %*% fixed_s)
    
    ch_Q <- chol(V_f_all)
    f_i_s <- backsolve(
      ch_Q,
      forwardsolve(t(ch_Q), mu_f_all) + rnorm(n * G)
    )
    
    ##------------------------------
    ## Block 3: sample variance of
    ## spline coefficients in f(t)
    ##------------------------------
    sigma_delta_s <- 1 / rgamma(
      n = 1,
      shape = a + L / 2,
      rate = b + 1 / 2 * sum(f_s[-(1:(K + 1))]^2)
    )
    
    ##------------------------------
    ## Block 4: sample subject-level
    ## random-effect variances
    ##------------------------------
    f_i_s_reshape <- matrix(f_i_s, byrow = TRUE, ncol = G)
    
    if(L == 1){
      sigma_delta_f_all_s <- 1 / rgamma(
        n = n,
        shape = a + L / 2,
        rate = b + 1 / 2 * f_i_s_reshape[, -(1:(K + 1))]^2
      )
    }else{
      sigma_delta_f_all_s <- 1 / rgamma(
        n = n,
        shape = a + L / 2,
        rate = b + 1 / 2 * rowSums(f_i_s_reshape[, -(1:(K + 1))]^2)
      )
    }
    
    sigma_phi_f_all_s <- 1 / rgamma(
      n = n,
      shape = a + (K + 1) / 2,
      rate = b + 1 / 2 * rowSums(f_i_s_reshape[, 1:(K + 1)]^2)
    )
    
    ##------------------------------
    ## Block 5: sample variances for
    ## tau(t)
    ##------------------------------
    sigma_delta_tau_s <- 1 / rgamma(
      n = 1,
      shape = a + L / 2,
      rate = b + 1 / 2 * sum(tau_s[-(1:K)]^2)
    )
    
    sigma_phi_tau_s <- 1 / rgamma(
      n = 1,
      shape = a + K / 2,
      rate = b + 1 / 2 * sum(tau_s[1:K]^2)
    )
    
    if(if_site){
      ##------------------------------
      ## Block 6: sample variances for
      ## gamma(t)
      ##------------------------------
      sigma_phi_gamma_s <- 1 / rgamma(
        n = 1,
        shape = a + K / 2,
        rate = b + 1 / 2 * sum(gamma_s[1:K]^2)
      )
      
      sigma_delta_gamma_s <- 1 / rgamma(
        n = 1,
        shape = a + L / 2,
        rate = b + 1 / 2 * sum(gamma_s[-(1:K)]^2)
      )
    }
    
    ##------------------------------
    ## Block 7: sample residual
    ## variances
    ##------------------------------
    mu_Y <- Y - Q_fixed %*% fixed_s - Q_f_all %*% f_i_s
    
    sigma_e_s <- 1 / rgamma(
      n = 2,
      shape = a + c(n_S0, n_S1) / 2,
      rate = b + 1 / 2 * c(sum(mu_Y[-Indx_S]^2), sum(mu_Y[Indx_S]^2))
    )
    
    ## Store posterior samples
    post_fixed[s, ] <- fixed_s
    post_f_i[s, ] <- f_i_s
    
    post_sigma_e[s, ] <- sigma_e_s
    post_sigma_delta[s] <- sigma_delta_s
    post_sigma_delta_f_all[s, ] <- sigma_delta_f_all_s
    post_sigma_phi_f_all[s, ] <- sigma_phi_f_all_s
    post_sigma_delta_tau[s] <- sigma_delta_tau_s
    post_sigma_phi_tau[s] <- sigma_phi_tau_s
    post_sigma_phi_gamma[s] <- sigma_phi_gamma_s
    post_sigma_delta_gamma[s] <- sigma_delta_gamma_s
  }
  
  post_f <- post_fixed[, 1:G]
  post_tau <- post_fixed[, (G + 1):(2 * G - 1)]
  
  if(if_site) post_gamma <- post_fixed[, (2 * G):(3 * G - 2)]
  
  if(if_covariates | if_site){
    post_beta <- matrix(post_fixed[, (3 * G - 1):(3 * G + n_p_all - 2)], nrow = nsamp)
  }else{
    post_beta <- array(0, c(nsamp, 1))
  }
  
  post_ATE <- Q_f[1:n_T, -1] %*% t(post_tau[-(1:nburn), ])
  Prob <- apply(post_ATE, 1, function(x) mean(x >= 0))
  
  res <- list(
    post_f = post_f[-(1:nburn), ],
    post_tau = post_tau[-(1:nburn), ],
    post_gamma = post_gamma[-(1:nburn), ],
    post_beta = post_beta[-(1:nburn), ],
    post_f_i = post_f_i[-(1:nburn), ],
    
    post_sigma_e = post_sigma_e[-(1:nburn), ],
    post_sigma_delta = post_sigma_delta[-(1:nburn)],
    post_sigma_delta_f_all = post_sigma_delta_f_all[-(1:nburn), ],
    post_sigma_phi_f_all = post_sigma_phi_f_all[-(1:nburn), ],
    post_sigma_delta_tau = post_sigma_delta_tau[-(1:nburn)],
    post_sigma_phi_tau = post_sigma_phi_tau[-(1:nburn)],
    post_sigma_phi_gamma = post_sigma_phi_gamma[-(1:nburn)],
    post_sigma_delta_gamma = post_sigma_delta_gamma[-(1:nburn)],
    post_ATE = rowMeans(post_ATE),
    Prob = Prob
  )
  
  if(if_model_diagnostic){
    
    pdf(paste0(plot_path, "_ATE.pdf"), height = 3, width = 6)
    plot(as.mcmc(res$post_ATE))
    dev.off()
    
    if(if_site){
      pdf(paste0(plot_path, "_Offsite_Bias_Shift.pdf"), height = 3, width = 6)
      if(is.null(ncol(res$post_beta))){
        plot(as.mcmc(res$post_beta))
      }else{
        plot(as.mcmc(res$post_beta[, 1]))
      }
      dev.off()
    }
    
    if(if_all_onsite){
      pdf(paste0(plot_path, "_error.pdf"), height = 3, width = 6)
      plot(as.mcmc(res$post_sigma_e[, 2]))
      dev.off()
    }else{
      pdf(paste0(plot_path, "_error.pdf"), height = 6, width = 6)
      plot(as.mcmc(res$post_sigma_e))
      dev.off()
    }
  }
  
  if(if_plot | if_model_diagnostic){
    
    t_new <- seq(0, 1, 0.01)
    Q_new <- array(0, c(length(t_new), G))
    
    for(k in 0:K){
      Q_new[, k + 1] <- t_new^k
    }
    
    for(l in 1:L){
      Q_new[, K + 1 + l] <- ifelse(t_new <= kappa[l], (t_new - kappa[l])^d, 0)
    }
    
    ## Trajectory of control from onsite
    f0_fitted_onsite <- res$post_f %*% t(Q_new)
    f0_ci_onsite <- HPDinterval(as.mcmc(f0_fitted_onsite))
    
    ## Trajectory of control from offsite
    f0_fitted_offsite <- res$post_f %*% t(Q_new) + 
      replicate(length(t_new), ifelse(is.null(ncol(res$post_beta)), res$post_beta, res$post_beta[, 1]))
    f0_ci_offsite <- HPDinterval(as.mcmc(f0_fitted_offsite))
    
    ## Trajectory of treatment effect onsite
    tau_fitted <- res$post_tau %*% t(Q_new[, -1])
    tau_ci <- HPDinterval(as.mcmc(tau_fitted))
    
    ## Trajectory of treatment from onsite
    f1_fitted_onsite <- tau_fitted + f0_fitted_onsite
    f1_ci_onsite <- HPDinterval(as.mcmc(f1_fitted_onsite))
    
    ## Trajectory of treatment from offsite
    if(if_site){
      gamma_fitted <- res$post_gamma %*% t(Q_new[, -1])
      gamma_ci <- HPDinterval(as.mcmc(gamma_fitted))
      f1_fitted_offsite <- tau_fitted + gamma_fitted + f0_fitted_offsite
    }else{
      f1_fitted_offsite <- tau_fitted + f0_fitted_offsite
    }
    
    f1_ci_offsite <- HPDinterval(as.mcmc(f1_fitted_offsite))
    
    if(if_model_diagnostic) pdf(paste0(plot_path, "_fitted_curves.pdf"), height = 8, width = 8)
    par(mfrow = c(2, 2))
    
    ## Trajectory of control onsite
    plot(t_new, colMeans(f0_fitted_onsite), type = "n", ylim = range(Y),
         xlab = "Time", ylab = "Outcome", main = "Onsite")
    lines(t_new, f0_ci_onsite[, 1], lty = 2)
    lines(t_new, f0_ci_onsite[, 2], lty = 2)
    polygon(c(t_new, rev(t_new)),
            c(f0_ci_onsite[, 2], rev(f0_ci_onsite[, 1])),
            col = "grey", border = NA)
    lines(t_new, colMeans(f0_fitted_onsite), lwd = 3, col = "orange")
    
    ## Trajectory of treatment onsite
    lines(t_new, f1_ci_onsite[, 1], lty = 2)
    lines(t_new, f1_ci_onsite[, 2], lty = 2)
    polygon(c(t_new, rev(t_new)),
            c(f1_ci_onsite[, 2], rev(f1_ci_onsite[, 1])),
            col = "grey", border = NA)
    lines(t_new, colMeans(f1_fitted_onsite), lwd = 3, col = "light blue")
    
    lines(t_all[X_all == 0], Y[X_all == 0], type = "p", cex = 0.5, col = "orange")
    lines(t_all[X_all == 1], Y[X_all == 1], type = "p", cex = 0.5, col = "light blue")
    
    for(t in unique(t_all)) lines(t, mean(data$Y[t_all == t & X_all == 1]),
                                  type = "p", pch = 4, cex = 1.5, col = "light blue", lwd = 2)
    for(t in unique(t_all)) lines(t, mean(data$Y[t_all == t & X_all == 0]),
                                  type = "p", pch = 3, cex = 1.5, col = "orange", lwd = 2)
    
    ## Trajectory of control offsite
    plot(t_new, colMeans(f0_fitted_offsite), type = "n", ylim = range(Y),
         xlab = "Time", ylab = "Outcome", main = "Offsite")
    lines(t_new, f0_ci_offsite[, 1], lty = 2)
    lines(t_new, f0_ci_offsite[, 2], lty = 2)
    polygon(c(t_new, rev(t_new)),
            c(f0_ci_offsite[, 2], rev(f0_ci_offsite[, 1])),
            col = "grey", border = NA)
    lines(t_new, colMeans(f0_fitted_offsite), lwd = 3, col = "orange")
    
    ## Trajectory of treatment offsite
    lines(t_new, f1_ci_offsite[, 1], lty = 2)
    lines(t_new, f1_ci_offsite[, 2], lty = 2)
    polygon(c(t_new, rev(t_new)),
            c(f1_ci_offsite[, 2], rev(f1_ci_offsite[, 1])),
            col = "grey", border = NA)
    lines(t_new, colMeans(f1_fitted_offsite), lwd = 3, col = "light blue")
    
    lines(t_all[X_all == 0], Y[X_all == 0], type = "p", cex = 0.5, col = "orange")
    lines(t_all[X_all == 1], Y[X_all == 1], type = "p", cex = 0.5, col = "light blue")
    
    for(t in unique(t_all)) lines(t, mean(data$Y[t_all == t & X_all == 1]),
                                  type = "p", pch = 4, cex = 1.5, col = "light blue", lwd = 2)
    for(t in unique(t_all)) lines(t, mean(data$Y[t_all == t & X_all == 0]),
                                  type = "p", pch = 3, cex = 1.5, col = "orange", lwd = 2)
    
    ## Trajectory of treatment effect
    plot(t_new, colMeans(tau_fitted), type = "n", ylim = range(tau_ci),
         xlab = "Time", ylab = "Treatment effect")
    lines(t_new, tau_ci[, 1], lty = 2)
    lines(t_new, tau_ci[, 2], lty = 2)
    polygon(c(t_new, rev(t_new)),
            c(tau_ci[, 2], rev(tau_ci[, 1])),
            col = "grey", border = NA)
    lines(t_new, colMeans(tau_fitted), lwd = 3, col = "black")
    lines(t_new, Delta(t_new), lwd = 3, lty = 2, col = "red")
    
    ## Trajectory of treatment-specific offsite bias
    if(if_site){
      plot(t_new, colMeans(gamma_fitted), type = "n", ylim = range(gamma_ci),
           xlab = "Time", ylab = "Treatment-offsite bias")
      lines(t_new, gamma_ci[, 1], lty = 2)
      lines(t_new, gamma_ci[, 2], lty = 2)
      polygon(c(t_new, rev(t_new)),
              c(gamma_ci[, 2], rev(gamma_ci[, 1])),
              col = "grey", border = NA)
      lines(t_new, colMeans(gamma_fitted), lwd = 3, col = "black")
      lines(t_new, Gamma_bias(t_new), lwd = 3, lty = 2, col = "red")
    }
    
    par(mfrow = c(1, 1))
    if(if_model_diagnostic) dev.off()
  }
  
  return(res)
}




