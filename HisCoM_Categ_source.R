mu_prob <- function(cumprob, nobs, ncat1) {
  mat <- matrix(cumprob, nobs, ncat1, TRUE)
  mat <- rbind(mat[, 1], diff(t(mat)), 1 - mat[, ncat1])
  mat <- c(mat)
  return(mat)
}

norm    <- function(u) sqrt(sum(u^2)) 
HisCoM_Categ_POM <- function(y, path, path_var, indx, data, maxiter, lambda1, lambda2, tol){
  #path: list of pathway
  #path_var: variable list in all Pathways
  ##################################################
  #########Response Variable (Phenotype)############
  ##################################################
  y <- as.numeric(factor(y))
  nobs <- length(y)
  ncat <- nlevels(factor(y))
  ncat1 <- ncat-1
  Y <- rep(y, each = ncat1)
  Intercept <- rep.int(seq(ncat1), length(y))
  y_mat <- as.numeric(Y == Intercept)
  ncase <- length(y_mat)
  ###
  X.all <- data[,match(path_var, colnames(data))]
  xnames <- colnames(X.all)
  X_mat_1 <- apply(X.all, 2, function(co) rep(co, each = ncat1))
  X_mat_1 <- matrix(X_mat_1, ncol = ncol(X_mat_1), dimnames = NULL)
  X_mat_1 <- scale(X_mat_1)*sqrt(nobs/(nobs-1))
  #diag(t(X_mat_1)%*%X_mat_1)
  X_mat_0 <- model.matrix(~factor(Intercept)-1 )
  X_mat_0 <- matrix(X_mat_0, ncol = ncol(X_mat_0), dimnames = NULL)
  
  X_mat <- cbind(X_mat_0, X_mat_1)
  
  ##############
  
  nvar <- c()    ############How many Variables Per Group 
  for (i in 1:length(unique(path))){
    nvar <- c(nvar, sum(path==unique(path)[i]))
  }
  ndset <- length(nvar)   ####Total Number of Group 
  sum_nvar <- sum(nvar)    #######Total number of metabolites in all pathways
  W1 = matrix(0, sum_nvar,ndset)
  kk = 0
  for (j in 1:ndset) {
    Nj            = nvar[j]
    k            = kk + 1
    kk            = kk + Nj
    W1[k:kk,j]        = 99 * ones(Nj, 1) ##library(pracma)
  }
  
  windex        = which(W1 == 99)                     ### For w* vector    #w_star_99
  num_windex <- length(windex)                       #w_star_99
  W <- W1
  W[windex] <- runif(num_windex) #rand(num_windex,1)
  W_new <- as.numeric(W[windex])
  ###
  I_mat <- diag(ncat1)
  W2 <- adiag(I_mat, W1)  ## library(magic)
  W2 <- as.matrix(W2)
  w_star_99        = which(W2 == 99)
  w_Kro_idx        = which(t(W2) == 99)
  F_mat <- X_mat %*% W2
  
  beta_new <- ginv(t(F_mat)%*% F_mat) %*% t(F_mat) %*% y_mat
  
  
  est_new <- c(W_new, beta_new)
  
  converge <- F
  iter <- 0
  
  family <- make.link("logit") 
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  exclude <- seq(ncat, nobs * ncat, ncat)
  ncat1 <- ncat-1
  
  Nsubs <- rep(1:nobs, each=ncat1)
  nbeta <- length(beta_new)
  
  ans <- diag(ncat1)
  ans[seq(2, ncat1^2, ncat1 + 1)] <- -1
  
  while(iter < maxiter){
    W_old            <- W_new
    W2[w_star_99]    <- W_old
    beta_old         <- beta_new
    est_old          <- c(W_old, beta_old)
    F_mat            <- X_mat %*% W2
    eta              <- drop(F_mat %*% beta_old)   
    
    fitproball       <- mu_prob(linkinv(eta), nobs, ncat1)  ###linkinv(eta) produce cumulative probability
    fitprob          <- fitproball[-exclude]
    dummy            <- mu.eta(eta)
    pi               <- fitprob
    resids           <- y_mat - pi
    
    ###### Update W
    kk = ncat1
    for (j in 1:ndset) {
      Nj              = nvar[j]
      k               = kk + 1
      kk              = kk + Nj
      X_j             = X_mat[,k:kk]
      B_mat <- as.matrix(X_mat[,k:kk])*beta_old[(ncat1+j),]
      beta_old1 <- beta_old
      beta_old1[(ncat1+j)] <- 0
      z2 <- F_mat %*% beta_old1
      ##
      Sw_mat <- matrix(0, Nj, 1, FALSE)         #gradient:S
      Hw_mat <- matrix(0, Nj, Nj, FALSE)  ##Second derivative
      for(i in 1:nobs){
        selector         <- Nsubs == i
        mueta_i1         <- dummy[selector]
        mueta_i          <- apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = ncat1))
        dpi_eta          <- ans*mueta_i
        B_mat_i          <- B_mat[selector,]
        Pi_vct           <- pi[selector]
        V_mat            <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
        V_inv            <- ginv(V_mat)  #inversematrix
        res_i            <- resids[selector]
        G_i              <- t(dpi_eta) %*% V_inv %*% dpi_eta
        Z_i              <- eta[selector] + (ginv(dpi_eta)%*%res_i)
        Z_1i             <- Z_i - z2[selector]
        Sw_300           <- (t(B_mat_i) %*% G_i %*% Z_1i)   
        Sw_mat           <- Sw_mat + Sw_300
        Hw_300           <- (t(B_mat_i) %*% G_i %*% B_mat_i)
        Hw_mat           <- Hw_mat + Hw_300
      }
      ##
      w_j                <-   ginv(Hw_mat + lambda1*diag(Nj)) %*% Sw_mat 
      w_j                <-   sqrt(nobs)*w_j/ norm(X_j%*%w_j) 
      W2[k:kk,(ncat1+j)] <- w_j
      F_mat[, (ncat1+j)] <- X_j %*% w_j
    }
    W_new <- W2[w_star_99]
    
    ###beta_update
    
    Sb_mat <-  matrix(0, nbeta, 1, FALSE)    #gradient:S
    Hb_mat <-  matrix(0, nbeta, nbeta, FALSE)     ##Second derivative
    for(i in 1:nobs){
      selector <- Nsubs == i
      mueta_i1 <- dummy[selector]
      mueta_i <- apply(t(as.matrix(mueta_i1)), 2, function(x) rep(x, each = ncat1))
      dpi_eta <- ans*mueta_i
      F_mat_i <- F_mat[selector,]
      Pi_vct <- pi[selector]
      V_mat <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
      V_inv <- ginv(V_mat)
      res_i <- resids[selector]
      G_i <- t(dpi_eta) %*% V_inv %*% dpi_eta
      Z_i <- eta[selector] + (ginv(dpi_eta)%*%res_i)   
      Sb_300 <- (t(F_mat_i) %*% G_i %*% Z_i)  
      Sb_mat <- Sb_mat + Sb_300
      Hb_300 <-  (t(F_mat_i) %*% G_i %*% F_mat_i) 
      Hb_mat <- Hb_mat + Hb_300
    }
    
    p_beta <- rep(lambda2,length(beta_old)) #lambda2*diag(length(beta_old))
    p_beta[c(1:ncat1, indx)] <- 0
    beta_new <- ginv(Hb_mat + diag(p_beta)) %*% Sb_mat
    
    est_new <- c(W_new, beta_new)
    crit <- sum(abs(est_new - est_old))
    iter <- iter + 1
    if(iter%%50==0){
      cat("iter = ", iter, " | diff = ", crit, "\n")
    }
    if (crit <= tol) {
      break
    }
  }
  beta_coef <- beta_new
  weight_coef <- W_new
  W2[w_star_99] <- weight_coef
  
  ## deviance calculate
  #eta_F <- drop(X_mat %*% W2 %*% beta_coef)
  #fitproball <- mu_prob(linkinv(eta_F), nobs, ncat1)  
  fitprob <- fitproball
  Y1 <- rep(y, each = ncat)
  Intercept <- rep.int(seq(ncat), length(y))
  Y1_mat <- as.numeric(Y1 == Intercept)
  like_Li <- sum(Y1_mat * log(fitprob))
  dev <- 2*sum(Y1_mat*log((Y1_mat/fitprob) + 0.00000001)) 
  ##
  var_naem <- c("beta01", "beta02", unique(path))
  beta_coef <- data.frame(var_naem, as.numeric(beta_coef))
  colnames(beta_coef) <- c("pathway", "coefficients")
  ##
  fit <- list()
  fit$xnames <- xnames
  fit$nvar <- nvar
  fit$beta_coef <- beta_coef
  fit$weight_coef <- weight_coef
  #fit$sd <- sdP
  fit$W <- W2
  fit$LikeLi <- like_Li
  fit$deviance <- dev
  fit$crit <- crit
  fit$Fm <- F_mat
  fit$Xm <- X_mat
  fit
}

#########################################HisCoM-Categ for Adjucent Category logit model##################
#####################
HisCoM_Categ_ACL <- function(y, path, path_var, indx, data, maxiter, lambda1, lambda2, tol){
  #path: list of pathway
  #path_var: variable list in all Pathways
  ##################################################
  #########Response Variable (Phenotype)############
  ##################################################
  y <- as.numeric(factor(y))
  nobs <- length(y)
  ncat <- nlevels(factor(y))
  ncat1 <- ncat-1
  Y <- rep(y, each = ncat1)
  Intercept <- rep.int(seq(ncat1), length(y))
  y_mat <- as.numeric(Y == Intercept)
  ncase <- length(y_mat)
  ###
  X.all <- data[,match(path_var, colnames(data))]
  xnames <- colnames(X.all)
  X_mat_1 <- apply(X.all, 2, function(co) rep(co, each = ncat1))
  X_mat_1 <- matrix(X_mat_1, ncol = ncol(X_mat_1), dimnames = NULL)
  X_mat_1 <- scale(X_mat_1)*sqrt(ncase/(ncase-1))
  #diag(t(X_mat_1)%*%X_mat_1)
  X_mat_0 <- model.matrix(~factor(Intercept)-1 )
  X_mat_0 <- matrix(X_mat_0, ncol = ncol(X_mat_0), dimnames = NULL)
  
  X_mat <- cbind(X_mat_0, X_mat_1)
  
  #####
  dummy <- ncat - 1
  dummy.matrix <- diag(rep.int(1, dummy))
  dummy.matrix[upper.tri(dummy.matrix)] <- 1
  X_mat[, 1:dummy] <- kronecker(rep.int(1, nrow(X_mat) / dummy), dummy.matrix)
  if (dummy != ncol(X_mat)) {
    X_mat[, -c(1:dummy)] <- X_mat[, -c(1:dummy)] * 
      rep(dummy:1,nrow(X_mat) / dummy )
  }
  
  ##############
  nvar <- c()    ############How many Variables Per Group 
  for (i in 1:length(unique(path))){
    nvar <- c(nvar, sum(path==unique(path)[i]))
  }
  ndset <- length(nvar)   ####Total Number of Group 
  sum_nvar <- sum(nvar)    #######Total number of metabolites in all pathways
  W1 = matrix(0, sum_nvar,ndset)
  kk = 0
  for (j in 1:ndset) {
    Nj            = nvar[j]
    k            = kk + 1
    kk            = kk + Nj
    W1[k:kk,j]        = 99 * ones(Nj, 1) ##library(pracma)
  }
  
  windex        = which(W1 == 99)                     ### For w* vector    #w_star_99
  num_windex <- length(windex)                       #w_star_99
  W <- W1
  W[windex] <- runif(num_windex) #rand(num_windex,1)
  W_new <- as.numeric(W[windex])
  ###
  I_mat <- diag(ncat1)
  W2 <- adiag(I_mat, W1)  ## library(magic)
  W2 <- as.matrix(W2)
  w_star_99        = which(W2 == 99)
  w_Kro_idx        = which(t(W2) == 99)
  F_mat <- X_mat %*% W2
  
  beta_new <- ginv(t(F_mat)%*% F_mat) %*% t(F_mat) %*% y_mat
  
  
  est_new <- c(W_new, beta_new)
  
  converge <- F
  iter <- 0
  
  ncat1 <- ncat-1
  
  Nsubs <- rep(1:nobs, each=ncat1)
  nbeta <- length(beta_new)
  
  
  while(iter < maxiter){
    W_old            <- W_new
    W2[w_star_99]    <- W_old
    beta_old         <- beta_new
    wb_old <- W2 %*% beta_old
    est_old          <- c(W_old, beta_old)
    F_mat            <- X_mat %*% W2
    F_mat2 <- F_mat[,-c(1:ncat1)]
    F_mat            <- cbind( F_mat[,c(1:ncat1)], scale(F_mat2)*sqrt(ncase/(ncase-1)))
    eta              <- drop(F_mat %*% beta_old)   
    
    fitprob <- exp(matrix(eta, length(eta) / ncat1, ncat1, TRUE))
    fitprob <- fitprob / (1 + .rowSums(fitprob, nobs, ncat1, FALSE))
    fitproball <- as.vector(t(cbind(fitprob, 1 - .rowSums(fitprob,nobs, ncat1, FALSE))))
    fitprob <- as.vector(t(fitprob))
    dummy            <- fitprob 
    pi               <- fitprob
    resids           <- y_mat - pi
    
    ###### Update W
    kk = ncat1
    for (j in 1:ndset) {
      Nj              = nvar[j]
      k               = kk + 1
      kk              = kk + Nj
      X_j             = X_mat[,k:kk]
      B_mat <- as.matrix(X_mat[,k:kk])*beta_old[(ncat1+j),]
      beta_old1 <- beta_old
      beta_old1[(ncat1+j)] <- 0
      z2 <- F_mat %*% beta_old1
      ##
      Sw_mat <- matrix(0, Nj, 1, FALSE)         #gradient:S
      Hw_mat <- matrix(0, Nj, Nj, FALSE)  ##Second derivative
      for(i in 1:nobs){
        selector         <- Nsubs == i
        B_mat_i          <- B_mat[selector,]
        Pi_vct           <- pi[selector]
        V_mat            <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
        V_inv            <- ginv(V_mat)  #inversematrix
        res_i            <- resids[selector]
        Z_i              <- eta[selector] + (V_inv %*% res_i)
        Z_1i             <- Z_i - z2[selector]
        Sw_300           <- (t(B_mat_i) %*% V_mat  %*% Z_1i)   
        Sw_mat           <- Sw_mat + Sw_300
        Hw_300           <- (t(B_mat_i) %*% V_mat  %*% B_mat_i)
        Hw_mat           <- Hw_mat + Hw_300
      }
      ##
      w_j                <-   ginv(Hw_mat + lambda1*diag(Nj)) %*% Sw_mat 
      w_j                <-   sqrt(nobs)*w_j/ norm(X_j%*%w_j)  #sqrt(ncase)*
      W2[k:kk,(ncat1+j)] <- w_j
      F_mat[, (ncat1+j)] <- X_j %*% w_j
    }
    W_new <- W2[w_star_99]
    
    ###beta_update
    
    Sb_mat <-  matrix(0, nbeta, 1, FALSE)    #gradient:S
    Hb_mat <-  matrix(0, nbeta, nbeta, FALSE)     ##Second derivative
    for(i in 1:nobs){
      selector <- Nsubs == i
      F_mat_i <- F_mat[selector,]
      Pi_vct <- pi[selector]
      V_mat <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
      V_inv <- ginv(V_mat)
      res_i <- resids[selector]
      Z_i <- eta[selector] + (V_inv %*% res_i)   
      Sb_300 <- (t(F_mat_i) %*% V_mat  %*%  Z_i)  
      Sb_mat <- Sb_mat + Sb_300
      Hb_300 <-  (t(F_mat_i) %*% V_mat  %*% F_mat_i) 
      Hb_mat <- Hb_mat + Hb_300
    }
    
    p_beta <- rep(lambda2,length(beta_old)) #lambda2*diag(length(beta_old))
    p_beta[c(1:ncat1, indx)] <- 0
    beta_new <- ginv(Hb_mat + diag(p_beta)) %*% Sb_mat
    
    wb_new <- W2 %*% beta_new
    
    est_new <- c(W_new, beta_new)
    #crit <- sum(abs(est_new - est_old))
    crit <- sum(abs(wb_new - wb_old))
    iter <- iter + 1
    if(iter%%50==0){
      cat("iter = ", iter, " | diff = ", crit, "\n")
    }
    if (crit <= tol) {
      break
    }
  }
  beta_coef <- beta_new
  weight_coef <- W_new
  W2[w_star_99] <- weight_coef
  ##
  fitprob <- fitproball
  Y1 <- rep(y, each = ncat)
  Intercept <- rep.int(seq(ncat), length(y))
  Y1_mat <- as.numeric(Y1 == Intercept)
  like_Li <- sum(Y1_mat * log(fitprob))
  dev <- 2*sum(Y1_mat*log((Y1_mat/fitprob) + 0.00000001)) 
  #
  var_naem <- c("beta01", "beta02", unique(path))
  beta_coef <- data.frame(var_naem, as.numeric(beta_coef))
  colnames(beta_coef) <- c("pathway", "coefficients")
  ##
  fit <- list()
  fit$xnames <- xnames
  fit$nvar <- nvar
  fit$beta_coef <- beta_coef
  fit$weight_coef <- weight_coef
  fit$W <- W2
  fit$LikeLi <- like_Li
  fit$deviance <- dev
  fit$crit <- crit
  fit$Fm <- F_mat
  fit$Xm <- X_mat
  fit
}
#########################################HisCoM-Categ for Baseline Category logit model##################
#####################
HisCoM_Categ_BCL <- function(y, path, path_var, indx, data, maxiter, lambda1, lambda2, tol){
  #path: list of pathway
  #path_var: variable list in all Pathways
  ##################################################
  #########Response Variable (Phenotype)############
  ##################################################
  y <- as.numeric(factor(y))
  nobs <- length(y)
  ncat <- nlevels(factor(y))
  ncat1 <- ncat-1
  Y <- rep(y, each = ncat1)
  Intercept <- rep.int(seq(ncat1), length(y))
  y_mat <- as.numeric(Y == Intercept)
  ncase <- length(y_mat)
  
  ###############################
  ###
  X.all <- data[,match(path_var, colnames(data))]
  Xinit_mat <- apply(X.all, 2, function(co) rep(co, each = ncat1))
  Xinit_mat <- matrix(Xinit_mat, ncol = ncol(Xinit_mat), dimnames = NULL)
  Xinit_mat <- scale(Xinit_mat)*sqrt(ncase/(ncase-1)) #nobs
  #diag(t(X_mat_1)%*%X_mat_1)
  X_mat <- model.matrix(~factor(Intercept)-1 )
  ######
  if (ncol(Xinit_mat) > 1) {
    Xinit_mat <- cbind(X_mat, Xinit_mat)
  } else {
    Xinit_mat <- X_mat
  }
  if (ncol(Xinit_mat) != (ncat - 1)) {
    X_inter <- X_mat
    for (i in ncat:ncol(Xinit_mat)) X_mat <- cbind(X_mat, X_inter * Xinit_mat[, i])
  }
  
  X_mat <- matrix(X_mat, ncol = ncol(X_mat), dimnames = NULL)
  X_mat <- X_mat[, c(matrix(seq(ncol(X_mat)), ncol = ncat - 1, byrow = TRUE))]
  #diag(t(X_mat)%*%X_mat)
  ##########  
  ##############
  
  nvar <- c()    ############How many Variables Per Group 
  for (i in 1:length(unique(path))){
    nvar <- c(nvar, sum(path==unique(path)[i]))
  }
  ndset <- length(nvar)   ####Total Number of Group 
  sum_nvar <- sum(nvar)    #######Total number of metabolites in all pathways
  W1 = matrix(0, sum_nvar,ndset)
  kk = 0
  for (j in 1:ndset) {
    Nj            = nvar[j]
    k            = kk + 1
    kk            = kk + Nj
    W1[k:kk,j]        = 99 * ones(Nj, 1) ##library(pracma)
  }
  W <- adiag(1, W1, 1,W1)
  windex        = which(W == 99)
  w_Kro_idx        = which(t(W) == 99)
  num_windex <- length(windex)                       #w_star_99
  #set.seed(2)
  W[windex] <- runif(num_windex) #rand(num_windex,1)
  W_new <- as.numeric(W[windex])
  ###
  

  
  
  
  #####
  F_mat <- X_mat %*% W
  beta_new <- ginv(t(F_mat)%*% F_mat) %*% t(F_mat) %*% y_mat
  est_new <- c(W_new, beta_new)
  
  converge <- F
  iter <- 0
  wb_new <- W %*% beta_new
  
  
  Nsubs <- rep(1:nobs, each=ncat1)
  nbeta <- length(beta_new)
  err_rate=0.0001
  
  while(iter < maxiter){
    W_old            <- W_new
    W[windex]       <- W_old
    wb_old          <- wb_new
    beta_old         <- beta_new
    est_old          <- c(W_old, beta_old)
    F_mat            <- X_mat %*% W
    eta              <- drop(F_mat %*% beta_old)
    ##
    fitprob <- exp(matrix(eta, length(eta) / ncat1, ncat1, TRUE))
    fitprob <- fitprob / (1 + .rowSums(fitprob, nobs, ncat1, FALSE))
    fitproball <- as.vector(t(cbind(fitprob, 1 - .rowSums(fitprob,nobs, ncat1, FALSE))))
    fitprob <- as.vector(t(fitprob))
    dummy <- fitprob
    
    pi <- as.vector(fitprob)  
    resids <- (y_mat - pi)
    
    A_mat <- kronecker(X_mat, t(beta_old))[,w_Kro_idx]
    ###### Update W
    kk = 0
    kk1 <- sum_nvar
    for (j in 1:ndset) {
      Nj              = nvar[j]
      ###
      k               = kk + 1
      kk              = kk + Nj
      k1               = kk1 + 1
      kk1              = kk1 + Nj
      ##
      A_j             = A_mat[,c(k:kk, k1:kk1)]
      X_j             = X_mat[,c((k:kk)+1, (k1:kk1)+1)]
      beta_old1 <- beta_old
      beta_old1[c(1+j, ndset +2 +j)] <- 0
      z2 <- F_mat %*% beta_old1
      ##
      Sw_mat <- matrix(0, 2*Nj, 1, FALSE)         #gradient:S
      Hw_mat <- matrix(0, 2*Nj, 2*Nj, FALSE)  ##Second derivative
      for(i in 1:nobs){
        selector         <- Nsubs == i
        ##
        A_mat_i          <- A_j[selector,]
        Pi_vct           <- pi[selector]
        V_mat            <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
        res_i            <- resids[selector]
        Z_i              <- eta[selector] + (ginv(V_mat) %*% res_i)
        Z_1i             <- Z_i - z2[selector]
        Sw_300           <- (t(A_mat_i) %*% V_mat %*% Z_1i)   
        Sw_mat           <- Sw_mat + Sw_300
        Hw_300           <- (t(A_mat_i) %*% V_mat %*% A_mat_i)
        Hw_mat           <- Hw_mat + Hw_300
      }
      ##
      w_j                <-   ginv(Hw_mat + lambda1*diag(2*Nj)) %*% Sw_mat 
      rrf                <-   X_j %*% t(adiag(t(w_j[1:Nj]), t(w_j[(Nj+1):(2*Nj)])))
      ###
      w11 <- w_j[1:Nj]*sqrt(nobs)/sqrt(sum(as.numeric(rrf[,1])**2))
      w22 <- w_j[(Nj+1):(2*Nj)]*sqrt(nobs)/sqrt(sum(as.numeric(rrf[,2])**2))
      ff1 <- (X_mat[,((k:kk))+1] %*%  t(t(w11)))
      sum(ff1**2)
      ff2 <- (X_mat[,(k1:kk1)+1] %*%  t(t(w22)))
      sum(ff2**2)
      W[(k:kk)+1,(1+j)] <- w11   
      W[(k1:kk1)+2,(ndset+2+j)] <-  w22 
      
      F_mat[, (1+j)] <- ff1 
      F_mat[, (ndset +2 +j)] <- ff2 
    }
    
    W_new <- W[windex]

    ###beta_update
    
    Sb_mat <-  matrix(0, nbeta, 1, FALSE)    #gradient:S
    Hb_mat <-  matrix(0, nbeta, nbeta, FALSE)     ##Second derivative
    for(i in 1:nobs){
      selector <- Nsubs == i
      F_mat_i <- F_mat[selector,]
      Pi_vct <- pi[selector]
      V_mat <- diag(Pi_vct) - Pi_vct %*% t(Pi_vct)
      V_inv <- ginv(V_mat)
      res_i <- resids[selector]
      Z_i <- eta[selector] + (V_inv %*% res_i)   
      Sb_300 <- (t(F_mat_i) %*% V_mat %*% Z_i)  
      Sb_mat <- Sb_mat + Sb_300
      Hb_300 <-  (t(F_mat_i) %*% V_mat %*% F_mat_i) 
      Hb_mat <- Hb_mat + Hb_300
    }
    
    p_beta <- rep(lambda2,length(beta_old)) 
    p_beta[c(1,indx+1, ndset+2,indx+ndset+2)] <- 0
    beta_new <- ginv(Hb_mat + diag(p_beta)) %*% Sb_mat
    
    est_new <- c(W_new, beta_new)

    wb_new <- W %*% beta_new
    crit <- max(abs(abs(as.numeric(wb_new) - as.numeric(wb_old))))
    
    iter <- iter + 1
    if(iter%%50==0){
      cat("iter = ", iter, " | diff = ", crit, "\n")
    }
    if (crit <= tol) {
      break
    }
  }
  
  
  
  beta_coef <- beta_new
  weight_coef <- W_new
  W[windex] <- weight_coef
  ####Log-likelihood & deviance####
  fitprob <- fitproball
  Y1 <- rep(y, each = ncat)
  Intercept <- rep.int(seq(ncat), length(y))
  Y1_mat <- as.numeric(Y1 == Intercept)
  like_Li <- sum(Y1_mat * log(fitprob))
  dev <- 2*sum(Y1_mat*log((Y1_mat/fitprob) + 0.00000001))
  
  ##
  var_naem <- c("beta01",  unique(path), "beta02",unique(path))
  beta_coef <- data.frame(var_naem, as.numeric(beta_coef))
  colnames(beta_coef) <- c("pathway", "coefficients")
  
  ##
  fit <- list()
  fit$nvar <- nvar
  fit$beta_coef <- beta_coef
  fit$weight_coef <- weight_coef
  fit$W <- W
  fit$LikeLi <- like_Li
  fit$deviance <- dev
  fit$crit <- crit
  fit$Fm <- F_mat
  fit$Xm <- X_mat
  fit
}
