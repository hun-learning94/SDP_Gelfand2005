SDPlearn = function(Y.obs, Xs, coords.train, 
                    phigrid_d = 0.1, 
                    B = 2e3, burnin = 5e2)
{
  Tn = ncol(Y.obs) # num of data replications
  N = nrow(Y.obs) # size of each replication (# of sites)
  p = dim(Xs)[2]
  #### priors
  library(Rcpp)
  library(RcppArmadillo)
  sourceCpp("rdist_cpp.cpp")
  ## Beta, tau for data model
  ## 1. Beta : g - prior
  ## beta ~ N(0, n * tau^2 * solve(nrow(Y) * crossprod(X)))
  ## 2. tau
  a_tau = 2; b_tau = 2
  ## nu for DP
  a_nu = 3; b_nu = 0.005
  ## G0 for DP
  ## 1. phi
  d = rdist_cpp(coords.train, base::crossprod)
  maxd = max(d)
  b_phi = 3/(phigrid_d * maxd); ngrid = 200
  phi_grid = seq(b_phi/ngrid, b_phi, length=ngrid)
  ## 2. sigma2
  a_sig = 2; b_sig = 2
  
  #### support functions
  ## for step 1. update theta
  H = function(d, phi) exp(-phi * d)
  library(mvtnorm) # for dmvnorm
  
  calc_Lambda = function(tau, sigma, Hphi_inv){
    n = nrow(Hphi)
    return(solve(diag(n)/tau^2 + Hphi_inv/sigma^2))
  }
  
  # initial NA imputation
  NA_IDX = matrix(F, nrow=N, ncol=Tn)
  for(t_ in 1:Tn) NA_IDX[,t_] = is.na(Y.obs[,t_])
  Y.obs.imp = Y.obs
  for(t_ in 1:Tn) Y.obs.imp[NA_IDX[,t_],t_] = 0
  
  #### placeholders
  # B = 2e3 # num of iterations
  # burnin = 5e2
  
  Y.comp.samp = array(rep(matrix(nrow = N, ncol = Tn), B),
                      dim = c(N, Tn, B))
  Y.comp.samp[,,1] = Y.obs.imp
  THETA.samp = array(rep(matrix(nrow = N, ncol = Tn), B),
                     dim = c(N, Tn, B)) # N X Tn X B 3D array of sampled latent processes
  Ws.samp = vector("list", length = B) # list of length B, cluster assignments
  BETA.samp = matrix(nrow = p, ncol = B)
  TAU.samp = rep(NA_real_, length = B)
  NU.samp = rep(NA_real_, length = B)
  PHI.samp = rep(NA_real_, length = B)
  SIG.samp = rep(NA_real_, length = B)
  
  #### initial values
  THETA.samp[,,1] =  matrix(Y.obs.imp[,1], nrow = N, ncol = Tn)
  BETA.samp[, 1] = solve(crossprod(X), t(X) %*% Y.obs.imp[,1])
  TAU.samp[1] = 1
  NU.samp[1] = 1
  PHI.samp[1] = 1
  SIG.samp[1] = 1
  Tjs = c(nrow(Y.obs))
  Ws = rep(1, length = Tn)
  Ws.samp[[1]] = Ws
  Ws_count = function(Ws){
    res = table(Ws)
    attributes(res) = NULL
    return(res)
  }
  Tjs = Ws_count(Ws)
  Tjstars = matrix(THETA.samp[,1,1], nrow = N)
  sumXtX = matrix(rowSums(matrix(apply(Xs,3, crossprod), ncol = dim(Xs)[3])), p, p)
  sumXtX_inv = solve(sumXtX)
  
  # site specific centering of Y
  # site.mean = rowMeans(Y.obs, na.rm=T)
  # Y.obs.raw = Y.obs
  # Y.obs = sweep(Y.obs ,1,site.mean, "-")
  
  for(b in 2:B){
    if((b %% 10)==0){
      cat(b, "\n")
      if(b > 100){
        traceidx = 100:b
        opar = par(no.readonly = T)
        par(mfrow = c(2, 4))
        par(mar = c(2.5, 3, 2, 1), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0))
        plot(BETA.samp[1, traceidx], type = "l")
        plot(BETA.samp[2, traceidx], type = "l")
        plot(BETA.samp[3, traceidx], type = "l")
        plot(TAU.samp[traceidx], type = "l")
        plot(NU.samp[traceidx], type = "l")
        plot(PHI.samp[traceidx], type = "l")
        plot(SIG.samp[traceidx], type = "l")
        plot(1:3)
        par(opar)
      }
    } 
    ## step 1. update theta
    Hphi_inv = solve(H(d, PHI.samp[b-1]))
    Lambda = solve(diag(N)/TAU.samp[b-1]^2 + Hphi_inv/SIG.samp[b-1]^2)
    ERR = matrix(nrow = N, ncol = Tn)
    
    Ws = Ws.samp[[b-1]]
    THETA = THETA.samp[, , b-1]
    for(t_ in 1:Tn){
      ## q0: prob of sampling from new cluster
      XBeta = matrix(Xs[, , t_], ncol = p) %*% BETA.samp[,b-1,drop=F]
      ERR[,t_] = Y.obs.imp[,t_] - XBeta
      log_q0 = log(NU.samp[b-1]) +
        0.5 * determinant(Lambda, logarithm = T)$modulus -
        0.5/TAU.samp[b-1]^2 * (sum(ERR[,t_]^2) - t(ERR[,t_])%*%Lambda%*%ERR[,t_]/TAU.samp[b-1]^2) -
        N * 0.5 * (log(2*pi)) - N*log(TAU.samp[b-1]) - N*log(SIG.samp[b-1]) +
        0.5 * determinant(Hphi_inv, logarithm = T)$modulus
      attributes(log_q0) = NULL
      
      ## qjs: prob of assigning one of existing clusters
      Ws_tmp = Ws[-t_] # exclude current cluster
      THETA_tmp = THETA[, -t_]
      Tjstars_tmp = unique(THETA_tmp, MARGIN = 2)
      Wsstars_tmp = unique(Ws_tmp)
      Tjs_tmp = rep(0, length = length(Wsstars_tmp))
      for(i in seq_along(Wsstars_tmp)){
        t_st = Wsstars_tmp[i]
        Tjs_tmp[i] = sum(Ws_tmp == t_st)
      }
      
      log_qjs = rep(0, length(Tjs_tmp))
      for(t_st in seq_along(log_qjs)){
        log_qjs[t_st] = log(Tjs_tmp[t_st]) + 
          dmvnorm(Y.obs.imp[,t_], XBeta + Tjstars_tmp[,t_st], TAU.samp[b-1]^2*diag(N), log = T)
      }
      log_qs = c(log_q0, log_qjs)
      qs = exp(log_qs - max(log_qs))
      qs = qs/sum(qs)
      
      if(runif(1) < qs[1]){                                            # sample new cluster
        THETA[, t_] = rmvnorm_Sig(1,
                                  Lambda %*% ERR[,t_] / TAU.samp[b-1]^2,
                                  Lambda)
        Ws[t_] = length(Tjs_tmp) + 1                               # update cluster assignment
      } else {                                                      # assign existing cluster
        fuck = log_qs[-1]
        fuck = exp(fuck - max(fuck))
        idx = sample.int(length(Tjs_tmp), 1, prob = fuck/sum(fuck))
        THETA[, t_] = Tjstars_tmp[, idx]
        Ws[t_] = Wsstars_tmp[idx]                                  # update cluster assignment
      }
      
    }
    
    Ws.samp[[b]] = Ws
    #print(length(unique(Ws)))
    THETA.samp[,,b] = THETA
    Tjstars = unique(THETA.samp[, , b], MARGIN = 2)
    Wsstars = unique(Ws)
    Tjs = rep(0, length = length(Wsstars))
    for(i in seq_along(Wsstars)){
      t_st = Wsstars[i]
      Tjs[i] = sum(Ws == t_st)
    }
    
    # print(length(Tjs))
    # print(ncol(Tjstars))
    # 
    ## step 2. update theta*
    for(i in seq_along(Tjs)){
      Sigma_tmp = solve(Tjs[i] * diag(N) / TAU.samp[b-1]^2 + Hphi_inv / SIG.samp[b-1]^2)
      mu_tmp = Sigma_tmp %*% rowSums(ERR[, Ws == Wsstars[i], drop = F]) / TAU.samp[b-1]^2
      Tjstars[, i] = rmvnorm_Sig(1, mu = mu_tmp, Sigma = Sigma_tmp)
      THETA.samp[, Ws == Wsstars[i], b] = rep(Tjstars[, i], sum(Ws == Wsstars[i]))
    }
    
    
    ## step 3. update beta, tau2
    ## note that we have used g-prior
    ## 1. Beta|tau,D,... ~ N(Beta; mu, Sig)
    ##    with g-prior where g = n*T
    sum_Y_XB = matrix(0, nrow = p, ncol = 1)
    for(t_ in 1:Tn) sum_Y_XB = sum_Y_XB + t(Xs[,,t_]) %*% matrix(Y.obs.imp[,t_] - THETA.samp[,t_,b], ncol = 1)
    Sig_tmp = sumXtX_inv * TAU.samp[b-1]^2
    mu_tmp = Sig_tmp %*% sum_Y_XB / TAU.samp[b-1]^2
    BETA.samp[,b] = rmvnorm_Sig(1, mu_tmp, Sig_tmp)
    
    ## 2. tau|Beta,D,... ~ IG(a, b) (rinvgamma(a,b) = 1/rgamma(a, rate = b))
    ##                        a = a_tau + nT/2
    ##                        b = b_tau + 0.5 * sum_{t=1}^T ||Yt - XtB - Thetat||^2
    sum_Y_XB_Theta = 0
    for(t_ in 1:Tn) sum_Y_XB_Theta = sum_Y_XB_Theta + 
      sum((Y.obs.imp[,t_] - matrix(Xs[, , t_], ncol = p) %*% BETA.samp[,b,drop=F] - THETA.samp[,t_,b])^2)
    TAU.samp[b] = 1/rgamma(1,
                           a_tau + N*Tn*0.5,
                           b_tau + sum_Y_XB_Theta*0.5)
    TAU.samp[b] = sqrt(TAU.samp[b])
    
    ## step 4. update nu, sigma2, phi
    ## 1. nu = p*G(a_nu+T*, b_nu-log(eta)) + (1-p)*G(a_nu+T*-1, b_nu-log(eta))
    ##         eta ~ Beta(nu+1, T)
    ##         p = (a_nu + T* - 1) / (T*(b_nu - log(eta)) + a_nu + T* - 1)
    l_eta_tmp = log(rbeta(1, NU.samp[b-1]+1, Tn))
    p_tmp = (a_nu + length(Tjs) - 1) / (Tn*(b_nu - l_eta_tmp) + a_nu + length(Tjs) - 1)
    NU.samp[b] = ifelse(runif(1) < p_tmp, 
                        rgamma(a_nu + length(Tjs), b_nu - l_eta_tmp), 
                        rgamma(a_nu + length(Tjs) - 1, b_nu - l_eta_tmp))
    ## 2. sigma2 ~ IG(a, b)
    ##                a = a_sig + nT/2
    ##                b = 0.5 * sum_{t=1}^T* t(theta*) %*% Hphi_inv %*% theta*
    sig_tmp = 0
    for(t_st in seq_along(Tjs)) sig_tmp = sig_tmp + t(Tjstars[, t_st]) %*% Hphi_inv %*% Tjstars[, t_st]
    SIG.samp[b] = 1/rgamma(1,
                           a_sig + N*Tn*0.5,
                           b_sig + sig_tmp*0.5)
    SIG.samp[b] = sqrt(SIG.samp[b])
    ## 3. phi ~ f(0, b_phi) \propto |Hphi_inv|^(T*/2) exp(- sum_{t=1}^T* t(theta*) %*% Hphi_inv %*% theta* / 2sig^2)
    ##          b_phi = set in advance
    phi_p = rep(-Inf, length(phi_grid))
    for(k in seq_along(phi_grid)){
      Hphi_inv = solve(H(d, phi_grid[k]))
      phi_tmp = 0
      for(t_st in seq_along(Tjs)) phi_tmp = phi_tmp + t(Tjstars[, t_st]) %*% Hphi_inv %*% Tjstars[, t_st]
      res = 0.5*length(Tjs) * determinant(Hphi_inv, logarithm = T)$modulus - 0.5 * phi_tmp / SIG.samp[b]^2
      attributes(res) = NULL
      phi_p[k] = res
    }
    phi_p = exp(phi_p - max(phi_p)) # to prevent numerical overflow
    phi_p = phi_p / sum(phi_p)
    PHI.samp[b] = sample(phi_grid, 1, prob = phi_p)
    
    # step 5. impute missing data
    for(t_ in 1:Tn){
      NA_idx = NA_IDX[,t_]
      if(sum(NA_idx) == 0) next
      Y.obs.imp[NA_idx, t_] = matrix(Xs[NA_idx, , t_],ncol = p) %*% BETA.samp[,b,drop=F] + THETA.samp[NA_idx, t_,b] +
        TAU.samp[b] * rnorm(sum(NA_idx))
    }
    Y.comp.samp[,,b] = Y.obs.imp
  }
  
  traceidx = (1+burnin):B
  return(list(
    Y.comp.samp = Y.comp.samp[,,(1+burnin):B],
    THETA.samp = THETA.samp[,,(1+burnin):B],
    Ws.samp = Ws.samp[(1+burnin):B],
    BETA.samp = BETA.samp[,(1+burnin):B],
    TAU.samp = TAU.samp[(1+burnin):B],
    NU.samp = NU.samp[(1+burnin):B],
    PHI.samp = PHI.samp[(1+burnin):B],
    SIG.samp = SIG.samp[(1+burnin):B]
  ))
}










