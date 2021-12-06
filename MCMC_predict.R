SDPpredict = function(RES, 
                      X_preds,
                      coords.train, coords.pred){
  
  Y.comp.samp = RES$Y.comp.samp
  Ws.samp = RES$Ws.samp
  THETA.samp = RES$THETA.samp
  BETA.samp = RES$BETA.samp
  TAU.samp = RES$TAU.samp
  NU.samp = RES$NU.samp
  PHI.samp = RES$PHI.samp
  SIG.samp = RES$SIG.samp
  
  Tn = dim(Y.comp.samp)[2] 
  N = dim(Y.comp.samp)[1] 
  Npred = nrow(coords.pred)
  p = dim(BETA.samp)[1]
  B = dim(Y.comp.samp)[3]
  
  library(Rcpp)
  library(RcppArmadillo)
  sourceCpp("rdist_cpp.cpp")
  
  
  # distance matrix with train, pred combined
  d = rdist_cpp(rbind(coords.train, coords.pred), base::crossprod)
  train_idx = 1:nrow(coords.train)
  pred_idx = nrow(coords.train)+ 1:nrow(coords.pred)
  H = function(d, phi) exp(-phi * d)
  library(mvtnorm) # for dmvnorm
  
  Y.pred.samp = array(rep(matrix(nrow = Npred, ncol = Tn), B),
                      dim = c(Npred, Tn, B))
  ETA.pred.samp = array(rep(matrix(nrow = Npred, ncol = Tn), B),
                      dim = c(Npred, Tn, B))
  
  for(b in 1:B){
    if((b %% 100)==0) cat(b, "\n")
    Hphi = H(d, PHI.samp[b])
    Hphi_11 = Hphi[train_idx, train_idx]
    Hphi_11_inv = solve(Hphi_11)
    Hphi_12 = Hphi[train_idx, pred_idx]
    Hphi_21 = Hphi[pred_idx, train_idx]
    Hphi_22 = Hphi[pred_idx, pred_idx]
    cond_Sig = Hphi_22 - Hphi_21 %*% Hphi_11_inv %*% Hphi_12
    cond_mu = Hphi_21 %*% Hphi_11_inv
    
    # sample new theta*
    Ws = Ws.samp[[b]]
    THETA = THETA.samp[,,b]
    Tjstars = unique(THETA, MARGIN = 2)
    Wsstars = unique(Ws)
    Tjs = rep(0, length = length(Wsstars))
    for(i in seq_along(Wsstars)){
      Tjs[i] = sum(Ws == Wsstars[i])
    }
    
    Tjstars_pred = matrix(0, nrow = nrow(coords.pred), ncol = length(Wsstars))
    for(i in seq_along(Wsstars)){
      Tjstars_pred[, i] = rmvnorm_Sig(1,
                                      cond_mu %*% Tjstars[,i],
                                      cond_Sig)
    }
    
    # sample new Y
    prob = c(Tjs, NU.samp[b])
    prob = prob/sum(prob)
    for(t_ in 1:Tn){
      idx = sample.int(length(prob), 1, prob= prob)
      mean_tmp = matrix(X_preds[, , t_],ncol = p) %*% BETA.samp[,b,drop=F]
      if(idx == length(prob)){ # sample from base measure
        sig_tmp = SIG.samp[b]^2 * Hphi_22
        ETA.pred.samp[, t_, b]= mean_tmp
        Y.pred.samp[, t_, b] = rmvnorm_Sig(1,
                                           mean_tmp,
                                           sig_tmp) + TAU.samp[b] * rnorm(Npred)
      } else { # from existing cluster
        ETA.pred.samp[, t_, b]= mean_tmp + Tjstars_pred[,idx]
        Y.pred.samp[, t_, b] = mean_tmp + Tjstars_pred[,idx] + TAU.samp[b] * rnorm(Npred)
      }
    }
  }
  
  return(list(Y.pred.samp = Y.pred.samp, ETA.pred.samp = ETA.pred.samp))
}