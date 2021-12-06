# response: Mixture of two GP
rmvnorm_Sig = function(n, mu, Sigma){
  p = length(mu)
  Z = matrix(rnorm(p*n), p, n)
  L = t(chol(Sigma))
  X = L %*% Z
  X = sweep(X, 1, mu, FUN = "+")
  return(X)
}

# nonstationary Gaussian
nonstGP = function(N, X, coords, Beta, sigma, rho, tau){
  # process model eta ~ MixGP(Xb, sig^2*H(rho))
  # beta: regression coefficients
  # sigma: GP variance
  # rho: GP correlation
  # H(rho): nonstationary correlation function
  # tau: observation error (nugget)
  
  X = cbind(1, X)
  
  # nonstationary covariance function
  H_nonstan = function(xy, sigma, rho){
    x = xy[,1]; x = x - median(x)
    y = xy[,2]; y = y - median(y)
    xy = cbind(x,y)
    res = matrix(NA_real_, nrow(xy), nrow(xy))
    for(i in 1:nrow(xy)){
      for(j in 1:i){ # fill lower triangle components, including diagonal
        res[i,j] = 
          max(1, sigma * sqrt(sum(xy[i,]^2))) * 
          max(1, sigma * sqrt(sum(xy[j,]^2))) *
          exp(-rho*sqrt(sum((xy[i,] - xy[j,])^2)))
      }
    }
    res[upper.tri(res)] = t(res)[upper.tri(res)]
    return(res)
  }
  
  Sigma = H_nonstan(coords, sigma, rho)
  
  Y = matrix(NA_real_, nrow = nrow(coords), ncol = N)
  #ETA = matrix(NA_real_, nrow = nrow(coords), ncol = N)
  for(n in 1:N){
    eta = rmvnorm_Sig(1, X %*% Beta, Sigma)
    #ETA[,n] = eta
    Y[,n] = eta + rnorm(nrow(eta))*tau # data model y ~ GP(eta, tau^2)
  }
  
  
  Y = round(Y, 2)
  return(Y)
  #return(list(Y=Y, ETA= ETA))
}


# plot covariances
plotcov = function(Y.pred, coords.pred,
                   Y.true,
                   Y.obs, coords.obs,
                   showplot = T
){
  norm_ = function(x, y) sqrt(sum((x-y)^2))
  # distance, empricial covariances, theoretical covariances of 
  # the specific test site vis-a-vis others in train data
  d = Cov.pred = Cov.true = rep(0, nrow(coords.obs))
  for(i in 1:nrow(coords.obs)){
    d[i] = norm_(coords.pred, coords.obs)
    Cov.pred[i] = cov(Y.pred, Y.obs[,i])
    Cov.true[i] = cov(Y.true, Y.obs[,i])
  }
  idx = base::order(d, descending = F)
  d = d[idx]
  Cov.pred = Cov.pred[idx]
  Cov.true = Cov.true[idx]
  if(showplot){
    plot(d, Cov.pred, type = "p", col = "black")
    points(d, Cov.true, col = "red")
    main("True vs Posterior covariances")
  }
  return(invisible(list(d = d, Cov.pred = Cov.pred, Cov.true = Cov.true)))
}