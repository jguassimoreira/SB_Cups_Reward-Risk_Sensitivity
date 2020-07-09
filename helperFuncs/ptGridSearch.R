ptGridSearch = function(data, bound, bin) {
  
  bounds = matrix(bound, length(bound)/2, 2, byrow = TRUE) #order of pairs is rho, lambda, mu
  nBin = matrix(bin, 1, 3) # same as above, alpha, beta, lambda self, lambda other
  nParam = dim(bounds)[1] #number of parameters (should always be three)
  p = matrix(list(),1,nParam)
  loglik = array(NaN, dim=c(nBin[1], nBin[2], nBin[3]))
  
  for (i in 1:nParam) {
    
    range = linspace(bounds[i,1],bounds[i,2],nBin[i]+1)
    p[[1,i]] = range[2:length(range)]
    
  }
  
  #Do the grid search 
  params = matrix(NaN, 1, 3)
  for (i in 1:nBin[1]) {
    params[1,1] = p[[1,1]][i] #rho
    for (o in 1:nBin[2]) {
      params[1,2] = p[[1,2]][o] #lambda
      for (j in 1:nBin[3]) {
        params[1,3] = p[[1,3]][j] #mu
        tmp = ptFitMod(params,data)
        loglik[i,o,j] = -tmp #negate the negative in the model fitting function since we're doing a grid search here
      }
    }
  }
  
  loglik = (loglik - mean(loglik)/std(loglik))
  lik = exp(loglik)
  
  tmpRho = rowSumArray3(lik) #sum across other parameter's for a given parameter's likelihood
  tmpLam = colSumArray3(lik) #needed to compute marginal likelihood below
  tmpMu = d3SumArray3(lik)

  marglik_Rho = tmpRho/sum(tmpRho) #compute marginal likelihoods of all parameters
  marglik_Lam = tmpLam/sum(tmpLam)
  marglik_Mu =  tmpMu/sum(tmpMu)  
  
  ML_Rho = p[[1]][max(marglik_Rho) == marglik_Rho] #Max likelihood estimates
  ML_Lam = p[[2]][max(marglik_Lam) == marglik_Lam]
  ML_Mu = p[[3]][max(marglik_Mu) == marglik_Mu]
  
  EV_Rho = sum(p[[1]] * marglik_Rho) #Expected value estimates
  EV_Lam = sum(p[[2]] * marglik_Lam)
  EV_Mu = sum(p[[3]] * marglik_Mu)
  
  start_Rho = mean(c(ML_Rho, EV_Rho)) #average the ML and EV estimates to create starting values for below
  start_Lam = mean(c(ML_Lam, EV_Lam))
  start_Mu = mean(c(ML_Mu, EV_Mu))
  
  return(list(start_Rho, start_Lam, start_Mu))
  
}