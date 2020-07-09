tmGridSearch = function(data, bound, bin) {
  
  bounds = matrix(bound, length(bound)/2, 2, byrow = TRUE) #order of pairs is rho, lambda, mu
  nBin = matrix(bin, 1, 2) # same as above, alpha, beta, lambda self, lambda other
  nParam = dim(bounds)[1] #number of parameters (should always be three)
  p = matrix(list(),1,nParam)
  loglik = array(NaN, dim=c(nBin[1], nBin[2]))
  
  for (i in 1:nParam) {
    
    range = linspace(bounds[i,1],bounds[i,2],nBin[i]+1)
    p[[1,i]] = range[2:length(range)]
    
  }
  
  #Do the grid search 
  params = matrix(NaN, 1, 2)
  for (i in 1:nBin[1]) {
    params[1,1] = p[[1,1]][i] #tau
    for (o in 1:nBin[2]) {
      params[1,2] = p[[1,2]][o] #beta
        tmp = tmFitMod(params,data)
        loglik[i,o] = -tmp #negate the negative in the model fitting function since we're doing a grid search here
      }
    }

  
  loglik = (loglik - mean(loglik)/std(loglik))
  lik = exp(loglik)
  
  tmpTau = rowSums(lik) #sum across other parameter's for a given parameter's likelihood
  tmpBeta = colSums(lik) #needed to compute marginal likelihood below

  marglik_Tau = tmpTau/sum(tmpTau) #compute marginal likelihoods of all parameters
  marglik_Beta = tmpBeta/sum(tmpBeta)
  
  ML_Tau = p[[1]][max(marglik_Tau) == marglik_Tau] #Max likelihood estimates
  ML_Beta = p[[2]][max(marglik_Beta) == marglik_Beta]
  
  EV_Tau = sum(p[[1]] * marglik_Tau) #Expected value estimates
  EV_Beta = sum(p[[2]] * marglik_Beta)
  
  start_Tau = mean(c(ML_Tau, EV_Tau)) #average the ML and EV estimates to create starting values for below
  start_Beta = mean(c(ML_Beta, EV_Beta))
  
  return(list(start_Tau, start_Beta))
  
}