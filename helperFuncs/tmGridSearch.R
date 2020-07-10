#################################################################
## function to estimate target model parameters using grid search
tmGridSearch = function(data, bound, bin) {
  
  bounds = matrix(bound, length(bound)/2, 2, byrow = TRUE) #order of pairs is tau and beta
  nBin = matrix(bin, 1, 2) # same as above, tau, beta
  nParam = dim(bounds)[1] #number of parameters (should always be two for this model)
  p = matrix(list(),1,nParam) #matrix describing probabilities of safe v risky decision for each trial given the parameter
  loglik = array(NaN, dim=c(nBin[1], nBin[2])) #set our grid
  
  for (i in 1:nParam) {
    #create linearly spaced intervals (equal to the number of bins) for potential parameter values
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

  
  loglik = (loglik - mean(loglik)/std(loglik)) #standardize
  lik = exp(loglik) #make them likelihoods
  
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