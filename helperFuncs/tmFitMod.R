tmFitMod = function(params, data) {
  
  tau = params[1]
  beta = params[2]
  
  #define current ratio
  i = data$EV / data$SD
  
  #define delta
  d = tau - i
  d[which(i < 0)] = (-tau) - i[which(i < 0)]
  
  
  #calculate risk probability
  r = 1 / (1 + exp(beta*d))
  
  r[r == 1] = .99 #prob can't be ==1 or 0, else it won't work. So change to smidge below/above.
  r[r == 0] = .01
  #transform to loglikelihood
  loglik = sum(data$Decision*log(r) + (1 - data$Decision)*log(1-r), na.rm = TRUE)
  #return negative loglikelihood, needed for optimization (since optim minimizes a function)
  
  if(params[1] < 0) return(99999)
  if(params[1] > 10) return(99999)
  if(params[2] < .01) return(99999)
  if(params[2] > 20) return(99999)
  
  return(-loglik)
  
  
  
}