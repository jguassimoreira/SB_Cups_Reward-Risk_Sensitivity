ptFitMod = function(params, data) {
  
  rho = params[1]
  lam = params[2]
  mu = params[3]
  
  #utility for nonzero outcome of risky decision
  uNonZero = ifelse(data$nonZeroAmount > 0, 
                    data$riskProb * (data$nonZeroAmount ^ rho),
                    (data$riskProb * -lam * (-data$nonZeroAmount) ^ rho))
  #loss utility for zero out of risky decision (ultimately just zero)
  uZero = (1-data$riskProb) * (data$zeroAmount^rho)
  #total utility for risky decision
  uTot = uNonZero + uZero
  #utility for certain decisiion
  uCer = ifelse(data$certainAmount > 0, 
                1 * (data$certainAmount ^ rho),
                (1 * -lam * (-data$certainAmount) ^ rho))
  #use softmax to find probability of taking a risk
  prob = (1+(exp(-mu*(uTot - uCer))))^-1
  
  prob[prob == 1] = .99 #prob can't be ==1 or 0, else it won't work. So change to smidge below/above.
  prob[prob == 0] = .01
  #transform to loglikelihood
  loglik = sum(data$dec*log(prob) + (1 - data$dec)*log(1-prob))
  #return negative loglikelihood, needed for optimization (since optim minimizes a function)
  
  if(params[1] < 0) return(99999)
  if(params[1] > 10) return(99999)
  if(params[2] < 0) return(99999)
  if(params[2] > 10) return(99999)
  if(params[3] < .01) return(99999)
  if(params[3] > 20) return(99999)
  
  #if(p[1] + p[2] < 0) return(99999)
  #if(p[1] + p[3] < 0) return(99999)
  #if(p[1] + p[2] + p[3] < 0) return(99999)
  #if(p[4] < 0) return(99999)
  #if(p[4] > 1) return(99999)
  #if(p[5] < 0) return(99999)
  
  return(-loglik)
  
  
  
}