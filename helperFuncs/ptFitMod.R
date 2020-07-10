#########################################
## prospect theory model fitting function
ptFitMod = function(params, data) {
  
  #input(s): parameters -- rho (risk aversion), lambda (loss aversion), mu (choice consistency) -- and data
  #output(s): parameter estimates for rho (risk aversion), lambda (loss aversion), mu (choice consistency), # of iterations it took to converge, the error message, and the value of the likelihood function
  
  rho = params[1] #set rho
  lam = params[2] #set lambda
  mu = params[3] #set mu
  
  #Implementing the model used by P. Sokol-Hessner's lab (see Sokol-Hessner et al., 2009 (PNAS), 2015 (Psych Science))
  
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
  #use softmax to find probability of taking a risks
  prob = (1+(exp(-mu*(uTot - uCer))))^-1
  
  prob[prob == 1] = .99 #prob can't be ==1 or 0, else it won't work. So change to smidge below/above.
  prob[prob == 0] = .01
  #transform to loglikelihood
  loglik = sum(data$dec*log(prob) + (1 - data$dec)*log(1-prob))
  #return negative loglikelihood, needed for optimization (since optim minimizes a function)

  #these next couple lines 'hack' optim so that we can perform 'constrained' optimization with the algorithms that don't allow for it
  if(params[1] < 0) return(99999)
  if(params[1] > 10) return(99999)
  if(params[2] < 0) return(99999)
  if(params[2] > 10) return(99999)
  if(params[3] < .01) return(99999)
  if(params[3] > 20) return(99999)

  return(-loglik) #returning the negative here to be congruent with optim
  
  
  
}