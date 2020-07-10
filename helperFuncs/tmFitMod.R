################################
## target model fitting function
tmFitMod = function(params, data) {
  
  #input(s): parameters -- tau (target threshold) and data
  #output(s): parameter estimates for tau, # of iterations it took to converge, the error message, and the value of the likelihood function
  
  tau = params[1]
  beta = params[2] #choice consistency here is notated with beta (consistent with origin of the model, see below)
  
  #Implementing the model used in Wallsten et al 2005 Psych Review (they used it for the BART, what I have here is an adaptation)
  
  #define current reward/risk ratio for each trial
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
  
  #these next couple lines 'hack' optim so that we can perform 'constrained' optimization with the algorithms that don't allow for it
  if(params[1] < 0) return(99999)
  if(params[1] > 10) return(99999)
  if(params[2] < .01) return(99999)
  if(params[2] > 20) return(99999)
  
  return(-loglik) #returning the negative here to be congruent with optim
  
  
  
}