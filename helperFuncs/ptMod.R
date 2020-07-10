############################################################
## custom function to fit prospect theory model to cups data
ptMod = function(inDat, gbound, gbin, obound, grid) {
  
  #input(s): a dataframe containing a subject's data, bounds for the gridsearch, number of bins for the grid search, and an option indicating whether we want to get starting values from a grid search
  #output(s): parameter estimates for rho (risk aversion), lambda (loss aversion), mu (choice consistency), # of iterations it took to converge, the error message, and the value of the likelihood function
  
  #grid is a number that indicates whether you need to conduct a grid search or not 
  #if it is zero, that means you want to use the grid search
  #if it is non zero, you don't want grid search, you want starting values to be randomly initialized
  
  dat = ptReformatCups(inDat) #reformat the cups data to get it just the way we need it (decisions, risk probability, nonzero outcome amount)
  
  fun <- function(pm) ptFitMod(pm, data=dat) #create the function handle 
  
  if (grid == 0) { 
  
    startVals = ptGridSearch(dat, gbound, gbin)   #run grid search to get the starting values
    
  } else {
    
    set.seed(grid) #randomly initialized startving values if you don't want to use grid search
    startVals = list(runif(1, 0, 10),
                     runif(1, 0, 10), 
                     runif(1, 0, 25))
    
  }
  
  
  #I was using the commented chunk of text below before and it wasn't converging as well. Then I 'hacked' optim to performed 'constrainted' optimzation with Nelder-Mead (optim's default)
  #conParams = optim(c(startVals[[1]], startVals[[2]], startVals[[3]]),
  #                  method="L-BFGS-B",
  #                  fn=fun,
  #                  lower=c(1e-7, 1e-7, .4),
  #                  upper = c(obound[1], obound[2], obound[3]))
  conParams = optim(c(startVals[[1]], startVals[[2]], startVals[[3]]),
                    fn=fun) #run optim
  
  return(list(conParams$par[1], conParams$par[2], conParams$par[3], conParams$counts, conParams$convergence, conParams$message, conParams$value))
  
}