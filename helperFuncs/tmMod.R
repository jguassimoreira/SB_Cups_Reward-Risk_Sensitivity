###################################################
## custom function to fit target model to cups data
tmMod = function(inDat, gbound, gbin, obound, grid) {
  
  #grid is a number that indicates whether you need to conduct a grid search or not
  #if it is zero, that means you want to use the grid search
  #if it is non zero, you don't want grid search, you want parameters to be randomly initialized
  
  fun = function(pm) tmFitMod(pm, data=inDat) #create the function handle 
  
  if (grid == 0) {
  
    startVals = tmGridSearch(inDat, gbound, gbin)   #run grid search to get the starting values
    
  } else {
    
    set.seed(grid) #randomly initialized startving values if you don't want to use grid search
    startVals = list(runif(1, 0, 10), 
                     runif(1, 0, 25))
    
  }
  
  
  #I was using the commented chunk of text below before and it wasn't converging as well. Then I 'hacked' optim to performed 'constrainted' optimzation with Nelder-Mead (optim's default)
  #conParams = optim(c(startVals[[1]], startVals[[2]]),
  #                  method="L-BFGS-B",
  #                  fn=fun,
  #                  lower=c(.05, 1),
  #                  upper = c(obound[1], obound[2]))
  conParams = optim(c(startVals[[1]], startVals[[2]]),
                    fn=fun)
  
  return(list(conParams$par[1], conParams$par[2], conParams$counts, conParams$convergence, conParams$message, conParams$value))
  
}