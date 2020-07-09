ptMod = function(inDat, gbound, gbin, obound, grid) {
  
  #grid is a number that indicates whether you need to conduct a grid search or not
  #if it is zero, that means you want to use the grid search
  #if it is non zero, you don't want grid search, you want parameters to be randomly initialized
  
  dat = ptReformatCups(inDat)
  
  fun <- function(pm) ptFitMod(pm, data=dat)
  
  if (grid == 0) {
  
    startVals = ptGridSearch(dat, gbound, gbin)   
    
  } else {
    
    set.seed(grid)
    startVals = list(runif(1, 0, 10),
                     runif(1, 0, 10), 
                     runif(1, 0, 25))
    
  }
  
  
  
  #conParams = optim(c(startVals[[1]], startVals[[2]], startVals[[3]]),
  #                  method="L-BFGS-B",
  #                  fn=fun,
  #                  lower=c(1e-7, 1e-7, .4),
  #                  upper = c(obound[1], obound[2], obound[3]))
  conParams = optim(c(startVals[[1]], startVals[[2]], startVals[[3]]),
                    fn=fun)
  
  return(list(conParams$par[1], conParams$par[2], conParams$par[3], conParams$counts, conParams$convergence, conParams$message, conParams$value))
  
}