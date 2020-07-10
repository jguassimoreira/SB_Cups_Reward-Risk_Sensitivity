#################################################
## custom K Nearest Neighbors imputation function
imputeKNN = function(dat, var, k) {
  
  #input(s): data (including target variable -- the one we want to impute -- and support variables -- the ones we're going to use to identify nearest neighbors), name of variable to impute, and number of neighbors
  #output(s): dataframe with a column for the target variable with complete cases (observed + imputed)
  
  outDat = dat[sprintf('%s', var)]
  
  #loop over subjects 
  for (c in 1:nrow(dat)) { 
    
    #if the case is complete, skip it
    if (!is.na(dat[c,sprintf('%s', var)])) next
    
    #create new dataframe that includes all variables and only complete cases 
    completeDat = dat[!is.na(dat[sprintf('%s',var)]),]
    
    #add the case in the current loop (remember it has to be missing)
    completeCDat = rbind(completeDat, dat[c,])
    
    #compute the distance between the current, missing case and the other cases based on the *support* variables
    dists = as.matrix(dist(completeCDat[,-which(names(completeCDat) == sprintf('%s', var))]))[nrow(completeCDat), ]
    #remove the last row, which just contains the distance of the missing case with itself
    validDists = as.vector(dists)[-length(dists)]
    
    #find the k closest cases
    neighbrs = order(validDists)[1:k]
    wts = c(k:1) #create ranks used for weights (closer neighbors --> higher weights)
    
    imputedDat = sum(completeDat[sprintf('%s', var)][neighbrs,] * wts) / sum(wts) #take weighted average of neighbors
    
    outDat[c,1] = imputedDat #place the imputed observation into the output dataframe
        
  }
  
  return(outDat[[1]])
  
}