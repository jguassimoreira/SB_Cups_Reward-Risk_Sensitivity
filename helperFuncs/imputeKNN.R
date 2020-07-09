imputeKNN = function(dat, var, k) {
  
  outDat = dat[sprintf('%s', var)]
  
  for (c in 1:nrow(dat)) {
    
    if (!is.na(dat[c,sprintf('%s', var)])) next
    
    completeDat = dat[!is.na(dat[sprintf('%s',var)]),]
    
    completeCDat = rbind(completeDat, dat[c,])
    
    dists = as.matrix(dist(completeCDat[,-which(names(completeCDat) == sprintf('%s', var))]))[nrow(completeCDat), ]
    validDists = as.vector(dists)[-length(dists)]
    
    neighbrs = order(validDists)[1:k]
    wts = c(k:1)
    
    imputedDat = sum(completeDat[sprintf('%s', var)][neighbrs,] * wts) / sum(wts)
    
    outDat[c,1] = imputedDat
        
  }
  
  return(outDat[[1]])
  
}