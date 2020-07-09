cleanPTParams = function(inDat) {
  
  names = names(inDat); n = strsplit(names, "_")
  likMask = lapply(n, `[[`, 1) == 'lik'
  rhoMask = lapply(n, `[[`, 1) == 'rho'
  lamMask = lapply(n, `[[`, 1) == 'lambda'
  
  outDat = data.frame(ID = inDat$ID,
                       rho = rep(NA, length(inDat$ID)),
                       lambda = rep(NA,length(inDat$ID)),
                       stringsAsFactors = FALSE)
  
  for (r in 1:dim(inDat)[1]) {
    
    lDat = round(inDat[r,likMask],5)
    paramIDs = which(lDat == min(lDat))
    rho = rowMeans(inDat[r,rhoMask][paramIDs])
    lambda = rowMeans(inDat[r,lamMask][paramIDs])
    outDat[r,'rho'] = rho
    outDat[r,'lambda'] = lambda
    
  }
  
  return(outDat)
  
}