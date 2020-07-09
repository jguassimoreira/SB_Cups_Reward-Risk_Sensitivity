getCupsDat = function(subList) {
  
  cupsDat = tibble(ID="Null", Context=999, Decision=999, EV=999, SD=999, TrialType=999, TrialNum=999)
  
  for (s in subList) {
    
    dat = read_csv(s)
    dat$TrialNum = c(1:length(dat[[1]]))
    
    cupsDat = full_join(cupsDat, dat)
    
  }
  
  return(cupsDat)
  
}