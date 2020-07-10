##########################################
## custom function to load in cups data ##
getCupsDat = function(subList) {
  
  #input(s): a vector of paths to each subject's cups data
  #output(s): a tibble containing trial-by-trial cups data for every subject in the study
  
  cupsDat = tibble(ID="Null", Context=999, Decision=999, EV=999, SD=999, TrialType=999, TrialNum=999) #create empty tibble with placeholder for first row
  
  for (s in subList) { #loop over vector of subject paths
    
    dat = read_csv(s) #read in the data
    dat$TrialNum = c(1:length(dat[[1]])) #add a trial number variable
    
    cupsDat = full_join(cupsDat, dat) #join the data from the current iteration with the aggregate dataset
    
  }
  
  return(cupsDat) #output the tibble
  
}