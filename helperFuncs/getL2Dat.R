################################################
## custom function to load subject-level data ##
getL2Dat = function(L2DatPath) {
  
  #input(s): path to csv containing subject level data (already compiled on SANDLab server)
  #output(s): a dataframe containing the data in said csv 
  
  L2Dat = read_csv(L2DatPath) #read in data, output the dataframe
  
}