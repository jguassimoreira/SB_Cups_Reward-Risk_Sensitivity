ptReformatCups = function(dat) {
  
  #start by bringing over the decisions from the raw data
  outDat = data.frame(dec = dat$cueResp.corr)
  
  #reformat the certain amount (+2 for gain trials, -2 for loss trials)
  outDat$certainAmount = case_when(dat$TrialType == 1 ~ 2, dat$TrialType == 0 ~ -2)
  
  #reformat the nonZero amounts for risky choices 
  outDat$nonZeroAmount = case_when(dat$imageDec == "CueLoss24.jpg" ~ -4,
                             dat$imageDec == "CueLoss26.jpg" ~ -6,
                             dat$imageDec == "CueLoss210.jpg" ~ -10,
                             dat$imageDec == "CueLoss34.jpg" ~ -4,
                             dat$imageDec == "CueLoss36.jpg" ~ -6,
                             dat$imageDec == "CueLoss310.jpg" ~ -10,
                             dat$imageDec == "CueLoss54.jpg" ~ -4,
                             dat$imageDec == "CueLoss56.jpg" ~ -6,
                             dat$imageDec == "CueLoss510.jpg" ~ -10,
                             dat$imageDec == "CueWin24.jpg" ~ 4,
                             dat$imageDec == "CueWin26.jpg" ~ 6,
                             dat$imageDec == "CueWin210.jpg" ~ 10,
                             dat$imageDec == "CueWin34.jpg" ~ 4,
                             dat$imageDec == "CueWin36.jpg" ~ 6,
                             dat$imageDec == "CueWin310.jpg" ~ 10,
                             dat$imageDec == "CueWin54.jpg" ~ 4,
                             dat$imageDec == "CueWin56.jpg" ~ 6,
                             dat$imageDec == "CueWin510.jpg" ~ 10)
  
  #create a column with the zero amounts for the risky choices
  outDat$zeroAmount = rep(0, length(outDat$certainAmount))
  
  #Probabilities of nonZero amounts for risky choices
  outDat$riskProb = case_when(dat$imageDec == "CueLoss24.jpg" ~ (1/2),
                              dat$imageDec == "CueLoss26.jpg" ~ (1/2),
                              dat$imageDec == "CueLoss210.jpg" ~ (1/2),
                              dat$imageDec == "CueLoss34.jpg" ~ (1/3),
                              dat$imageDec == "CueLoss36.jpg" ~ (1/3),
                              dat$imageDec == "CueLoss310.jpg" ~ (1/3),
                              dat$imageDec == "CueLoss54.jpg" ~ (1/5),
                              dat$imageDec == "CueLoss56.jpg" ~ (1/5),
                              dat$imageDec == "CueLoss510.jpg" ~ (1/5),
                              dat$imageDec == "CueWin24.jpg" ~ (1/2),
                              dat$imageDec == "CueWin26.jpg" ~ (1/2),
                              dat$imageDec == "CueWin210.jpg" ~ (1/2),
                              dat$imageDec == "CueWin34.jpg" ~ (1/3),
                              dat$imageDec == "CueWin36.jpg" ~ (1/3),
                              dat$imageDec == "CueWin310.jpg" ~ (1/3),
                              dat$imageDec == "CueWin54.jpg" ~ (1/5),
                              dat$imageDec == "CueWin56.jpg" ~ (1/5),
                              dat$imageDec == "CueWin510.jpg" ~ (1/5))
  
  return(outDat)
  
}