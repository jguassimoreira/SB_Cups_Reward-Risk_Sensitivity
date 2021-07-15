###########################################################################
######################## Analysis code #################################### 
########### JF Guassi Moreira (jguassimoreira[at]ucla[dot]edu) ############
############################ Spring 2021 ##################################
setwd("P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/")


#########
## Load the libraries we'll need
#########
library(readr)
library(dplyr)
library(tidyr)
library(ggmcmc)
library(lattice)
library(rjags)
library(R2jags)
library(ggplot2)
library(hrbrthemes)
library(ggpubr)

#########
## Import helper functions, set path to data
#########
scriptPath = file.path("P:", "Social_Behavior_(SB_study)", "SB_Scripts", "wave1",
                       "Behavioral_tasks", "SB_Cups_Self", "SB_Cups_Reward-Risk_Sensitivity", "helperFuncs") #path to helper funcs
file.sources = Sys.glob(file.path(sprintf("%s", scriptPath), "*.R")) #Construct paths to helper funcs, save in vector
sapply(file.sources,FUN=source) #use sapply to source each file path in the helper func vector
sbPath = file.path("P:", "Social_Behavior_(SB_study)", "Data", "Behavioral_data") #set path to data on SANDLab server


#########
## Load in the data
#########
subList = Sys.glob("P:/Social_Behavior_(SB_study)/Data/Behavioral_data/SB*/wave1/Lab_session/Raw/*_CupsSelf_Raw.csv") #Get path to each subjects' cups data

cupsDat = getCupsDatCompMod(subList, mod = "prospect") #Compile Cups data (Level 1, i.e., trial level), uses helper func
cupsDat = cupsDat[-1,] #The first row is a bunch of 999 placeholders because I coded this like a doofus the first time around

L2Dat = getL2Dat("P:/Social_Behavior_(SB_study)/Data/Study_logs/SB_group_list.csv") #Load the level two (i.e., between subjects) data, uses helper func (though it's not really necessary)
L2Dat = L2Dat[L2Dat$ID %in% unique(cupsDat$ID),]

#L2Dat = L2Dat[L2Dat$PI==0,]

L2Dat = rbind(L2Dat[L2Dat$PI==0,], L2Dat[L2Dat$PI==1,])

L2Dat$ID_num = 1:length(L2Dat$ID)

L2Dat = select(L2Dat, ID = ID, ID_num, Age, PI, Sex); L2Dat = arrange(L2Dat, ID)

allDat = inner_join(cupsDat,L2Dat) #merge the trial- and subject-level data

allDat = allDat[order(allDat$ID_num),]

allDat = rbind(allDat[allDat$trialType==1,], allDat[allDat$trialType==0,])




#########
## Prospect theory modeling
#########

y = allDat$dec
certain = allDat$certainAmount
nonZero = allDat$nonZeroAmount
riskProb = allDat$riskProb
trialType = allDat$trialType
piStat = L2Dat$PI
id = allDat$ID_num

sink("compMod_ptMod1_priorSet1-groupDiff.txt")
cat("
# PT model 1  
model
{

  for (j in 1:83) {
  
  #begin by defining participant j has unique parameter values for
  #4 parameters: alpha (risk attitude, gain domain), beta (risk attitude, loss domain)
  #lambda (loss aversion), and mu (choice consistency)
  #alpha and beta are constrainted [0,1], so they're transformed to probit scale
  #lambda and mu are also constrained [0,inf), so we assume they come from lognormal dist
  
  #these are just the initial transformations
  rho[j] <- exp(lrho[j])
  lambda[j] <- exp(llambda[j])
	mu[j]   <- exp(lmu[j]) 
  
  #now we'll say that each subject's parameters were drawn from normal distributions
  #centered around the group mean 
  
  lrho[j] ~ dnorm(lmu.rho,ltau.rho)
  llambda[j] ~ dnorm(lmu.lambda, ltau.lambda)
  lmu[j] ~ dnorm(lmu.mu, ltau.mu)
  
  }
  
  for (j in 84:145) {
  
  rho[j] <- exp(lrho.pi[j])
  lambda[j] <- exp(llambda.pi[j])
	mu[j]   <- exp(lmu.pi[j]) 
	
	lrho.pi[j] ~ dnorm(lmu.rho.pi,ltau.rho.pi)
  llambda.pi[j] ~ dnorm(lmu.lambda.pi, ltau.lambda.pi)
  lmu.pi[j] ~ dnorm(lmu.mu.pi, ltau.mu.pi)
  
  }

  ##Priors for hyper distributions 
  
  #for comparison youths
  
  #for rho (risk attitudes)
  lmu.rho ~ dunif(-5,5) 
  ltau.rho = pow(lsigma.rho,-2)
  lsigma.rho ~ dunif(0,1)

  #for lambda (loss aversion)
  lmu.lambda ~ dunif(-5,5)
  ltau.lambda = pow(lsigma.lambda,-2)
  lsigma.lambda ~ dunif(0,1.5)
  
  #for mu (choice consistency)
  lmu.mu ~ dunif(0, 7.5)
  ltau.mu = pow(lsigma.mu,-2)
  lsigma.mu ~ dunif(0,2)
  
  #for pi youths
  
  #for rho (risk attitudes)
  lmu.rho.pi ~ dunif(-5,5) #1/900
  ltau.rho.pi = pow(lsigma.rho.pi,-2)
  lsigma.rho.pi ~ dunif(0,1)

  #for lambda (loss aversion)
  lmu.lambda.pi ~ dunif(-5,5)
  ltau.lambda.pi = pow(lsigma.lambda.pi,-2)
  lsigma.lambda.pi ~ dunif(0,1.5)
  
  #for mu (choice consistency)
  lmu.mu.pi ~ dunif(0, 7.5)
  ltau.mu.pi = pow(lsigma.mu.pi,-2)
  lsigma.mu.pi ~ dunif(0,2)
  
  # To obtain the mean of the hyper distribution on the wanted scale:
  
  #comparison
	mu.rho <- exp(lmu.rho) 
	mu.lambda <- exp(lmu.lambda)
  mu.mu   <- exp(lmu.mu)
  
  #pi and difference
  mu.rho.pi <- exp(lmu.rho.pi) 
	mu.lambda.pi <- exp(lmu.lambda.pi)
  mu.mu.pi   <- exp(lmu.mu.pi)

  mu.rho.diff = mu.rho.pi - mu.rho
  mu.lambda.diff = mu.lambda.pi - mu.lambda
  mu.mu.diff = mu.mu.pi - mu.mu

  #define model fitting
  
  for (i in 1:3945) {
  
  #So, jags gets angry because it can't evaluate a negative prospect with an exponent without me making the
  #negative a positive; this is a problem if I'm looping over positive and negative trials
  #so I'm splitting them into two loops. This first loop is for the positive trials
  
  #compute value of certain and uncertain prospects
  valCertain[i] = 1 * (pow(certain[i], rho[id[i]]))
  valRisk[i] = riskProb[i] * (pow(nonZero[i], rho[id[i]]))
  
  #compute the model implied probability of making the risky choice, relate it back to the observed data
  probRisk[i] <- (1)/(1+exp((mu[id[i]] * -1)*(valRisk[i] - valCertain[i])))
  
  y[i] ~ dbern(probRisk[i])
  
  
  }
  
  for (i in 3946:7886) {
  
  #now loop over the loss trials
  
  valCertain[i] = 1 * (-1 * lambda[id[i]]) * (pow((-1 * certain[i]), rho[id[i]]))
  valRisk[i] = riskProb[i] * (-1 * lambda[id[i]]) * (pow((-1 * nonZero[i]), rho[id[i]]))

  probRisk[i] <- (1)/(1+exp((mu[id[i]] * -1)*(valRisk[i] - valCertain[i])))
  
  y[i] ~ dbern(probRisk[i])
  
  }


}
    ",fill = TRUE)
sink()


#inits = rep(list(list(
#  lmu.alpha = 0.7, lsigma.alpha = 1,
#  lmu.beta = 0.7, lsigma.beta = 1,
#  lmu.lambda = 0, lsigma.lambda = 0.5,
#  lmu.mu = 0, lsigma.mu = 0.5)),5)

inits = rep(list(list(
  lmu.rho = 0 ,lsigma.rho = 0.15,
  lmu.lambda = 0, lsigma.lambda = 0.25,
  lmu.mu = 0, lsigma.mu = 0.5,
  lmu.rho.pi = 0 ,lsigma.rho.pi = 0.15,
  lmu.lambda.pi = 0, lsigma.lambda.pi = 0.25,
  lmu.mu.pi = 0, lsigma.mu.pi = 0.5)),10)


ptData = list(
  y = y,
  certain = certain,
  nonZero = nonZero,
  riskProb = riskProb,
  #piStat = piStat,
  #trialType = trialType,
  id = id
)

parameters = c(
  "lmu.rho", "mu.rho", "lsigma.rho",
  "lmu.lambda", "mu.lambda", "lsigma.lambda", 
  "lmu.mu", "mu.mu", "lsigma.mu",
  "lmu.rho.pi", "mu.rho.pi", "lsigma.rho.pi",
  "lmu.lambda.pi", "mu.lambda.pi", "lsigma.lambda.pi", 
  "lmu.mu.pi", "mu.mu.pi", "lsigma.mu.pi",
  "mu.rho.diff", "mu.lambda.diff", "mu.mu.diff"
  
)

#parameters = c(
#  "alpha", "lmu.alpha", "mu.alpha", "lsigma.alpha",
#  "beta", "lmu.beta", "mu.beta", "lsigma.beta",
#  "lambda", "lmu.lambda", "mu.lambda", "lsigma.lambda", 
#  "mu", "lmu.mu", "mu.mu", "lsigma.mu"
 # 
#)


#parameters = c(
#  "alpha", "mu.phi.alpha", "mu.alpha", "sigma.phi.alpha",
#  "beta", "mu.phi.beta", "mu.beta", "sigma.phi.beta",
#  "lambda", "lmu.lambda", "mu.lambda", "lsigma.lambda", 
#  "mu", "lmu.mu", "mu.mu", "lsigma.mu"
#  
#)


proc.time()
run1 = jags(data = ptData, inits = inits, parameters.to.save = parameters, model.file="compMod_ptMod1_priorSet1-groupDiff.txt", 
            n.chains=10, n.iter=60000, n.burnin=30000, n.thin=3)
proc.time()


saveRDS(run1, file="P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/output-priorSet1-groupDiff")

run1.mcmc <- as.mcmc(run1$BUGSoutput$sims.matrix)
summary(run1.mcmc)
# ggmcmc needs to reformat the data. 
# restructure the data as a ggmcmc object
run1.ggs <- ggs(run1.mcmc)
ggmcmc(run1.ggs, file = "P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/test.pdf")
#ggmcmc(run1.ggs, file = "P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/ggmcmc-output-priorSet2-groupDiff.pdf") 


#make some nice plots
pt_priorSet1 = readRDS("P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/output-priorSet1-groupDiff")
pt_priorSet2 = readRDS("P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/output-priorSet2-groupDiff")
pt_priorSet3 = readRDS("P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/output-priorSet3-groupDiff")
pt_priorSet4 = readRDS("P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/output-priorSet4-groupDiff")

pt_priorSet1_forDensity = data.frame(samps = as.vector(pt_priorSet1$BUGSoutput$sims.matrix[,c("mu.rho.diff", "mu.lambda.diff", "mu.mu.diff")]),
                                     param = rep(c("Rho", "Lambda", "Mu"), each = 100000))

pt_priorSet2_forDensity = data.frame(samps = as.vector(pt_priorSet2$BUGSoutput$sims.matrix[,c("mu.rho.diff", "mu.lambda.diff", "mu.mu.diff")]),
                                     param = rep(c("Rho", "Lambda", "Mu"), each = 100000))

pt_priorSet3_forDensity = data.frame(samps = as.vector(pt_priorSet3$BUGSoutput$sims.matrix[,c("mu.rho.diff", "mu.lambda.diff", "mu.mu.diff")]),
                                     param = rep(c("Rho", "Lambda", "Mu"), each = 100000))

pt_priorSet4_forDensity = data.frame(samps = as.vector(pt_priorSet4$BUGSoutput$sims.matrix[,c("mu.rho.diff", "mu.lambda.diff", "mu.mu.diff")]),
                                     param = rep(c("Rho", "Lambda", "Mu"), each = 100000))

library(ggplot2); library(hrbrthemes)

ggplot(data=pt_priorSet1_forDensity, aes(x=samps, group=param, fill=param)) +
  geom_density(adjust=1.5) +
  theme_ipsum(base_size=18) +
  facet_wrap(~param) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.2, "lines"),
    axis.ticks.x=element_blank()
  ) + geom_line(data=data.frame(x=-.1:1, y = 1.5))


ggplot(data=pt_priorSet2_forDensity, aes(x=samps, group=param, fill=param)) +
  geom_density(adjust=1.5) +
  theme_ipsum(base_size=18) +
  facet_wrap(~param) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.2, "lines"),
    axis.ticks.x=element_blank()
  )

ggplot(data=pt_priorSet3_forDensity, aes(x=samps, group=param, fill=param)) +
  geom_density(adjust=1.5) +
  theme_ipsum(base_size=18) +
  facet_wrap(~param) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.2, "lines"),
    axis.ticks.x=element_blank()
  )

ggplot(data=pt_priorSet4_forDensity, aes(x=samps, group=param, fill=param)) +
  geom_density(adjust=1.5) +
  theme_ipsum(base_size=18) +
  facet_wrap(~param) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.2, "lines"),
    axis.ticks.x=element_blank()
  )

compute_FXsize_from_post = function(diff_vec, sigma1_vec, sigma2_vec) {
  
  out = diff_vec / sqrt((sigma1_vec+sigma2_vec)/2)
  
  return(out)
  
}



lambdaSet1_d = compute_FXsize_from_post(pt_priorSet1$BUGSoutput$sims.matrix[,"mu.lambda.diff"], exp(pt_priorSet1$BUGSoutput$sims.matrix[,"lsigma.lambda"]), exp(pt_priorSet1$BUGSoutput$sims.matrix[,"lsigma.lambda.pi"]))
lambdaSet2_d = compute_FXsize_from_post(pt_priorSet2$BUGSoutput$sims.matrix[,"mu.lambda.diff"], exp(pt_priorSet2$BUGSoutput$sims.matrix[,"lsigma.lambda"]), exp(pt_priorSet2$BUGSoutput$sims.matrix[,"lsigma.lambda.pi"]))
lambdaSet3_d = compute_FXsize_from_post(pt_priorSet3$BUGSoutput$sims.matrix[,"mu.lambda.diff"], exp(pt_priorSet3$BUGSoutput$sims.matrix[,"lsigma.lambda"]), exp(pt_priorSet3$BUGSoutput$sims.matrix[,"lsigma.lambda.pi"]))
lambdaSet4_d = compute_FXsize_from_post(pt_priorSet4$BUGSoutput$sims.matrix[,"mu.lambda.diff"], exp(pt_priorSet4$BUGSoutput$sims.matrix[,"lsigma.lambda"]), exp(pt_priorSet4$BUGSoutput$sims.matrix[,"lsigma.lambda.pi"]))

rhoSet1_d = compute_FXsize_from_post(pt_priorSet1$BUGSoutput$sims.matrix[,"mu.rho.diff"], exp(pt_priorSet1$BUGSoutput$sims.matrix[,"lsigma.rho"]), exp(pt_priorSet1$BUGSoutput$sims.matrix[,"lsigma.rho.pi"]))
rhoSet2_d = compute_FXsize_from_post(pt_priorSet2$BUGSoutput$sims.matrix[,"mu.rho.diff"], exp(pt_priorSet2$BUGSoutput$sims.matrix[,"lsigma.rho"]), exp(pt_priorSet2$BUGSoutput$sims.matrix[,"lsigma.rho.pi"]))
rhoSet3_d = compute_FXsize_from_post(pt_priorSet3$BUGSoutput$sims.matrix[,"mu.rho.diff"], exp(pt_priorSet3$BUGSoutput$sims.matrix[,"lsigma.rho"]), exp(pt_priorSet3$BUGSoutput$sims.matrix[,"lsigma.rho.pi"]))
rhoSet4_d = compute_FXsize_from_post(pt_priorSet4$BUGSoutput$sims.matrix[,"mu.rho.diff"], exp(pt_priorSet4$BUGSoutput$sims.matrix[,"lsigma.rho"]), exp(pt_priorSet4$BUGSoutput$sims.matrix[,"lsigma.rho.pi"]))

muSet1_d = compute_FXsize_from_post(pt_priorSet1$BUGSoutput$sims.matrix[,"mu.mu.diff"], exp(pt_priorSet1$BUGSoutput$sims.matrix[,"lsigma.mu"]), exp(pt_priorSet1$BUGSoutput$sims.matrix[,"lsigma.mu.pi"]))
muSet2_d = compute_FXsize_from_post(pt_priorSet2$BUGSoutput$sims.matrix[,"mu.mu.diff"], exp(pt_priorSet2$BUGSoutput$sims.matrix[,"lsigma.mu"]), exp(pt_priorSet2$BUGSoutput$sims.matrix[,"lsigma.mu.pi"]))
muSet3_d = compute_FXsize_from_post(pt_priorSet3$BUGSoutput$sims.matrix[,"mu.mu.diff"], exp(pt_priorSet3$BUGSoutput$sims.matrix[,"lsigma.mu"]), exp(pt_priorSet3$BUGSoutput$sims.matrix[,"lsigma.mu.pi"]))
muSet4_d = compute_FXsize_from_post(pt_priorSet4$BUGSoutput$sims.matrix[,"mu.mu.diff"], exp(pt_priorSet4$BUGSoutput$sims.matrix[,"lsigma.mu"]), exp(pt_priorSet4$BUGSoutput$sims.matrix[,"lsigma.mu.pi"]))

library(bayestestR)

median(lambdaSet1_d)
median(lambdaSet2_d)
median(lambdaSet3_d)
median(lambdaSet4_d)

hdi(lambdaSet1_d, ci = 0.89, verbose = TRUE)
hdi(lambdaSet2_d, ci = 0.89, verbose = TRUE)
hdi(lambdaSet3_d, ci = 0.89, verbose = TRUE)
hdi(lambdaSet4_d, ci = 0.89, verbose = TRUE)

median(rhoSet1_d)
median(rhoSet2_d)
median(rhoSet3_d)
median(rhoSet4_d)

hdi(rhoSet1_d, ci = 0.89, verbose = TRUE)
hdi(rhoSet2_d, ci = 0.89, verbose = TRUE)
hdi(rhoSet3_d, ci = 0.89, verbose = TRUE)
hdi(rhoSet4_d, ci = 0.89, verbose = TRUE)

median(muSet1_d)
median(muSet2_d)
median(muSet3_d)
median(muSet4_d)

hdi(muSet1_d, ci = 0.89, verbose = TRUE)
hdi(muSet2_d, ci = 0.89, verbose = TRUE)
hdi(muSet3_d, ci = 0.89, verbose = TRUE)
hdi(muSet4_d, ci = 0.89, verbose = TRUE)

pt_d_df = data.frame(lambdaSet1_d, lambdaSet2_d, lambdaSet3_d, lambdaSet4_d,
                     rhoSet1_d, rhoSet2_d, rhoSet3_d, rhoSet4_d,
                     muSet1_d, muSet2_d, muSet3_d, muSet4_d)

posterior_ggdensity_plot = function(df, posterior, color, hdi, rope) {
  
  out = ggplot(df, aes(x=posterior)) + 
    geom_density(fill=color, alpha = 0.5) + theme_classic(base_size=21) +
    labs(x="Effect Size", y = "Density") + 
    geom_segment(x = hdi[1], y = .15, xend = hdi[2], yend = .15, color = "black", size=1.5) +
    geom_vline(xintercept = rope[1], color = "black", size=.75, linetype="dashed") +
    geom_vline(xintercept = rope[2], color = "black", size=.75, linetype="dashed")
  
  return(out)
  
}

lambdaSet1_d_hist = ggplot(pt_d_df, aes(x=lambdaSet1_d)) + 
  geom_density(fill="#5c3b52", alpha = 0.5) + theme_classic(base_size=18) +
  labs(x="Effect Size", y = "Density") + 
  geom_segment(x = -.17, y = .15, xend = .08, yend = .15, color = "black", size=1.5) +
  geom_vline(xintercept = -0.1, color = "black", size=.75, linetype="dashed") +
  geom_vline(xintercept = 0.1, color = "black", size=.75, linetype="dashed")

rhoSet1_d_hist = ggplot(pt_d_df, aes(x=rhoSet1_d)) + 
  geom_density(fill="#d07d6b", alpha = 0.5) + theme_classic(base_size=18) +
  labs(x="Effect Size", y = "Density") + 
  geom_segment(x = -.25, y = .15, xend = .08, yend = .15, color = "black", size=1.5) +
  geom_vline(xintercept = -0.1, color = "black", size=.75, linetype="dashed") +
  geom_vline(xintercept = 0.1, color = "black", size=.75, linetype="dashed")

muSet1_d_hist = ggplot(pt_d_df, aes(x=muSet1_d)) + 
  geom_density(fill="#a66d44", alpha = 0.5) + theme_classic(base_size=18) +
  labs(x="Effect Size", y = "Density") + 
  geom_segment(x = -.03, y = .15, xend = .03, yend = .15, color = "black", size=1.5) +
  geom_vline(xintercept = -0.1, color = "black", size=.75, linetype="dashed") +
  geom_vline(xintercept = 0.1, color = "black", size=.75, linetype="dashed")

lambdaSet1_d_hist = posterior_ggdensity_plot(df = pt_d_df, 
                                             posterior = lambdaSet1_d,
                                             color = "#5c3b52",
                                             hdi = c(-.17, .08),
                                             rope = c(-.1, .1))
rhoSet1_d_hist = posterior_ggdensity_plot(df = pt_d_df, 
                                             posterior = rhoSet1_d,
                                             color = "#d07d6b",
                                             hdi = c(-.25, .08),
                                             rope = c(-.1, .1))
muSet1_d_hist = posterior_ggdensity_plot(df = pt_d_df, 
                                             posterior = muSet1_d,
                                             color = "#a66d44",
                                             hdi = c(-.03, .03),
                                             rope = c(-.1, .1))


lambdaSet2_d_hist = posterior_ggdensity_plot(df = pt_d_df, 
                                             posterior = lambdaSet2_d,
                                             color = "#5c3b52",
                                             hdi = c(-.23, .06),
                                             rope = c(-.1, .1))
rhoSet2_d_hist = posterior_ggdensity_plot(df = pt_d_df, 
                                          posterior = rhoSet2_d,
                                          color = "#d07d6b",
                                          hdi = c(0.14, .47),
                                          rope = c(-.1, .1))
muSet2_d_hist = posterior_ggdensity_plot(df = pt_d_df, 
                                         posterior = muSet2_d,
                                         color = "#a66d44",
                                         hdi = c(-.11, -.03),
                                         rope = c(-.1, .1))


lambdaSet3_d_hist = posterior_ggdensity_plot(df = pt_d_df, 
                                             posterior = lambdaSet3_d,
                                             color = "#5c3b52",
                                             hdi = c(-.19, .08),
                                             rope = c(-.1, .1))
rhoSet3_d_hist = posterior_ggdensity_plot(df = pt_d_df, 
                                          posterior = rhoSet3_d,
                                          color = "#d07d6b",
                                          hdi = c(-0.23, .11),
                                          rope = c(-.1, .1))
muSet3_d_hist = posterior_ggdensity_plot(df = pt_d_df, 
                                         posterior = muSet3_d,
                                         color = "#a66d44",
                                         hdi = c(-.03, .03),
                                         rope = c(-.1, .1))


lambdaSet4_d_hist = posterior_ggdensity_plot(df = pt_d_df, 
                                             posterior = lambdaSet4_d,
                                             color = "#5c3b52",
                                             hdi = c(-.14, .05),
                                             rope = c(-.1, .1))
rhoSet4_d_hist = posterior_ggdensity_plot(df = pt_d_df, 
                                          posterior = rhoSet4_d,
                                          color = "#d07d6b",
                                          hdi = c(-0.28, .08),
                                          rope = c(-.1, .1))
muSet4_d_hist = posterior_ggdensity_plot(df = pt_d_df, 
                                         posterior = muSet4_d,
                                         color = "#a66d44",
                                         hdi = c(-.03, .03),
                                         rope = c(-.1, .1))


png("pt_post_plots.png", width = 20, height = 16, units = "in", res = 900)
ggarrange(lambdaSet1_d_hist, rhoSet1_d_hist, muSet1_d_hist,
          lambdaSet2_d_hist, rhoSet2_d_hist, muSet2_d_hist,
          lambdaSet3_d_hist, rhoSet3_d_hist, muSet3_d_hist,
          lambdaSet4_d_hist, rhoSet4_d_hist, muSet4_d_hist,
          ncol = 3, nrow = 4)
dev.off()


run1.mcmc <- as.mcmc(run1$BUGSoutput$sims.matrix)
summary(run1.mcmc)
# ggmcmc needs to reformat the data. 
# restructure the data as a ggmcmc object
run1.ggs <- ggs(run1.mcmc)
ggmcmc(run1.ggs, file = "P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/test.pdf")

ggmcmc(ggs(as.mcmc(pt_priorSet1)), file = "P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/pt_priorSet1_diag.pdf")
ggmcmc(ggs(as.mcmc(pt_priorSet2)), file = "P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/pt_priorSet2_diag.pdf")
ggmcmc(ggs(as.mcmc(pt_priorSet3)), file = "P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/pt_priorSet3_diag.pdf")
ggmcmc(ggs(as.mcmc(pt_priorSet4)), file = "P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/pt_priorSet4_diag.pdf")


#########
## Target Model
#########

#load data for target model
subList = Sys.glob("P:/Social_Behavior_(SB_study)/Data/Behavioral_data/SB*/wave1/Lab_session/Clean/*CupsSelf.csv") #Get path to each subjects' cups data

cupsDat = getCupsDatCompMod(subList, mod = "target") #Compile Cups data (Level 1, i.e., trial level), uses helper func
cupsDat = cupsDat[-1,] #The first row is a bunch of 999 placeholders because I coded this like a doofus the first time around

L2Dat = getL2Dat("P:/Social_Behavior_(SB_study)/Data/Study_logs/SB_group_list.csv") #Load the level two (i.e., between subjects) data, uses helper func (though it's not really necessary)
L2Dat = L2Dat[L2Dat$ID %in% unique(cupsDat$ID),]

L2Dat = rbind(L2Dat[L2Dat$PI==0,], L2Dat[L2Dat$PI==1,])

L2Dat$ID_num = 1:length(L2Dat$ID)

L2Dat = select(L2Dat, ID = ID, ID_num, Age, PI, Sex); L2Dat = arrange(L2Dat, ID)

allDat = inner_join(cupsDat,L2Dat) #merge the trial- and subject-level data

allDat = allDat[order(allDat$ID_num),]

allDat = rbind(allDat[allDat$TrialType==1,], allDat[allDat$TrialType==0,])


y = allDat$Decision
evsdRatio = allDat$EV / allDat$SD
#riskProb = allDat$riskProb
trialType = allDat$TrialType
#piStat = L2Dat$PI
id = allDat$ID_num

sink("compMod_tgMod1_priorSet2-groupDiff.txt")
cat("
# tg model 1  
model
{

  for (j in 1:82) {
  
  #begin by defining participant j has unique parameter values for
  #2 parameters: tau (target threshold) and mu (choice consistency)
  
  #these are just the initial transformations
  taug[j] <- exp(ltaug[j])
  taul[j] <- exp(ltaul[j])
	mu[j]   <- exp(lmu[j]) 
  
  #now we'll say that each subject's parameters were drawn from normal distributions
  #centered around the group mean 
  
  ltaug[j] ~ dnorm(lmu.taug,ltau.taug)
  ltaul[j] ~ dnorm(lmu.taul,ltau.taul)
  lmu[j] ~ dnorm(lmu.mu, ltau.mu)
  
  }
  
  for (j in 83:144) {
  
  #these are just the initial transformations
  taug[j] <- exp(ltaug.pi[j])
  taul[j] <- exp(ltaul.pi[j])
	mu[j]   <- exp(lmu.pi[j]) 
  
  #now we'll say that each subject's parameters were drawn from normal distributions
  #centered around the group mean 
  
  ltaug.pi[j] ~ dnorm(lmu.taug,ltau.taug.pi)
  ltaul.pi[j] ~ dnorm(lmu.taul,ltau.taul.pi)
  lmu.pi[j] ~ dnorm(lmu.mu, ltau.mu.pi)
  
  }
  

  ##Priors for hyper distributions 
  
  #for comparison youths
  
  #for tau (target threshold)
  lmu.taug ~ dunif(-4,1) 
  ltau.taug = pow(lsigma.taug,-2)
  lsigma.taug ~ dunif(0,1)
  
  lmu.taul ~ dunif(-4,1) 
  ltau.taul = pow(lsigma.taul,-2)
  lsigma.taul ~ dunif(0,1)
  
  #for mu (choice consistency)
  lmu.mu ~ dunif(-1, 1.61)
  ltau.mu = pow(lsigma.mu,-2)
  lsigma.mu ~ dunif(0,2)
  
  #for PI youths
  
  #for rho (target threshold)
  lmu.taug.pi ~ dunif(-4,1) 
  ltau.taug.pi = pow(lsigma.taug.pi,-2)
  lsigma.taug.pi ~ dunif(0,1)
  
  lmu.taul.pi ~ dunif(-4,1) 
  ltau.taul.pi = pow(lsigma.taul.pi,-2)
  lsigma.taul.pi ~ dunif(0,1)
  
  #for mu (choice consistency)
  lmu.mu.pi ~ dunif(-1, 1.61)
  ltau.mu.pi = pow(lsigma.mu.pi,-2)
  lsigma.mu.pi ~ dunif(0,2)
  
  # To obtain the mean of the hyper distribution on the wanted scale:
  
  #comparison
	mu.taug <- exp(lmu.taug) 
	mu.taul <- exp(lmu.taul) 
  mu.mu   <- exp(lmu.mu)
  
  #PI and difference
  mu.taug.pi <- exp(lmu.taug.pi) 
	mu.taul.pi <- exp(lmu.taul.pi) 
  mu.mu.pi   <- exp(lmu.mu.pi)
  
  mu.taug.diff = mu.taug.pi - mu.taug
  mu.taul.diff = mu.taul.pi - mu.taul
  mu.mu.diff = mu.mu.pi - mu.mu


  #define model fitting
  
  for (i in 1:3889) {
  
  #So, jags gets angry because it can't evaluate a negative prospect with an exponent without me making the
  #negative a positive; this is a problem if I'm looping over positive and negative trials
  #so I'm splitting them into two loops. This first loop is for the positive trials
  
  #compute the difference between observed evsdratio and tau
  delta[i] = evsdRatio[i] - taug[id[i]]
  
  #compute the model implied probability of making the risky choice, relate it back to the observed data
  probRisk[i] <- (1)/(1+exp((mu[id[i]] * -1)*(delta[i])))
  
  y[i] ~ dbern(probRisk[i])
  
  
  }
  
  for (i in 3890:7740) {
  
  #now loop over the loss trials
  
  #compute the difference between observed evsdratio and tau
  delta[i] = evsdRatio[i] - taul[id[i]]
  
  #compute the model implied probability of making the risky choice, relate it back to the observed data
  probRisk[i] <- (1)/(1+exp((mu[id[i]] * -1)*(delta[i])))
  
  y[i] ~ dbern(probRisk[i])
  
  }


}
    ",fill = TRUE)
sink()


#inits = rep(list(list(
#  lmu.alpha = 0.7, lsigma.alpha = 1,
#  lmu.beta = 0.7, lsigma.beta = 1,
#  lmu.lambda = 0, lsigma.lambda = 0.5,
#  lmu.mu = 0, lsigma.mu = 0.5)),5)

inits = rep(list(list(
  lmu.taug = 0 ,lsigma.taug = 0.15,
  lmu.taul = 0, lsigma.taul = 0.25,
  lmu.mu = 0, lsigma.mu = 0.5,
  lmu.taug.pi = 0 ,lsigma.taug.pi = 0.15,
  lmu.taul.pi = 0, lsigma.taul.pi = 0.25,
  lmu.mu.pi = 0, lsigma.mu.pi = 0.5)),10)


tgData = list(
  y = y,
  evsdRatio = evsdRatio,
  #piStat = piStat,
  #trialType = trialType,
  id = id
)

parameters = c(
  "lmu.taug", "mu.taug", "lsigma.taug",
  "lmu.taul", "mu.taul", "lsigma.taul", 
  "lmu.mu", "mu.mu", "lsigma.mu",
  "lmu.taug.pi", "mu.taug.pi", "lsigma.taug.pi",
  "lmu.taul.pi", "mu.taul.pi", "lsigma.taul.pi", 
  "lmu.mu.pi", "mu.mu.pi", "lsigma.mu.pi",
  "mu.taug.diff", "mu.taul.diff", "mu.mu.diff"
  
)

#parameters = c(
#  "alpha", "lmu.alpha", "mu.alpha", "lsigma.alpha",
#  "beta", "lmu.beta", "mu.beta", "lsigma.beta",
#  "lambda", "lmu.lambda", "mu.lambda", "lsigma.lambda", 
#  "mu", "lmu.mu", "mu.mu", "lsigma.mu"
# 
#)


#parameters = c(
#  "alpha", "mu.phi.alpha", "mu.alpha", "sigma.phi.alpha",
#  "beta", "mu.phi.beta", "mu.beta", "sigma.phi.beta",
#  "lambda", "lmu.lambda", "mu.lambda", "lsigma.lambda", 
#  "mu", "lmu.mu", "mu.mu", "lsigma.mu"
#  
#)


proc.time()
run1 = jags(data = tgData, inits = inits, parameters.to.save = parameters, model.file="compMod_tgMod1_priorSet2-groupDiff.txt", 
            n.chains=10, n.iter=60000, n.burnin=30000, n.thin=3)
proc.time()

saveRDS(run1, file="P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/output-tg-priorSet2-groupDiff")


run1.mcmc <- as.mcmc(run1)
summary(run1.mcmc)
# ggmcmc needs to reformat the data. 
# restructure the data as a ggmcmc object
run1.ggs <- ggs(run1.mcmc)
ggmcmc(run1.ggs, file = "P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/ggmcmc-output-tg-priorSet6-groupDiff.pdf") 


#make some plots and extract information for tables
tg_priorSet1 = readRDS("P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/output-tg-priorSet1-groupDiff")
tg_priorSet2 = readRDS("P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/output-tg-priorSet2-groupDiff")

taugSet1_d = compute_FXsize_from_post(tg_priorSet1$BUGSoutput$sims.matrix[,"mu.taug.diff"], exp(tg_priorSet1$BUGSoutput$sims.matrix[,"lsigma.taug"]), exp(tg_priorSet1$BUGSoutput$sims.matrix[,"lsigma.taug.pi"]))
taugSet2_d = compute_FXsize_from_post(tg_priorSet2$BUGSoutput$sims.matrix[,"mu.taug.diff"], exp(tg_priorSet2$BUGSoutput$sims.matrix[,"lsigma.taug"]), exp(tg_priorSet2$BUGSoutput$sims.matrix[,"lsigma.taug.pi"]))

taulSet1_d = compute_FXsize_from_post(tg_priorSet1$BUGSoutput$sims.matrix[,"mu.taul.diff"], exp(tg_priorSet1$BUGSoutput$sims.matrix[,"lsigma.taul"]), exp(tg_priorSet1$BUGSoutput$sims.matrix[,"lsigma.taul.pi"]))
taulSet2_d = compute_FXsize_from_post(tg_priorSet2$BUGSoutput$sims.matrix[,"mu.taul.diff"], exp(tg_priorSet2$BUGSoutput$sims.matrix[,"lsigma.taul"]), exp(tg_priorSet2$BUGSoutput$sims.matrix[,"lsigma.taul.pi"]))

muSet1_d = compute_FXsize_from_post(tg_priorSet1$BUGSoutput$sims.matrix[,"mu.mu.diff"], exp(tg_priorSet1$BUGSoutput$sims.matrix[,"lsigma.mu"]), exp(tg_priorSet1$BUGSoutput$sims.matrix[,"lsigma.mu.pi"]))
muSet2_d = compute_FXsize_from_post(tg_priorSet2$BUGSoutput$sims.matrix[,"mu.mu.diff"], exp(tg_priorSet2$BUGSoutput$sims.matrix[,"lsigma.mu"]), exp(tg_priorSet2$BUGSoutput$sims.matrix[,"lsigma.mu.pi"]))

library(bayestestR)

median(taugSet1_d)
median(taugSet2_d)

hdi(taugSet1_d, ci = 0.89, verbose = TRUE)
hdi(taugSet2_d, ci = 0.89, verbose = TRUE)

median(taulSet1_d)
median(taulSet2_d)

hdi(taulSet1_d, ci = 0.89, verbose = TRUE)
hdi(taulSet2_d, ci = 0.89, verbose = TRUE)

median(muSet1_d)
median(muSet2_d)

hdi(muSet1_d, ci = 0.89, verbose = TRUE)
hdi(muSet2_d, ci = 0.89, verbose = TRUE)

##code for posterior densities

#put posteriors into dataframe
tm_d_df = data.frame(taugSet1_d, taugSet2_d, 
                     taulSet1_d, taulSet2_d, 
                     muSet1_d, muSet2_d)

#create histograms/density plots
taugSet1_d_hist = posterior_ggdensity_plot(df = tm_d_df, 
                                             posterior = taugSet1_d,
                                             color = "#ceaab7",
                                             hdi = c(-.17, 1.59),
                                             rope = c(-.1, .1))
taulSet1_d_hist = posterior_ggdensity_plot(df = tm_d_df, 
                                          posterior = taulSet1_d,
                                          color = "#edddd4",
                                          hdi = c(-.02, 2.12),
                                          rope = c(-.1, .1))
muSet1_d_hist = posterior_ggdensity_plot(df = tm_d_df, 
                                         posterior = muSet1_d,
                                         color = "#3c748c",
                                         hdi = c(-.02, 1.98),
                                         rope = c(-.1, .1))



taugSet2_d_hist = posterior_ggdensity_plot(df = tm_d_df, 
                                           posterior = taugSet2_d,
                                           color = "#ceaab7",
                                           hdi = c(-.17, 1.59),
                                           rope = c(-.1, .1))
taulSet2_d_hist = posterior_ggdensity_plot(df = tm_d_df, 
                                           posterior = taulSet2_d,
                                           color = "#edddd4",
                                           hdi = c(-.02, 2.12),
                                           rope = c(-.1, .1))
muSet2_d_hist = posterior_ggdensity_plot(df = tm_d_df, 
                                         posterior = muSet2_d,
                                         color = "#3c748c",
                                         hdi = c(-.02, 1.98),
                                         rope = c(-.1, .1))

#arrange and save
png("tm_post_plots.png", width = 18, height = 14, units = "in", res = 900)
ggarrange(taugSet1_d_hist, taulSet1_d_hist, muSet1_d_hist,
          taugSet2_d_hist, taulSet2_d_hist, muSet2_d_hist,
          ncol = 3, nrow = 2)
dev.off()

#now to plot diagnostics

#first put the data from BUGS output array into a dataframe
tg_priorSet1_df = data.frame(sample = 1:10000, tg_priorSet1$BUGSoutput$sims.array[,,c('mu.taug', 'mu.taug.pi', 'mu.taug.diff',
                                                                    'mu.taul', 'mu.taul.pi', 'mu.taul.diff',
                                                                    'mu.mu', 'mu.mu.pi', 'mu.mu.diff')])

tg_priorSet1_df_long = reshape(tg_priorSet1_df, direction='long', 
        varying=names(tg_priorSet1_df)[-1], 
        timevar='chain',
        times=c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10'),
        v.names=c('mu.taug', 'mu.taug.pi', 'mu.taug.diff',
                  'mu.taul', 'mu.taul.pi', 'mu.taul.diff',
                  'mu.mu', 'mu.mu.pi', 'mu.mu.diff'),
        idvar='sample')

ggplot(data = tg_priorSet1_df_long) + 
  geom_line(mapping = aes(x = sample, y = mu.taul))#, 
                          #group = chain, 
                          #color = chain)
  #)

ggmcmc(ggs(as.mcmc(tg_priorSet1)), file = "P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/tg_priorSet1_diag.pdf")
ggmcmc(ggs(as.mcmc(tg_priorSet2)), file = "P:/Social_Behavior_(SB_study)/SB_Scripts/wave1/Behavioral_tasks/SB_Cups_Self/SB_Cups_Reward-Risk_Sensitivity/tg_priorSet2_diag.pdf")


