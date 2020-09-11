###########################################################################
######################## Analysis code #################################### 
########### JF Guassi Moreira (jguassimoreira[at]ucla[dot]edu) ############
####################### Spring/Summer 2020 ################################

#########
## Load the libraries we'll need
#########
library(readr)
library(dplyr)
library(tidyr)
library(lmerTest)
library(pracma)
library(wesanderson)
library(RColorBrewer)
library(ggplot2)
library(plyr)
library(gridExtra)

#########
## Import helper functions, set path to data
#########
scriptPath = file.path("~", "Documents", "Scripts", "R", "helperFuncs") #path to helper funcs
file.sources = Sys.glob(file.path(sprintf("%s", scriptPath), "*.R")) #Construct paths to helper funcs, save in vector
sapply(file.sources,FUN=source) #use sapply to source each file path in the helper func vector
sbPath = file.path("P:", "Social_Behavior_(SB_study)", "Data", "Behavioral_data") #set path to data on SANDLab server


#########
## Load in the data
#########

## Note: If you are downloading this script from github/osf, the data have already compiled and placed into a .csv on the osf project page, so you could skip this chunk and load the data in yourself

subList = Sys.glob("P:/Social_Behavior_(SB_study)/Data/Behavioral_data/SB*/wave1/Lab_session/Clean/*_CupsSelf.csv") #Get path to each subjects' cups data

cupsDat = getCupsDat(subList) #Compile Cups data (Level 1, i.e., trial level), uses helper func
cupsDat = cupsDat[-1,] #The first row is a bunch of 999 placeholders because I coded this like a doofus the first time around

L2Dat = getL2Dat("P:/Social_Behavior_(SB_study)/Data/Study_logs/SB_group_list.csv") #Load the level two (i.e., between subjects) data, uses helper func (though it's not really necessary)

allDat = inner_join(L2Dat, cupsDat) #merge the trial- and subject-level data

#########
## Part 1 of analysis strategy -- multilevel logistic regression
#########

## First we need to do some grand mean centering
allDat$EV = allDat$EV - mean(allDat$EV); allDat$SD = allDat$SD - mean(allDat$SD) #grand mean center EV (reward) and then SD (risk)
allDat$Age = allDat$Age - mean(allDat$Age) #grand mean center age
allDat$IQ_percentile = allDat$IQ_percentile - mean(allDat$IQ_percentile) #grand mean center IQ scores (obtained from the verbal and matrix reasoning subtests of the WASI-II)
allDat$DOSPERT = allDat$DOSPERT - mean(allDat$DOSPERT) #grand mean center DOSPERT (domain specific risk taking scale), i.e., self-reports of 'real-world' risk-taking
allDat$TrialNum = allDat$TrialNum - mean(allDat$TrialNum) #grand mean center trial number for an additional ancillary/supplementary analysis


## run the main models reported in Table 2
mod1 = glmer(Decision ~ EV + SD + TrialType + (1+EV+SD+TrialType|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa")) #Model 1 from Table 2
mod2 = glmer(Decision ~ EV + SD + TrialType + PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa")) #Model 2 from Table 2
mod3 = glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa")) #Model 3 from Table 2
mod4 = glmer(Decision ~ EV*PI + EV*Age + SD*PI + SD*Age + TrialType*PI + TrialType*Age + Sex + IQ_percentile + (1+EV|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa")) #Model 4 from Table 2

## Additional/ancillary models 
mod5 = glmer(Decision ~ EV*PI + SD*PI*TrialType  + Age + Sex + IQ_percentile + (1+SD|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa")) #additional model testing 3 way interaction between group(PI), age, and risk (SD) -- as noted in the manuscript, it does not converge
mod3dos = glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + DOSPERT + (1+EV+SD|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa")) #additional model testing whether DOSPERT scores were related to behavior on the Cups Task
mod3tn = glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + TrialNum + (1+EV|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa")) #Follow up model 3 by adding trial number (model won't converge with both EV and SD as random effects)


## Recentering EV and re-running model 3 to get a sense of where along the range of EV values there were significant group differences in risk-taking likelihoods
rcDat_m5 = allDat; rcDat_m5$EV = rcDat_m5$EV - (-5) #recenter at minus 5 EV
rcDat_m2 = allDat; rcDat_m2$EV = rcDat_m2$EV - (-2) #recenter at minus 2 EV
rcDat_p2 = allDat; rcDat_p2$EV = rcDat_p2$EV - 2 #recenter at plus 2 EV
rcDat_p5 = allDat; rcDat_p5$EV = rcDat_p5$EV - 5 #recenter at plus 5 EV

mod3i_rcM5 = glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = rcDat_m5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 150000))) #EV centered at minus 5 
mod3i_rcM2 = glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = rcDat_m2, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 150000))) #EV centered at minus 2 
mod3i_rcP2 = glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = rcDat_p2, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 150000))) #EV centered at plus 2 
mod3i_rcP5 = glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = rcDat_p5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 150000))) #EV centered at plus 5 

## Running gain and loss trials separately to unpack the main effect of SD
gainDat = inner_join(L2Dat, cupsDat); gainDat = gainDat[gainDat$EV > 0,] #gain trials
lossDat = inner_join(L2Dat, cupsDat); lossDat = lossDat[lossDat$EV < 0,] #loss trials

#repeat grand mean centering within trial type
gainDat$SD = gainDat$SD - mean(gainDat$SD); gainDat$EV = gainDat$EV - mean(gainDat$EV)
gainDat$Age = gainDat$Age - mean(gainDat$Age)
gainDat$IQ_percentile = gainDat$IQ_percentile - mean(gainDat$IQ_percentile)

lossDat$SD = lossDat$SD - mean(lossDat$SD); lossDat$EV = lossDat$EV - mean(lossDat$EV)
lossDat$Age = lossDat$Age - mean(lossDat$Age)
lossDat$IQ_percentile = lossDat$IQ_percentile - mean(lossDat$IQ_percentile)

#ran up to mod 3 separately, discovered the apparent differences in SD*PI interactions between gain and loss trials (don't over interpret this, since we couldn't formally evaluate for 3 way interaction in model 5)
gainMod1 = glmer(Decision ~ EV + SD + (1+EV+SD|ID), data = gainDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
lossMod1 = glmer(Decision ~ EV + SD + (1+EV+SD|ID), data = lossDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
gainMod2 = glmer(Decision ~ EV + SD + PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = gainDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
lossMod2 = glmer(Decision ~ EV + SD + PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = lossDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
gainMod3 = glmer(Decision ~ EV*PI + SD*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = gainDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
lossMod3 = glmer(Decision ~ EV*PI + SD*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = lossDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Remove outlier (became apparent when plotting random effects in Figure 3)
allDat_noOutlier = inner_join(L2Dat[-(L2Dat$ID=='8366515152'),], cupsDat)

#grand mean center again
allDat_noOutlier$SD = allDat_noOutlier$SD - mean(allDat_noOutlier$SD); allDat_noOutlier$EV = allDat_noOutlier$EV - mean(allDat_noOutlier$EV)
allDat_noOutlier$Age = allDat_noOutlier$Age - mean(allDat_noOutlier$Age)
allDat_noOutlier$IQ_percentile = allDat_noOutlier$IQ_percentile - mean(allDat_noOutlier$IQ_percentile)

#run the models
mod1 = glmer(Decision ~ EV + SD + TrialType + (1+EV+SD+TrialType|ID), data = allDat_noOutlier, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mod2 = glmer(Decision ~ EV + SD + TrialType + PI + Age + Sex + IQ_percentile + (1+EV+SD+TrialType|ID), data = allDat_noOutlier, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mod3 = glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + (1+EV+SD+TrialType|ID), data = allDat_noOutlier, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mod4 = glmer(Decision ~ EV*PI + EV*Age + SD*PI + SD*Age + TrialType*PI + TrialType*Age + Sex + IQ_percentile +  (1+EV+SD|ID), data = allDat_noOutlier, family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Do a few supplementary/ancillary analyses with only the PI subjects
piDat = inner_join(L2Dat, cupsDat); piDat = piDat[piDat$PI == 1,]
piDat$EV = piDat$EV - mean(piDat$EV); piDat$SD = piDat$SD - mean(piDat$SD) #grand mean center EV (reward) and then SD (risk)
piDat$Age = piDat$Age - mean(piDat$Age) #grand mean center age
piDat$IQ_percentile = piDat$IQ_percentile - mean(piDat$IQ_percentile) #grand mean center IQ scores (obtained from the verbal and matrix reasoning subtests of the WASI-II)
piDat$IPPA_Parent = piDat$IPPA_Parent - mean(piDat$IPPA_Parent, na.rm = T) #grand mean center IPPA parent variable, also needed for ancillary/supplementary analysis
piDat$timePI = piDat$timePI - mean(piDat$timePI, na.rm = T) #grand mean center time spent in institutionalization variable for PI participants

piModIPPA = glmer(Decision ~ EV*IPPA_Parent + SD*IPPA_Parent + TrialType + Age + Sex + IQ_percentile + (1+EV|ID), data = piDat, family = binomial, control = glmerControl(optimizer = "bobyqa")) #does relationship quality buffer effect of EV and SD? 
piModtimePI = glmer(Decision ~ EV*timePI + SD*timePI + TrialType + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = piDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))  #does time spent in institutionalization buffer effect of EV and SD? 

#########
## Plot results from this part of the analysis strategy
#########

## Start with model implied probabilities in a Cleveland plot (i.e., dotchart) broken down by trial type, EV, and group 
posEV = c(rep(linspace(.8, 5, 5), each = 2)) #positive EVs
negEV = c(rep(linspace(-.8, -5, 5), each = 2)) #negative EVs
age=0 #create the plot for someone at the mean age
grpPI = c(rep(c(1,0), 5)) #create a vector for PI and comparison individuals

#create model implied logits -- these coefficients are taken from model 3 -- effect of gender is averaged, and WASI-II score is assumed to be zero (i.e., at grand mean)
posLogits = 0.751 + (0.674*posEV) + (0.141*0) + (-1.290*1) + (-0.656*grpPI) + (-0.043*age) + (0*-0.005) + (-.335*posEV*grpPI) + (-.077*0*grpPI) + (1.266*1*grpPI) + 0.0993085 #this last bit here is the avg of gender effect
posProb = exp(posLogits) / (1 + exp(posLogits)) #convert to probabilities

#create model implied logits -- these coefficients are taken from model 3 -- effect of gender is averaged, and WASI-II score is assumed to be zero (i.e., at grand mean)
negLogits = 0.751 + (0.674*negEV) + (0.141*0) + (-1.290*0) + (-0.656*grpPI) + (-0.043*age) + (0*-0.005) + (-.335*negEV*grpPI) + (-.077*0*grpPI) + (1.266*0*grpPI) + 0.0993085 #this last bit here is the avg of gender effect
negProb = exp(negLogits) / (1 + exp(negLogits)) #convert to probabilities

posGroups = (rep(c("EV: +0.80", "EV: +1.85", "EV: +2.90", "EV: +3.95", "EV: +5.00"), each = 2)) #Define positive EVs (these are somewhat arbitrary, but I wanted evenly spaced bins from .8 to 5)
negGroups = (rep(c("EV: -0.80", "EV: -1.85", "EV: -2.90", "EV: -3.95", "EV: -5.00"), each = 2)) #Define negative EVs
labs = (rep(c("PI", "Comp"), 5)) #labels for the plot - PI : previously institutionalized, comp : comparison

## Save the plots
png(file.path("C:", "Users", "jguas", "Documents", "Publications & Posters", "Papers", "2020", "SB_Cups", "NHB", "bothEV_dotchart.png"), height = 7, width = 12.5, units = "in", res = 750)
par(mfrow=c(1,2)) #save the plots side by side
dotchart(posProb, labels = labs, 
         groups = as.factor(posGroups), gcolor = 'black',
         color = 'black',
         cex = 1, cex.lab = 1.25, pch = 20, xlab = "Probability of Risky Decision", main = "Gain Trials")
dotchart(negProb, labels = labs, 
         groups = as.factor(negGroups), gcolor = 'black',
         color = 'black',
         cex = 1, cex.lab = 1.25, pch = 20, xlab = "Probability of Risky Decision", main = "Loss Trials")
dev.off()

png(file.path("C:", "Users", "jguas", "Documents", "Publications & Posters", "Papers", "2020", "SB_Cups", "NHB", "bothEV_dotchart_stacked.png"), height = 13, width = 7.5, units = "in", res = 750)
par(mfrow=c(2,1)) #save the plots stacked (what is currently show in the manuscript as of 07.09.2020)
dotchart(posProb, labels = labs, 
         groups = as.factor(posGroups), gcolor = 'black',
         color = 'black',
         cex = 1, cex.lab = 1.25, pch = 20, xlab = "Probability of Risky Decision", main = "Gain Trials")
dotchart(negProb, labels = labs, 
         groups = as.factor(negGroups), gcolor = 'black',
         color = 'black',
         cex = 1, cex.lab = 1.25, pch = 20, xlab = "Probability of Risky Decision", main = "Loss Trials")
dev.off()



## Now we're going to plot fixed effects and random effects of reward on prob(risky choice) by group (i.e., visualize the reward x group interaction)

ranefEV = aggregate(. ~ ID, data = allDat, mean)[,c("ID", "PI")] #From long file, get everyone's IDs and group status in wide format (i.e., 1 row with sub id and group status) 
ranefEV$ranef = ranef(mod3)$ID[,2] #get random effects from model 3
ranefEV$b = ranefEV$ranef + 0.674091 + (-0.335033*ranefEV$PI) #add the random effect to the fixed effect of reward and the adjustment for PI individuals to recover each subject specific association between reward and p(risky choice)
#the above are in logits!

fixefEV = data.frame(EV = rep(seq(-5,5,.05), 2), #create a dataframe with fixed effects for plotting
                     logit = c(rep(seq(-5,5,.05), 2)*0.674091, rep(seq(-5,5,.05), 2)*(0.672888-0.335033)),
                     Group = rep(c("Comp", "PI"), each = 402)) 

fixefEV$prob = exp(fixefEV$logit) / (1 + exp(fixefEV$logit)) #create a column that has probabilities (in addition to the logits)


se = 0.050996 * 2 #SE * Z value - standard error

#plot 1 fixed effect of reward on log odds(risky choice)
p1 = ggplot() + geom_line(aes(y = logit, x = EV, colour = Group), size=1.5,
                           data = fixefEV, stat="identity") +
  theme_bw() +
  labs(y = "Log Odds of Risky Decision") +
  scale_color_manual(values=brewer.pal(n = 8, name = "Dark2")[c(3,5)]) +
  geom_ribbon(aes(x = fixefEV$EV, ymin=fixefEV$logit-se, ymax=fixefEV$logit+se, fill=fixefEV$Group),alpha=0.3, show.legend=F) +
  scale_fill_manual(values=brewer.pal(n = 8, name = "Dark2")[c(3,5)]) +
  ylim(-4,4)

#plot 2 fixed effect of reward on p(risky choice)
p2 = ggplot() + geom_line(aes(y = prob, x = EV, colour = Group), size=1.5,
                          data = fixefEV, stat="identity") +
  theme_bw() +
  labs(y = "Probability of Risky Decision") +
  scale_color_manual(values=brewer.pal(n = 8, name = "Dark2")[c(3,5)]) +
  #geom_ribbon(aes(x = fixefEV$EV, ymin=fixefEV$prob-seProb, ymax=fixefEV$prob+seProb, fill=as.factor(fixefEV$Group)),alpha=0.3) +
  scale_fill_manual(values=brewer.pal(n = 8, name = "Dark2")[c(3,5)])

#Ok here's my dirty secret -- I wasn't sure how to get many individual lines in one plot, and was too tired/impatient to dig deeper, so I went mega-doofus and used a loop with ggplot. Don't tell anyone, else I'll be exiled
#plot individual random effects in logits
p1loop = ggplot() + geom_line(aes(y = logit, x = EV, colour = Group), size=1.5,
                             data = fixefEV, stat="identity") +
  theme_bw() +
  labs(y = "Log Odds of Risky Decision") +
  scale_color_manual(values=brewer.pal(n = 8, name = "Dark2")[c(3,5)]) +
  ylim(-4,4)

for (r in 1:dim(ranefEV)[1]) { 
  
  print(r)
  if (ranefEV$PI[r] == 0 ) {
    
    rcol = brewer.pal(n = 8, name = "Set2")[3]
    
  } else {
    
    rcol = brewer.pal(n = 8, name = "Set2")[5]
    
  }
  
  p1loop = p1loop + geom_line(aes_string(x=seq(-5,5,.05), y=(seq(-5,5,.05) * ranefEV$b[r])), color = rcol)#geom_segment(aes(x = -5, xend = 5, y = ranefEV$b[r]*-5, yend = ranefEV$b[r]*5), color = rcol)
  
}

#Plot individual random effects in probability metric
p2loop = ggplot() + geom_line(aes(y = prob, x = EV, colour = Group), size=1.5,
                              data = fixefEV, stat="identity") +
  theme_bw() +
  labs(y = "Probability of Risky Decision") +
  scale_color_manual(values=brewer.pal(n = 8, name = "Dark2")[c(3,5)])


for (r in 1:dim(ranefEV)[1]) { 
  
  print(r)
  if (ranefEV$PI[r] == 0 ) {
    
    rcol = brewer.pal(n = 8, name = "Set2")[3]
    
  } else {
    
    rcol = brewer.pal(n = 8, name = "Set2")[5]
    
  }
  
  subProb = seq(-5,5,.05) * ranefEV$b[r]
  subProb = exp(subProb) / (1+exp(subProb))
  
  p2loop = p2loop + geom_line(aes_string(x=seq(-5,5,.05), y=subProb), color = rcol)
  
}

#save all of these plots in 1 figure
grid.arrange(p1, p1loop, p2, p2loop, nrow = 2, ncol=2)
t = arrangeGrob(p1, p1loop, p2, p2loop, nrow = 2, ncol=2)
ggsave(file.path("C:", "Users", "jguas", "Documents", "Publications & Posters", "Papers", "2020", "SB_Cups", "NHB", "ran_fix_logit_prob.png"), t, dpi = 900, width = 10, height = 6.67, units = c("in"))

#########
## Part 2 of analysis strategy -- computational modeling
#########

## Start with prospect theory (PT)
## The way this workflow is organized is that starting values for likelihood function optimization are determined via grid search, and then they are piped into optim
## Because starting values can sometimes affect the optimization routine, I'm going to use three different starting values
## The way I did this was to run the grid search 3 times, defining the grid differently each time and ending up with a different value
## I could have chosen the values randomly, but I figured repeating the grid search three times was more often going to leave me with starting values that were in the correct neighborhood

#define list of the three different grids - first vector are pairs of bounds for each of the three PT parameters (rho, lambda, mu), second vector is n bins for the grid
initList = list(list(c(0,2,0,4.5,1,8),c(8,8,10)),
                list(c(0,2.5,0,5,.5,10),c(12,12,15)),
                list(c(0,2.5,0,5,.5,15),c(15,15,20)))

#define a dataframe that will ultimately hold the parameter estimates for each optimization (3 optimizations, each with a diff start val), the convergence code, how many iterations it took, and the final likelihood estimate
outParamsPT = data.frame(ID = L2Dat$ID, 
                         rho_init1 = rep(NA,length(L2Dat$ID)),
                         lambda_init1 = rep(NA,length(L2Dat$ID)),
                         mu_init1 = rep(NA,length(L2Dat$ID)),
                         conv_init1 = rep(NA, length(L2Dat$ID)), #convergence code for the first optimization (i.e., with first set of starting values)
                         count_init1 = rep(NA, length(L2Dat$ID)), #number of iterations it took optim to converge with this set of starting values
                         lik_init1 = rep(NA, length(L2Dat$ID)), #value of the likelihood function with this set of starting values
                         rho_init2 = rep(NA,length(L2Dat$ID)),
                         lambda_init2 = rep(NA,length(L2Dat$ID)),
                         mu_init2 = rep(NA, length(L2Dat$ID)),
                         conv_init2 = rep(NA, length(L2Dat$ID)), 
                         count_init2 = rep(NA, length(L2Dat$ID)), 
                         lik_init2 = rep(NA, length(L2Dat$ID)),
                         rho_init3 = rep(NA,length(L2Dat$ID)),
                         lambda_init3 = rep(NA,length(L2Dat$ID)),
                         mu_init3 = rep(NA, length(L2Dat$ID)),
                         conv_init3 = rep(NA, length(L2Dat$ID)), 
                         count_init3 = rep(NA, length(L2Dat$ID)), 
                         lik_init3 = rep(NA, length(L2Dat$ID)),
                         stringsAsFactors = FALSE)

#obound = c(2.5,5,15) #specify the upper bound we're going to use in optim -- not needed after I tweaked to apply a constraint a different way

#Loop over subjects to read in their raw cups data
for (sub in outParamsPT$ID) {
  
  rawDat = read_csv(sprintf("%s/%s/wave1/Lab_session/Raw/%s_CupsSelf_Raw.csv", sbPath, sub, sub)) #read in data
  rawDat = rawDat[!is.na(rawDat$cueResp.corr),] #remove NA rows from the raw data
  
  for (i in 1:length(initList)) { #loop through the initList (i.e., run the grid search and then optimization for each set of grids)
    
      gbound = initList[[i]][[1]];gbin = initList[[i]][[2]] #assign the bounds and bin numbers to separate variables
    
      ptOut = ptMod(rawDat, gbound, gbin, obound, 0) #fit the prospect theory model
      
      #assign output to dataframe (probably a more efficient way to do this)
      outParamsPT[outParamsPT$ID == sub, paste('rho_init', as.character(i), sep="")] = ptOut[[1]]
      outParamsPT[outParamsPT$ID == sub, paste('lambda_init', as.character(i), sep="")] = ptOut[[2]]
      outParamsPT[outParamsPT$ID == sub, paste('mu_init', as.character(i), sep="")] = ptOut[[3]]
      outParamsPT[outParamsPT$ID == sub, paste('conv_init', as.character(i), sep="")] = ptOut[[5]]
      outParamsPT[outParamsPT$ID == sub, paste('count_init', as.character(i), sep="")] = ptOut[[4]][1]
      outParamsPT[outParamsPT$ID == sub, paste('lik_init', as.character(i), sep="")] = ptOut[[7]]
    
  }
  
  
}

#Set any values that didn't converge to be missing
outParamsPT[outParamsPT$conv_init1>0,c('rho_init1', 'lambda_init1', 'mu_init1')] = NA
outParamsPT[outParamsPT$conv_init2>0,c('rho_init2', 'lambda_init2', 'mu_init2')] = NA
outParamsPT[outParamsPT$conv_init3>0,c('rho_init3', 'lambda_init3', 'mu_init3')] = NA

#Even if optim said it converged I'm going to count values right at the boundary as not converging
outParamsPT$lambda_init1[which(outParamsPT$lambda_init1==10)] = NA
outParamsPT$lambda_init2[which(outParamsPT$lambda_init2==10)] = NA
outParamsPT$lambda_init3[which(outParamsPT$lambda_init3==10)] = NA

outParamsPT$lambda_init1[which(outParamsPT$lambda_init1<.05)] = NA
outParamsPT$lambda_init2[which(outParamsPT$lambda_init2<.05)] = NA
outParamsPT$lambda_init3[which(outParamsPT$lambda_init3<.05)] = NA

outParamsPT$rho_init1[which(outParamsPT$rho_init1<.05)] = NA
outParamsPT$rho_init2[which(outParamsPT$rho_init2<.05)] = NA
outParamsPT$rho_init3[which(outParamsPT$rho_init3<.05)] = NA

outParamsPT$mu_init1[which(outParamsPT$mu_init1<.01)] = NA
outParamsPT$mu_init2[which(outParamsPT$mu_init2<.01)] = NA
outParamsPT$mu_init3[which(outParamsPT$mu_init3<.01)] = NA

outParamsPT$mu_init1[which(outParamsPT$mu_init1>9.9)] = NA
outParamsPT$mu_init2[which(outParamsPT$mu_init2>9.9)] = NA
outParamsPT$mu_init3[which(outParamsPT$mu_init3>9.9)] = NA

#Take only the data we need (subject ID, parameter estimates (averaged over the three sets))
outConv = data.frame(ID = outParamsPT$ID,
                     rho = rowMeans(as.matrix(outParamsPT[,c('rho_init1', 'rho_init2', 'rho_init3')]), na.rm = TRUE),
                     lambda = rowMeans(as.matrix(outParamsPT[,c('lambda_init1', 'lambda_init2', 'lambda_init3')]), na.rm = TRUE), 
                     mu = rowMeans(as.matrix(outParamsPT[,c("mu_init1", "mu_init2", "mu_init3")]), na.rm = TRUE),
                     stringsAsFactors = FALSE)

## Now analyze prospect theory parameters
compDat = inner_join(outConv, L2Dat) #join the computational modeling data with the other data
t.test(compDat$rho ~ compDat$PI) #independent samples t-tests for all computational parameters
t.test(log(compDat$lambda) ~ compDat$PI) #earlier versions of this script also peeped differences based on age or sex for exploratory puposes, removed here for clarity
t.test(log(compDat$mu) ~ compDat$PI)

#testing effect of ELA on these parameters while adjusting for covariates
rhoMod1 = lm(rho ~ PI + Age + Sex + IQ_percentile, dat = compDat)
lambdaMod1 = lm(log(lambda) ~ PI + Age + Sex + IQ_percentile, dat = compDat)
muMod1 = lm(log(mu) ~ PI + Age + Sex + IQ_percentile, dat = compDat)

## Now impute the missing values and re-analyze
#Here we're using group, gender (labeled as sex here due to an earlier confusion), age, dospert scores and WASI-II scores
impRho = imputeKNN(compDat[,c('rho', 'PI', 'Sex', 'Age', 'DOSPERT', 'IQ_percentile')], 'rho', 5)
compDat$impRho = impRho
impLambda = imputeKNN(compDat[,c('lambda', 'PI', 'Sex', 'Age', 'DOSPERT', 'IQ_percentile')], 'lambda', 5)
compDat$impLambda = impLambda
impMu = imputeKNN(compDat[,c('mu', 'PI', 'Sex', 'Age', 'DOSPERT', 'IQ_percentile')], 'mu', 5)
compDat$impMu = impMu

#rerun the analyses (t-tests and regression analyses)
t.test(compDat$impRho ~ compDat$PI)
t.test(log(compDat$impLambda) ~ compDat$PI)
t.test(log(compDat$impMu) ~ compDat$PI)

impRhoMod = lm(impRho ~ PI + Age + Sex + IQ_percentile, dat = compDat)
impLambdaMod = lm(log(impLambda) ~ PI + Age + Sex + IQ_percentile, dat = compDat)
impMuMod = lm(log(impMu) ~ PI + Age + Sex + IQ_percentile, dat = compDat)

## Now we're going to do the target model (TM)
## Same workflow as before

#list of conditions for grid searches for initial parameters
initList = list(list(c(0,2.5,1,8),c(8,10)),
                list(c(0,3,.5,10),c(12,15)),
                list(c(0,3,.5,12),c(15,20)))

#notably, we need to re-estimate a mu (choice consistency) parameter for this model as well. I suppose we could have used the one from PT, but that seemed inappropriate here
#but we only care about tau right now. 
#dataframe to hold our output. Again, we're only including tau here (i.e., we're estimating mu, but not doing anything more with it)
outParamsTM = data.frame(ID = L2Dat$ID, 
                         tau_init1 = rep(NA,length(L2Dat$ID)),
                         conv_init1 = rep(NA, length(L2Dat$ID)), 
                         count_init1 = rep(NA, length(L2Dat$ID)), 
                         lik_init1 = rep(NA, length(L2Dat$ID)), 
                         tau_init2 = rep(NA,length(L2Dat$ID)),
                         conv_init2 = rep(NA, length(L2Dat$ID)), 
                         count_init2 = rep(NA, length(L2Dat$ID)), 
                         lik_init2 = rep(NA, length(L2Dat$ID)),
                         tau_init3 = rep(NA,length(L2Dat$ID)),
                         conv_init3 = rep(NA, length(L2Dat$ID)), 
                         count_init3 = rep(NA, length(L2Dat$ID)), 
                         lik_init3 = rep(NA, length(L2Dat$ID)),
                         stringsAsFactors = FALSE)

obound = c(3,15) #upper bounds for optim, but again, as before, not needed after the 'hack' I mentioned above

#Loop over subjects, fit the TM, and save the parameter estimate of tau
for (sub in outParamsTM$ID) {
  
  inDat = read_csv(sprintf("%s/%s/wave1/Lab_session/Clean/%s_CupsSelf.csv", sbPath, sub, sub))
  
  for (i in 1:length(initList)) {
    
    #same deal as before, set up conditions for our grid searches
    gbound = initList[[i]][[1]];gbin = initList[[i]][[2]]
    
    #fit the model
    tmOut = tmMod(inDat, gbound, gbin, obound, 0)
    
    #save what we need (similar values as before)
    outParamsTM[outParamsTM$ID == sub, paste('tau_init', as.character(i), sep="")] = tmOut[[1]]
    outParamsTM[outParamsTM$ID == sub, paste('conv_init', as.character(i), sep="")] = tmOut[[4]]
    outParamsTM[outParamsTM$ID == sub, paste('count_init', as.character(i), sep="")] = tmOut[[3]][1]
    outParamsTM[outParamsTM$ID == sub, paste('lik_init', as.character(i), sep="")] = tmOut[[6]]
    
  }
  
  
}

#Set any values that didn't converge to be missing
outParamsTM[outParamsTM$conv_init1>0,c('tau_init1')] = NA
outParamsTM[outParamsTM$conv_init2>0,c('rho_init2')] = NA
outParamsTM[outParamsTM$conv_init3>0,c('rho_init3')] = NA

#same deal as before -- any implausible low tau values are treated as missing (non-converged)
outParamsTM$tau_init1[outParamsTM$tau_init1 < .01] = NA
outParamsTM$tau_init2[outParamsTM$tau_init2 < .01] = NA
outParamsTM$tau_init3[outParamsTM$tau_init3 < .01] = NA

#take mean of the three parameter estimates
outConv = data.frame(ID = outParamsTM$ID,
                     tau = rowMeans(as.matrix(outParamsTM[,c('tau_init1', 'tau_init2', 'tau_init3')]), na.rm = TRUE), 
                     stringsAsFactors = FALSE)

## Run statistical tests
compDat = inner_join(outConv, L2Dat)
t.test(log(compDat$tau) ~ compDat$PI)


tauMod1 = lm(log(tau) ~ PI+Age+Sex+IQ_percentile, dat = compDat)

## Impute the data and re-run analyses
impTau = imputeKNN(compDat[,c('tau', 'PI', 'Sex', 'Age', 'DOSPERT', 'IQ_percentile')], 'tau', 5)
compDat$impTau = impTau
t.test(log(compDat$impTau) ~ compDat$PI)

impTauMod1 = lm(log(impTau) ~ PI+Age+Sex+IQ_percentile, dat = compDat)


#########
## Analyze DOSPERT data (self-reported 'real world' risk taking)
#########

#independent samples t test
t.test(compDat$DOSPERT ~ compDat$PI)

#adjust based on covariates
dosMod = lm(DOSPERT ~ Age + PI + Sex + IQ_percentile, data = compDat)

#########
## Additional supplemental analyses suggested by reviewers
#########
