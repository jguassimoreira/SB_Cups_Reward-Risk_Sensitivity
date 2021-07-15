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
mod3iq = glmer(Decision ~ EV + SD + TrialType + Age + Sex + IQ_percentile*PI + (1+EV+SD|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa")) #Model 3 from Table 2


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

