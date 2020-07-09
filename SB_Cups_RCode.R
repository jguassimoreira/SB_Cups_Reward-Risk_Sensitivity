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

##Import helper functions
scriptPath = file.path("~", "Documents", "Scripts", "R", "helperFuncs")
sbPath = file.path("P:", "Social_Behavior_(SB_study)", "Data", "Behavioral_data")
file.sources = Sys.glob(file.path(sprintf("%s", scriptPath), "*.R"))
sapply(file.sources,FUN=source)


### Model 1 - Modeling risky behaviors on the cups task as a function of trial-level variables
### and between subject predictors

#Get path to subjects' data
subList = Sys.glob("P:/Social_Behavior_(SB_study)/Data/Behavioral_data/SB*/wave1/Lab_session/Clean/*_CupsSelf.csv")

#Compile Cups Dat (Level 1, trial level)
cupsDat = getCupsDat(subList)
cupsDat = cupsDat[-1,]

#Get Level 2 (individual level) data
L2Dat = getL2Dat("P:/Social_Behavior_(SB_study)/Data/Study_logs/SB_group_list.csv")

#merge
allDat = inner_join(L2Dat, cupsDat)

#grand mean center EV(return), SD(risk), age, and IQ
allDat$SD = allDat$SD - mean(allDat$SD); allDat$EV = allDat$EV - mean(allDat$EV)
allDat$Age = allDat$Age - mean(allDat$Age)
allDat$IQ_percentile = allDat$IQ_percentile - mean(allDat$IQ_percentile)
allDat$DOSPERT = allDat$DOSPERT - mean(allDat$DOSPERT)

#run the models
mod1 <- glmer(Decision ~ EV + SD + TrialType + (1+EV+SD+TrialType|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mod2 <- glmer(Decision ~ EV + SD + TrialType + PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mod3 <- glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mod4 <- glmer(Decision ~ EV*PI + EV*Age + SD*PI + SD*Age + TrialType*PI + TrialType*Age + Sex + IQ_percentile + (1+EV|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mod5 <- glmer(Decision ~ EV*PI + SD*PI*TrialType  + Age + Sex + IQ_percentile + (1+SD|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))

#Add DOSPERT into model 3
mod3dos <- glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + DOSPERT + (1+EV+SD|ID), data = allDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))


#Let's do some recentering the get some regions of significance
rcDat_m5 = allDat; rcDat_m5$EV = rcDat_m5$EV - (-5) #recenter at minus 5 EV
rcDat_m2 = allDat; rcDat_m2$EV = rcDat_m2$EV - (-2) #recenter at minus 2 EV
rcDat_p2 = allDat; rcDat_p2$EV = rcDat_p2$EV - 2 #recenter at plus 2 EV
rcDat_p5 = allDat; rcDat_p5$EV = rcDat_p5$EV - 5 #recenter at plust 5 EV

mod3i_rcM5 <- glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = rcDat_m5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 150000)))
mod3i_rcM2 <- glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = rcDat_m2, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 150000)))
mod3i_rcP2 <- glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = rcDat_p2, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 150000)))
mod3i_rcP5 <- glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = rcDat_p5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 150000)))


#Run with gain and loss trials separate

gainDat = inner_join(L2Dat, cupsDat); gainDat = gainDat[gainDat$EV > 0,]
lossDat = inner_join(L2Dat, cupsDat); lossDat = lossDat[lossDat$EV < 0,]

gainDat$SD = gainDat$SD - mean(gainDat$SD); gainDat$EV = gainDat$EV - mean(gainDat$EV)
gainDat$Age = gainDat$Age - mean(gainDat$Age)
gainDat$IQ_percentile = gainDat$IQ_percentile - mean(gainDat$IQ_percentile)

lossDat$SD = lossDat$SD - mean(lossDat$SD); lossDat$EV = lossDat$EV - mean(lossDat$EV)
lossDat$Age = lossDat$Age - mean(lossDat$Age)
lossDat$IQ_percentile = lossDat$IQ_percentile - mean(lossDat$IQ_percentile)

gainMod1 <- glmer(Decision ~ EV + SD + (1+EV+SD|ID), data = gainDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
lossMod1 <- glmer(Decision ~ EV + SD + (1+EV+SD|ID), data = lossDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
gainMod2 <- glmer(Decision ~ EV + SD + PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = gainDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
lossMod2 <- glmer(Decision ~ EV + SD + PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = lossDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
gainMod3 <- glmer(Decision ~ EV*PI + SD*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = gainDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
lossMod3 <- glmer(Decision ~ EV*PI + SD*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = lossDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))

#Run without SB334 (outlier)
allDat_no334 = inner_join(L2Dat[-(L2Dat$ID=='SB334'),], cupsDat)

allDat_no334$SD = allDat_no334$SD - mean(allDat_no334$SD); allDat_no334$EV = allDat_no334$EV - mean(allDat_no334$EV)
allDat_no334$Age = allDat_no334$Age - mean(allDat_no334$Age)
allDat_no334$IQ_percentile = allDat_no334$IQ_percentile - mean(allDat_no334$IQ_percentile)

#run the models
mod1 <- glmer(Decision ~ EV + SD + TrialType + (1+EV+SD+TrialType|ID), data = allDat_no334, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mod2 <- glmer(Decision ~ EV + SD + TrialType + PI + Age + Sex + (1+EV+SD+TrialType|ID), data = allDat_no334, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mod3 <- glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + (1+EV+SD+TrialType|ID), data = allDat_no334, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mod4 <- glmer(Decision ~ EV*PI + EV*Age + SD*PI + SD*Age + TrialType*PI + TrialType*Age + Sex + (1+EV+SD|ID), data = allDat_no334, family = binomial, control = glmerControl(optimizer = "bobyqa"))

mod2i <- glmer(Decision ~ EV + SD + TrialType + PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = allDat_no334, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mod3i <- glmer(Decision ~ EV*PI + SD*PI + TrialType*PI + Age + Sex + IQ_percentile + (1+EV+SD|ID), data = allDat_no334, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mod4i <- glmer(Decision ~ EV*PI + EV*Age + SD*PI + SD*Age + TrialType*PI + TrialType*Age + Sex + IQ_percentile + (1+EV|ID), data = allDat_no334, family = binomial, control = glmerControl(optimizer = "bobyqa"))



##Make some plots
#Model implied probabilities
posEV = c(rep(linspace(.8, 5, 5), each = 2))
negEV = c(rep(linspace(-.8, -5, 5), each = 2))
age=0
#loAge = -(sd(allDat$Age) - (.5*sd(allDat$Age))); meanAge = 0; hiAge = (sd(allDat$Age) + (.5*sd(allDat$Age)))
#age = c(rep(c(loAge, loAge, meanAge, meanAge, hiAge, hiAge), 5))
grpPI = c(rep(c(1,0), 5))

posLogits = 0.751 + (0.674*posEV) + (0.141*0) + (-1.290*1) + (-0.656*grpPI) + (-0.043*age) + (0*-0.005) + (-.335*posEV*grpPI) + (-.077*0*grpPI) + (1.266*1*grpPI) + 0.0993085 #this last bit here is to take avg of gender effect
posProb = exp(posLogits) / (1 + exp(posLogits))

negLogits = 0.751 + (0.674*negEV) + (0.141*0) + (-1.290*0) + (-0.656*grpPI) + (-0.043*age) + (0*-0.005) + (-.335*negEV*grpPI) + (-.077*0*grpPI) + (1.266*0*grpPI) + 0.0993085 #this last bit here is to take avg of gender effect
negProb = exp(negLogits) / (1 + exp(negLogits))

posGroups = (rep(c("EV: +0.80", "EV: +1.85", "EV: +2.90", "EV: +3.95", "EV: +5.00"), each = 2))
negGroups = (rep(c("EV: -0.80", "EV: -1.85", "EV: -2.90", "EV: -3.95", "EV: -5.00"), each = 2))
labs = (rep(c("PI", "Comp"), 5))

#my_cols = wes_palette("FantasticFox1", n = 5) #c("gray", "#E69F00", "#56B4E9", "purple", "green")
#my_cols = brewer.pal(n = 8, name = "Dark2")[c(1,2,3,4,8)]


#png(file.path("C:", "Users", "jguas", "Documents", "Publications & Posters", "Papers", "2020", "SB_Cups", "CD_SpecialIssue", "posEV_dotchart.png"), height = 5, width = 5, units = "in", res = 750)
png(file.path("C:", "Users", "jguas", "Documents", "Publications & Posters", "Papers", "2020", "SB_Cups", "CD_SpecialIssue", "posEV_dotchart.png"), height = 5, width = 5, units = "in", res = 750)
dcp = dotchart(posProb, labels = labs, 
         groups = as.factor(posGroups), gcolor = 'black',
         color = 'black',
         cex = 1, cex.lab = 1.25, pch = 20, xlab = "Probability of Risky Decision", main = "Gain Trials")
dev.off()

png(file.path("C:", "Users", "jguas", "Documents", "Publications & Posters", "Papers", "2020", "SB_Cups", "CD_SpecialIssue", "negEV_dotchart.png"), height = 5, width = 5, units = "in", res = 750)
dcn = dotchart(negProb, labels = labs, 
         groups = as.factor(negGroups), gcolor = 'black',
         color = 'black',
         cex = 1, cex.lab = 1.25, pch = 20, xlab = "Probability of Risky Decision", main = "Loss Trials")
dev.off()

png(file.path("C:", "Users", "jguas", "Documents", "Publications & Posters", "Papers", "2020", "SB_Cups", "NHB", "bothEV_dotchart.png"), height = 7, width = 12.5, units = "in", res = 750)
par(mfrow=c(1,2))
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
par(mfrow=c(2,1))
dotchart(posProb, labels = labs, 
         groups = as.factor(posGroups), gcolor = 'black',
         color = 'black',
         cex = 1, cex.lab = 1.25, pch = 20, xlab = "Probability of Risky Decision", main = "Gain Trials")
dotchart(negProb, labels = labs, 
         groups = as.factor(negGroups), gcolor = 'black',
         color = 'black',
         cex = 1, cex.lab = 1.25, pch = 20, xlab = "Probability of Risky Decision", main = "Loss Trials")
dev.off()

#dev.off()
#dotchart(posProb, labels = labs, 
#         groups = as.factor(posGroups), gcolor = my_cols,
#         color = my_cols[as.factor(posGroups)],
#         cex = .5, cex.lab = 1.2, pch = 19, xlab = "Probability of Risky Decision")
#dev.off()

#png(file.path("C:", "Users", "jguas", "Documents", "Publications & Posters", "Papers", "2020", "SB_Cups", "CD_SpecialIssue", "negEV_dotchart.png"), height = 5, width = 5, units = "in", res = 750)
#dotchart(negProb, labels = labs, 
#         groups = as.factor(negGroups), gcolor = my_cols,
#         color = my_cols[as.factor(posGroups)],
#         cex = .5, cex.lab = 1.2, pch = 19, xlab = "Probability of Risky Decision")
#dev.off()



#Break down EV by group

ranefEV = aggregate(. ~ ID, data = allDat, mean)[,c("ID", "PI")]
ranefEV$ranef = ranef(mod3)$ID[,2]
ranefEV$b = ranefEV$ranef + 0.674091 + (-0.335033*ranefEV$PI)

#fixefEV = data.frame(EV = rep(sort(unique(allDat$EV))[-c(2,6,12,16)], 2),
#                     logit = c(sort(unique(allDat$EV))[-c(2,6,12,16)]*0.674091, sort(unique(allDat$EV))[-c(2,6,12,16)]*(0.672888-0.335033)),
#                     Group = rep(c("Comp", "PI"), each = 14))

fixefEV = data.frame(EV = rep(seq(-5,5,.05), 2),
                     logit = c(rep(seq(-5,5,.05), 2)*0.674091, rep(seq(-5,5,.05), 2)*(0.672888-0.335033)),
                     Group = rep(c("Comp", "PI"), each = 402))

fixefEV$prob = exp(fixefEV$logit) / (1 + exp(fixefEV$logit))


se = 0.050996 * 2 #SE * Z value

p1 = ggplot() + geom_line(aes(y = logit, x = EV, colour = Group), size=1.5,
                           data = fixefEV, stat="identity") +
  theme_bw() +
  labs(y = "Log Odds of Risky Decision") +
  scale_color_manual(values=brewer.pal(n = 8, name = "Dark2")[c(3,5)]) +
  geom_ribbon(aes(x = fixefEV$EV, ymin=fixefEV$logit-se, ymax=fixefEV$logit+se, fill=fixefEV$Group),alpha=0.3, show.legend=F) +
  scale_fill_manual(values=brewer.pal(n = 8, name = "Dark2")[c(3,5)]) +
  ylim(-4,4)
png(file.path("C:", "Users", "jguas", "Documents", "Publications & Posters", "Papers", "2020", "SB_Cups", "CD_SpecialIssue", "fixed_effect_ev.png"), height = 4, width = 5.5, units = "in", res = 600)
p1
dev.off()

p2 = ggplot() + geom_line(aes(y = prob, x = EV, colour = Group), size=1.5,
                          data = fixefEV, stat="identity") +
  theme_bw() +
  labs(y = "Probability of Risky Decision") +
  scale_color_manual(values=brewer.pal(n = 8, name = "Dark2")[c(3,5)]) +
  #geom_ribbon(aes(x = fixefEV$EV, ymin=fixefEV$prob-seProb, ymax=fixefEV$prob+seProb, fill=as.factor(fixefEV$Group)),alpha=0.3) +
  scale_fill_manual(values=brewer.pal(n = 8, name = "Dark2")[c(3,5)])

p1loop = ggplot() + geom_line(aes(y = logit, x = EV, colour = Group), size=1.5,
                             data = fixefEV, stat="identity") +
  theme_bw() +
  labs(y = "Log Odds of Risky Decision") +
  scale_color_manual(values=brewer.pal(n = 8, name = "Dark2")[c(3,5)]) +
  ylim(-4,4)
#p2 = p1 + geom_abline(aes(slope = ranefEV$b[1], intercept = 0)) + xlim(-5,5)
#p3 = p1 + geom_segment(aes(x = -5, xend = 5, y = ranefEV$b[1]*-5, yend = ranefEV$b[1]*5))


for (r in 1:dim(ranefEV)[1]) { #
  
  print(r)
  if (ranefEV$PI[r] == 0 ) {
    
    rcol = brewer.pal(n = 8, name = "Set2")[3]
    
  } else {
    
    rcol = brewer.pal(n = 8, name = "Set2")[5]
    
  }
  
  p1loop = p1loop + geom_line(aes_string(x=seq(-5,5,.05), y=(seq(-5,5,.05) * ranefEV$b[r])), color = rcol)#geom_segment(aes(x = -5, xend = 5, y = ranefEV$b[r]*-5, yend = ranefEV$b[r]*5), color = rcol)
  
}


p2loop = ggplot() + geom_line(aes(y = prob, x = EV, colour = Group), size=1.5,
                              data = fixefEV, stat="identity") +
  theme_bw() +
  labs(y = "Probability of Risky Decision") +
  scale_color_manual(values=brewer.pal(n = 8, name = "Dark2")[c(3,5)])
#p2 = p1 + geom_abline(aes(slope = ranefEV$b[1], intercept = 0)) + xlim(-5,5)
#p3 = p1 + geom_segment(aes(x = -5, xend = 5, y = ranefEV$b[1]*-5, yend = ranefEV$b[1]*5))


for (r in 1:dim(ranefEV)[1]) { #
  
  print(r)
  if (ranefEV$PI[r] == 0 ) {
    
    rcol = brewer.pal(n = 8, name = "Set2")[3]
    
  } else {
    
    rcol = brewer.pal(n = 8, name = "Set2")[5]
    
  }
  
  subProb = seq(-5,5,.05) * ranefEV$b[r]
  subProb = exp(subProb) / (1+exp(subProb))
  
  p2loop = p2loop + geom_line(aes_string(x=seq(-5,5,.05), y=subProb), color = rcol)#geom_segment(aes(x = -5, xend = 5, y = ranefEV$b[r]*-5, yend = ranefEV$b[r]*5), color = rcol)
  
}

grid.arrange(p1, p1loop, p2, p2loop, nrow = 2, ncol=2)
t = arrangeGrob(p1, p1loop, p2, p2loop, nrow = 2, ncol=2)
ggsave(file.path("C:", "Users", "jguas", "Documents", "Publications & Posters", "Papers", "2020", "SB_Cups", "NHB", "ran_fix_logit_prob.png"), t, dpi = 900, width = 10, height = 6.67, units = c("in"))

png(file.path("C:", "Users", "jguas", "Documents", "Publications & Posters", "Papers", "2020", "SB_Cups", "CD_SpecialIssue", "random_effect_ev.png"), height = 4, width = 5, units = "in", res = 600)
ploop
dev.off()

##Plot EV/SD effects across group seperately for gain and loss

#gain
pEV_gain = ggplot() + geom_line(aes_string(x = -5:5, y = (1.226671+(c(-5:5)*0.878650))), color = brewer.pal(n = 8, name = "Set2")[3], size=1.5) +
  geom_line(aes_string(x = -5:5, y = ((1.226671-.169318)+(c(-5:5)*(.878650-.516078)))), color = brewer.pal(n = 8, name = "Set2")[5], size=1.5) +
  theme_bw() +
  labs(y = "Log Odds of Risky Decision") +
  ylim(-6,6)

pEV_loss = ggplot() + geom_line(aes_string(x = -5:5, y = (-.731316+(c(-5:5)*0.748608))), color = brewer.pal(n = 8, name = "Dark2")[3], size=1.5) +
  geom_line(aes_string(x = -5:5, y = ((-.731316+.148945)+(c(-5:5)*(.748608-.522223)))), color = brewer.pal(n = 8, name = "Dark2")[5], size=1.5) +
  theme_bw() +
  labs(y = "Log Odds of Risky Decision") +
  ylim(-6,6)

pEV_gain_loss = ggplot() + geom_line(aes_string(x = 0:5, y = (1.226671+(c(0:5)*0.878650))), color = brewer.pal(n = 8, name = "Set2")[3], size=1.5) +
  geom_line(aes_string(x = 0:5, y = ((1.226671-.169318)+(c(0:5)*(.878650-.516078)))), color = brewer.pal(n = 8, name = "Set2")[5], size=1.5) +
  geom_line(aes_string(x = -5:0, y = (-.731316+(c(-5:0)*0.748608))), color = brewer.pal(n = 8, name = "Dark2")[3], size=1.5) +
  geom_line(aes_string(x = -5:0, y = ((-.731316+.148945)+(c(-5:0)*(.748608-.522223)))), color = brewer.pal(n = 8, name = "Dark2")[5], size=1.5) +
  theme_bw() +
  labs(y = "Log Odds of Risky Decision", x = "EV") +
  ylim(-6,6)

### Model 2 - Interrogating mechanism by following up with prospect theory

#initList = list(list(c(0,5,0,5,0,12),c(8,8,10)),
#                list(c(0,10,0,10,0,15),c(12,12,15)),
#                list(c(0,10,0,10,0,25),c(15,15,20)))

initList = list(list(c(0,2,0,4.5,1,8),c(8,8,10)),
                list(c(0,2.5,0,5,.5,10),c(12,12,15)),
                list(c(0,2.5,0,5,.5,15),c(15,15,20)))


outParamsPT = data.frame(ID = L2Dat$ID, 
                         rho_init1 = rep(NA,length(L2Dat$ID)),
                         lambda_init1 = rep(NA,length(L2Dat$ID)),
                         mu_init1 = rep(NA,length(L2Dat$ID)),
                         conv_init1 = rep(NA, length(L2Dat$ID)), 
                         count_init1 = rep(NA, length(L2Dat$ID)), 
                         lik_init1 = rep(NA, length(L2Dat$ID)), 
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
#gbound = c(0,10,0,10,0,15)# c(0, 5, 0, 5, 0, 15)
#gbin = c(12,12,15)# c(8, 8, 10)
obound = c(2.5,5,15) #c(10,10,35)


for (sub in outParamsPT$ID) {
  
  rawDat = read_csv(sprintf("%s/%s/wave1/Lab_session/Raw/%s_CupsSelf_Raw.csv", sbPath, sub, sub))
  rawDat = rawDat[!is.na(rawDat$cueResp.corr),]
  
  for (i in 1:length(initList)) {
    
      gbound = initList[[i]][[1]];gbin = initList[[i]][[2]]
    
      ptOut = ptMod(rawDat, gbound, gbin, obound, 0)
      
      outParamsPT[outParamsPT$ID == sub, paste('rho_init', as.character(i), sep="")] = ptOut[[1]]
      outParamsPT[outParamsPT$ID == sub, paste('lambda_init', as.character(i), sep="")] = ptOut[[2]]
      outParamsPT[outParamsPT$ID == sub, paste('mu_init', as.character(i), sep="")] = ptOut[[3]]
      outParamsPT[outParamsPT$ID == sub, paste('conv_init', as.character(i), sep="")] = ptOut[[5]]
      outParamsPT[outParamsPT$ID == sub, paste('count_init', as.character(i), sep="")] = ptOut[[4]][1]
      outParamsPT[outParamsPT$ID == sub, paste('lik_init', as.character(i), sep="")] = ptOut[[7]]
    
  }
  
  
}

outParamsPT[outParamsPT$conv_init1>0,c('rho_init1', 'lambda_init1', 'mu_init1')] = NA
outParamsPT[outParamsPT$conv_init2>0,c('rho_init2', 'lambda_init2', 'mu_init2')] = NA
outParamsPT[outParamsPT$conv_init3>0,c('rho_init3', 'lambda_init3', 'mu_init3')] = NA

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


outConv = data.frame(ID = outParamsPT$ID,
                     rho = rowMeans(as.matrix(outParamsPT[,c('rho_init1', 'rho_init2', 'rho_init3')]), na.rm = TRUE),
                     lambda = rowMeans(as.matrix(outParamsPT[,c('lambda_init1', 'lambda_init2', 'lambda_init3')]), na.rm = TRUE), 
                     mu = rowMeans(as.matrix(outParamsPT[,c("mu_init1", "mu_init2", "mu_init3")]), na.rm = TRUE),
                     stringsAsFactors = FALSE)



## Need to run test on params (PI v Comp, correlate with age, etc.)
compDat = inner_join(outConv, L2Dat)
t.test(compDat$rho ~ compDat$PI)
t.test(compDat$lambda ~ compDat$PI)
t.test(log(compDat$lambda) ~ compDat$PI)
t.test(log(compDat$mu) ~ compDat$PI)
t.test(compDat$rho ~ compDat$Sex)
t.test(compDat$lambda ~ compDat$Sex)
cor.test(compDat$rho, compDat$Age)
cor.test(compDat$lambda, compDat$Age)

rhoMod1 = lm(rho ~ PI + Age + Sex + IQ_percentile, dat = compDat)
lambdaMod1 = lm(log(lambda) ~ PI + Age + Sex + IQ_percentile, dat = compDat)
muMod1 = lm(log(mu) ~ PI + Age + Sex + IQ_percentile, dat = compDat)


plot(compDat$Age, log(compDat$lambda), pch=21, 
     bg=c("red","blue")[unclass(as.factor(compDat$PI))])

#use imputed data
impRho = imputeKNN(compDat[,c('rho', 'PI', 'Sex', 'Age', 'DOSPERT', 'IQ_percentile')], 'rho', 5)
compDat$impRho = impRho
impLambda = imputeKNN(compDat[,c('lambda', 'PI', 'Sex', 'Age', 'DOSPERT', 'IQ_percentile')], 'lambda', 5)
compDat$impLambda = impLambda
impMu = imputeKNN(compDat[,c('mu', 'PI', 'Sex', 'Age', 'DOSPERT', 'IQ_percentile')], 'mu', 5)
compDat$impMu = impMu

t.test(compDat$impRho ~ compDat$PI)
t.test(log(compDat$impLambda) ~ compDat$PI)
t.test(log(compDat$impMu) ~ compDat$PI)

impRhoMod = lm(impRho ~ PI + Age + Sex + IQ_percentile, dat = compDat)
impLambdaMod = lm(log(impLambda) ~ PI + Age + Sex + IQ_percentile, dat = compDat)
impMuMod = lm(log(impMu) ~ PI + Age + Sex + IQ_percentile, dat = compDat)




plot(compDat$Age, log(compDat$impLambda), pch=21, 
     bg=c("red","blue")[unclass(as.factor(compDat$PI))])
plot(compDat$Age, compDat$impLambda, pch=21, 
     bg=c("red","blue")[unclass(as.factor(compDat$PI))])
plot(compDat$Age, compDat$impRho, pch=21, 
     bg=c("red","blue")[unclass(as.factor(compDat$PI))])
plot(compDat$Age, log(compDat$impMu), pch=21, 
     bg=c("red","blue")[unclass(as.factor(compDat$PI))])



### Model 3 - Interrogating mechanism by following up with a target model
initList = list(list(c(0,2.5,1,8),c(8,10)),
                list(c(0,3,.5,10),c(12,15)),
                list(c(0,3,.5,12),c(15,20)))


outParamsTM = data.frame(ID = L2Dat$ID, 
                         tau_init1 = rep(NA,length(L2Dat$ID)),
                         conv_init1 = rep(NA, length(L2Dat$ID)), 
                         count_init1 = rep(NA, length(L2Dat$ID)), 
                         lik_init1 = rep(NA, length(L2Dat$ID)), 
                         tau_init2 = rep(NA,length(L2Dat$ID)),
                         conv_init2 = rep(NA, length(L2Dat$ID)), 
                         count_init2 = rep(NA, length(L2Dat$ID)), 
                         lik_init2 = rep(NA, length(L2Dat$ID)),
                         tu_init3 = rep(NA,length(L2Dat$ID)),
                         conv_init3 = rep(NA, length(L2Dat$ID)), 
                         count_init3 = rep(NA, length(L2Dat$ID)), 
                         lik_init3 = rep(NA, length(L2Dat$ID)),
                         stringsAsFactors = FALSE)
#gbound = c(0,10,0,10,0,15)# c(0, 5, 0, 5, 0, 15)
#gbin = c(12,12,15)# c(8, 8, 10)
obound = c(3,15) #c(10,10,35)


for (sub in outParamsTM$ID) {
  
  inDat = read_csv(sprintf("%s/%s/wave1/Lab_session/Clean/%s_CupsSelf.csv", sbPath, sub, sub))
  
  for (i in 1:length(initList)) {
    
    gbound = initList[[i]][[1]];gbin = initList[[i]][[2]]
    
    tmOut = tmMod(inDat, gbound, gbin, obound, 0)
    
    outParamsTM[outParamsTM$ID == sub, paste('tau_init', as.character(i), sep="")] = tmOut[[1]]
    outParamsTM[outParamsTM$ID == sub, paste('conv_init', as.character(i), sep="")] = tmOut[[4]]
    outParamsTM[outParamsTM$ID == sub, paste('count_init', as.character(i), sep="")] = tmOut[[3]][1]
    outParamsTM[outParamsTM$ID == sub, paste('lik_init', as.character(i), sep="")] = tmOut[[6]]
    
  }
  
  
}

outParamsTM$tau_init1[outParamsTM$tau_init1 < .01] = NA
outParamsTM$tau_init2[outParamsTM$tau_init2 < .01] = NA
outParamsTM$tau_init3[outParamsTM$tau_init3 < .01] = NA

outConv = data.frame(ID = outParamsTM$ID,
                     tau = rowMeans(as.matrix(outParamsTM[,c('tau_init1', 'tau_init2', 'tau_init3')]), na.rm = TRUE), 
                     stringsAsFactors = FALSE)

compDat = inner_join(outConv, L2Dat)
t.test(log(compDat$tau) ~ compDat$PI)


tauMod1 = lm(log(tau) ~ PI+Age+Sex+IQ_percentile, dat = compDat)


plot(compDat$Age, compDat$tau, pch=21, 
     bg=c("red","blue")[unclass(as.factor(compDat$PI))])

##impute the data

impTau = imputeKNN(compDat[,c('tau', 'PI', 'Sex', 'Age', 'DOSPERT', 'IQ_percentile')], 'tau', 5)
compDat$impTau = impTau
t.test(log(compDat$impTau) ~ compDat$PI)

impTauMod1 = lm(log(impTau) ~ PI+Age+Sex+IQ_percentile, dat = compDat)


plot(compDat$Age, log(compDat$impTau), pch=21, 
     bg=c("red","blue")[unclass(as.factor(compDat$PI))])

## Model 4 -- examining self-reported risk-taking
plot(compDat$Age, compDat$DOSPERT, pch=21, 
          bg=c("red","blue")[unclass(as.factor(compDat$PI))])
compDat$cAge = compDat$Age - mean(compDat$Age)
compDat$cAgeSq = compDat$cAge * compDat$cAge

hist(compDat$DOSPERT)
t.test(compDat$DOSPERT ~ compDat$PI)
t.test(compDat$DOSPERT ~ compDat$Sex)
cor.test(compDat$DOSPERT, compDat$Age)
mod1 = lm(DOSPERT ~ Age + PI, data = compDat)
mod2 = lm(DOSPERT ~ Age + PI + Sex, data = compDat)
mod3 = lm(DOSPERT ~ Age*PI, data = compDat)

mod5 = lm(DOSPERT ~ cAge + cAgeSq, data = compDat)
mod6 = lm(DOSPERT ~ cAge + cAgeSq + PI, data = compDat)
mod7 = lm(DOSPERT ~ Age + IQ_percentile, data = L2Dat)

dosMod = lm(DOSPERT ~ Age + PI + Sex + IQ_percentile, data = compDat)
