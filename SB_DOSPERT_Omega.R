library(readr)
library(psych)
library(dplyr)

##This script obtains model-based reliability estimates for the SB wave 1 DOSPERT data (i.e., omega and omega-hierarchical)
##Written by JF Guassi Moreira, 04.13.2020

sbPath = file.path("P:", "Social_Behavior_(SB_study)", "Data") #set the path to the SB data

dosDat = read_csv(sprintf("%s/Questionnaire_data/wave1/SB_questionnaires_all_ages_CLEANMasterCopy.csv", sbPath)) #load in the wave 1, clean quesitonnaire data

#Helper function to get the data we want from the big master data spreadsheet
gather_dat = function(labs, dat) {
  
  out = dat[,labs]
  out = drop_na(out)
  return(out)
  
}

#Get labels for each version of the administered quesitonnaire

#Likelihood of engaging
childDosLikeLabs = paste('DOSPERTc_L', 1:40, sep = '') #40 child items. To do: get rid of hard coded item number
adolDosLikeLabs = paste('DOSPERTa_L', 1:39, sep = '') #39 teen/adult items. To do: get rid of hard coded item number

#Perceived riskiness of engaging
childDosRiskLabs = paste('DOSPERTc_R', 1:40, sep = '') #40 child items. To do: get rid of hard coded item number
adolDosRiskLabs = paste('DOSPERTa_R', 1:39, sep = '') #39 teen/adult items. To do: get rid of hard coded item number

#Perceived benefits of engaging 
childDosBenefitLabs = paste('DOSPERTc_B', 1:40, sep = '') #40 child items. To do: get rid of hard coded item number
adolDosBenefitLabs = paste('DOSPERTa_B', 1:39, sep = '') #39 teen/adult items. To do: get rid of hard coded item number

#Compile a list of subsets we want to get omega values for
#Child:Teen crossed with Likelihood:Riskiness:Benefits
subSets = c("childDosLike", "childDosRisk", "childDosBenefit",
               "adolDosLike", "adolDosRisk", "adolDosBenefit")

#Loop over each term in the list above, extract the relevant data, and save to a dataframe
for (set in subSets) {
  
  assign(set, gather_dat(get(paste(set, "Labs", sep="")), dosDat))
  
}

#Loop each set again, this time submitting the data to the omega function, which computes model-based reliability
for (set in subSets) {
  
  assign(paste(set, "omega", sep=""), omega(get(set), 5, plot=FALSE, fm = "ml"))
  #run the bifactor EFA in omega. pull out 5 factors (corresponding to the DOSPERT subscales), use maximum likelihood estimation
  
}

#put our results into a nice matrix
omegaMat = matrix(c(rep(subSets, each = 2),
                    rep(c("O_T", "O_H"), length(subSets)),
                    childDosLikeomega$omega.tot,
                    childDosLikeomega$omega_h,
                    childDosRiskomega$omega.tot,
                    childDosRiskomega$omega_h,
                    childDosBenefitomega$omega.tot,
                    childDosBenefitomega$omega_h,
                    adolDosLikeomega$omega.tot,
                    adolDosLikeomega$omega_h,
                    adolDosRiskomega$omega.tot,
                    adolDosRiskomega$omega_h,
                    adolDosBenefitomega$omega.tot,
                    adolDosBenefitomega$omega_h), nrow = 12, ncol = 3, byrow = FALSE)
print(omegaMat)
