
#########################################################
#########################################################
#                    Analyse data                       #
#########################################################
#########################################################

# this is an R script that can be used to to run confirmatory analysis planned in 
# the research project titled: "Expectancy of the effectiveness of unconventional hypnosis techniques"
# for more details, see the research protocol posted at the project's OSF site: https://osf.io/7khc6/

# Created by Zoltan Kekecs, PhD, Researchr, Lund University, Department of Psychology

#########################################################
#                                                       #
#                      Packages                         #
#                                                       #
#########################################################

library(BayesFactor) # for Bayesian analyses
library(psych) # for descriptives and to determine skew

#########################################################
#                    Load data                          #
#########################################################

# you can either load the data from an external file here, which should be stored in the R object called data
# or data can be simulated with the simulation script for this study, for example by running to following code:

# to simulate no difference between conventional and unconventional techniques
# source("https://raw.githubusercontent.com/kekecsz/Sham_hypnosis_Project_Stage1/master/Data_simulation_PLB_HYP_nodif.R")

# to simulate a clear difference between conventional and unconventional techniques (conventional evoking higher expectancy)
# source("https://raw.githubusercontent.com/kekecsz/Sham_hypnosis_Project_Stage1/master/Data_simulation_PLB_HYP_dif.R")

#########################################################
#           Set current stopping point                  #
#########################################################

# Stopping points for analysis expressed in eligiple participants in each of the 3 groups (150 means that the analysis will be performed when N = 150 eligible participants in each of the 3 grups)
# this is either 150, 200, or 240

analyze_at_N = 150

#########################################################
#     Compute Expectancy difference scores              #
#########################################################

# Compute difference scores between expectancy of conventional and unconventional technique
data$diff_effect_pain = data[,"effect_pain_con"] - data[,"effect_pain_unc"]
data$diff_effect_induce = data[,"effect_induce_con"] - data[,"effect_induce_unc"]
data$diff_choose_over_pharmanalg = data[,"choose_over_pharmanalg_con"] - data[,"choose_over_pharmanalg_unc"]
data$diff_benefit_perc = data[,"benefit_perc_con"] - data[,"benefit_perc_unc"]



#########################################################
#         Define eligible participants                  #
#########################################################

# All anayses will be performed on
# - adults
# create a vector with row numbers containing adults only
adult = data[,"age_range"] != "under 18"
data_expl = data[adult, ]

# The confirmatory analysis will be conducted on the sample of 
# - adult 
# - university students 
# - who did not use hypnosis before.
eligible_for_confanal = data[,"age_range"] != "under 18" & data[, "university_student"] == "yes" & data[, "tried_hypnosis"] == "no"
data_conf = data[eligible_for_confanal,]



#########################################################
#  Create subsamples by group for confifmatory analysis #
#########################################################

# create subgroups by unconventional technique group
data_embed_conf = data_conf[data_conf[,"unc_type"] == "unc_embed",]
data_subliminal_conf = data_conf[data_conf[,"unc_type"] == "unc_subliminal",]
data_whitenoise_conf = data_conf[data_conf[,"unc_type"] == "unc_whitenoise",]

# only use the first analyze_at_N observations from each group as per pre-registration
# if there are analyze_at_N observations or more
if(nrow(data_embed_conf) > analyze_at_N){data_embed_conf = data_embed_conf[1:analyze_at_N,]}
if(nrow(data_subliminal_conf) > analyze_at_N){data_subliminal_conf = data_subliminal_conf[1:analyze_at_N,]}
if(nrow(data_whitenoise_conf) > analyze_at_N){data_whitenoise_conf = data_whitenoise_conf[1:analyze_at_N,]}

# Descriptives
# percentage of eligible participants for confirmatory analysis
nrow(data[eligible_for_confanal,])/nrow(data)  


#########################################################
#         Descriptives and model diagnostics            #
#########################################################

### Explore normality of the difference scores
describe(data_embed_conf[, c("diff_effect_pain", "diff_effect_induce")])
describe(data_subliminal_conf[, c("diff_effect_pain", "diff_effect_induce")])
describe(data_whitenoise_conf[, c("diff_effect_pain", "diff_effect_induce")])

show_graphs = T
if(show_graphs == T){
  hist(data_embed_conf[,"diff_effect_pain"], breaks = 20)
  qqnorm(data_embed_conf[,"diff_effect_pain"])
  qqline(data_embed_conf[,"diff_effect_pain"])
  
  hist(data_embed_conf[,"diff_effect_induce"])
  qqnorm(data_embed_conf[,"diff_effect_induce"])
  qqline(data_embed_conf[,"diff_effect_induce"])
  
  hist(data_subliminal_conf[,"diff_effect_pain"], breaks = 20)
  qqnorm(data_subliminal_conf[,"diff_effect_pain"])
  qqline(data_subliminal_conf[,"diff_effect_pain"])
  
  hist(data_subliminal_conf[,"diff_effect_induce"])
  qqnorm(data_subliminal_conf[,"diff_effect_induce"])
  qqline(data_subliminal_conf[,"diff_effect_induce"])
  
  hist(data_whitenoise_conf[,"diff_effect_pain"], breaks = 20)
  qqnorm(data_whitenoise_conf[,"diff_effect_pain"])
  qqline(data_whitenoise_conf[,"diff_effect_pain"])
  
  hist(data_whitenoise_conf[,"diff_effect_induce"])
  qqnorm(data_whitenoise_conf[,"diff_effect_induce"])
  qqline(data_whitenoise_conf[,"diff_effect_induce"])
}




#########################################################
#             Transformation for analysis               #
#########################################################

# As described in the pre-registration, if during the model diagnostics it is 
# suspected that skewness might violate the assumptions of the statistical model,
# we will use an appropriate transformation to normalize data. 

# As suggested by Tabachnick and Fidell (2007) and Howell (2007), the following guidelines
# will be used when transforming data:
# If your data distribution is...

# Moderately positive skewness: Square-Root
# NEWX = SQRT(X)

# Substantially positive skewness: Logarithmic (Log 10)
# NEWX = LG10(X)

# Substantially positive skewness (with zero values): Logarithmic (Log 10)
# NEWX = LG10(X + C)

# Moderately negative skewness: Square-Root
# NEWX = SQRT(K - X)

# Substantially negative skewness Logarithmic (Log 10)
# NEWX = LG10(K - X)

# where C = a constant added to each score so that the smallest score is 1.
# K = a constant from which each score is subtracted so that the smallest score is 1; usually equal to the largest score + 1.

# References
# Howell, D. C. (2007). Statistical methods for psychology (6th ed.). Belmont, CA: Thomson Wadsworth.
# Tabachnick, B. G., & Fidell, L. S. (2007). Using multivariate statistics (5th ed.). Boston: Allyn and Bacon.


# for example

# data_embed_conf$diff_effect_pain_transf = log(max(data_embed_conf[, "diff_effect_pain"])+1-data_embed_conf[, "diff_effect_pain"])
# data_embed_conf$diff_effect_induce_transf = log(max(data_embed_conf[, "diff_effect_induce"])+1-data_embed_conf[, "diff_effect_induce"])
# data_subliminal_conf$diff_effect_pain_transf = log(max(data_subliminal_conf[, "diff_effect_pain"])+1-data_subliminal_conf[, "diff_effect_pain"])
# data_subliminal_conf$diff_effect_induce_transf = log(max(data_subliminal_conf[, "diff_effect_induce"])+1-data_subliminal_conf[, "diff_effect_induce"])
# data_whitenoise_conf$diff_effect_pain_transf = log(max(data_whitenoise_conf[, "diff_effect_pain"])+1-data_whitenoise_conf[, "diff_effect_pain"])
# data_whitenoise_conf$diff_effect_induce_transf = log(max(data_whitenoise_conf[, "diff_effect_induce"])+1-data_whitenoise_conf[, "diff_effect_induce"])

#########################################################
#                 Confirmatory analysis                 #
#########################################################

# table that will be filled with confirmatory results

conf_results = as.data.frame(matrix(NA, nrow = 3, ncol = 2))
names(conf_results) = c("diff_effect_pain", "diff_effect_induce")
row.names(conf_results) = c("embed", "subliminal", "whitenoise")

# confirmatory analysis on non-transformed data
# rscale = 1 corresponds to the Cauchy prior distribution with a scaling constant rscale = "wide" in the BayesFactor package
# in case of large differences, approximation is used, this warning message is suppressed by suppressMessages()
conf_results["embed", "diff_effect_pain"] = suppressMessages(1/matrix(ttestBF(data_embed_conf[,"diff_effect_pain"], mu = 0, rscale = 1, nullInterval = c(0, Inf)))[1])
conf_results["embed", "diff_effect_induce"] = suppressMessages(1/matrix(ttestBF(data_embed_conf[,"diff_effect_induce"], mu = 0, rscale = 1, nullInterval = c(0, Inf)))[1])

conf_results["subliminal", "diff_effect_pain"] = suppressMessages(1/matrix(ttestBF(data_subliminal_conf[,"diff_effect_pain"], mu = 0, rscale = 1, nullInterval = c(0, Inf)))[1])
conf_results["subliminal", "diff_effect_induce"] = suppressMessages(1/matrix(ttestBF(data_subliminal_conf[,"diff_effect_induce"], mu = 0, rscale = 1, nullInterval = c(0, Inf)))[1])

conf_results["whitenoise", "diff_effect_pain"] = suppressMessages(1/matrix(ttestBF(data_whitenoise_conf[,"diff_effect_pain"], mu = 0, rscale = 1, nullInterval = c(0, Inf)))[1])
conf_results["whitenoise", "diff_effect_induce"] = suppressMessages(1/matrix(ttestBF(data_whitenoise_conf[,"diff_effect_induce"], mu = 0, rscale = 1, nullInterval = c(0, Inf)))[1])


# confirmatory analysis on transformed data if necessary
# If transformation is necessary, analysis results will be reported both with 
# and without transformation, but if skewness is considered problematic, 
# the analysis of the transformed data will be used for final statistical inference 

# conf_results["embed", "diff_effect_pain"] = 1/matrix(ttestBF(data_embed_conf[,"diff_effect_pain_transf"], mu = 0, rscale = 1, nullInterval = c(0, Inf)))[1]
# conf_results["embed", "diff_effect_induce"] = 1/matrix(ttestBF(data_embed_conf[,"diff_effect_induce_transf"], mu = 0, rscale = 1, nullInterval = c(0, Inf)))[1]

# conf_results["subliminal", "diff_effect_pain"] = 1/matrix(ttestBF(data_subliminal_conf[,"diff_effect_pain_transf"], mu = 0, rscale = 1, nullInterval = c(0, Inf)))[1]
# conf_results["subliminal", "diff_effect_induce"] = 1/matrix(ttestBF(data_subliminal_conf[,"diff_effect_induce_transf"], mu = 0, rscale = 1, nullInterval = c(0, Inf)))[1]

# conf_results["whitenoise", "diff_effect_pain"] = 1/matrix(ttestBF(data_whitenoise_conf[,"diff_effect_pain_transf"], mu = 0, rscale = 1, nullInterval = c(0, Inf)))[1]
# conf_results["whitenoise", "diff_effect_induce"] = 1/matrix(ttestBF(data_whitenoise_conf[,"diff_effect_induce_transf"], mu = 0, rscale = 1, nullInterval = c(0, Inf)))[1]

conf_results

#########################################################
#               Statistical inference                   #
#########################################################


conf_test_H = if(
  # We will conclude that our primary hypothesis is supported if there is 
  # **at least one** unconventional technique where M0 is supported (Bayes Factor (0,1) is higher
  # than 3) for both the expectancy of pain reduction and hypnosis effectiveness.
  sum(apply(conf_results, 1, min)>3)>0){"TRUE"} else if(
    
    # We will conclude that our primary hypothesis is rejected if for **all three** 
    # unconventional techniques, M1 is supported (Bayes Factor (0,1) is lower than 
    # 1/3) for both the expectancy of pain reduction and hypnosis effectiveness.
    sum(apply(conf_results, 1, max)<1/3)>2){"FALSE"} else {
      # In all other cases, we will conclude that the study did not yield conclusive 
      # evidence to support or reject the main hypothesis.
      "INCONCLUSIVE"
    }


### final statistical inference abut the study hypothesis.

# We hypothesize that at least one of the unconventional hypnosis techniques tested in this study will have a 
# comparable (or higher) expected effectiveness, in both pain reduction and hypnosis induction, to the 
# conventional hypnosis technique tested in this study among people who did not try hypnosis before.

# if conf_test_H is "TRUE", we conclude that this primary hypothesis is supported
# if conf_test_H is "FALSE", we conclude that this primary hypothesis is rejected
# In all other cases, we will conclude that the study did not yield conclusive 
# evidence to support or reject the main hypothesis.

print(conf_test_H)

