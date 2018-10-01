#########################################################
#########################################################
#              Data simulation script                   #
#########################################################
#########################################################

# this is an R script that can be used to simulate a dataset that will have similar characteristics to the dataset
# expected in the research project titled: "Expectancy of the effectiveness of unconventional hypnosis techniques"
# for more details, see the research protocol posted at the project's OSF site: https://osf.io/7khc6/

# Created by Zoltan Kekecs, PhD, Researchr, Lund University, Department of Psychology

#########################################################
#                                                       #
#                      Packages                         #
#                                                       #
#########################################################

library(MASS) # for mvrnorm



#########################################################
#                                                       #
#                   Custom functions                    #
#                                                       #
#########################################################

# function to ganarate random timestamp
# adapted from https://stackoverflow.com/questions/14720983/efficiently-generate-a-random-sample-of-times-and-dates-between-two-dates
generate_timestamp <- function(N, st="2012/01/01", et="2012/12/31") {
  st <- as.POSIXct(as.Date(st))
  et <- as.POSIXct(as.Date(et))
  dt <- as.numeric(difftime(et,st,unit="sec"))
  ev <- sort(runif(N, 0, dt))
  rt <- st + ev
  return(rt)
}


#########################################################
#                                                       #
#          Parameters for simulating data               #
#                                                       #
#########################################################


### set parameters to simulate data

## study parameters

N_per_group = 240
num_conditions_per_person = 2
num_conf_questions_per_condition = 2


## outcome parameters

# variable names of the continuous outcome measures
outcome_names = c("effect_pain_con", "effect_induce_con", "choose_over_pharmanalg_con", "benefit_perc_con", "effect_pain_unc", "effect_induce_unc", "choose_over_pharmanalg_unc", "benefit_perc_unc", "effect_pain_hyp_in_general")


# correlation of the outcome variables
cor_Q1_Q2 = 0.5 #correlation between responses given to the two main questions (within the same subject)
cor_RH_PH = 0.3 #correlation between responses given to the same questions between the real hypnosis and the placebo hypnosis conditions (within the same subject)
cor_RH_Q1_PH_Q2 = 0.1 #correlation between responses given to the two main questions between the real hypnosis and the placebo hypnosis conditions (within the same subject)

# correlation matrix of the outcomes
correlation_matrix = matrix(c(1, cor_Q1_Q2, cor_Q1_Q2, cor_Q1_Q2, cor_RH_PH, cor_RH_Q1_PH_Q2, cor_RH_Q1_PH_Q2, cor_RH_Q1_PH_Q2, cor_Q1_Q2, # the meaning of the columns: effec_pain_con	effec_induce_con	choose_over_pharmanalg_con	benefit_perc_con	effec_pain_unc	effec_induce_unc	choose_over_pharmanalg_unc	benefit_perc_unc, effect_pain_hyp_in_general.
                              cor_Q1_Q2, 1, cor_Q1_Q2, cor_Q1_Q2, cor_RH_Q1_PH_Q2, cor_RH_PH, cor_RH_Q1_PH_Q2, cor_RH_Q1_PH_Q2, cor_Q1_Q2, 
                              cor_Q1_Q2, cor_Q1_Q2, 1, cor_Q1_Q2, cor_RH_Q1_PH_Q2, cor_RH_Q1_PH_Q2, cor_RH_PH, cor_RH_Q1_PH_Q2, cor_Q1_Q2, 
                              cor_Q1_Q2, cor_Q1_Q2, cor_Q1_Q2, 1, cor_RH_Q1_PH_Q2, cor_RH_Q1_PH_Q2, cor_RH_Q1_PH_Q2, cor_RH_PH, cor_Q1_Q2, 
                              cor_RH_PH, cor_RH_Q1_PH_Q2, cor_RH_Q1_PH_Q2, cor_RH_Q1_PH_Q2, 1, cor_Q1_Q2, cor_Q1_Q2, cor_Q1_Q2, cor_Q1_Q2, 
                              cor_RH_Q1_PH_Q2, cor_RH_PH, cor_RH_Q1_PH_Q2, cor_RH_Q1_PH_Q2, cor_Q1_Q2, 1, cor_Q1_Q2, cor_Q1_Q2, cor_Q1_Q2, 
                              cor_RH_Q1_PH_Q2, cor_RH_Q1_PH_Q2, cor_RH_PH, cor_RH_Q1_PH_Q2, cor_Q1_Q2, cor_Q1_Q2, 1, cor_Q1_Q2, cor_Q1_Q2, 
                              cor_RH_Q1_PH_Q2, cor_RH_Q1_PH_Q2, cor_RH_Q1_PH_Q2, cor_RH_PH, cor_Q1_Q2, cor_Q1_Q2, cor_Q1_Q2, 1, cor_Q1_Q2,
                              cor_Q1_Q2, cor_Q1_Q2, cor_Q1_Q2, cor_Q1_Q2, cor_Q1_Q2, cor_Q1_Q2, cor_Q1_Q2, cor_Q1_Q2, 1),
                            ncol = 9, byrow = TRUE)

# means and SDs of the outcomes
# based on Kendrick, C., Koep, L., Johnson, A., Fisher, W., & Elkins, G. (2013). Feasibility 
# of a sham hypnosis: Empirical data and implications for randomized trials of hypnosis. 
# Contemporary Hypnosis and Integrative Therapy, 29(4), 317-331.
# table 3, showing the ratings of people in the hypnosis and sham groups,
# extracting the mean ratings from columns "Beneficial" and "Receipt of hypnosis".
# the following ratings are expected
mean_expectancy = mean(c(4.6, 4.24, 4.68, 4.56)) # mean expectancy = mean(c(4.6, 4.24, 4.68, 4.56))  = 4.52
sd_expectancy = mean(c(0.5, 0.831, 0.476, 0.651)) # sd of expectancy = mean(c(0.5, 0.831, 0.476, 0.651)) = 0.6145
SDs = c(sd_expectancy*15, sd_expectancy*1.5, sd_expectancy*1.5, sd_expectancy*15,sd_expectancy*15, sd_expectancy*1.5, sd_expectancy*1.5, sd_expectancy*15, sd_expectancy*1.5) # the order of the data: effec_pain_con	effec_induce_con	choose_over_pharmanalg_con	benefit_perc_con	effec_pain_unc	effec_induce_unc	choose_over_pharmanalg_unc	benefit_perc_unc, effect_pain_hyp_in_general.

effect_size_difference = 0.5 # the effect size (in cohen's d) of the difference between the expectancy of conventioal and unconventional hypnosis (positive number means that conventional hypnosis evokes higher expectancy), set this to 0 to produce equal expectancy.
means = c(mean_expectancy*15, mean_expectancy*1.5, mean_expectancy*1.5, mean_expectancy*15, mean_expectancy*15-(SDs[5]*effect_size_difference), mean_expectancy*1.5-(SDs[6]*effect_size_difference), mean_expectancy*1.5-(SDs[7]*effect_size_difference), mean_expectancy*15-(SDs[8]*effect_size_difference), mean_expectancy*1.5)

p_invalid_cases = 0.4 #(because of p is 0.1 for university_student = no, 0.1 for tried_hypnosis = yes, and 0.01 for age_range = "under 18")


# variables required for analysis and eligibility check

var_list = list(data.frame(factor_level = c("unc_whitenoise", "unc_subliminal", "unc_embed"), 
                           p = c(0.333, 0.333, 0.334),
                           effect_con = c(0, 0, 0),
                           effect_unc = c(0.1, 0, -0.1))) # the effect on the mean expressed in units of standard deviation (cohen's d)
names(var_list) = "unc_type"

var_list[["age_range"]] = data.frame(factor_level = c("under 18", "18 - 24", "25 - 34", "35 - 44", "45 - 54", "55 - 64", "65 - 74", "75 - 84", "85 or older"), 
                                     p = c(0.01, 0.7, 0.2, 0.05, 0.01, 0.01, 0.01, 0.005, 0.005),
                                     effect_con = c(0, 0, 0, 0, 0, 0, 0, 0, 0),
                                     effect_unc = c(0, 0, 0, 0, 0, 0, 0, 0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)

var_list[["university_student"]] = data.frame(factor_level = c("yes", "no"), 
                                              p = c(0.9, 0.1),
                                              effect_con = c(0, 0),
                                              effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)


var_list[["tried_hypnosis"]] = data.frame(factor_level = c("yes", "no"), 
                                          p = c(0.1, 0.9),
                                          effect_con = c(0, 0),
                                          effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)



sim_mod_dem_vars = T # whether to simulate moderator and demographic variables

if(sim_mod_dem_vars == T){
  ## Set parameters of moderating and demographic variables and their effect on the expectancy of conventional and unconventional techniques
  # factor_level lists the possible levels of the variable
  # p is the probability of any of the factor levels occuring
  # effect_con is the effect of the particular factor level on the expectancy of the conventional technique expressed in standardized mean difference (cohen's d)
  # effect_unc is the effect of the particular factor level on the expectancy of the unconventional technique expressed in standardized mean difference (cohen's d)
  
  var_list[["order"]] = data.frame(factor_level = c("unc_first", "unc_second"), 
                                   p = c(0.5, 0.5),
                                   effect_con = c(-0.05, 0.05),
                                   effect_unc = c(0.05, -0.05)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["country_of_residence"]] = data.frame(factor_level = c("Sweden", "USA", "Hungary"), 
                                                  p = c(0.43, 0.14, 0.43),
                                                  effect_con = c(0, 0, 0),
                                                  effect_unc = c(0, 0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["gender"]] = data.frame(factor_level = c("male", "female", "other", "prefer no to specify"), 
                                    p = c(0.48, 0.48, 0.01, 0.01),
                                    effect_con = c(-0.05, 0.05, 0, 0),
                                    effect_unc = c(0.05, -0.05, 0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["painful_dental_experience"]] = data.frame(factor_level = c("no", "few", "many"), 
                                                       p = c(0.2, 0.6, 0.2),
                                                       effect_con = c(0, 0, 0),
                                                       effect_unc = c(0, 0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  
  var_list[["tried_hypnosis_for_pain"]] = data.frame(factor_level = c("yes", "no"), 
                                                     p = c(0.1, 0.9), # this probability will be used later to generate yes-no askwers within those who answere "yes" to the tried_hypnosis question. Anyone who answered "no" to the tried_hypnosis question will have "no" for this question as well.
                                                     effect_con = c(0, 0),
                                                     effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  
  var_list[["hyp_knowledge"]] = data.frame(factor_level = 0:10, 
                                           p = c(0.3, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05),
                                           effect_con = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                           effect_unc = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["hyp_knowledge_source_lecture"]] = data.frame(factor_level = c("yes", "no"), 
                                                          p = c(0.1, 0.9), 
                                                          effect_con = c(0, 0),
                                                          effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["hyp_knowledge_source_scilit"]] = data.frame(factor_level = c("yes", "no"), 
                                                         p = c(0.1, 0.9), 
                                                         effect_con = c(0, 0),
                                                         effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["hyp_knowledge_source_book"]] = data.frame(factor_level = c("yes", "no"), 
                                                       p = c(0.1, 0.9), 
                                                       effect_con = c(0, 0),
                                                       effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["hyp_knowledge_source_internet"]] = data.frame(factor_level = c("yes", "no"), 
                                                           p = c(0.1, 0.9), 
                                                           effect_con = c(0, 0),
                                                           effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["hyp_knowledge_source_tvmovie"]] = data.frame(factor_level = c("yes", "no"), 
                                                          p = c(0.1, 0.9), 
                                                          effect_con = c(0, 0),
                                                          effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["hyp_knowledge_source_radiopodcast"]] = data.frame(factor_level = c("yes", "no"), 
                                                               p = c(0.1, 0.9), 
                                                               effect_con = c(0, 0),
                                                               effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["hyp_knowledge_source_exp_someoneelse"]] = data.frame(factor_level = c("yes", "no"), 
                                                                  p = c(0.1, 0.9), 
                                                                  effect_con = c(0, 0),
                                                                  effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["hyp_knowledge_source_exp_onme"]] = data.frame(factor_level = c("yes", "no"), 
                                                           p = c(0.1, 0.9), # this will be replaced by the data of tried_hypnosis
                                                           effect_con = c(0, 0),
                                                           effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["hyp_knowledge_source_exp_iused"]] = data.frame(factor_level = c("yes", "no"), 
                                                            p = c(0.1, 0.9), 
                                                            effect_con = c(0, 0),
                                                            effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
  var_list[["hyp_knowledge_source_other"]] = data.frame(factor_level = c("yes", "no"), 
                                                        p = c(0.1, 0.9), 
                                                        effect_con = c(0, 0),
                                                        effect_unc = c(0, 0)) # the effect on the mean expressed in units of standard deviation (cohen's d)
  
}





#########################################################
#                                                       #
#                    Simulate data                      #
#                                                       #
#########################################################

#########################################################
#                Simulate outcome data                  #
#########################################################

# calculate total sample size
total_N = round(N_per_group * (1 + p_invalid_cases))* length(var_list[["unc_type"]][,"factor_level"])

# generate outcome data based on above specified correlation matrix
data_pre = mvrnorm(total_N, mu = rep(0, ncol(correlation_matrix)),
                   Sigma = correlation_matrix,
                   empirical = F)
for(i in 1:ncol(correlation_matrix)){
  data_pre[,i] = data_pre[,i] * SDs[i] + means[i]
}
data_pre = as.data.frame(data_pre)
names(data_pre) = outcome_names

# we suppose that those who rated the effectiveness of the conventional techniques effectiveness for pain higher than the unconvetional one's, will choose the conventional technique
data_pre$choose_unc_over_con = apply(data_pre, 1, function(x) if(as.numeric(x["effect_pain_con"]) > as.numeric(x["effect_pain_unc"])){"no"}else{"yes"})
# here we add some noise into choose_unc_over_con so it is not 1:1 correlated with the highest of effect_pain_con or effect_pain_unc
data_pre$choose_unc_over_con = apply(data_pre, 1, function(x) if(sample(1:10, 1) > 9){if(x["choose_unc_over_con"] == "yes"){"no"} else {"yes"}} else {x["choose_unc_over_con"]})




#########################################################
#              Simulate moderating variables            #
#########################################################

# generating moderating and demographic variables based on above specified parameters
other_variables = data.frame(matrix(NA, nrow = total_N, ncol = length(names(var_list))))
names(other_variables) = names(var_list)
for(i in 1:length(names(var_list))){
  other_variables[,i] = sample(var_list[[i]][,"factor_level"], replace = T, total_N, prob=var_list[[i]][,"p"])
}

if(sim_mod_dem_vars == T){
  # adjusting tried_hypnosis_for_pain so it can only have "yes" answer if tried_hypnosis is also "yes"
  other_variables[,"tried_hypnosis_for_pain"] = apply(other_variables, 1, function(x) if(x["tried_hypnosis"] == "yes"){as.character(sample(var_list[["tried_hypnosis_for_pain"]][,"factor_level"], 1, replace = T, prob=var_list[["tried_hypnosis_for_pain"]][,"p"]))} else {"no"})
  # adjusting hyp_knowledge_source_exp_onme so it matches tried_hypnosis
  other_variables[,"hyp_knowledge_source_exp_onme"] = other_variables[,"tried_hypnosis"]
}

# combining outcomes and moderating variables
data_before_effects = cbind(data_pre, other_variables)


data = data_before_effects

if(sim_mod_dem_vars == T){
  
  #############################################################
  # Adjusting the outcomes with the effects of the moderators #
  #############################################################
  
  for(i in 1:length(names(var_list))){
    print(paste("adjusting for the effect of: [", i, "] ", names(var_list)[i], sep = ""))
    for(j in 1:length(var_list[[i]][,"effect_con"])){
      effect_pain_con = SDs[1] * var_list[[i]][j,"effect_con"]
      data[,"effect_pain_con"] = apply(data, 1, function(x) if(var_list[[i]][j, "factor_level"] == x[names(var_list)[i]]){as.numeric(x["effect_pain_con"]) + effect_pain_con} else {as.numeric(x["effect_pain_con"])})
      effect_induce_con = SDs[2] * var_list[[i]][j,"effect_con"]
      data[,"effect_induce_con"] = apply(data, 1, function(x) if(var_list[[i]][j, "factor_level"] == x[names(var_list)[i]]){as.numeric(x["effect_induce_con"]) + effect_induce_con} else {as.numeric(x["effect_induce_con"])})
      effect_choose_over_con = SDs[3] * var_list[[i]][j,"effect_con"]
      data[,"choose_over_pharmanalg_con"] = apply(data, 1, function(x) if(var_list[[i]][j, "factor_level"] == x[names(var_list)[i]]){as.numeric(x["choose_over_pharmanalg_con"]) + effect_choose_over_con} else {as.numeric(x["choose_over_pharmanalg_con"])})
      effect_benefit_perc_con = SDs[4] * var_list[[i]][j,"effect_con"]
      data[,"benefit_perc_con"] = apply(data, 1, function(x) if(var_list[[i]][j, "factor_level"] == x[names(var_list)[i]]){as.numeric(x["benefit_perc_con"]) + effect_benefit_perc_con} else {as.numeric(x["benefit_perc_con"])})
      
      effect_pain_unc = SDs[5] * var_list[[i]][j,"effect_unc"]
      data[,"effect_pain_unc"] = apply(data, 1, function(x) if(var_list[[i]][j, "factor_level"] == x[names(var_list)[i]]){as.numeric(x["effect_pain_unc"]) + effect_pain_unc} else {as.numeric(x["effect_pain_unc"])})
      effect_induce_unc = SDs[6] * var_list[[i]][j,"effect_unc"]
      data[,"effect_induce_unc"] = apply(data, 1, function(x) if(var_list[[i]][j, "factor_level"] == x[names(var_list)[i]]){as.numeric(x["effect_induce_unc"]) + effect_induce_unc} else {as.numeric(x["effect_induce_unc"])})
      effect_choose_over_unc = SDs[7] * var_list[[i]][j,"effect_unc"]
      data[,"choose_over_pharmanalg_unc"] = apply(data, 1, function(x) if(var_list[[i]][j, "factor_level"] == x[names(var_list)[i]]){as.numeric(x["choose_over_pharmanalg_unc"]) + effect_choose_over_unc} else {as.numeric(x["choose_over_pharmanalg_unc"])})
      effect_benefit_perc_unc = SDs[8] * var_list[[i]][j,"effect_unc"]
      data[,"benefit_perc_unc"] = apply(data, 1, function(x) if(var_list[[i]][j, "factor_level"] == x[names(var_list)[i]]){as.numeric(x["benefit_perc_unc"]) + effect_benefit_perc_unc} else {as.numeric(x["benefit_perc_unc"])})
      
      effect_pain_hyp_in_general = SDs[9] * var_list[[i]][j,"effect_con"]
      data[,"effect_pain_hyp_in_general"] = apply(data, 1, function(x) if(var_list[[i]][j, "factor_level"] == x[names(var_list)[i]]){as.numeric(x["effect_pain_hyp_in_general"]) + effect_pain_hyp_in_general} else {as.numeric(x["effect_pain_hyp_in_general"])})
      
    }
  }
  
  
  
  #########################################################
  #                Adding final variables                 #
  #########################################################
  
  data$time_of_data_collection = generate_timestamp(N = total_N, st="2018/10/10", et="2018/12/31")
  
  data$participant_ID = paste("ID_", 1:total_N, sep = "")
  
  data$why_choose = "some text here"
  
  data$comment = "other text here"
  
  
}

#########################################################
#                    Rounding                           #
#########################################################


# rounding data to better resemble actual raw data gathered in the survey
data[,outcome_names] = round(data[,outcome_names])
