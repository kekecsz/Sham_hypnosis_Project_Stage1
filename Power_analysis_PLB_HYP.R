#########################################################
#########################################################
#                  Power analysis script                #
#########################################################
#########################################################

# this is an R script that can be used to conduct a Monte-Carlo-simulation-based analysis of
# the study's operational characteristics in the research project 
# titled: "Expectancy of the effectiveness of unconventional hypnosis techniques"
# for more details, see the research protocol posted at the project's OSF site: https://osf.io/7khc6/

# Created by Zoltan Kekecs, PhD, Researchr, Lund University, Department of Psychology

#########################################################
#                                                       #
#                      Packages                         #
#                                                       #
#########################################################

library(MASS) # for mvrnorm
library(progress) # for progress bar

#########################################################
#                                                       #
#            Set parameters for the simulation          #
#                                                       #
#########################################################


# set parameters of the simulation power analysis

# number of iterations
iterations = 10000 

# the effect sizes to try in this stimulation (in cohen's d) of the difference 
# between the expectancy of conventioal and unconventional hypnosis 
# (positive number means that conventional hypnosis evokes higher expectancy), 
# set this to 0 to produce equal expectancy.
effect_size_difference_to_try = c(0, 0.3, 0.4, 0.5) 




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



simul_PLB_study <- function(N_per_group = 240,
                            cor_Q1_Q2 = 0.5, #correlation between responses given to the two main questions (within the same subject)
                            cor_RH_PH = 0.3, #correlation between responses given to the same questions between the real hypnosis and the placebo hypnosis conditions (within the same subject)
                            cor_RH_Q1_PH_Q2 = 0.1, #correlation between responses given to the two main questions between the real hypnosis and the placebo hypnosis conditions (within the same subject)
                            effect_size_difference = 0, # the effect size (in cohen's d) of the difference between the expectancy of conventioal and unconventional hypnosis (positive number means that conventional hypnosis evokes higher expectancy), set this to 0 to produce equal expectancy.
                            p_invalid_cases = 0.4, #(because of p is 0.1 for university_student = no, 0.1 for tried_hypnosis = yes, and 0.01 for age_range = "under 18")
                            sim_mod_dem_vars = F, # whether to simulate moderator and demographic variables
                            show_graphs = F){
 
  pb$tick() # progress bar tick
   
  #########################################################
  #                                                       #
  #          Parameters for simulating data               #
  #                                                       #
  #########################################################
  
  
  ### set parameters to simulate data
  
  ## study parameters
  
  N_per_group = N_per_group
  num_conditions_per_person = 2
  num_conf_questions_per_condition = 2
  
  
  ## outcome parameters
  
  # variable names of the continuous outcome measures
  outcome_names = c("effect_pain_con", "effect_induce_con", "choose_over_pharmanalg_con", "benefit_perc_con", "effect_pain_unc", "effect_induce_unc", "choose_over_pharmanalg_unc", "benefit_perc_unc", "effect_pain_hyp_in_general")
  
  
  # correlation of the outcome variables
  cor_Q1_Q2 = cor_Q1_Q2 #correlation between responses given to the two main questions (within the same subject)
  cor_RH_PH = cor_RH_PH #correlation between responses given to the same questions between the real hypnosis and the placebo hypnosis conditions (within the same subject)
  cor_RH_Q1_PH_Q2 = cor_RH_Q1_PH_Q2 #correlation between responses given to the two main questions between the real hypnosis and the placebo hypnosis conditions (within the same subject)
  
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
  
  effect_size_difference = effect_size_difference # the effect size (in cohen's d) of the difference between the expectancy of conventioal and unconventional hypnosis (positive number means that conventional hypnosis evokes higher expectancy), set this to 0 to produce equal expectancy.
  means = c(mean_expectancy*15, mean_expectancy*1.5, mean_expectancy*1.5, mean_expectancy*15, mean_expectancy*15-(SDs[5]*effect_size_difference), mean_expectancy*1.5-(SDs[6]*effect_size_difference), mean_expectancy*1.5-(SDs[7]*effect_size_difference), mean_expectancy*15-(SDs[8]*effect_size_difference), mean_expectancy*1.5)
  
  p_invalid_cases = p_invalid_cases
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
  
  
  
  sim_mod_dem_vars = sim_mod_dem_vars # whether to simulate moderator and demographic variables
  
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
  
  
  
  #########################################################
  #########################################################
  #                    Analyse data                       #
  #########################################################
  #########################################################
  
  #########################################################
  #                                                       #
  #                      Packages                         #
  #                                                       #
  #########################################################
  
  library(BayesFactor) # for Bayesian analyses
  library(psych) # for descriptives and to determine skew
  
  
  
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
  
  # only use the first N_per_group observations from each group as per pre-registration
  # if there are N_per_group observations or more
  if(nrow(data_embed_conf) > N_per_group){data_embed_conf = data_embed_conf[1:N_per_group,]}
  if(nrow(data_subliminal_conf) > N_per_group){data_subliminal_conf = data_subliminal_conf[1:N_per_group,]}
  if(nrow(data_whitenoise_conf) > N_per_group){data_whitenoise_conf = data_whitenoise_conf[1:N_per_group,]}
  
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
  
  show_graphs = show_graphs
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
  
  return(c(conf_test_H, effect_size_difference))
}








#########################################################
#########################################################
#              Run study simulations                    #
#########################################################
#########################################################


# set up table where results of the simulation power analysis will be stored
operational_characteristics = as.data.frame(matrix(NA, nrow = 4, ncol = 4))
names(operational_characteristics) = c("effect_size_difference", "p_accept_hip", "p_decline_hip", "p_inconclusive")






for(i in 1:length(effect_size_difference_to_try)){
  print(paste("Now simulating studies with effect_size_difference =", effect_size_difference_to_try[i]))
  
  pb <- progress_bar$new(
    format = " simulation progress [:bar] :percent eta: :eta",
    total = iterations, clear = FALSE, width= 60)
  
  
  
  out = replicate(iterations, simul_PLB_study(N_per_group = 240,
                                                    cor_Q1_Q2 = 0.5, #correlation between responses given to the two main questions (within the same subject)
                                                    cor_RH_PH = 0.3, #correlation between responses given to the same questions between the real hypnosis and the placebo hypnosis conditions (within the same subject)
                                                    cor_RH_Q1_PH_Q2 = 0.1, #correlation between responses given to the two main questions between the real hypnosis and the placebo hypnosis conditions (within the same subject)
                                                    effect_size_difference = effect_size_difference_to_try[i], # the effect size (in cohen's d) of the difference between the expectancy of conventioal and unconventional hypnosis (positive number means that conventional hypnosis evokes higher expectancy), set this to 0 to produce equal expectancy.
                                                    p_invalid_cases = 0.4, #(because of p is 0.1 for university_student = no, 0.1 for tried_hypnosis = yes, and 0.01 for age_range = "under 18")
                                                    sim_mod_dem_vars = F, # whether to simulate moderator and demographic variables
                                                    show_graphs = F
  ))
  
  operational_characteristics[i, "effect_size_difference"] = out[2,1]
  # chance for concluding the hypothesis is false
  operational_characteristics[i, "p_accept_hip"] = table(out[1,])["TRUE"]/sum(table(out[1,]))
  # chance for concluding the hypothesis is false
  operational_characteristics[i, "p_decline_hip"] = table(out[1,])["FALSE"]/sum(table(out[1,]))
  # chance for inconclusive result
  operational_characteristics[i, "p_inconclusive"] = table(out[1,])["INCONCLUSIVE"]/sum(table(out[1,]))
}



operational_characteristics
