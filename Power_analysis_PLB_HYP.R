library(BayesFactor)
library(progress)
library(MASS) # for mvrnorm


simulate_PLC_study <- function(N_per_group,
                               num_conditions,
                               num_questions_per_condition,
                               within,
                               correlation_matrix,
                               means,
                               SDs,
                               rscale,
                               nullInterval,
                               bf_for_all_or_atleastone){ # in case of multiple question, this parameter defines 
                                                          # whether we would like power to be estimated for finding
                                                          # that answers to at least one question is significantly 
                                                          # different between conditions or significantly the same 
                                                          # (set paramater to "one"); or for finding that answers to
                                                          # all of the questions are significantly different or the same
                                                          # between conditions.

  pb$tick() # progress bar tick
  
  data_pre = mvrnorm(N_per_group, mu = rep(0, ncol(correlation_matrix)),
                                   Sigma = correlation_matrix,
                                   empirical = F)
  for(i in 1:ncol(correlation_matrix)){
    data_pre[,i] = data_pre[,i] * SDs[i] + means[i]
  }
  
  
  if(within == F){   
    
    data = as.matrix(data_pre[,1:num_questions_per_condition])
    for(i in 1:(num_conditions-1)){
      data = rbind(data, 
                   as.matrix(data_pre[,((i)*num_questions_per_condition + 1):((i)*num_questions_per_condition + num_questions_per_condition)]))
    }  
    
    data = as.data.frame(data)
    
    question_names = NA
    for(i in 1:num_questions_per_condition){
      question_names[i] = paste("question", i, sep ="")
    }
    names(data) = c(question_names)
    

    group_names = NA
    for(i in 1:num_conditions){
      group_names[i] = paste("group", i, sep ="")
    }
    
    group = factor(rep(group_names, each = N_per_group))
    data = cbind(data, group)
    

  } else {


    same_questions = matrix(NA, ncol = num_conditions, nrow = num_questions_per_condition)
    for(i in 1:num_questions_per_condition){
      same_questions[i,] = i-1 + seq(1, ncol(data_pre), num_questions_per_condition)
    }
    
    sep_by_questions = list(NA)
    for(i in 1:num_questions_per_condition){
      sep_by_questions[[i]] = data_pre[,same_questions[i,]]
    }
    
    
    combs <- combn(num_conditions,2) # create all possible pairings of conditions
    
    
    condition_pairs = NA
    for(i in 1:ncol(combs)){
      condition_pairs[i] = paste("condition_pair_", combs[1,i], "_", combs[2,i], sep = "")
    }

    
    difference_list = list(NA)
    for(i in 1:num_questions_per_condition){
    
      difference_list[[i]] = sep_by_questions[[i]][,combs[1,1]] - sep_by_questions[[i]][,combs[2,1]]
      if(ncol(combs) > 1){
        for(j in 2:ncol(combs)){
          difference_list[[i]][(length(difference_list[[i]])+1):(length(difference_list[[i]])+N_per_group)] = sep_by_questions[[i]][,combs[1,j]] - sep_by_questions[[i]][,combs[2,j]]
        }
      }
    }

    
  
  
  data = as.data.frame(matrix(unlist(difference_list), ncol = num_questions_per_condition))
  
  question_names = NA
  for(i in 1:num_questions_per_condition){
    question_names[i] = paste("question", i, sep ="")
  }
  names(data) = c(question_names)
  
  group = factor(rep(condition_pairs, each = N_per_group))
  data = cbind(data, group)
  
  }
    

  # extract BF 
  
  ps = NA

  if(nullInterval[1] == -Inf){ alternative = "less" } else if(nullInterval[2] == Inf){
    alternative = "greater" } else {alternative = "two.sided"}
  
  if(length(levels(data$group)) == 1){
    for(i in 1:num_questions_per_condition){
      ps[i] = t.test(data[,i], mu = 0,
                     alternative = alternative)$p.value
    }
  } else if(length(levels(data$group)) == 2){
    for(i in 1:num_questions_per_condition){
      ps[i] = t.test(formula = as.formula(paste(question_names[i], " ~ ", "group")),
                     alternative = alternative, data = data)$p.value
    }
  } else {
      for(i in 1:num_questions_per_condition){
      ps[i] = summary(aov(formula = as.formula(paste(question_names[i], " ~ ", "group")), data = data))[[1]][1,5]
    }
  }
  
  if(bf_for_all_or_atleastone == "all"){
    t_pval = max(ps)
  } else {
    t_pval = min(ps)
  }

  
  bfs = NA

  if(length(levels(data$group)) == 1){
    for(i in 1:num_questions_per_condition){
      bfs[i] = suppressMessages(matrix(ttestBF(data[,i], mu = 0, rscale = rscale, nullInterval = nullInterval))[1])
    }
  } else if(length(levels(data$group)) == 2){
    for(i in 1:num_questions_per_condition){
      bfs[i] = suppressMessages(matrix(ttestBF(formula = as.formula(paste(question_names[i], " ~ ", "group")), rscale = rscale, nullInterval = nullInterval, data = data))[1])
    }
  } else {
    for(i in 1:num_questions_per_condition){
      bfs[i] = matrix(anovaBF(formula = as.formula(paste(question_names[i], " ~ ", "group")), rscaleFixed = rscale, rscaleRandom = "nuisance", data = data))[1]
    }
  }

  
  if(bf_for_all_or_atleastone == "all"){
    if(mean(log(bfs)) > 0){
      bf = min(bfs)
    } else {
      bf = max(bfs)
    }
  } else {
    if(mean(log(bfs)) > 0){
      bf = max(bfs)
    } else {
      bf = min(bfs)
    }
  }
  
  return(c(1/bf, t_pval))
}



iterations = 5000


# based on Kendrick, C., Koep, L., Johnson, A., Fisher, W., & Elkins, G. (2013). Feasibility 
# of a sham hypnosis: Empirical data and implications for randomized trials of hypnosis. 
# Contemporary Hypnosis and Integrative Therapy, 29(4), 317-331.
# table 3, showing the ratings of people in the hypnosis and sham groups,
# extracting the mean ratings from columns "Beneficial" and "Receipt of hypnosis".
# the following ratings are expected
# mean expectancy = mean(c(4.6, 4.24, 4.68, 4.56))  = 4.52
# sd of expectancy = mean(c(0.5, 0.831, 0.476, 0.651)) = 0.6145

mean_expectancy = mean(c(4.6, 4.24, 4.68, 4.56))
sd_expectancy = mean(c(0.5, 0.831, 0.476, 0.651))


# ratings are equal (M0 is correct)


pb <- progress_bar$new(
  format = " simulation progress [:bar] :percent eta: :eta",
  total = iterations, clear = FALSE, width= 60)



out_pre = replicate(iterations, 
                    simulate_PLC_study(N_per_group = 240,
                                       num_conditions = 2, # if more than 2 conditions, anova will be conducted
                                       num_questions_per_condition = 2, # if more than 2 questions, separate test will be run and the least favorable or the most favorable will be saved of the bfs and the ps based on parameter bf_for_all_or_atleastone
                                       within = T,
                                       correlation_matrix = matrix(c(1, 0.5,0.3,  0.1,
                                                                     0.5, 1, 0.1, 0.3,
                                                                     0.3, 0.1, 1, 0.5,
                                                                     0.1, 0.3, 0.5, 1),
                                                                   ncol = 4, byrow = TRUE),
                                       means = c(mean_expectancy, mean_expectancy, mean_expectancy, mean_expectancy),
                                       SDs = c(sd_expectancy, sd_expectancy, sd_expectancy, sd_expectancy),
                                       rscale = "wide",
                                       nullInterval = c(0, Inf), # if two tailed, set this to NULL, if one tailed, set this to c(-Inf,0) or c(0, Inf). c(0, Inf) means we expect higher values in the first condition compared to the second. This parameter is not used if 3 or more conditions are present and an anova is run. 
                                       bf_for_all_or_atleastone = "all")) # in case of multiple question, this parameter defines 
                                                                          # whether we would like power to be estimated for finding
                                                                          # that answers to at least one question is significantly 
                                                                          # different between conditions or significantly the same 
                                                                          # (set paramater to "one"); or for finding that answers to
                                                                          # all of the questions are significantly different or the same
                                                                          # between conditions.


out = as.data.frame(t(out_pre))
names(out) = c("BF", "t_pval")

# number of placebo conditions tested in parallel:
num_plac_cond = 3

#if M0 is correct
power_M0_ture = round(mean(out[,"BF"] > 3)^num_plac_cond, 3) # power to detect similarity of ALL placebo conditions with the real hypnosis condition
falsepos_M0_ture = round(1-((1-mean(out[,"BF"] < 0.3333))^num_plac_cond), 3) # chance for falsely saying that there is at least one difference

power_M0_ture
falsepos_M0_ture



# Real hypnosis evokes higher expectancy (real difference is 0.3*SD) (M1 is correct)


pb <- progress_bar$new(
  format = " simulation progress [:bar] :percent eta: :eta",
  total = iterations, clear = FALSE, width= 60)



out_pre = replicate(iterations, 
                    simulate_PLC_study(N_per_group = 240,
                                       num_conditions = 2, # if more than 2 conditions, anova will be conducted
                                       num_questions_per_condition = 2, # if more than 2 questions, separate test will be run and the least favorable or the most favorable will be saved of the bfs and the ps based on parameter bf_for_all_or_atleastone
                                       within = T,
                                       correlation_matrix = matrix(c(1, 0.5,0.3,  0.1,
                                                                     0.5, 1, 0.1, 0.3,
                                                                     0.3, 0.1, 1, 0.5,
                                                                     0.1, 0.3, 0.5, 1),
                                                                   ncol = 4, byrow = TRUE),
                                       means = c(mean_expectancy, mean_expectancy, mean_expectancy-(sd_expectancy*0.3), mean_expectancy-(sd_expectancy*0.3)),
                                       SDs = c(sd_expectancy, sd_expectancy, sd_expectancy, sd_expectancy),
                                       rscale = "wide",
                                       nullInterval = c(0, Inf), # if two tailed, set this to NULL, if one tailed, set this to c(-Inf,0) or c(0, Inf). c(0, Inf) means we expect higher values in the first condition compared to the second. This parameter is not used if 3 or more conditions are present and an anova is run. 
                                       bf_for_all_or_atleastone = "all")) # in case of multiple question, this parameter defines 
# whether we would like power to be estimated for finding
# that answers to at least one question is significantly 
# different between conditions or significantly the same 
# (set paramater to "one"); or for finding that answers to
# all of the questions are significantly different or the same
# between conditions.


out = as.data.frame(t(out_pre))
names(out) = c("BF", "t_pval")

# number of placebo conditions tested in parallel:
num_plac_cond = 3

#if the populations DO differ
power_M1_0.3SD_ture = round(mean(out[,"BF"] < 0.3333)^num_plac_cond, 3) # power to detect difference of ALL placebo conditions compard to the real hypnosis condition
falsepos_M1_0.3SD_ture = round(1-((1-mean(out[,"BF"] > 3))^num_plac_cond), 3) # chance for falsely saying that there is AT LEAST ONE similarity

power_M1_0.3SD_ture
falsepos_M1_0.3SD_ture


# Real hypnosis evokes higher expectancy (real difference is 0.5*SD) (M1 is correct)


pb <- progress_bar$new(
  format = " simulation progress [:bar] :percent eta: :eta",
  total = iterations, clear = FALSE, width= 60)



out_pre = replicate(iterations, 
                    simulate_PLC_study(N_per_group = 240,
                                       num_conditions = 2, # if more than 2 conditions, anova will be conducted
                                       num_questions_per_condition = 2, # if more than 2 questions, separate test will be run and the least favorable or the most favorable will be saved of the bfs and the ps based on parameter bf_for_all_or_atleastone
                                       within = T,
                                       correlation_matrix = matrix(c(1, 0.5,0.3,  0.1,
                                                                     0.5, 1, 0.1, 0.3,
                                                                     0.3, 0.1, 1, 0.5,
                                                                     0.1, 0.3, 0.5, 1),
                                                                   ncol = 4, byrow = TRUE),
                                       means = c(mean_expectancy, mean_expectancy, mean_expectancy-(sd_expectancy*0.5), mean_expectancy-(sd_expectancy*0.5)),
                                       SDs = c(sd_expectancy, sd_expectancy, sd_expectancy, sd_expectancy),
                                       rscale = "wide",
                                       nullInterval = c(0, Inf), # if two tailed, set this to NULL, if one tailed, set this to c(-Inf,0) or c(0, Inf). c(0, Inf) means we expect higher values in the first condition compared to the second. This parameter is not used if 3 or more conditions are present and an anova is run. 
                                       bf_for_all_or_atleastone = "all")) # in case of multiple question, this parameter defines 
# whether we would like power to be estimated for finding
# that answers to at least one question is significantly 
# different between conditions or significantly the same 
# (set paramater to "one"); or for finding that answers to
# all of the questions are significantly different or the same
# between conditions.


out = as.data.frame(t(out_pre))
names(out) = c("BF", "t_pval")

# number of placebo conditions tested in parallel:
num_plac_cond = 3

#if the populations DO differ
power_M1_0.5SD_ture = round(mean(out[,"BF"] < 0.3333)^num_plac_cond, 3) # power to detect difference of ALL placebo conditions compard to the real hypnosis condition
falsepos_M1_0.5SD_ture = round(1-((1-mean(out[,"BF"] > 3))^num_plac_cond), 3) # chance for falsely saying that there is AT LEAST ONE similarity

power_M1_0.5SD_ture
falsepos_M1_0.5SD_ture



