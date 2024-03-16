# function to format gbmt results to a dataframe `model_stats_df`
# `gbmt_model_list` - list of gbmt results
# `model_stats_df`` - dataframe containing all model results arranged by ascending BIC
get_gbmt_stats <- function(gbmt_model_list) {
  
  # initiate dataframe for results
  model_stats_df <- data.frame('model_list_number' = c(),
                               'model_param' = c(), 
                               'logLik' = c(),
                               'aic' = c(), 
                               'bic' = c(),
                               'caic' = c(),
                               'ssbic' = c(),
                               'hqic' = c())
  
  # for each individual gbmt model, we will extract statistics
  for(i in 1:length(gbmt_model_list)) {
    
    model_number <- i
    model_ll <- gbmt_model_list[[i]][["logLik"]] 
    model_ic <- gbmt_model_list[[i]][["ic"]] %>% t() %>% data.frame()
    model_param <- gbmt_model_list[[i]][["group_num"]]
    
    model_stats <- data.frame('model_list_number' = model_number,
                              'model_param' = model_param,
                              'logLik' = model_ll)
    
    model_stats_row <- cbind(model_stats, model_ic)
    
    model_stats_df <- rbind(model_stats_df, model_stats_row)
    
  }
  
  return(model_stats_df %>% 
           arrange(bic))
  
  
}

# function to assign group number to patient clinical data
# `day_results` - dataframe of clinical results used when running model
# `gbmt_model_list` - list of gbmt results
# `group_clinical_data` - dataframe merging the group assignments to the clinical dataframe
assign_group <- function(day_results, gbmt_model_list) {
  
  # initiate empty dataframe 
  group_df <- data.frame('anonymized_id' = c(), 'group' = c())
  
  # for each individual model in our list, we will extract the group assignments (`assign.list`)
  for (i in 1:length(gbmt_model_list[["assign.list"]])) {
    
    group_df_row <- data.frame(
      anonymized_id = gbmt_model_list[["assign.list"]][[i]]) %>%
      mutate(group = paste('Group', i, sep = ' ')
             )
    
    group_df <- rbind(group_df, group_df_row) %>% 
      mutate(anonymized_id = as.integer(anonymized_id))
  }
  
  # join dataframe of group assignments to clinical outcomes by `Patient_id`
  group_clinical_data <- left_join(group_df, 
                                   day_results, 
                                   by = "anonymized_id")
  
  return(group_clinical_data)
  
}

# rename each gbmt model result with the name of the group and polynomial degree tested
# `gbmt_model_list` - list of gbmt results
# `model_list` - new list containing renamed gbmt results
rename_gbmt_results <- function(gbmt_model_list) {
  
  for (i in seq(1:length(gbmt_model_list))) {
    
    group_num <- gbmt_model_list[[i]][["group_num"]]
    
    gbmt_model_list[[group_num]] <- gbmt_model_list[[i]]
    
  }
  
  # remove original unnamed lists
  models_to_remove <- which(names(gbmt_model_list)=="")
  
  
  # Remove lists with specified indices
  gbmt_model_list <- gbmt_model_list[-models_to_remove]
  
  return(gbmt_model_list)
  
}

# calculate the odds of correct classification for each group
# https://andrewpwheeler.com/2016/10/06/group-based-trajectory-models-in-stata-some-graphs-and-fit-statistics/
calc_occ <- function(gbmt_model) {
  
  appa_results <- gbmt_model[["appa"]]
  
  # initiate empty dataframe
  occ_results <- data.frame("group" = c(), "occ" = c())
  
  for(i in seq(1:length(appa_results))) {
    
    # calculate the estimated proportion of patients for each group
    # this is calculated using the sum of the posterior probabilities
    pi_j <- sum(gbmt_model[["posterior"]][,i]) /length(gbmt_model[["posterior"]][,i])
    
    # retrieve average posterior probability of assignment
    appa_j <- appa_results[i] %>% as.numeric()
    
    # calculate numerator and denominator for occ 
    n <- appa_j / (1 - appa_j) 
    d <- pi_j / (1 - pi_j)
    
    # calculate occ
    occ <- n/d
    
    # save to individual dataframe
    occ_df <- data.frame("group" = i, "occ" = occ)
    
    # combine result to overall dataframe
    occ_results <- rbind(occ_results, occ_df)
    
  }
  
  return(occ_results)
}

# From Klijn 2017: 
# "The mismatch is the difference between the estimated group probability pi_j 
# and the real proportion of the sample assigned to group (Pj = Nj /N). 
# In cases of a perfectly fitting model, the mismatch (pi_j – Pj) score is zero"
calc_mismatch <- function(gbmt_model) {
  
  # initiate empty dataframe
  mismatch_results <- data.frame("group" = c(), "mismatch" = c())
  
  for(i in seq(1:length(appa_results))) {
    
    # calculate the estimated proportion of patients for each group
    pi_j <- sum(gbmt_model[["posterior"]][,i]) /length(gbmt_model[["posterior"]][,i])
    
    # calculate total number of patients in cohort
    tol = length(gbmt_model[["assign"]])
  
    # real number assigned to group
    N <- length(gbmt_model[["assign.list"]][[i]]) 
    
    # real proportion assigned to group
    Pj = N / tol
    
    # calculate mismatch
    mismatch = pi_j - Pj 
    
    # save to individual dataframe
    mismatch_df <- data.frame("group" = i, "mismatch" = mismatch)
    
    # combine result to overall dataframe
    mismatch_results <- rbind(mismatch_results, mismatch_df)
    
  }
  
  return(mismatch_results)
  
  
}

# calculate (SD-GMP - standard deviation of assignment probability per group, based on the indiviudals assigned to that group)
calc_sdGmp <- function(gmbt_model) {
  
  # extra dataframe of posterior probabilities
  pp_df <- gbmt_model[["posterior"]] %>% 
    data.frame() %>% 
    rownames_to_column("anonymized_id") 
  
  # select the columns representing groups; these will be used next for referencing which columns to calculate the max value
  group_names <- pp_df %>% 
    select(contains("X")) %>% 
    names() %>% 
    as.vector()
  
  # determine max probability and remove extra columns
  pp_df <- pp_df %>% 
    mutate(max_prob = select(., group_names) %>% do.call(pmax, .)) %>% 
    select(anonymized_id, max_prob)
  
  
  # initiate empty dataframe
  sdGmp_results <- data.frame("group" = c(), "sdGmp" = c())
  
  for(i in seq(1:length(appa_results))) {
  
    # create dataframe for patients assigned to each group
    group_assignments <- data.frame(
      anonymized_id = gbmt_model[["assign.list"]][[i]]) %>%
      mutate(group = paste('group', i, sep = ' ')
      ) 
    
    # filter `pp_df` to include those of the current group_assignment
    group_pp <- merge(pp_df,
                      group_assignments, by = "anonymized_id")
    
    # calculate sd
    sdGmp <- sd(group_pp$max_prob)
    
    # save to individual dataframe
    sdGmp_df <- data.frame("group" = i, "sdGmp" = sdGmp)
    
    # combine result to overall dataframe
    sdGmp_results <- rbind(sdGmp_results, sdGmp_df)
    
  }
  
  return(sdGmp_results)
}
