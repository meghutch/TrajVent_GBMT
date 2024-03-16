### March 16, 2024
### Updated version of model with cleaned up dataset 
### Change parameters (num_start, model_num, scaling, etc) as needed

# packages
library(tidyverse)
library(gbmt)

# load in config file
source("/projects/b1196/project_5/TrajVent_GBMT/R/config.R")

# set working directory
setwd(set_wd)

# start time
start_time = Sys.time()

# define number of random starts
num_start = 100

# define model_num as the iterative number for each "replicate" of the analysis
# we plan to rerun this script 10x
model_num = 10

# import previous list of random_seeds
random_seeds = read.csv('random_seeds.csv')

# generate random_seed 
rand_seed = runif(n = 1, min = 0, max = 10000000) %>% round()
while(rand_seed %in% random_seeds$seed) {
  rand_seed = runif(n = 1, min = 0, max = 10000000) %>% round()
}
print(paste('random seed:', rand_seed))

#for (day in c(3, 5)) {

for (day in c(5)) {
  
  # s = scaling method; 2 = standardization; 4 = logarithmic ratio to the mean
  for (s in c(4)) {
    
    ## specify time window for analysis & load data
    if(day == 3) {
      
      ## load 3 day window
      dataset <- read.csv(paste0(set_data_dir, "04_3day_intub.csv"))
      day_window <- "3day"
      
    }
    
    if(day == 5) {
      
      ### load 5 day window
      dataset <- read.csv(paste0(set_data_dir, "04_5day_intub.csv"))
      day_window <- "5day"
      
    }
    
    if(s == 4) {
      
      # add a 1 to lung compliance (when we take log during scaling, this will convert to 0)
      dataset <- dataset %>% 
        mutate(Lung_Compliance = Lung_Compliance + 1)
    }
    
    # pre-process vent parameters
    day_vent <- dataset %>% 
      select(anonymized_id, new_intubation_day, is_prolonged_14days, 
             PEEP, FiO2, Plateau_Pressure, Lung_Compliance, Minute_Ventilation)
    
    # select vent params
    varNames <- c("PEEP", "FiO2", "Plateau_Pressure", "Lung_Compliance", "Minute_Ventilation")
    
    # vent parameters
    vent_params <- day_vent %>% 
      select(anonymized_id, new_intubation_day, varNames) %>% 
      data.frame()
    
    # gbmt model
    # `d` - Positive integer value indicating the polynomial degree of group trajectories. Default is 2.
    # `ng` - Positive integer value indicating the number of groups to create. Default is 1.
    # `scaling` - normalization method; 2 = standardization
    
    
    ### Run Models
    set.seed(rand_seed)
    model_list = list()
    
    # first for loop changes the group parameter (`ng`)
    for (i in seq(from = 2, to = 8, by = 1)) {
      
      # second for loop changes the polynomial degree (`d`)
      #for (j in c(1, 2, 3)) {
      for (j in c(3)) {
        
        print(paste('Groups=', i, j, sep = '_'))
        
        model_x <- gbmt(x.names = varNames, unit = "anonymized_id", time = "new_intubation_day", 
                        d = j, ng = i, data = vent_params, pruning = TRUE, 
                        scaling = s, maxit = 200, nstart = num_start)
        
        # add model name variable with number of groups
        # model_[ng]
        model_x$group_num <- paste('model', i, j, sep = '_')
        
        # append model to list
        model_list <- append(model_list, list(model_x))
      }
    }
    
    # create date variable
    date_var = Sys.time() %>% 
      format(., "%Y%m%d_%H%M")
    
    if(s == 0) {
      
      scale_var <- "no_normalization"
    }
    
    if(s == 2) {
      
      scale_var <- "standardization"
    }
    
    if(s == 4) {
      
      scale_var <- "logarithmic ratio to the mean"
    }
    
    # run model
    save(model_list, file = paste0('results/Models_logRatio_100_restarts/', date_var, '_Model_', model_num, '_', 'gbmt.rda'))
    
  }
}
# end time
end_time = Sys.time()

# save random seed to list
new_rand_seed_df <- data.frame("model_run" = model_num,
                               "seed" = rand_seed)
# create new data frame
random_seeds_new <- rbind(random_seeds, new_rand_seed_df)

write.csv(random_seeds_new, "random_seeds.csv", row.names = FALSE)

# total run time
print('run time')
print(end_time - start_time)

# print model information
print(paste('day_window:', day_window))

print(paste('random seed:', rand_seed))

print(paste('nstart:', num_start))

print(paste('standardization:', scale_var))

print(paste('model_num:', model_num))