### testing
group_number = 4


# calculate cluster stability index
calc_csi <- function(model_list, group_number) {
  
  
  # extract all models from the same group
  group_list = list()
  
  for(i in names(model_list)) {
    
    print(i)
    
    # group_numbers start at 2 so we will subtract 1
    model_run <- model_list[[i]][[group_number-1]]
    
    group_list[[i]] <- model_run
  }
  
  # assign patients to group
  
  # initialize empty dataframe
  assigned_groups_df <- data.frame()
  
  for(i in names(model_list)) {
    
    print(i)
    
    assign_group_df <- assign_group(day5, 
                                    gbmt_model_list = group_list[[i]]) %>% 
      distinct(anonymized_id, group) %>% 
      mutate(model = paste(i))
    
    assigned_groups_df <- rbind(assigned_groups_df, assign_group_df)
    
  }
  
  # determine the proportion of patients who were assigned to the same group across all models
  assigned_groups_df_counts <- assigned_groups_df %>% 
    group_by(anonymized_id, group) %>% 
    summarize(., n_distinct(model))
  
  hist(assigned_groups_df_counts$`n_distinct(model)`)
  
  
  
  
  
  # for each pair of models, calculate the proportion of patients who were in the same group
  eval_groups <- assigned_groups_df %>% 
    filter(model %in% c("Model_1", "Model_2"))

  csi <- assigned_groups_df %>% 
    group_by(anonymized_id, group) %>% 
    summarize(., n_distinct(model))
  
  # #####
  # c = 0
  # d = 0
  # 
  # for(i in 1:(length(names(model_list))-1)) {
  #   for(j in (i + 1):length(model_list)) {
  #     
  #     print("####### reviewing concordance ######")
  #     print(paste0(names(model_list)[i], "+", names(model_list)[j]))
  # 
  #     eval_groups <- assigned_groups_df %>% 
  #       filter(model %in% c(names(model_list)[i], names(model_list)[j])) %>% 
  #       arrange(model)
  #   
  #   m1 <- assigned_groups_df %>% 
  #     filter(model == c(names(model_list)[i]))
  #   
  #   #print(m1)
  #   
  #   m2 <- assigned_groups_df %>% 
  #     filter(model == c(names(model_list)[j]))
  #            
  #   #print(m2)
  #   
  #   ids <- unique(assigned_groups_df$anonymized_id)
  #   
  #   # compare each patient to see if they are in the same group
  #   for(id1 in 1:(length(ids)-1)) {
  #     for(id2 in (i + 1):length(ids)){
  #       
  #       id_1 = ids[id1]
  #       id_2 = ids[id2]
  #       
  #       #print(paste("#### eval ids:", id_1, "vs", id_2, "#####"))
  #       
  #       m1_compare <- m1 %>% 
  #         filter(anonymized_id %in% c(id_1, id_2)) %>% 
  #         distinct(group) %>% 
  #         n_distinct()
  #       
  #       m2_compare <- m2 %>% 
  #         filter(anonymized_id %in% c(id_1, id_2)) %>% 
  #         distinct(group) %>% 
  #         n_distinct()
  #       
  #       if(m1_compare == 1 & m2_compare == 1) {
  #         
  #         c = c + 1 # concordant
  #       } else {
  #         
  #         d = d + 1 # concordant
  #       }
  #       
  #       #print(paste('concordance', c))
  #       #print(paste('discordance', d))
  #     }
  #   }
  #   }
  #   print(paste('concordance', c))
  #   print(paste('discordance', d))
  # }
  
  
  ### with help from chatGPT
  c <- 0
  d <- 0
  
  # Initialize counts
  concordance_count <- 0
  discordance_count <- 0
  
  for (i in 1:(length(names(model_list))-1)) {
    for (j in (i + 1):length(model_list)) {
      
      print("####### reviewing concordance ######")
      print(paste0(names(model_list)[i], "+", names(model_list)[j]))
      
      eval_groups <- assigned_groups_df %>% 
        filter(model %in% c(names(model_list)[i], names(model_list)[j])) %>% 
        arrange(model)
      
      m1 <- assigned_groups_df %>% 
        filter(model == c(names(model_list)[i]))
      
      m2 <- assigned_groups_df %>% 
        filter(model == c(names(model_list)[j]))
      
      ids <- unique(assigned_groups_df$anonymized_id)
      
      # compare each patient to see if they are in the same group
      for (id1 in 1:(length(ids)-1)) {
        for (id2 in (id1 + 1):length(ids)) {  # Corrected the loop iteration
          
          id_1 <- ids[id1]
          id_2 <- ids[id2]
          
          m1_compare <- m1 %>% 
            filter(anonymized_id %in% c(id_1, id_2)) %>% 
            distinct(group) %>% 
            n_distinct()
          
          m2_compare <- m2 %>% 
            filter(anonymized_id %in% c(id_1, id_2)) %>% 
            distinct(group) %>% 
            n_distinct()
          
          if (m1_compare == 1 & m2_compare == 1) {
            c <- c + 1 # concordant
          } else {
            d <- d + 1 # discordant
          }
        }
      }
      # Update counts for each pair of models
      concordance_count <- concordance_count + c
      discordance_count <- discordance_count + d
      
      # Reset counts for the next pair of models
      c <- 0
      d <- 0
    }
  }
  # Output final counts
  print(paste('Total Concordance:', concordance_count))
  print(paste('Total Discordance:', discordance_count))
  

  
  # concordance
  c_stat <- concordance_count / (concordance_count+discordance_count)
  
  # Number of models
  total_models <- 10
  
  # Number of patients
  total_patients <- 291
  
  # Calculate combinations of models
  combinations_models <- choose(total_models, 2)
  
  # Calculate combinations of patients
  combinations_patients <- choose(total_patients, 2)
  
  # Total combinations possible
  total_combinations <- combinations_models * combinations_patients
  
  # Print the result
  print(total_combinations)
  
  total_combinations == (concordance_count+discordance_count)
  
  
  
  
  for (i in length(names(model_list)-1)[-length(model_list)]) {
    
    for (j in i+1:length(model_list)names(model_list)[length(model_list)) {
      
      print(paste(i, "+", j))
    }
  }
  
  
  for (i in 1:(length(model_list) - 1)) {
    
    model_num = model_list[i]
    
    print(model_num)
  }
    
    for (j in (i + 1):length(combs)) {
      
      print(paste(i, "+", j))
    }
  }
  
 
  
  #proportion of patients assigned to the same cluster in 2 models
  
  #akerlund
  # CSI was calculated between all pairs of the models of the same number of clusters, 
  # and median and interquartile range (IQR) was calculated. As the CSI, when numbers of 
  # clusters < < number of patients, by nature is higher for lower number of clusters, 
  # a penalty for the number of clusters was added by subtracting 1/n clusters from all median CSI. 
  # The optimal clustering was defined as number of clusters with the highest median CSI 
  # representing the most stable number of clusters). When describing the clusters, 
  # the model of the optimal number of clusters with the highest log-likelihood was chosen to represent 
  # our model.
  
  # Median CSI = 1 indicates perfect match, 0 indicates no matches between different models
  
  # import all models within results/Models_logRatio_100_restarts to calculate CSI
  
}
