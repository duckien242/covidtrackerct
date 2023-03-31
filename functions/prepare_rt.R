
prepare_rt <- function(variant_file,case_file) {

    var_merge = variant_file %>%
      left_join(case_file, by = c("week","state")) %>% #merges covidestim data with our data and keeps only 
      filter(!is.na(infections)) %>%
      select(-c(n.Other,freq.Other)) %>%
      select_if(function(col)max(col) != 0) #drops variants that dont have any infections otherwise estimate_R will throw an error
    
    #create 7 replicates of week to create days
    var_merge2 <- var_merge[rep(seq_len(nrow(var_merge)), 7), ]
    
    var_merge2_main <- var_merge2 %>%
      arrange(week) %>% #arrange by week for creating the days
      mutate(days = seq.Date(min(week), min(week)+nrow(.)-1, by = "days")) %>% #create variable for days
      mutate_if(is.numeric,
                .funs = ~./7) %>% #divide by 7 to convert weekly total to avg daily
      mutate_at(vars(starts_with("freq")), 
                .funs = ~.*7) %>% #revert freq back to original because freq is for the entire week
      mutate_if(is.numeric, 
                .funs = ~zoo::rollmean(., k = 7, fill = 0)) %>% #divide by 7 to get avg daily, roll mean to smooth
      mutate_at(vars(starts_with("freq")), 
                .funs = ~ifelse(. <= 0.02, 0, .)) %>% #remove 1 off sequences of variants to prevent wide Rt may need to be edited in the future to include
      mutate(across(starts_with("freq"), 
                    .fns = list("infections" = ~.*infections))) %>% #calculate variant specific number of infections
      mutate_at(vars(ends_with("infections")), 
                .funs = ~round(.,0)) #round infections to whole individuals 
    
    var_merge2_lowci <- var_merge2 %>%
      arrange(week) %>% #arrange by week for creating the days
      mutate(days = seq.Date(min(week), min(week)+nrow(.)-1, by = "days")) %>% #create variable for days
      mutate_if(is.numeric,
                .funs = ~./7) %>% #divide by 7 to convert weekly total to avg daily
      mutate_at(vars(starts_with("freq")), 
                .funs = ~.*7) %>% #revert freq back to original because freq is for the entire week
      mutate_if(is.numeric, 
                .funs = ~zoo::rollmean(., k = 7, fill = 0)) %>% #divide by 7 to get avg daily, roll mean to smooth
      mutate_at(vars(starts_with("freq")), 
                .funs = ~ifelse(. <= 0.02, 0, .)) %>% #remove 1 off sequences of variants to prevent wide Rt may need to be edited in the future to include
      mutate(across(starts_with("freq"), 
                    .fns = list("025" = ~.*infections_p2_5))) %>% #calculate variant specific number of infections (95lci)
      mutate_at(vars(ends_with("infections")), 
                .funs = ~round(.,0)) #round infections to whole individuals 
    
    var_merge2_upci <- var_merge2 %>%
      arrange(week) %>% #arrange by week for creating the days
      mutate(days = seq.Date(min(week), min(week)+nrow(.)-1, by = "days")) %>% #create variable for days
      mutate_if(is.numeric,
                .funs = ~./7) %>% #divide by 7 to convert weekly total to avg daily
      mutate_at(vars(starts_with("freq")), 
                .funs = ~.*7) %>% #revert freq back to original because freq is for the entire week
      mutate_if(is.numeric, 
                .funs = ~zoo::rollmean(., k = 7, fill = 0)) %>% #divide by 7 to get avg daily, roll mean to smooth
      mutate_at(vars(starts_with("freq")), 
                .funs = ~ifelse(. <= 0.02, 0, .)) %>% #remove 1 off sequences of variants to prevent wide Rt may need to be edited in the future to include
      mutate(across(starts_with("freq"), 
                    .fns = list("975" = ~.*infections_p97_5))) %>% #calculate variant specific number of infections (95uci)
      mutate_at(vars(ends_with("infections")), 
                .funs = ~round(.,0)) #round infections to whole individuals 
    
    var_merge3_main = var_merge2_main %>%
      select(days, ends_with(c("infections"))) %>% #remove unnecessary columns
      select_if(function(col) max(col) != 0) %>% #removes variants with no infections (after freq <= -.0.02) and thus no frequency
      pivot_longer(cols = c(starts_with("freq")), values_to = "I", names_to = "variant") %>% #pivot longer to split by variant for estimate_R
      mutate(variant = str_remove(variant, "freq\\.")) %>% #removed because it is annoying
      arrange(variant)
    
    var_merge3_lowci = var_merge2_lowci %>%
      select(days, ends_with(c("025"))) %>% #remove unnecessary columns
      select_if(function(col) max(col) != 0) %>% #removes variants with no infections (after freq <= -.0.02) and thus no frequency
      pivot_longer(cols = c(starts_with("freq")), values_to = "I", names_to = "variant") %>% #pivot longer to split by variant for estimate_R
      mutate(variant = str_remove(variant, "freq\\.")) %>% #removed because it is annoying
      mutate(variant = str_remove(variant, ".025")) %>% #removed because it is annoying
      arrange(variant)
    
    var_merge3_upci = var_merge2_upci %>%
      select(days, ends_with(c("975"))) %>% #remove unnecessary columns
      select_if(function(col) max(col) != 0) %>% #removes variants with no infections (after freq <= -.0.02) and thus no frequency
      pivot_longer(cols = c(starts_with("freq")), values_to = "I", names_to = "variant") %>% #pivot longer to split by variant for estimate_R
      mutate(variant = str_remove(variant, "freq\\.")) %>% #removed because it is annoying
      mutate(variant = str_remove(variant, ".975")) %>% #removed because it is annoying
      arrange(variant)
    
    colnames(var_merge3_lowci)[which(colnames(var_merge3_lowci)  == "I" )] <- 'I_025'
    colnames(var_merge3_upci)[which(colnames(var_merge3_upci)  == "I" )] <- 'I_975'
    
    var_merge3 = list(var_merge3_main,var_merge3_lowci,var_merge3_upci)
    return(var_merge3)
    
}