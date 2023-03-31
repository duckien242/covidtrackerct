#generates Rt calculation, smooths the line then merges with 
#the other variant data in a dataframe
#link to paper (Cori et al, 2013) : https://academic.oup.com/aje/article/178/9/1505/89262

#run for everything -Other
rt_fun= function(df){
  
  non0 <- min(which(df$I > 0)) #1st day with infections of variant to start the R estimate otherwise R estimate artificially high
  end0 <- max(which(df$I > 0)) #last day with infections of variant to start the R estimate otherwise R estimate artificially high
  df2 = df[non0:end0,] #dataframe filtered where there is the first case of variant
  
  
  #input of interval for R estimate
  #estimate burn in period for R estimate: t_end - t_start
  t_start<-seq(2, nrow(df2)-15) 
  t_end<- t_start + 15
  
  # configuration of input data for R estimate
  # essentially we are estimating the serial intervals of SARS-CoV-2 (Omicron variant specific) by drawing from two ( truncated normal ) distributions for the mean and standard deviation of the serial interval
  config <- make_config(list(mean_si = 3.5, std_mean_si = 1, min_mean_si = 1, max_mean_si = 6, # estimates for SARS-CoV-2 serial interval
                             std_si = 1, std_std_si = 0.5, min_std_si = 0.5, max_std_si = 1.5,
                             n1= 80, n2=20, t_start=t_start, t_end=t_end)
  )
  
  # main R estimate function:
  mean_Rt = estimate_R(df2$I, #will search for column named I which was created in the ci_fun but explicitly named here
                       method="uncertain_si",
                       config = config)
  
  #adds back in days that were filtered out to match the days in the main dataframe
  mean_Rt$R$t_start = mean_Rt$R$t_start + non0 
  mean_Rt$R$t_end = mean_Rt$R$t_end + non0
  
  # #binds them into a dataframe
  rt_df<-cbind.data.frame(day = mean_Rt$R$t_start,
                          Rt = mean_Rt$R$`Mean(R)`,
                          rtlowci = mean_Rt$R$`Quantile.0.05(R)`,
                          rtupci = mean_Rt$R$`Quantile.0.95(R)`)
  
  
  
  #merges the Rt value with the other variant data and renames Rt to have variant suffix
  merge = df %>%
    arrange(days)%>% #keep in week so that the day variable lines up with the first week
    mutate(day = 1:nrow(df)) %>% #used to merge with the estimate_R variable output for the day
    left_join(rt_df) # %>%
  #rename_with(.fn = ~paste0(name,"_",.), .cols = c("Rt", "rtlowci", "rtupci")) #renames the smooth_spline output to have variant prefix
}