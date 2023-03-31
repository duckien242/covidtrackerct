# function for preparing for sequence proportion per case plot
# takes SARS-CoV-2 case data file & sequence file and generate a dataframe file for plotting

perc_seq <- function(case_file, seq_file) {
  
  # create week column & daily case count column
  cases_file_filt = case_file %>% 
    mutate(case_count = cases - lag(cases)) %>%
    group_by(`week` = cut(`date`, "week")) 
  cases_file_filt = data.table(cases_file_filt)
  
  # sum rows by week to calculate weekly case count
  cases_file_filt_count = cases_file_filt[,list(case_count=sum(case_count)), by='week']
  # find cumulative case count by week
  cases_file_filt_cumu = cases_file_filt[,list(cases=max(cases)), by='week']
  # merge weekly count & cumulative count
  cases_file_filt.dt = merge(cases_file_filt_count,cases_file_filt_cumu, by = 'week')
  # filter out counts before 2021
  cases_file_filt.dt = cases_file_filt.dt[-c(1:which(cases_file_filt.dt$week == "2020-12-28")),]
  cases_file_filt.dt$week = as.Date(cases_file_filt.dt$week)
  
  # process sequence count file
  seq_file_filt.dt = data.table(seq_file)
  # sum rows by week to calculate weekly sequence counts
  seq_file_filt.dt = seq_file_filt.dt[,list(seq_count=sum(n)), by='week']
  # find cumulative sequence count
  seq_file_filt.dt = seq_file_filt.dt %>%
    mutate(seq_cumulative = cumsum(`seq_count`))
  # filter out counts before 2021
  seq_file_filt.dt = seq_file_filt.dt[-c(1:which(seq_file_filt.dt$week == "2020-12-28")),]
  
  # merge case & sequence count files
  seq_prop = merge(cases_file_filt.dt,seq_file_filt.dt, by = 'week',all.x = TRUE)
  seq_prop[is.na(seq_prop)] <- 0
  # calculate % sequenced & remove sequence counts
  seq_prop = seq_prop %>%
    mutate(percent = round(seq_count * 100/ case_count,2)) 
  # find log count
  seq_prop = seq_prop %>%
    mutate(log_cases = log(cases)) %>%
    mutate(log_seq = log(seq_cumulative))
  
}