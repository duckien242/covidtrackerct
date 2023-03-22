#### Variant Report for SARS-CoV-2 ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#initial commit date: 06/25/22
#original author: Kien Pham
#email k.pham@yale.edu or duckien242@gmail.com
#to be updated Tuesday afternoons

########################################################

# I. PACKAGE INSTALL ####
packages = c("dplyr","tidyr","ggplot2","RColorBrewer","shiny",
             "shinythemes","plotly","DT","data.table",
             "stats","ggpubr","stringr")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

rm(packages, package.check)

library(purrr)
library(readr)
library(DescTools)
library(lubridate)
library(zoo)
library(EpiEstim)
library(ggpubr)
library(httr)
library(jsonlite)

set.seed(1234)


########################################################
##### SET RECENT DATES (CHANGE FOR WEEKLY BUILD) #######

date_3_months = as.Date(lubridate::today() - 92) # CHANGE TO MOST RECENT 3 MONTHS FOR WHO VARIANT PLOT
date_3_weeks = as.Date(lubridate::today() - 22) # CHANGE TO MOST RECENT 3 WEEKS FOR LINEAGE PLOT
date_gisaid = as.Date("2023-03-15") # FOR FILTERING GISAID DOWNLOAD, REFERENCE THIS DATE

# I.1. MULTIPLOT FUNCTION ####
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

################# DATA CHECK LIST ######################

# 1. full metadata file (in data folder, or download from gisaid & change to correct template in data folder, CODE TBA)
# 2. cumulative case data (download from John Hopkins github: https://raw.githubusercontent.com/CSSEGISandData/SARS-CoV-2/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv)
### NOTE: AFTER MAR 9, 2023, JHU STOPPED UPDATING THEIR COVID DATA. We are now using NYTimes' github, which combines JHU data with federal case count. https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv summary table: variant surveillance for Connecticut (template in data or backup folder)
# 4. lab metadata file (for Ct value graph)
# 5. variant name call file: to be updated weekly following pango releases: https://github.com/cov-lineages/pango-designation/

################# STEPS #########################
# 1. update GISAID date in server.R & ui.R
# 2. download latest week's Connecticut metadata from GISAID based on GISAID date, search Location = Connecticut & Collection date = GISAID date, then download Patient metadata
# 3. download lineage designation from pangolin github, manually add latest (Omicron) lineages to variant_names.tsv file, with proper WHO variant designation
# 4. download latest GLab metadata sheet from Google Drive
# 5. change Table_new.csv to Table.csv, delete old Table.csv
# 6. Run App, check for any errors

########################################################
######### START SHINY APP SERVER #########
shinyServer(function(input,output){
  
# II. Read metadata files ########
  
  # we want to read metadata for all Connecticut sequences, format and filter it to add WHO variant designation and limit to sequences in the past 3 months
  
  # read master file
  metadata_CT_old = read.table("data/metadata_CT_all.tsv", sep = "\t", header = TRUE) 
  metadata_CT_old$`Collection.date` <- as.Date(metadata_CT_old$`Collection.date`)
  
  # read recent gisaid submission file.
  new_genomes_title<-list.files('data',pattern = "gisaid",full.names = TRUE,recursive=TRUE,include.dirs=TRUE) 
  new_genomes = read.table(new_genomes_title, sep = "\t", header = TRUE)
  new_genomes$`Collection.date` <- as.Date(new_genomes$`Collection.date`)
  new_genomes = new_genomes %>%
    filter(Collection.date > "2021-01-01")
  
  # add new submissions to master file
  metadata_CT_all = rbind(metadata_CT_old,new_genomes)
  metadata_CT_all = metadata_CT_all %>%
    filter(Lineage != "") %>%
    drop_na(Collection.date)
  
  metadata_CT_all = subset(metadata_CT_all, !duplicated(subset(metadata_CT_all, select=c(`Accession.ID`)))) # remove duplicate sequences
  
  write_tsv(metadata_CT_all,"data/metadata_CT_all.tsv") # write master tsv file for next week
  
  # filter badly named lineages by GISAID
  metadata_CT_all$pango_lineage <- metadata_CT_all$Lineage
  metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "BA.2.75 (marker override based on Emerging Variants AA substitutions)"] = "BA.2.75"
  metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "BF.6 (marker override based on Emerging Variants AA substitutions)"] = "BF.6"
  metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "XBB (marker override based on Emerging Variants AA substitutions)"] = "XBB"
  metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "BA.2" & metadata_CT_all["Collection.date"] > "2022-12-16"] = "XBB.1.5"
  
  # create week column
  #metadata_CT_all <- metadata_CT_all %>% mutate(CopyDate = `Collection.date`)
  metadata_CT_all <- metadata_CT_all %>% group_by(`week` = cut(`Collection.date`, "week")) 
  metadata_CT_all$week <- as.Date(metadata_CT_all$week)
  metadata_CT_all <- metadata_CT_all[!is.na(metadata_CT_all$`week`),]
  
  # filter metadata within last 3 months
  metadata_CT_recent = metadata_CT_all[metadata_CT_all$`Collection.date` >= date_3_months,] 
  metadata_CT_recent =  metadata_CT_recent[!is.na(metadata_CT_recent$`pango_lineage`),]  
  
  # II.1. Create WHO variant name file ######
  
  # read lineage designation csv from pangolin github
  names = read.csv("data/variant_names.csv")
  metadata_CT_who = merge(metadata_CT_all,names[,c("pango_lineage","who_variants")], by = "pango_lineage", all.x = TRUE)
  metadata_CT_who = metadata_CT_who %>%
    mutate(who_variants = replace(who_variants,is.na(metadata_CT_who$who_variants) & week >= "2022-05-05", "Other Omicrons")) # replace NAs to other Omicrons, may change if new variant appears
  
  # create cumulative & most recent 3 week variant files for summary table
  cumulative <- metadata_CT_who %>% 
    select(`Collection.date`,`who_variants`,`pango_lineage`) %>%
    group_by() %>% count(`who_variants`)
  colnames(cumulative) = c("WHO.label","Cumulative.sequenced.cases.")
  cumulative[which(is.na(cumulative$WHO.label)),1] <- "Other"
  
  last_3_weeks <- metadata_CT_who %>% 
    select(`Collection.date`,`who_variants`,`pango_lineage`) %>%
    subset(`Collection.date` >= date_3_weeks) %>% # update date every 3 weeks
    count(`who_variants`) %>%
    mutate(freq = round(100*n/sum(n),2)) 
  colnames(last_3_weeks) = c("WHO.label","Total.sequenced.from.past.3.weeks..","Percent.sequenced.from.past.3.weeks..")
  last_3_weeks[which(is.na(last_3_weeks$WHO.label)),1] <- "Other"
  
  # recategorize minor variants (in Connecticut) into Other
  metadata_CT_who[which(is.na(metadata_CT_who$who_variants)),length(colnames(metadata_CT_who))] <- 'Other'
  metadata_CT_who[which(metadata_CT_who$who_variants == "Eta"),length(colnames(metadata_CT_who))] <- 'Other'
  metadata_CT_who[which(metadata_CT_who$who_variants == "Kappa"),length(colnames(metadata_CT_who))] <- 'Other'
  metadata_CT_who[which(metadata_CT_who$who_variants == "Zeta"),length(colnames(metadata_CT_who))] <- 'Other'
  metadata_CT_who[which(metadata_CT_who$who_variants == "Lambda"),length(colnames(metadata_CT_who))] <- 'Other'

  
# III. Daily lineage freq in last 3 months #######
  
  # we will now filter data to create plot for SARS-CoV-2lineage frequency in Connecticut within the last 3 months
  
  # create dummy column for lineage counting
  metadata_CT_recent$`CopyDate` = as.numeric(metadata_CT_recent$`Collection.date`)
  
  # remove most recent 3 collection dates (inadequate data)
  metadata_CT_recent = metadata_CT_recent[-c(which(metadata_CT_recent$`Collection.date` == max(metadata_CT_recent$`Collection.date`))),] # omit latest 3 days due to reporting bias
  
  # have all lineages be displayed every day, ensuring no breaks when plotting
  lineages_daily_draft <- metadata_CT_recent %>% complete(`Collection.date`, nesting(`pango_lineage`), fill = list(CopyDate = 0))
  lineages_daily_draft$`CopyDate` = ifelse(lineages_daily_draft$`Collection.date` == 0,0,1) # for counting lineages
  
  # find cumulative lineage scores & filter out minor lineages & assign colors
  lineages_sum <- metadata_CT_recent %>% group_by(`pango_lineage`) %>% count(pango_lineage) 
  lineages_sum = lineages_sum %>%
    mutate(ShortLin = pango_lineage) 
  lineages_sum = lineages_sum %>%
    mutate(ShortLin = case_when(
      endsWith(ShortLin, "2.12.1") ~ "BA.2.12.1",
      startsWith(ShortLin, "BG") ~ "BA.2.12.1",
      startsWith(ShortLin, "BA.2.75") ~ "BA.2.75",
      startsWith(ShortLin, "BA.1") ~ "BA.1",
      startsWith(ShortLin, "BA.2") ~ "BA.2",
      startsWith(ShortLin, "BK") ~ "BA.2",
      startsWith(ShortLin, "BP") ~ "BA.2",
      startsWith(ShortLin, "BL") ~ "BA.2.75",
      startsWith(ShortLin, "BH") ~ "BA.2",
      startsWith(ShortLin, "BM") ~ "BA.2",
      startsWith(ShortLin, "CH") ~ "BA.2",
      startsWith(ShortLin, "BJ") ~ "BA.2",
      startsWith(ShortLin, "BR") ~ "BA.2",
      startsWith(ShortLin, "BA.4") ~ "BA.4",
      startsWith(ShortLin, "BA.4.6") ~ "BA.4.6",
      startsWith(ShortLin, "BA.5") ~ "BA.5",
      startsWith(ShortLin, "BU") ~ "BA.5",
      startsWith(ShortLin, "BV") ~ "BA.5",
      startsWith(ShortLin, "BE") ~ "BA.5",
      startsWith(ShortLin, "BT") ~ "BA.5",
      startsWith(ShortLin, "BF") ~ "BA.5",
      startsWith(ShortLin, "BQ.1") ~ "BA.5",
      startsWith(ShortLin, "BQ.2") ~ "BA.5",
      startsWith(ShortLin, "BQ.3") ~ "BA.5",
      startsWith(ShortLin, "CN") ~ "BA.5",
      startsWith(ShortLin, "CD") ~ "BA.5",
      startsWith(ShortLin, "CE") ~ "BA.5",
      startsWith(ShortLin, "CF") ~ "BA.5",
      startsWith(ShortLin, "CL") ~ "BA.5",
      startsWith(ShortLin, "CG") ~ "BA.5",
      startsWith(ShortLin, "CK") ~ "BA.5",
      startsWith(ShortLin, "BZ") ~ "BA.5",
      startsWith(ShortLin, "XBB") ~ "XBB",
      TRUE ~ "Other"
    ))
  
  # add html color code to lineages
  lineages_sum = lineages_sum %>%
    mutate(color =  case_when(
      ShortLin == "BA.1" ~ "#2b931a",
      (ShortLin == "BA.2" & n >= quantile(lineages_sum$n)[5]) ~ "#c5f367", 
      (ShortLin == "BA.2" & n < quantile(lineages_sum$n)[5] & n > quantile(lineages_sum$n)[3]) ~ "#faf994", 
      (ShortLin == "BA.2" & n < quantile(lineages_sum$n)[3] ) ~ "#f8ffda", 
      ShortLin == "BA.2.12.1" ~ "#78f0a1",
      ShortLin == "BA.2.75" ~ "#3eb480",
      ShortLin == "BA.4" ~ "#225fd6", 
      ShortLin == "BA.4.6" ~ "#76e4f5", 
      (ShortLin == "BA.5" & n >= quantile(lineages_sum$n)[5]) ~ "#f3152f", 
      (ShortLin == "BA.5" & n < quantile(lineages_sum$n)[5] & n > quantile(lineages_sum$n)[4]) ~ "#fcaa99",
      (ShortLin == "BA.2" & n < quantile(lineages_sum$n)[4] & n > quantile(lineages_sum$n)[3]) ~ "#faf994", 
      (ShortLin == "BA.5" & n < quantile(lineages_sum$n)[3] ) ~ "#ffe4ee", 
      (ShortLin == "XBB" ~ "#333333"),
      TRUE ~ "#fef9f3"
    )
    )
  row.names(lineages_sum) = lineages_sum$pango_lineage

  
  # reclassify lineages, placing minor lineages as Other
  lineages_daily_draft$`pango_lineage` = ifelse(lineages_daily_draft$`pango_lineage`%in% rownames(lineages_sum),lineages_daily_draft$`pango_lineage`,"Other")
  
  lineages_daily <- lineages_daily_draft %>% group_by(`Collection.date`,`pango_lineage`) %>% 
    summarise_at(vars(`CopyDate`),list(n = sum)) %>%
    mutate(freq = round(100*n/sum(n),2)) 
  
  colnames(lineages_daily)[2] <- "Lineage names"
  
  # insert 3-day rolling average
  lineages_daily = lineages_daily %>%
    group_by(`Lineage names`) %>%
    mutate(freq = as.numeric(freq)) %>%
    mutate(movavg = round(rollmean(freq, 3, na.pad = TRUE),2)) 
  
  # sort by cumulative counts
  lineages_daily_sorted = lineages_daily[order(-lineages_daily$freq),]

  lineages_daily_sorted = lineages_daily_sorted %>%
    mutate(ShortLin = `Lineage names`) 
  lineages_daily_sorted = lineages_daily_sorted %>%
    mutate(ShortLin = case_when(
      endsWith(ShortLin, "2.12.1") ~ "BA.2.12.1",
      startsWith(ShortLin, "BG") ~ "BA.2.12.1",
      startsWith(ShortLin, "BA.2.75") ~ "BA.2.75",
      startsWith(ShortLin, "BA.1") ~ "BA.1",
      startsWith(ShortLin, "BA.2") ~ "BA.2",
      startsWith(ShortLin, "BK") ~ "BA.2",
      startsWith(ShortLin, "BP") ~ "BA.2",
      startsWith(ShortLin, "BL") ~ "BA.2.75",
      startsWith(ShortLin, "BH") ~ "BA.2",
      startsWith(ShortLin, "BM") ~ "BA.2",
      startsWith(ShortLin, "CH") ~ "BA.2",
      startsWith(ShortLin, "BJ") ~ "BA.2",
      startsWith(ShortLin, "BR") ~ "BA.2",
      startsWith(ShortLin, "BA.4") ~ "BA.4",
      startsWith(ShortLin, "BA.4.6") ~ "BA.4.6",
      startsWith(ShortLin, "BA.5") ~ "BA.5",
      startsWith(ShortLin, "BU") ~ "BA.5",
      startsWith(ShortLin, "BV") ~ "BA.5",
      startsWith(ShortLin, "BE") ~ "BA.5",
      startsWith(ShortLin, "BT") ~ "BA.5",
      startsWith(ShortLin, "BF") ~ "BA.5",
      startsWith(ShortLin, "BQ.1") ~ "BA.5",
      startsWith(ShortLin, "BQ.2") ~ "BA.5",
      startsWith(ShortLin, "BQ.3") ~ "BA.5",
      startsWith(ShortLin, "CN") ~ "BA.5",
      startsWith(ShortLin, "CD") ~ "BA.5",
      startsWith(ShortLin, "CE") ~ "BA.5",
      startsWith(ShortLin, "CF") ~ "BA.5",
      startsWith(ShortLin, "CL") ~ "BA.5",
      startsWith(ShortLin, "CG") ~ "BA.5",
      startsWith(ShortLin, "CK") ~ "BA.5",
      startsWith(ShortLin, "BZ") ~ "BA.5",
      startsWith(ShortLin, "XBB") ~ "XBB",
      TRUE ~ "Other"
    ))
    rm(lineages_daily_draft)
    
  # write csv file of lineage frequency in last 3 months
  write.csv(lineages_daily_sorted,"outputs/SARS-CoV-2 lineage frequency in Connecticut in the last 3 months.csv")


# IV. Weekly variant freq since Jan 2021 ###########
  
  # We will create a SARS-CoV-2 WHO variant frequency plot for Connecticut after January 2021
  
  who_weekly <- metadata_CT_who %>% group_by(`week`) %>% count(who_variants) %>%
    mutate(freq = round(n/sum(n),2)) 
  
  # filter out sequences before 2021
  who_weekly = who_weekly[c(which(who_weekly$week > '2020-12-31')),] 
  
  # omit latest week (inadequate data)
  who_weekly = who_weekly[-c(which(who_weekly$week == max(who_weekly$week))),] 
  who_weekly$`Variant names` <- factor(who_weekly$who_variants, 
                                       levels = c("Omicron (BA.5)",
                                                  "Omicron (BA.4)",
                                                  "Omicron (BA.2)",
                                                  "Omicron (BA.1)",
                                                  "Omicron (BA.3)",
                                                  "Omicron (XBB)",
                                                  "Delta",
                                                  "Alpha",
                                                  "Iota",
                                                  "Beta",
                                                  "Gamma",
                                                  "Lambda",
                                                  "Epsilon",
                                                  "Other"))
  # write csv file of WHO variant frequency in Connecticut after Jan 2021
  write.csv(who_weekly,"outputs/SARS-CoV-2 variant frequency in Connecticut.csv")
  
  # set html colors for WHO variants
  who_colors = c("Other" = "#c4c1c0",
                 "Alpha" = "#1160f6",
                 "Delta" = "#11adf6",
                 "Beta" = "#c81ac8",
                 "Gamma" = "#c81ac8",
                 "Lambda" = "#aef516",
                 "Epsilon" = "#13d84f",
                 "Iota" = "#05d0ab",
                 "Mu" = "#9e5220",
                 "Omicron (BA.1)" = "#2b931a",
                 "Omicron (BA.2)" = "#fdb902",
                 "Omicron (BA.3)" = "#e58b00",
                 "Omicron (BA.4)" = "#f47110",
                 "Omicron (BA.5)" = "#a90505",
                 "Omicron (XBB)" = "#333333"
  )

  
# V. Summary table ########
  
  # Table of the composition of SARS-CoV-2 sequence submissions (to GISAID) in Connecticut
  
  # read last week's summary table
  table1_old = read.csv("data/summary_table.csv")
  table1_old = table1_old %>%
    select(-`X`)
  
  # merge old summary table with cumulative & recent variant data
  tab = merge(table1_old,cumulative[,c("WHO.label","Cumulative.sequenced.cases.")], by = "WHO.label", all.x = TRUE)
  table1_new = merge(tab,last_3_weeks[,c("WHO.label","Total.sequenced.from.past.3.weeks..","Percent.sequenced.from.past.3.weeks..")], by = "WHO.label", all.x = TRUE)
  table1_new[is.na(table1_new)] <- 0 # assign 0 to NA values
  table1_new_draft <- table1_new 
  
  # y is new report, x is old
  # find changes since last week's report
  table1_new$`Percent.change.from.previous.report.y` = round((table1_new_draft$`Percent.sequenced.from.past.3.weeks...y` - table1_new_draft$`Percent.sequenced.from.past.3.weeks...x`),2)
  table1_new = table1_new[,-c(4:7)] # remove old entries
  table1_new[which(table1_new$`WHO label` == 0),1] <- "Other"
  colnames(table1_new) = c("WHO label",
                           "Pango lineage",
                           "CDC classification",
                           "Cumulative sequenced cases*",
                           "Total sequenced from past 3 weeks**",
                           "Percent sequenced from past 3 weeks**",
                           "Percent change from previous report")
  # write csv for summary table
  write.csv(table1_new,"data/summary_table_new.csv")
  write.csv(table1_new,"outputs/Summary table of SARS-CoV-2 variant frequency.csv")
  
  recent_variants_count = length(which(table1_new$`Total sequenced from past 3 weeks**` != 0)) # find non-zero variants in past 3 weeks
  
# VI. Sequence proportion ######
  
  # Create plot tracking SARS-CoV-2 weekly sequence proportion compared to case count in Connecticut 
  
  # read cumulative case file
  cases = fread("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv", sep = ",")
  
  # filter for Connecticut & remove columns
  cases_CT_t = cases %>%
    filter(Province_State == "Connecticut") %>%
    select(-c(UID,iso2,iso3,code3,FIPS,Admin2,Country_Region,Lat,Long_,Combined_Key,Province_State)) 
  cases_CT = tibble(apply(cases_CT_t,2,sum),mdy(colnames(cases_CT_t)))
  colnames(cases_CT) = c("case_cumulative","date")
  
  
  # read cumulative case file
  cases = fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv", sep = ",")
  cases_CT = cases %>%
    filter(state == "Connecticut") 
  
  # create week column & daily case count column
  cases_CT_filt = cases_CT %>% 
    mutate(case_count = cases - lag(cases)) %>%
    group_by(`week` = cut(`date`, "week")) 
  cases_CT_filt = data.table(cases_CT_filt)

  # sum rows by week to calculate weekly case count
  cases_CT_filt_count = cases_CT_filt[,list(case_count=sum(case_count)), by='week']
  # find cumulative case count by week
  cases_CT_filt_cumu = cases_CT_filt[,list(cases=max(cases)), by='week']
  # merge weekly count & cumulative count
  cases_CT_filt.dt = merge(cases_CT_filt_count,cases_CT_filt_cumu, by = 'week')
  # filter out counts before 2021
  cases_CT_filt.dt = cases_CT_filt.dt[-c(1:which(cases_CT_filt.dt$week == "2020-12-28")),]
  cases_CT_filt.dt$week = as.Date(cases_CT_filt.dt$week)
  
  # find sequence counts
  seq_CT = metadata_CT_all %>%
    group_by(`week`,`Collection.date`) %>%
    count()
  seq_CT_filt.dt = data.table(seq_CT)
  # sum rows by week to calculate weekly sequence counts
  seq_CT_filt.dt = seq_CT_filt.dt[,list(seq_count=sum(n)), by='week']
  # find cumulative sequence count
  seq_CT_filt.dt = seq_CT_filt.dt %>%
    mutate(seq_cumulative = cumsum(`seq_count`))
  # filter out counts before 2021
  seq_CT_filt.dt = seq_CT_filt.dt[-c(1:which(seq_CT_filt.dt$week == "2020-12-28")),]
  
  # merge case & sequence count files
  subsampler = merge(cases_CT_filt.dt,seq_CT_filt.dt, by = 'week',all.x = TRUE)
  subsampler[is.na(subsampler)] <- 0
  # calculate % sequenced & remove sequence counts
  subsampler = subsampler %>%
    mutate(percent = round(seq_count * 100/ case_count,2)) 
  
  rm(cases,cases_CT,cases_CT_filt,seq_CT,cases_CT_filt_cumu,cases_CT_filt_count,cases_CT_t, cases_CT_filt.dt, metadata_CT_old)
  gc()
  
  # find log count
  subsampler = subsampler %>%
    mutate(log_cases = log(cases)) %>%
    mutate(log_seq = log(seq_cumulative))
  
  # write sequence proportion csv file
  write.csv(subsampler,"outputs/SARS-CoV-2 case count and sequenced proportion in Connecticut.csv")
  
  # manually create y-axis for sequence proportion plot (we will need 2 y-axis)
  ay2 <- list(
    tickfont = list(size=11, color = "orange"),
    titlefont=list(size=14.6, color = "orange"),
    min = 0,
    max = 100,
    overlaying = "y",
    nticks = 5,
    side = "right",
    title = "proportion of cases sequenced"
  )
  
  ay1 <- list(
    tickfont = list(size=11, color = "blue"),
    titlefont=list(size=14.6, color = "blue"),
    min = 0,
    max = 100,
    nticks = 5,
    side = "left",
    title = "case count"
  )
  
  
# VII. Ct Values by Lineages #############
  
  # Create plot to compare PCR cycle threshold (Ct value) of SARS-CoV-2 samples that were sequenced and submitted by Grubaugh Lab
  
  glab = read.csv("data/GLab_SC2_sequencing_data - Sample metadata.csv")
  gc()

  glab_filt <- glab %>% mutate(Ctvals = Lineage) 
  
  #include only samples with lineage assignments
  glab_filt = glab_filt %>%
    mutate(Ctvals = case_when(
      endsWith(Ctvals, "2.12.1") ~ "BA.2.12.1",
      startsWith(Ctvals, "BG") ~ "BA.2.12.1",
      startsWith(Ctvals, "BA.2.75") ~ "BA.2.75",
      startsWith(Ctvals, "BA.1") ~ "BA.1",
      startsWith(Ctvals, "BA.2") ~ "BA.2",
      startsWith(Ctvals, "BK") ~ "BA.2",
      startsWith(Ctvals, "BP") ~ "BA.2",
      startsWith(Ctvals, "BL") ~ "BA.2.75",
      startsWith(Ctvals, "BH") ~ "BA.2",
      startsWith(Ctvals, "BM") ~ "BA.2",
      startsWith(Ctvals, "CH") ~ "BA.2",
      startsWith(Ctvals, "BJ") ~ "BA.2",
      startsWith(Ctvals, "BR") ~ "BA.2",
      startsWith(Ctvals, "BA.4") ~ "BA.4",
      startsWith(Ctvals, "BA.4.6") ~ "BA.4.6",
      startsWith(Ctvals, "BA.5") ~ "BA.5",
      startsWith(Ctvals, "BU") ~ "BA.5",
      startsWith(Ctvals, "BV") ~ "BA.5",
      startsWith(Ctvals, "BE") ~ "BA.5",
      startsWith(Ctvals, "BT") ~ "BA.5",
      startsWith(Ctvals, "BF") ~ "BA.5",
      startsWith(Ctvals, "BQ.1") ~ "BA.5",
      startsWith(Ctvals, "BQ.2") ~ "BA.5",
      startsWith(Ctvals, "BQ.3") ~ "BA.5",
      startsWith(Ctvals, "CN") ~ "BA.5",
      startsWith(Ctvals, "CD") ~ "BA.5",
      startsWith(Ctvals, "CE") ~ "BA.5",
      startsWith(Ctvals, "CF") ~ "BA.5",
      startsWith(Ctvals, "CL") ~ "BA.5",
      startsWith(Ctvals, "CG") ~ "BA.5",
      startsWith(Ctvals, "CK") ~ "BA.5",
      startsWith(Ctvals, "BZ") ~ "BA.5",
      startsWith(Ctvals, "XBB") ~ "XBB",
      TRUE ~ "Other"
    ))
  Ct_values <- glab_filt %>% dplyr::filter(`Ctvals` %in% c("BA.2","BA.5","XBB"))
  rm(glab,glab_filt)
  gc()
  
  
# VIII. Rt Values by Variant #############
  
  # Calculate and plot effective reproduction number (Rt) for SARS-CoV-2 variant in Connecticut, based on serial interval estimation
  
  # Making API request using the GET() function and specifying the API’s URL:
  url = paste("api2.covidestim.org/runs?geo_name=eq.Connecticut&run_date=gt.",lubridate::today() - 21,"&select=*%2Ctimeseries(*)",sep = "")
  res = GET(url)
  
  # convert the raw Unicode into a character vector that resembles the JSON format to obtain case data
  data = fromJSON(rawToChar(res$content))
  infect_import = data[[8]][[1]]
  rm(data,res,url)
  
  # filter columns for use
  infect_import$run_id = "Connecticut"
  colnames(infect_import)[1] <- "state"
  infect = infect_import %>%
    select(state,
           date,
           infections,
           infections_p2_5,
           infections_p97_5) %>%
    #*****************************************
    rename(week = date) %>% #rename for merging with var_data
    mutate(week = as.Date(week))
  
  # variant frequency data (who_weekly obtained in shinyapp pipeline)
  who_weekly_rt = who_weekly %>%
    select(-c(`Variant names`)) %>%
    mutate(week = week +3) %>% # create week to match 
    filter(!(who_variants == "Omicron (BA.2)" & week >= "2022-08-15")) # remove one-off BA.2 infections
    
  who_weekly_rt = data.frame(who_weekly_rt)
  var_data = reshape(who_weekly_rt, idvar = "week", timevar = "who_variants", direction = "wide")
  var_data[is.na(var_data)] <- 0
  var_data$state = "Connecticut"
  
  #*******************************************************************************
  #####MERGE DATA#####
  #*#*******************************************************************************
  var_merge = var_data %>%
    left_join(infect, by = c("week","state")) %>% #merges covidestim data with our data and keeps only 
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

  
  #creates list of dataframes by variant for calculating estimate_R
  var_list = var_merge3_main %>%
    group_by(variant) %>%
    arrange(variant) %>%
    group_split()
  
  
  var_name = unique(var_merge3_main$variant)
  
  rm(var_merge2_lowci,var_merge2_main,var_merge2_upci)
  gc()
  
  #*******************************************************************************
  #RT FUNCTION ####
  #*#******************************************************************************
  
  #Rt Calculation
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
  
  
  #run estimate_R function on 
  rt = list()
  
  for(i in seq_along(var_list)){
    rt[[i]] = rt_fun(var_list[[i]])
  }
  
  rt_comb = bind_rows(rt)
  
  rt_comb = rt_comb %>%
    mutate_all(~gsub("_infections", "", .)) 
  
  rt_comb = rt_comb %>%
    mutate(days = as.Date(days)) %>%
    left_join(var_merge3_lowci, by = c("days","variant")) %>% 
    left_join(var_merge3_upci, by = c("days","variant")) 
  
  rt_comb$variant <- factor(rt_comb$variant, 
                            levels = c("Omicron (BA.5)",
                                       "Omicron (XBB)",
                                       "Omicron (BA.4)",
                                       "Omicron (BA.2)",
                                       "Omicron (BA.1)",
                                       "Omicron (BA.3)",
                                       "Delta",
                                       "Other"))
  write.csv(rt_comb, "outputs/SARS-CoV-2 effective reproduction number in Connecticut since November 2021.csv")

######## IX. TOTAL SEQUENCES ########
  # estimate total sequences by Grubaugh Lab, may differ from submitted sequences on GISAID database
  
  old_total = read.csv("data/total_sequences.csv")
  old_total = old_total %>%
    select(-c(X))
   # mutate(date = mdy(date))
  
  total_sequences = old_total$total[nrow(old_total)] + nrow(new_genomes[grep("Yale", new_genomes$Virus.name),]) # recalculate total
  
  today = as.character(today())
  updated_total = c(today,total_sequences)
  
  new_total = old_total %>%
    rbind(updated_total) %>%
    mutate(date = as.Date(date)) 
  new_total = new_total %>%
    subset(!duplicated(subset(new_total, select=c(date)))) # delete duplicated rows by date (in case of multiple test runs)
  
  write.csv(new_total,"data/total_sequences.csv")
  
  new_total_fixed = new_total
  colnames(new_total_fixed) = c("Date of website update","Total number of sequences")
  write.csv(new_total_fixed,"outputs/Total number of SARS-CoV-2 sequences by Grubaugh Lab")
  
######## X. PLOTTING ##############
  
  # 1. Daily lineage freq in last 3 months
  
  output$lin_freq <- renderPlotly({
    
    if (input$button == "3-day rolling average") {
      lineage_rollavg = ggplot(lineages_daily_sorted,aes(x=`Collection.date`,y=`movavg`, group=`Lineage names`, color = `Lineage names`)) +
        geom_line() +
        labs(
          title = "Lineage frequencies in Connecticut (last 3 months): 3-day rolling average",
          y = "Lineage frequency (%)",
          x = "Time (days)",
          subtitle = "% for lattest week shown subject to change"
        )+
        theme_light() + 
        scale_color_manual(values=lineages_sum$color) +
        theme(
          plot.title = element_text(size = 15) ,
          axis.text.x = element_text(hjust = 1, size = 10, color = "black")
        ) 
      
      ggplotly(lineage_rollavg)
      
    } else {
      
      # daily count plot
      lineage_plot = ggplot(lineages_daily_sorted,aes(x=`Collection.date`,y=freq, group=`Lineage names`, color = `Lineage names`)) +
        geom_line() +
        labs(
          title = "Daily lineage frequencies in Connecticut (last 3 months)",
          y = "Lineage frequency (%)",
          x = "Time (days)",
          subtitle = "% for lattest week shown subject to change"
        )+
        theme_light() + 
        scale_color_manual(values=lineages_sum$color) +
        theme(
          plot.title = element_text(size = 15) ,
          axis.text.x = element_text(hjust = 1, size = 10, color = "black")
        ) 
      
      ggplotly(lineage_plot)
    }
    
  })
  
  # 2. Weekly variant freq since 1-2021
  
  output$var_freq <- renderPlotly({

      if (input$button_who == "Line chart display") {
        who_plot = ggplot(who_weekly,aes(x=`week`,y=freq, group = `Variant names`, color = `Variant names`)) +
          geom_line() +
          labs(
            title = "Weekly variant frequencies in Connecticut (since January 2021)",
            y = "Variant frequency (%)",
            x = "Time (weeks)"
          ) + 
          scale_color_manual(values = c(who_colors)) +
          scale_x_date(date_breaks="3 month", date_labels="%m-%Y") +
          theme_light() +
          theme(
            plot.title = element_text(size = 15,hjust = 0.5)
          )
        ggplotly(who_plot)
      } else {
        
        who_plot_bar = ggplot(who_weekly,aes(x=`week`,y=freq, group = `Variant names`, fill = `Variant names`)) +
          geom_bar(position = "fill", stat="identity") +
          labs(
            title = "Weekly variant frequencies in Connecticut (since January 2021)",
            y = "Variant frequency"
          ) + 
          scale_fill_manual(values = c(who_colors)) +
          scale_x_date(date_breaks="3 month", date_labels="%m-%Y") +
          theme_light() +
          theme(
            plot.title = element_text(size = 15,hjust = 0.5)
          )
        ggplotly(who_plot_bar)
      }
    
  })
  
  # 2.a. Rt values by variant\
  output$rt <- renderPlot({
    
    who_colors_rt = c(
                   "Delta" = "#11adf6",
                   "Omicron (BA.1)" = "#2b931a",
                   "Omicron (BA.2)" = "#fdb902",
                   "Omicron (BA.4)" = "#f47110", 
                   "Omicron (XBB)" = "#333333", 
                   "Omicron (BA.5)" = "#a90505"
    )

    if (input$button_rt == "Effective reproduction number") {
      
      rt_plot = ggplot(rt_comb, aes(x = `days`, y = as.numeric(Rt), group = variant, color = variant)) +
        geom_line() + 
        geom_ribbon(aes(ymin = as.numeric(rtlowci), 
                        ymax = as.numeric(rtupci), fill = variant), 
                    alpha=0.1, 
                    linetype="blank",
                    color="grey") +
        geom_hline(yintercept = 1) +
        labs(
          title = "Variant effective reproduction number (Rt) in Connecticut (since Nov 2021)",
          y = "Rt value",
          x = "Time (weeks)"
        ) + 
        scale_fill_manual(values = c(who_colors_rt)) +
        scale_color_manual(values = c(who_colors_rt)) +
        scale_x_date(date_breaks="3 month", date_labels="%m-%Y") +
        theme_light()+
        theme(
          plot.title = element_text(size = 24,hjust = 0.5),
          axis.title.x = element_text(size = 17, hjust = 0.5),
          axis.text.x = element_text(size = 14, hjust = 0.5),
          axis.text.y = element_text(size = 14, hjust = 0.5),
          legend.text = element_text(size = 14, hjust = 0.5),
          legend.title = element_text(size = 18, hjust = 0.5),
          axis.title.y = element_text(size = 17, hjust = 0.5)
        )
      rt_plot
      
    } else {
      rt_comb_filt = rt_comb %>%
        filter(days < max(days)-2)
      I_plot = ggplot(rt_comb_filt, aes(x = `days`, y = as.numeric(I), group = variant, color = variant)) +
        geom_line() + 
        labs(
          title = "Variant estimated infections in Connecticut (since Nov 2021)",
          y = "Estimated infections (cases)",
          x = "Time (weeks)"
        ) + 
        geom_ribbon(aes(ymin = as.numeric(I_025), 
                        ymax = as.numeric(I_975), fill = variant), 
                    alpha=0.1, 
                    linetype="blank",
                    color="grey") +
        scale_fill_manual(values = c(who_colors_rt)) +
        scale_color_manual(values = c(who_colors_rt)) +
        scale_x_date(date_breaks="3 month", date_labels="%m-%Y") +
        theme_light() +
        theme(
          plot.title = element_text(size = 22,hjust = 0.5)
        )
      I_plot
    }
    
    
  })
  
  
  # 3. summary table
  
  output$table1 <- DT::renderDataTable(
    datatable(table1_new,
              filter = "top",
              caption = htmltools::tags$caption(
                style = 'caption-side: top; text-align: center;font-size:20px',

                'Summary table: SARS-CoV-2 variants in Connecticut.'),
              options = list(pageLength = recent_variants_count, 
                             dom = 'Bfrtip',
                             selection="multiple",
                             order = list(list(5,'desc'),list(4,'desc')),
                             initComplete = JS(
                               "function(settings, json) {",
                               "$(this.api().table().header()).css({'background-color': '#504d49', 'color': '#fff'});",
                               "}"))
              )
  )
  
  
  # 4. Ct values by lineage
  
  output$CT <- renderPlot({

    intercept_char <- Ct_values %>% dplyr::filter(`Ctvals` == c("BA.2"))
    intercept <- as.numeric(intercept_char$`Yale.N1.FAM.`)
    intercept <- mean(intercept,na.rm = TRUE)
    rm(intercept_char)
    
    Ct_graph <- ggplot(data = Ct_values, aes(x= `Ctvals`, y = as.numeric(`Yale.N1.FAM.`)))+
      geom_hline(yintercept = intercept,colour = "lightgrey")+
      geom_boxplot(width=0.4,outlier.shape = NA,colour = "#666666")+
      geom_jitter(aes(col =`Ctvals`), alpha = 0.2,size = 0.2, stroke = 2,shape = 21, width = 0.15)+
      labs(title="PCR cycle threshold (CT) values by lineage",x=NULL, y = "Ct (N)")+
      scale_color_brewer(palette="Dark2")+
      coord_fixed(ratio = 0.06)+
      scale_y_reverse(breaks = seq(10, 40, by = 5))+
      theme_light()+
      theme(
        plot.title=element_text(size=20),
        legend.position="none",axis.text = element_text(size = 10, color = "black",face ="bold"), 
                            axis.title=element_text(size=12, face ="bold"))
    Ct_graph
    
  })
  

  
  # 5. Subsampler
  
  subsampler_short = subsampler %>% # last 10 months for better viewing of recent trends
    filter(week > (as.Date(lubridate::today() - 92) - 300))

  subsampler_plot = plot_ly(subsampler) %>%
    add_bars(x = ~week, y = ~case_count, colors = "blue", 
             name = "weekly cases",data = subsampler, showlegend=TRUE, inherit=FALSE) %>%
    add_lines(x = ~week, y = ~percent, colors ="red", yaxis = "y2", 
              name = "%seq",data = subsampler, showlegend=TRUE, inherit=FALSE) %>%
    layout(title = "Proportion of SARS-CoV-2 cases sequenced in Connecticut",
           yaxis = ay1,
           yaxis2 = ay2)
  
  subsampler_log = plot_ly(subsampler) %>%
    add_trace(x = ~week, y = ~log_cases, colors = "blue", 
              name = "log cumulative cases",data = subsampler, showlegend=TRUE, inherit=FALSE, mode = 'lines+markers') %>%
    add_trace(x = ~week, y = ~log_seq, colors ="red", 
              name = "log cumulative sequences",data = subsampler, showlegend=TRUE, inherit=FALSE, mode = 'lines+markers') %>%
    layout(title = "Log cumulative SARS-CoV-2 cases and sequences in Connecticut",
           yaxis = ay1)
  
  subsampler_plot_short = plot_ly(subsampler_short) %>%
    add_bars(x = ~week, y = ~case_count, colors = "blue", 
             name = "weekly cases",data = subsampler_short, showlegend=TRUE, inherit=FALSE) %>%
    add_lines(x = ~week, y = ~percent, colors ="red", yaxis = "y2", 
              name = "%seq",data = subsampler_short, showlegend=TRUE, inherit=FALSE) %>%
    layout(title = "Proportion of SARS-CoV-2 cases sequenced in Connecticut",
           yaxis = ay1,
           yaxis2 = ay2)
  
  subsampler_log_short = plot_ly(subsampler_short) %>%
    add_trace(x = ~week, y = ~log_cases, colors = "blue", 
              name = "log cumulative cases",data = subsampler_short, showlegend=TRUE, inherit=FALSE, mode = 'lines+markers') %>%
    add_trace(x = ~week, y = ~log_seq, colors ="red", 
              name = "log cumulative sequences",data = subsampler_short, showlegend=TRUE, inherit=FALSE, mode = 'lines+markers') %>%
    layout(title = "Log cumulative SARS-CoV-2 cases and sequences in Connecticut",
           yaxis = ay1)
  
  output$subsampler <- renderPlotly({
    
    if (input$button_shorten == TRUE & 
        input$button_subsampler == "Proportion of cases sequenced") {
          ggplotly(subsampler_plot)
      } else {
        if (input$button_shorten == FALSE & 
            input$button_subsampler == "Proportion of cases sequenced") {
        ggplotly(subsampler_plot_short)
      } else {
        if (input$button_shorten == TRUE & 
            input$button_subsampler == "Log cumulative cases-sequences") {
          ggplotly(subsampler_log)
      } else {
        ggplotly(subsampler_log_short)
      }
      }
      }
  
  })
  
  
  # 6. Total sequences 
  
  output$total_seq <- renderValueBox({
    
    total_sequences = new_total$total[nrow(new_total)] 

    valueBox(
            value = HTML("<p style='font-size:40px', text-align: center>",
                         total_sequences,"</p>"),
            subtitle = HTML("<p text-align: center>",
                            "Yale SARS-CoV-2 Genomic Surveillance Initiative: <br>Genomes Published<br>"),
            color = "orange",
            width = 4
    )
  })
  
  
  # 7. Nextstrain 
  
  observe({ 
    
    test <<- paste0("https://nextstrain.org/groups/grubaughlab-public/CT-SARS-CoV-2/connecticut?c=clade_membership") 
  })
  
  output$nextstrain <- renderUI({
    input$Member
    my_test <- tags$iframe(src=test, height = 800, width = 700)
    print(my_test)
    my_test
  })
  
  output$nextstrain_caption <- renderInfoBox({
    
    nextstrain_txt = "We created a custom Nextstrain page to visualize the relatedness of sequenced SARS-CoV-2 cases in the state. Below shows a phylogenetic tree of historic sequences in Connecticut, colored by pango lineage assignments."
    
    infoBox(title = HTML("Phylogenetics:"),
            value = HTML("<p style='font-size:18px; align: center'>",
                         nextstrain_txt,"</p>"),
            icon = icon("creative-commons-sampling"),
            fill = TRUE,
            width = 4
    )
  })
  
  # 8. Twitter 

  output$twitter_caption <- renderInfoBox({
    
    twitter_txt = ""
    
    infoBox(title = HTML("Twitter:"),
            value = HTML("<p style='font-size:18px'>",
                         twitter_txt,"</p>"),
            color = "orange",
            icon = icon("twitter"),
            fill = TRUE,
            width = 4
    )
  })
  
  # 9. About
  output$about_tracker <- renderInfoBox({
    
    infoBox(title = HTML("CovidTrackerCT:"),
            value = HTML("<p style='font-size:18px'>",
                        "CovidTrackerCT was built and is maintained by members of the <a href = 'https://grubaughlab.com/'>Grubaugh Lab </a> at the Yale School of Public Health, New Haven, CT. 
                        We started this website in order to share our weekly reports on SARS-CoV-2 genomic epidemiology and communicate the results of our research in a clear, up-to-date way to the people for whom they matter most; residents of Connecticut. "),
            color = "orange",
            fill = TRUE,
            width = 4
    )
  })
  
  output$about_grubaughlab <- renderInfoBox({
    
    infoBox(title = HTML("Grubaugh Lab:"),
            value = HTML("<p style='font-size:18px'>",
                         "The Grubaugh lab at the Yale School of Public Health does research in viral sequencing, evolution, and transmission. As it became clear that SARS-CoV-2 is a threat, we focused all of our efforts on helping in any way we could with the response, and are now applying our experience and knowledge to understand its spread in Connecticut."),
            color = "orange",
            fill = TRUE,
            width = 4
    )
  })
  
  output$about_impact <- renderInfoBox({
    
    infoBox(title = HTML("Impact Yale:"),
            value = HTML("<p style='font-size:18px'>",
                         "Our effort is part of a larger, multidisciplinary consortium of Yale laboratories called 
                         IMPACT (<strong>I</strong>mplementing <strong>M</strong>edical and <strong>P</strong>ublic health <strong>A</strong>ction against Coronavirus in <strong>CT</strong>). 
                         Together, these laboratories conduct a number of studies, including the Yale SARS-CoV-2 Genomic Surveillance Initiative. Samples being sequenced in this project are compiled in a biorepository curated by IMPACT. Many of the protocols we use benefitted from the input of many researchers at Yale and around the world, for which we are grateful."),
            color = "orange",
            fill = TRUE,
            width = 4
    )
  })
  
  # 10. Download buttons
  output$downloadLineages <- downloadHandler(
    filename = function() {
      paste("SARS-CoV-2 lineage frequency in Connecticut since ",date_3_months,".csv",sep ='')
    },
    content = function(con) {
      write.csv(lineages_daily_sorted,con)
    }
   )
  
  output$downloadVariants <- downloadHandler(
    filename = function() {
      paste("SARS-CoV-2 variant frequency in Connecticut.csv",sep ='')
    },
    content = function(con) {
      write.csv(who_weekly,con)
    }
  )
  
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("Summary table of SARS-CoV-2 variant frequency.csv",sep ='')
    },
    content = function(con) {
      write.csv(table1_new,con)
    }
  )
  
  output$downloadRt <- downloadHandler(
    filename = function() {
      paste("SARS-CoV-2 effective reproduction number in Connecticut since November 2021.csv",sep ='')
    },
    content = function(con) {
      write.csv(rt_comb,con)
    }
  )
  
  output$downloadSubsampler <- downloadHandler(
    filename = function() {
      paste("SARS-CoV-2 case count and sequenced proportion in Connecticut.csv",sep ='')
    },
    content = function(con) {
      write.csv(subsampler,con)
    }
  )


})
