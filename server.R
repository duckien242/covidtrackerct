#### Variant Report for SARS-CoV-2 ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#initial commit date: 06/25/22
#original author: Kien Pham
#email k.pham@yale.edu or duckien242@gmail.com
#to be updated Tuesday afternoons


########################################################
##### SET RECENT DATES (CHANGE FOR WEEKLY BUILD) #######

date_3_months = as.Date("2022-04-13") # CHANGE TO MOST RECENT 3 MONTHS
date_3_weeks = as.Date("2022-06-22") # CHANGE TO MOST RECENT 3 WEEKS
date_gisaid = as.Date("2022-06-24") # FOR FILTERING GISAID DOWNLOAD
total_sequences = 17292 # SEARCH "YALE" ON GISAID

########################################################

# I. PACKAGE INSTALL ####
packages = c("dplyr","tidyr","ggplot2","RColorBrewer","shiny","shinythemes","plotly","DT","data.table")
library(readr)
library(lubridate)
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

set.seed(1234)

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
# 2. cumulative case data (download from John Hopkins github: https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv)
# 3. Table 1: variant surveillance (template in data or backup folder)
# 4. lab metadata file (for Ct value graph)
# 5. variant name call file: to be updated weekly following pango releases: https://github.com/cov-lineages/pango-designation/

########################################################
######### START SHINY APP SERVER #########
shinyServer(function(input,output){
  
# II. Read metadata files ########
  
  # read master file
  metadata_CT_old = read.table("data/metadata_CT_all.tsv", sep = "\t", header = TRUE) 
  metadata_CT_old$`Collection.date` <- as.Date(metadata_CT_old$`Collection.date`)
  # read recent gisaid submission file.
  new_genomes_title<-list.files('data',pattern = "gisaid",full.names = TRUE,recursive=TRUE,include.dirs=TRUE) 
  new_genomes = read.table(new_genomes_title, sep = "\t", header = TRUE)
  new_genomes$`Collection.date` <- as.Date(new_genomes$`Collection.date`)
  # concat new submissions to masterfile
  metadata_CT_all = rbind(metadata_CT_old,new_genomes)
  metadata_CT_all = subset(metadata_CT_all, !duplicated(subset(metadata_CT_all, select=c(`Accession.ID`)))) # remove duplicate sequences
  write_tsv(metadata_CT_all,"data/metadata_CT_all.tsv") # write tsv file for next week
  
  metadata_CT_all$pango_lineage <- metadata_CT_all$Lineage
  metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "BA.2.12.1 (marker override: BA.2.12 + Spike_L452Q => BA.2.12.1)"] = "BA.2.12.1"
  
  # create week date column
  metadata_CT_all <- metadata_CT_all %>% mutate(CopyDate = `Collection.date`)
  metadata_CT_all <- metadata_CT_all %>% group_by(`week` = cut(`CopyDate`, "week")) 
  metadata_CT_all$week <- as.Date(metadata_CT_all$week)
  metadata_CT_all <- metadata_CT_all[!is.na(metadata_CT_all$`week`),]
  # filter by last 3 months
  metadata_CT_recent = metadata_CT_all[metadata_CT_all$`Collection.date` >= date_3_months,] 
  metadata_CT_recent =  metadata_CT_recent[!is.na(metadata_CT_recent$`pango_lineage`),]  
  
  # II.1. Create WHO variant name file ######
  names = read.csv("data/variant_names.csv")
  metadata_CT_who = merge(metadata_CT_all,names[,c("pango_lineage","who_variants")], by = "pango_lineage", all.x = TRUE)
  
  # create cumulative & most recent 3 week variant files for Table 1
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
  
  # recategorize minor variants into Other
  metadata_CT_who[which(is.na(metadata_CT_who$who_variants)),length(colnames(metadata_CT_who))] <- 'Other'
  metadata_CT_who[which(metadata_CT_who$who_variants == "Eta"),length(colnames(metadata_CT_who))] <- 'Other'
  metadata_CT_who[which(metadata_CT_who$who_variants == "Kappa"),length(colnames(metadata_CT_who))] <- 'Other'
  metadata_CT_who[which(metadata_CT_who$who_variants == "Zeta"),length(colnames(metadata_CT_who))] <- 'Other'
  metadata_CT_who[which(metadata_CT_who$who_variants == "Lambda"),length(colnames(metadata_CT_who))] <- 'Other'

  
# III. Daily lineage freq in last 3 months #######

  metadata_CT_recent$`CopyDate` = as.character(metadata_CT_recent$`CopyDate`)
  # remove most recent 3 collection dates (inadequate data)
  metadata_CT_recent = metadata_CT_recent[-c(which(metadata_CT_recent$`Collection.date` == max(metadata_CT_recent$`Collection.date`))),] # omit latest 3 days
  metadata_CT_recent = metadata_CT_recent[-c(which(metadata_CT_recent$`Collection.date` == max(metadata_CT_recent$`Collection.date`))),] # omit latest 3 days
  metadata_CT_recent = metadata_CT_recent[-c(which(metadata_CT_recent$`Collection.date` == max(metadata_CT_recent$`Collection.date`))),] # omit latest 3 days
  
  # have all lineages be displayed every day, ensuring no breaks when plotting
  lineages_daily_draft <- metadata_CT_recent %>% complete(`Collection.date`, nesting(`pango_lineage`), fill = list(CopyDate = 0))
  lineages_daily_draft$`CopyDate` = ifelse(lineages_daily_draft$`CopyDate` == 0,0,1) # for counting lineages
  
  # find cumulative lineage scores & filter out minor lineages
  lineages_sum <- metadata_CT_recent %>% group_by(`pango_lineage`) %>% count(pango_lineage)  %>%
    filter(n > 5)
  rownames(lineages_sum) = lineages_sum$pango_lineage
  
  # reclassify lineages, placing minor lineages as Other
  lineages_daily_draft$`pango_lineage` = ifelse(lineages_daily_draft$`pango_lineage`%in% rownames(lineages_sum),lineages_daily_draft$`pango_lineage`,"Other")
  
  lineages_daily <- lineages_daily_draft %>% group_by(`Collection.date`,`pango_lineage`) %>% 
    summarise_at(vars(`CopyDate`),list(n = sum)) %>%
    mutate(freq = round(100*n/sum(n),2)) 
  
  colnames(lineages_daily)[2] <- "Lineage names"
  
  # assign color palettes to lineages
  lineage_col = c("#000000","#f90808","#f96708","#f9a508",
                  "#f1d326","#ffff0d","#daff0d","#ebe7e7",
                  "#ffe0d8","#ffe7d8","#fff2d8","#fff9d8","#d3d7fe","#f6fed3","#e3fed3",
                  "#fffed8","#f9ffd8","#edffd8","#d8ffdc", "#f5e5fa","#faf5e5","#faede5","#d3fcfe",
                  "#d8fff4","#d8fdff","#d8f3ff","#d8e2ff","#d0cdff","#f0cdff","#fce2fc","#eaede9",
                  "#f1f1e4","#f1fffd","#f1faff","#f5f1ff","#fff1f1","#fffbf1")
  lineages_daily_sorted = lineages_daily[order(-lineages_daily$freq),]
  lineages_names = unique(lineages_daily_sorted$`Lineage names`)
  names(lineage_col) <- lineages_names
  rm(lineages_names,lineages_daily_draft,lineages_sum)
  write.csv(lineages_daily_sorted,"outputs/lineages_frequency_3_months.csv")

# IV. Weekly variant freq since Jan 2021 ###########
  
  who_weekly <- metadata_CT_who %>% group_by(`week`) %>% count(who_variants) %>%
    mutate(freq = round(n/sum(n),2)) 
  # filter out sequences before 2021
  who_weekly = who_weekly[c(which(who_weekly$week > '2020-12-31')),] 
  # omit latest week (inadequate data)
  who_weekly = who_weekly[-c(which(who_weekly$week == max(who_weekly$week))),] 
  who_weekly$`Variant names` = who_weekly$who_variants
  write.csv(who_weekly,"outputs/variant_frequency.csv")
  
# V. Table 1 ########
  # read last week's Table 1
  table1_old = read.csv("data/Table1.csv")
  table1_old = table1_old[,-c(1)]
  
  # merge old Table 1 with cumulative & recent variant data
  tab = merge(table1_old,cumulative[,c("WHO.label","Cumulative.sequenced.cases.")], by = "WHO.label", all.x = TRUE)
  table1_new = merge(tab,last_3_weeks[,c("WHO.label","Total.sequenced.from.past.3.weeks..","Percent.sequenced.from.past.3.weeks..")], by = "WHO.label", all.x = TRUE)
  table1_new[is.na(table1_new)] <- 0 # assign 0 to NA values
  # find changes since last week's report
  table1_new$`Percent.change.from.previous.report.y` = table1_new$`Percent.sequenced.from.past.3.weeks...y` - table1_new$`Percent.sequenced.from.past.3.weeks...x`
  table1_new = table1_new[,-c(4:7)]
  table1_new[which(table1_new$`WHO label` == 0),1] <- "Other"
  colnames(table1_new) = c("WHO label",
                           "Pango lineage",
                           "CDC classification",
                           "Cumulative sequenced cases*",
                           "Total sequenced from past 3 weeks**",
                           "Percent sequenced from past 3 weeks**",
                           "Percent change from previous report")
  write.csv(table1_new,"data/Table1_new.csv")
  write.csv(table1_new,"outputs/Table1_new.csv")
  
# VI. Subsampler ######
  
  # read cumulative case file
  cases_path = list.files(path = "data/",pattern = "time_series*",full.names = TRUE,recursive=TRUE,include.dirs=TRUE)
  cases = fread(cases_path)
  # filter for Connecticut & remove columns
  cases_CT_t = cases %>%
    filter(Province_State == "Connecticut") %>%
    select(-c(UID,iso2,iso3,code3,FIPS,Admin2,Country_Region,Lat,Long_,Combined_Key,Province_State)) 
  cases_CT = tibble(apply(cases_CT_t,2,sum),mdy(colnames(cases_CT_t)))
  colnames(cases_CT) = c("case_cumulative","date")
  # create week column & daily case count column
  cases_CT_filt = cases_CT %>% 
    mutate(CopyDate = `date`) %>%
    mutate(case_count = case_cumulative - lag(case_cumulative)) %>%
    group_by(`week` = cut(`CopyDate`, "week")) %>%
    select(-c(`CopyDate`,`case_cumulative`)) 
  cases_CT_filt.dt = data.table(cases_CT_filt)
  # sum rows by week to calculate weekly case count
  cases_CT_filt.dt = cases_CT_filt.dt[,list(case_count=sum(case_count)), by='week']
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
  # filter out counts before 2021
  seq_CT_filt.dt = seq_CT_filt.dt[-c(1:which(seq_CT_filt.dt$week == "2020-12-28")),]
  
  # merge case & sequence count files
  subsampler = merge(cases_CT_filt.dt,seq_CT_filt.dt, by = 'week',all.x = TRUE)
  subsampler[is.na(subsampler)] <- 0
  # calculate % sequenced & remove sequence counts
  subsampler = subsampler %>%
    mutate(percent = round(seq_count * 100/ case_count,2)) %>%
    select(-c(seq_count))
  
  rm(cases,cases_CT,cases_CT_filt,seq_CT,seq_CT_filt.dt,cases_CT_filt.dt,cases_path)
  write.csv("outputs/subsampler.csv")
  
  ylim.prim <- c(0, 70000)   
  ylim.sec <- c(0, 100)    
  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- ylim.prim[1] - b*ylim.sec[1]
  
  subsampler_plot <- ggplot(subsampler, aes(x=`week`)) +
    geom_col(aes(y=case_count), color = "darkblue", fill = "#e0f4ff") + 
    geom_line(aes(y=percent*b+a), color = "red", size = 1.5) +
    labs(title = "Percentage (%) of COVID-19 cases sequenced in Connecticut",
         x = "Time (weeks)")+
    xlab("Week") +
    scale_y_continuous("weekly cases", sec.axis = sec_axis(~ (. - a)/b, name = "%seq")) +
    theme_light() +
    theme(plot.title=element_text(size=30, face ="bold"),
          axis.text.y.right = element_text(color = "darkred", size = 10), 
          axis.title.y.right = element_text(color = "darkred", size = 10),
          axis.text.y.left = element_text(color = "darkblue", size = 10), 
          axis.title.y.left = element_text(color = "darkblue", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1)) 
  
# VII. Ct Values by Lineages #############
  
  glab = read.csv("data/GLab_SC2_sequencing_data - Sample metadata.csv")
  ######Assign BA.1 vs BA.2 ######
  glab_filt <- glab %>% mutate(Ctvals = Lineage) 
  glab_filt$Ctvals <- substr(glab_filt$Ctvals, 1,4)
  
  
  #include only samples with lineage assignments
  Ct_values <- glab_filt %>% dplyr::filter(`Ctvals` %in% c("BA.1", "BA.2","BA.4","BA.5"))
  rm(glab,glab_filt)
  
######## IX. PLOTTING ##############
  
  # 1. Daily lineage freq in last 3 months
  
  output$lin_freq <- renderPlotly({
    
    lineage_plot = ggplot(lineages_daily,aes(x=`Collection.date`,y=freq, group=`Lineage names`, color = `Lineage names`)) +
      geom_line() +
      labs(
        title = "Daily lineage frequencies in Connecticut (last 3 months)",
        y = "Lineage frequency (%)",
        x = "Time (days)",
        subtitle = "% for lattest week shown subject to change"
      )+
      theme_light() + 
      scale_color_manual(values=lineage_col) +
      theme(
        plot.title = element_text(face = "bold", size = 22) 
      ) 
    ggplotly(lineage_plot)
  
  })
  
  # 2. Weekly variant freq since 1-2021
  
  output$var_freq <- renderPlotly({
    
    who_colors = c("Other" = "#c4c1c0",
                   "Alpha" = "#1160f6",
                   "Delta" = "#11adf6",
                   "Beta" = "#c81ac8",
                   "Gamma" = "#c81ac8",
                   "Lambda" = "#aef516",
                   "Epsilon" = "#13d84f",
                   "Iota" = "#05d0ab",
                   "Mu" = "#9e5220",
                   "Omicron (BA.1)" = "#f1fa30",
                   "Omicron (BA.2)" = "#fdb902",
                   "Omicron (BA.3)" = "#e58b00",
                   "Omicron (BA.4)" = "#f47110", 
                   "Omicron (BA.5)" = "#d04805"
    )
    who_plot = ggplot(who_weekly,aes(x=`week`,y=freq, group = `Variant names`, color = `Variant names`)) +
      geom_line() +
      labs(
        title = "Weekly variant frequencies in Connecticut (since January 2021)",
        y = "Variant frequency (%)",
        x = "Time (weeks)"
      ) + 
      scale_color_manual(values = c(who_colors)) +
      theme_light() +
      theme(
        plot.title = element_text(face = "bold", size = 22)
      ) 
    ggplotly(who_plot)
    
  })
  
  
  # 3. Table 1
  
  output$table1 <- DT::renderDataTable(
    datatable(table1_new,
              filter = "top",
              caption = htmltools::tags$caption(
                style = 'caption-side: bottom; text-align: center',

                'Table 1: ', htmltools::em('SARS-CoV-2 variants in Connecticut.')),
              options = list(pageLength = 20, 
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

    intercept_char <- Ct_values %>% dplyr::filter(`Ctvals` == c("BA.1"))
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
        plot.title=element_text(size=30, face ="bold"),
        legend.position="none",axis.text = element_text(size = 10, color = "black",face ="bold"), 
                            axis.title=element_text(size=12, face ="bold"))
    Ct_graph
    
  })
  
  output$Ct_caption <- renderInfoBox({
    
    Ct_caption_txt = "A Ct value is defined as the number of amplification cycles required at which point the diagnostic result of the real-time PCR changes from negative 
                      (not detectable) to positive (detectable). 
                      The Ct value reflects the amount of virus copies in a sample, with a lower Ct value indicating more virus. See diagram to the right."
    infoBox(title = HTML("Note:"),
            value = HTML("<p style='font-size:18px'>",
                         Ct_caption_txt,"</p>"),
            icon = icon("creative-commons-sa"),
            fill = TRUE
    )
  })
  
  # 5. Subsampler
  
  output$subsampler <- renderPlot({
  
    subsampler_plot
    
  })
  
  
  # 6. Total sequences info box
  
  output$total_seq <- renderInfoBox({
    

    infoBox(title = HTML("Yale SARS-CoV-2 Genomic Surveillance Initiative: Genomes Published<br>"),
            value = HTML("<p style='font-size:40px', text-align: center>",
                         total_sequences,"</p>"),
            color = "orange",
            fill = TRUE,
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
    
    nextstrain_txt = "We created a custom Nextstrain page to visualize the relatedness of sequenced COVID-19 cases in the state. Below shows a phylogenetic tree of historic sequences in Connecticut, colored by pango lineage assignments."
    
    infoBox(title = HTML("Phylogenetics:"),
            value = HTML("<p style='font-size:18px'>",
                         nextstrain_txt,"</p>"),
            color = "orange",
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
  
})