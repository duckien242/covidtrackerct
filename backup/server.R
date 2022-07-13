#CREATE Report on Lineages & Ct values (N) FROM GLAB SEQUENCING DATA GOOGLE SHEETS
#to be updated Tuesday afternoons

# I. PACKAGE INSTALL ####
packages = c("dplyr","ggplot2","RColorBrewer","shiny","shinythemes","plotly","DT","htmltools")
library(readr)

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

######### START SHINY APP SERVER #########
shinyServer(function(input,output){
  
# II. Read metadata files ########
  
  metadata_CT_all = read.csv("data/metadata_CT_all.csv")
  metadata_CT_all = metadata_CT_all[,-c(1)]
  metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "BA.2.12.1 (marker override: BA.2.12 + Spike_L452Q => BA.2.12.1)"] = "BA.2.12.1"
  metadata_CT_recent = metadata_CT_all[metadata_CT_all$`date` >= "2022-03-22",] # CHANGE DATE IN WEEKLY UPDATE
  metadata_CT_recent =  metadata_CT_recent[!is.na(metadata_CT_recent$`pango_lineage`),]
  
  metadata_CT_all$`date` <- as.Date(metadata_CT_all$`date`)
  metadata_CT_all <- metadata_CT_all %>% mutate(CopyDate = `date`)
  metadata_CT_all <- metadata_CT_all %>% group_by(`week` = cut(`CopyDate`, "week"))
  metadata_CT_all$week <- as.Date(metadata_CT_all$week)
  metadata_CT_all$`date` <- as.Date(metadata_CT_all$`date`)
  #delete data without week assigned
  metadata_CT_all <- metadata_CT_all[!is.na(metadata_CT_all$`week`),]
  
  # II.1. Create WHO variant name file ######
  names = read.csv("data/variant_names.csv")
  metadata_CT_who = merge(metadata_CT_all,names[,c("pango_lineage","who_variants")], by = "pango_lineage", all.x = TRUE)
  
  cumulative <- metadata_CT_who %>% 
    select(`date`,`who_variants`,`pango_lineage`) %>%
    group_by() %>% count(`who_variants`)
  colnames(cumulative) = c("WHO.label","Cumulative.sequenced.cases.")
  
  last_3_weeks <- metadata_CT_who %>% 
    select(`date`,`who_variants`,`pango_lineage`) %>%
    subset(`date` >= "2022-06-02") %>% # update date every 3 weeks
    count(`who_variants`) %>%
    mutate(freq = round(100*n/sum(n),2)) 
  colnames(last_3_weeks) = c("WHO.label","Total.sequenced.from.past.3.weeks..","Percent.sequenced.from.past.3.weeks..")
  
  metadata_CT_who[which(is.na(metadata_CT_who$who_variants)),31] <- 'Other'
  metadata_CT_who[which(metadata_CT_who$who_variants == "Eta"),31] <- 'Other'
  metadata_CT_who[which(metadata_CT_who$who_variants == "Kappa"),31] <- 'Other'
  metadata_CT_who[which(metadata_CT_who$who_variants == "Zeta"),31] <- 'Other'
  metadata_CT_who = metadata_CT_who[c(which(metadata_CT_who$`date` >= "2021-01-01")),]
  
# III. Daily lineage freq in last 3 months #######

  lineages_daily <- metadata_CT_recent %>% group_by(`date`) %>% count(pango_lineage) %>%
    mutate(freq = round(100*n/sum(n),2)) %>%
    subset(freq > 5)
  colnames(lineages_daily)[2] <- "Lineage names"

# IV. Weekly variant freq since Jan 2021 ###########
  
  who_weekly <- metadata_CT_who %>% group_by(`week`) %>% count(who_variants) %>%
    mutate(freq = round(100*n/sum(n),2)) 
  colnames(who_weekly)[2] <- "Variant names"
  
# V. Table 1 ########
  table1_old = read.csv("data/Table1.csv")
  
  tab = merge(table1_old,cumulative[,c("WHO.label","Cumulative.sequenced.cases.")], by = "WHO.label", all.x = TRUE)
  
  table1_new = merge(tab,last_3_weeks[,c("WHO.label","Total.sequenced.from.past.3.weeks..","Percent.sequenced.from.past.3.weeks..")], by = "WHO.label", all.x = TRUE)
  table1_new[is.na(table1_new)] <- 0
  table1_new$`Percent.change.from.previous.report.y` = table1_new$`Percent.sequenced.from.past.3.weeks...y` - table1_new$`Percent.sequenced.from.past.3.weeks...x`
  table1_new = table1_new[,-c(4:7)]
  table1_new[which(table1_new$`WHO.label` == 0),1] <- "Other"
  colnames(table1_new) = c("WHO label",
                           "Pango lineage",
                           "CDC classification",
                           "Cumulative sequenced cases*",
                           "Total sequenced from past 3 weeks**",
                           "Percent sequenced from past 3 weeks**",
                           "Percent change from previous report")
  
  write.csv(table1_new,"data/Table1_new.csv")
  
# VI. Subsampler ######
  
  cases = read.table("data/matrix_cases_epiweeks.tsv", sep = "\t", header = TRUE)
  cases_CT = t(cases[which(cases$code == "CT"),])
  cases_CT_filt = cases_CT[-c(1:3)]
  
  sampling = read.table("data/weekly_sampling_proportions.tsv", sep = "\t", header = TRUE)
  sampling_CT = t(sampling[which(sampling$code == "CT"),])
  sampling_CT_filt = sampling_CT[-c(1:3)]
  
  ### UPDATE THIS EVERY WEEK
  week_CT = sort(c(unique(metadata_CT_all$`week`),"2022-06-13","2022-06-20","2022-06-27")) 
  
  subsampler_draft = data.frame(week_CT,cases_CT_filt,sampling_CT_filt)
  subsampler = subsampler_draft[-c(1:41),]
  subsampler$cases_CT_filt = as.numeric(subsampler$cases_CT_filt)
  subsampler$sampling_CT_filt = 100*as.numeric(subsampler$sampling_CT_filt)
  
  ylim.prim <- c(0, 70000)   
  ylim.sec <- c(0, 50)    
  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- ylim.prim[1] - b*ylim.sec[1]
  
  rm(cases,cases_CT,cases_CT_filt,sampling,sampling_CT,sampling_CT_filt,week_CT)
  
  subsampler_plot <- ggplot(subsampler, aes(x=`week_CT`)) +
    geom_col(aes(y=cases_CT_filt), color = "darkblue", fill = "#e0f4ff") + 
    geom_line(aes(y=sampling_CT_filt*b+a), color = "red", size = 1.5) +
    #geom_line(aes(y=sampling_weekly_filtered*b+a), color = "darkred")+
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
    
    # plot
    lineage_plot = ggplot(lineages_daily,aes(x=`date`,y=freq, group=`Lineage names`, color = `Lineage names`)) +
      geom_line() +
      labs(
        title = "Daily lineage frequencies in Connecticut (last 3 months)",
        y = "Lineage frequency (%)",
        x = "Time (days)",
        caption = "% for lattest week shown subject to change"
      ) + 
      scale_color_brewer(type = "seq", palette = "Spectral") +
      scale_x_discrete(breaks = c("2022-04-01","2022-05-01","2022-06-01")) +
      theme_light() +
      theme(
        plot.title = element_text(face = "bold", size = 22),
        axis.text.x = element_text(hjust = 1, size = 10, color = "black")
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
    
    total_sequences = 16739
    
    infoBox(title = HTML("Yale SARS-CoV-2 Genomic Surveillance Initiative: Genomes Published<br>"),
            value = HTML("<p style='font-size:40px'>",
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
    my_test <- tags$iframe(src=test, height = 800, width = 500)
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
  
  output$tweet <-   renderUI({
    tagList(
      tags$blockquote(class = "twitter-tweet",
                      tags$a(href = "https://twitter.com/NathanGrubaugh/status/1514688634317856770")),
      tags$script('twttr.widgets.load(document.getElementById("tweet"));')
    )
  })

  output$twitter_caption <- renderInfoBox({
    
    twitter_txt = "Follow us on Twitter @CovidCT"
    
    infoBox(title = HTML("Twitter:"),
            value = HTML("<p style='font-size:18px'>",
                         twitter_txt,"</p>"),
            color = "orange",
            icon = icon("twitter"),
            fill = TRUE,
            width = 4
    )
  })
  
})