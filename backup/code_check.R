#CREATE Report on Lineages & Ct values (N) FROM GLAB SEQUENCING DATA GOOGLE SHEETS
#to be updated Tuesday afternoons

#PACKAGE INSTALL ####
packages = c("readxl", "ggalluvial", "readr", "dplyr","tidyverse","ggplot2","RColorBrewer","cowplot","shiny","shinythemes","packcircles","plotly")

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


#shinyServer(function(input,output){
########Weekly Lineages
Lineages = read.csv("data/lineages.csv")
Lineages_Sum = read.csv("data/lineages_freq.csv")
dat.sum = read.csv("data/dat_sum.csv")
summary_plot = read.csv("data/pcr_trends.csv")

cases_weekly = read.table("data/matrix_cases_epiweeks.tsv", sep = "\t", header = TRUE)
cases_weekly_filtered = cases_weekly[-c(1:4),1]
sampling_weekly = read.table("data/weekly_sampling_proportions.tsv", sep = "\t", header = TRUE)
sampling_weekly_filtered = sampling_weekly[-c(1:4),3]
time = cases_weekly[-c(1:4),2]
subsampler = data.frame(time,cases_weekly_filtered,sampling_weekly_filtered)
subsampler$cases_weekly_filtered = as.numeric(subsampler$cases_weekly_filtered)
subsampler$sampling_weekly_filtered = as.numeric(subsampler$sampling_weekly_filtered)
#create color palette
#if new lineages are assigned, palette needs to be updated

OraBlu <- c("AY.103" = "#A233FF",
            "AY.25.1" = "#FF33AE",
            "B.1.517" = "#CC0066", 
            "BA.1" = "#006064",
            "BA.1.1" = "#0097A7",
            "BA.1.1.1"= "#138D75",
            "BA.1.1.14" = "#1F618D", 
            "BA.1.1.16" = "#b8fed9",
            "BA.1.15" = "#2ECC71",
            "BA.1.20" = "#7DD3AB",
            "BA.2" = "#F57C00", 
            "BA.2.1" = "#FFC107",
            "BA.2.10" = "#FF8C3F", 
            "BA.2.10.1" = "#FF6600",
            "BA.2.12.1" = "#FF6666",
            "BA.2.3.4" = "#cc1719", 
            "BA.2.6" = "#ffb370",
            "BA.2.7" = "#df1e32", 
            "BA.2.9" = "#7a0310",
            "BA.2.12" = "#FD8D3C", 
            "BA.2.13" = "#E31A1C",
            "BA.2.14" = "#ea4648",
            "BA.2.18" = "#ea4648",
            "BA.2.20" = "#ea4651",
            "BA.2.22" = "#5d1c1c",
            "BA.2.23" = "#f29091", 
            "BA.2.26" = "#ea6d46",
            "BA.4" = "#663399",
            "BA.5" = "#0d9c08")

#####Visualize Top 100 clonotype representation
#fetch legend for plotting together with bubble plots
legendL <- get_legend(ggplot(Lineages_Sum, aes(x = as.Date(`week`), y = `freq`))+ 
                        geom_col(aes(fill = `Lineage` ), width = 3)+
                        guides(fill=guide_legend(ncol=2))+
                        scale_fill_manual(values=c(OraBlu))+
                        theme(legend.text = element_text(size = 10),legend.title = element_blank())
)

table_1 = read.csv("data/Table1.csv",sheet = 3)

output$table1 <- DT::renderDataTable(
  datatable(table_1)
)

output$lin_freq <- renderPlot({
  LineagePlot <- ggplot(Lineages_Sum, aes(x = as.Date(`week`), y = `freq`, alluvium = `Lineage`)) + 
    geom_alluvium(aes(col = `Lineage`, fill = `Lineage`), alpha = 0.2,show.legend=FALSE) +
    geom_col(aes(fill = `Lineage` ), width = 3,show.legend=FALSE) +
    theme_classic() + 
    scale_x_date(breaks = "1 week") + 
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black",face ="bold"))+
    ggtitle("Lineages per week") + ylab("% of all samples") +
    coord_fixed(ratio = 0.2)+
    scale_fill_manual(values=c(OraBlu)) +
    scale_color_manual(values=c(OraBlu))
  
  LineagePlot
  
})

output$weekly_total <- renderPlot({
  
  WeeklyTotal_Seq <- Lineages %>% group_by(week) %>% count()
  samples_seq <- ggplot(WeeklyTotal_Seq, aes(x = as.Date(`week`), y=`n`, fill = `n`)) +
    geom_bar(width = 3, stat = 'identity') +
    theme_classic() + scale_x_date(breaks = "1 week") + xlab("Week") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
    ggtitle("Sequenced Samples") + ylab("Samples") +
    scale_fill_gradient2(low='red', mid = 'yellow', high='darkgreen', midpoint = 90)+
    coord_fixed(ratio = 0.2) 
  samples_seq
})

output$bubble <- renderPlot({
  
  WeeklyBubbles_1L <- ggplot(dat.sum)+
    geom_polygon(aes(x,y, group = id, fill = factor(Fill)), colour = "black", alpha = 1,show.legend=FALSE)+
    scale_y_reverse()+ theme_void()+
    scale_fill_manual(values = OraBlu)+
    guides(NULL)+
    coord_fixed(ratio = 1/1)+
    facet_wrap(.~sample_name, nrow = 1)+
    coord_fixed(ratio = 1/1)+
    theme(strip.text = element_text(face="bold", size=14,lineheight=5.0))+
    labs(caption= "bubble size represents lineage frequency")
  WeeklyBubbles_1L
})

output$trendplot <- renderPlot({
  trendPlot <- ggplot(summary_plot, aes(as.Date(`Week`))) + 
    geom_point(aes(y=`X.SGTF`), color = "darkred") +
    geom_line(aes(y=`X.SGTF`), color = "darkred")+
    geom_point(aes(y=`X.ORFTF`), color = "darkblue") + 
    geom_line(aes(y=`X.ORFTF`), color = "darkblue")+
    scale_x_date(breaks = "1 week") + xlab("Week") +
    scale_y_continuous("%SGTF", sec.axis = sec_axis(trans=~., name = "%ORFTF")) +
    theme_classic() +
    theme(axis.text.y.right = element_text(color = "darkblue"), 
          axis.title.y.right = element_text(color = "darkblue"),
          axis.text.y.left = element_text(color = "darkred"), 
          axis.title.y.left = element_text(color = "darkred"),
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle("Variant qPCR trends") +
    coord_fixed(ratio = 0.35)
  trendPlot 
})

output$sample_plot <- renderPlot({
  samples <- ggplot(summary_plot, aes(x = as.Date(`Week`), y=`n`, fill = `n`)) +
    geom_bar(width = 3, stat = 'identity') +
    theme_classic() + scale_x_date(breaks = "1 week") + xlab("Week") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
    ggtitle("Processed samples") + ylab("Samples") +
    scale_fill_gradient2(low='lightblue', mid = 'blue', high='darkblue', midpoint = 120)+
    coord_fixed(ratio = 0.2) 
  samples 
})

output$subsampler <- renderPlot({
  
  ylim.prim <- c(0, 70000)   
  ylim.sec <- c(0, 100)    
  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- ylim.prim[1] - b*ylim.sec[1]
  
  subsampler_plot <- ggplot(subsampler, aes(x=`time`)) +
    geom_col(aes(y=cases_weekly_filtered), color = "darkblue") + 
    geom_point(aes(y=sampling_weekly_filtered*b+a), color = "darkred") +
    #geom_line(aes(y=sampling_weekly_filtered*b+a), color = "darkred")+
    xlab("Week") +
    scale_y_continuous("weekly cases", sec.axis = sec_axis(~ (. - a)/b, name = "%seq")) +
    theme_classic() +
    theme(axis.text.y.right = element_text(color = "darkred"), 
          axis.title.y.right = element_text(color = "darkred"),
          axis.text.y.left = element_text(color = "darkblue"), 
          axis.title.y.left = element_text(color = "darkblue"),
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle("Total COVID-19 cases in Connecticut") 
  subsampler_plot
  
})
#})