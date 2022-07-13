#CREATE Report on PCR trends & sample numbers FROM GLAB SEQUENCING DATA GOOGLE SHEETS
#to be updated Friday afternoons

#PACKAGE INSTALL ####
packages = c("readxl", "googlesheets4", "ggalluvial", "readr", "dplyr","tidyverse","ggplot2","RColorBrewer","cowplot","plotly")

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


########PCR Trends#######################
#####Load Latest GSheet data############
#get colnames
Gsheetnames = read_sheet("1dUm-OtDOxvbS9OdCfnBrjBF9H5AKphauzzWxY-h4oUQ", 
                         sheet = "Sample metadata", 
                         range = cell_rows(1)
)
#select relevant columns only
Gsheetnames <- Gsheetnames[,c(1:15,22,23)]
#read data (data range has to be adapted!-> row names are Yale-ID-1)
#earliest Yale-ID
EY = 16965
#latest Yale-ID
LY = 20100
#fetches data from Gsheet, takes a while to run
Gsheet =read_sheet("1dUm-OtDOxvbS9OdCfnBrjBF9H5AKphauzzWxY-h4oUQ", 
                   sheet = "Sample metadata", range = cell_rows(EY:LY), 
                   col_names = FALSE,
                   col_types = "ccncDccDccnnnnc------cc-"
)
#assing correct colnames
colnames(Gsheet) <- colnames(Gsheetnames)

#only use data without "do not report"
Gsheet <- Gsheet %>% filter(!grepl('do not report', Filter))
#only use data without "not sequencing"
Gsheet <- Gsheet %>% filter(!grepl('not sequencing', `Tube Label`))

######Assign Week######
trends <- Gsheet %>% mutate(CopyDate = `Collection-date`)
#only use data for which there is a N value 
trends <- trends[!is.na(trends$`Yale-N1(FAM)`),]
trends <- trends %>% group_by(`week` = cut(`CopyDate`, "week"))
trends$week <- as.Date(trends$week)
#delete data without week assigned
trends <- trends[!is.na(trends$`week`),]
#assign SGTF or ORFTF
#SGTF
SGTF <- as.data.frame(trends$`Yale-69/70(HEX)`)
colnames(SGTF) = c("SGTF")
SGTF$SGTF[!is.na(SGTF$SGTF)]<-"FALSE"
SGTF$SGTF[is.na(SGTF$SGTF)]<-"TRUE"
#ORFTF
ORFTF <- as.data.frame(trends$`Yale-ORF1a(Cy5)`)
colnames(ORFTF) = c("ORFTF")
ORFTF$ORFTF[!is.na(ORFTF$ORFTF)]<-"FALSE"
ORFTF$ORFTF[is.na(ORFTF$ORFTF)]<-"TRUE"
#Join
TF <- cbind(trends,SGTF,ORFTF)
#Count
SGTF_week <- TF %>% count(week, as.factor(SGTF),.drop = FALSE)
colnames(SGTF_week) = c("week","SGTF","SGTF_n")
ORFTF_week <- TF %>% count(week, as.factor(ORFTF),.drop = FALSE)
colnames(ORFTF_week) = c("week","ORFTF","ORFTF_n")
summary <- cbind(SGTF_week,ORFTF_week)
summary <- summary[,c(1,2,3,5,6)]
colnames(summary) <- c("Week", "SGTF","SGTF(n)","ORFTF","ORFTF(n)")
#filter for plotting
summary_plot <- summary %>% dplyr::filter(`SGTF` %in% c("TRUE"))
WeeklyTotal <- trends %>% group_by(week) %>% count()
summary_plot <- cbind(summary_plot, WeeklyTotal)
summary_plot <- summary_plot[,c(-6)]
summary_plot <- summary_plot %>% mutate(`%SGTF` = (`SGTF(n)`/`n`)*100)
summary_plot <- summary_plot %>% mutate(`%ORFTF` = (`ORFTF(n)`/`n`)*100)
####Plotting
#preparing for two y-axis
ylim.prim <- c(0, 100)   
ylim.sec <- c(0, 100)    
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]
#Plot PCR Trends
trendPlot <- ggplot(summary_plot, aes(as.Date(`Week`))) + 
  geom_point(aes(y=`X.SGTF`), color = "darkred") +
  geom_line(aes(y=`X.SGTF`), color = "darkred")+
    geom_point(aes(y=`X.ORFTF`), color = "darkblue") + 
  geom_line(aes(y=`X.ORFTF`), color = "darkblue")+
    scale_x_date(breaks = "1 week") + xlab("Week") +
  scale_y_continuous("%SGTF", sec.axis = sec_axis(~ (. - a)/b, name = "%ORFTF")) +
  theme_classic() +
  theme(axis.text.y.right = element_text(color = "darkblue"), 
        axis.title.y.right = element_text(color = "darkblue"),
        axis.text.y.left = element_text(color = "darkred"), 
        axis.title.y.left = element_text(color = "darkred"),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Variant qPCR trends") +
  coord_fixed(ratio = 0.35)
ggplotly(trendPlot) 

#Plot Sample numbers (Variant PCR)
samples <- ggplot(summary_plot, aes(x = as.Date(`Week`), y=`n`, fill = `n`)) +
  geom_bar(width = 3, stat = 'identity') +
  theme_classic() + scale_x_date(breaks = "1 week") + xlab("Week") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  ggtitle("Processed samples") + ylab("Samples") +
  scale_fill_gradient2(low='lightblue', mid = 'blue', high='darkblue', midpoint = 120)+
  coord_fixed(ratio = 0.2) 
samples 

P <- plot_grid(trendPlot,samples, rel_widths = c(1, 1.1), ncol = 2, align = "h", axis = "l")
P
ggsave(filename="~/Documents/PostDoc/Projects/SC2_Surveillance/Trends.png", units = "cm", height = 10, scale = 1)

cases_weekly = read.table("data/matrix_cases_epiweeks.tsv", sep = "\t", header = TRUE)
cases_weekly_filtered = cases_weekly[-c(1:4),1]
sampling_weekly = read.table("data/weekly_sampling_proportions.tsv", sep = "\t", header = TRUE)
sampling_weekly_filtered = sampling_weekly[-c(1:4),3]
time = cases_weekly[-c(1:4),2]
subsampler = data.frame(time,cases_weekly_filtered,sampling_weekly_filtered)
subsampler$cases_weekly_filtered = as.numeric(subsampler$cases_weekly_filtered)
subsampler$sampling_weekly_filtered = as.numeric(subsampler$sampling_weekly_filtered)

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
ggplotly(subsampler_plot, dynamicTicks = TRUE) %>%
  rangeslider() %>%
  layout(hovermode = "x")

