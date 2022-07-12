#CREATE Report on Lineages & Ct values (N) FROM GLAB SEQUENCING DATA GOOGLE SHEETS
#to be updated Tuesday afternoons

#PACKAGE INSTALL ####
packages = c("readxl", "googlesheets4", "ggalluvial", "readr", "dplyr","tidyverse","ggplot2","RColorBrewer","cowplot","openxlsx")

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

# Set authentication token to be stored in a folder called `.secrets`
#options(gargle_oauth_cache = ".secrets")

# Authenticate manually
#gs4_auth()


# If successful, the previous step stores a token file.
# Check that a file has been created with:
list.files(".secrets/")

# Check that the non-interactive authentication works by first deauthorizing:
gs4_deauth()

# Authenticate using token. If no browser opens, the authentication works.
gs4_auth(cache = ".secrets", email = "k.pham@yale.edu")


#ss <- gs4_get("https://docs.google.com/spreadsheets/d/1cHCdfy1WiIOFFJDEpo3U97KaiPep0Tv_dvtIEEcnE3c")
#sheet_append(ss, data.frame(time=Sys.time()))
#####Load Latest GSheet data############
#get colnames
#Gsheetnames = read_sheet("1dUm-OtDOxvbS9OdCfnBrjBF9H5AKphauzzWxY-h4oUQ", 
#                         sheet = "New data for Shiny", 
#                         range = cell_rows(1)
#)


#drop columns that are not needed
#Gsheetnames <- Gsheetnames[,c(1:15,22,23)]
#read data (data range has to be adapted!-> row names are Yale-ID-1)

#Gsheet =read_sheet("1dUm-OtDOxvbS9OdCfnBrjBF9H5AKphauzzWxY-h4oUQ", 
#                   sheet = "New data for Shiny", 
##                   col_names = TRUE
#)
#####Load Latest GSheet data############
#get colnames
Gsheetnames = read_sheet("1cHCdfy1WiIOFFJDEpo3U97KaiPep0Tv_dvtIEEcnE3c", 
                         sheet = "Sample metadata", 
                         range = cell_rows(1)
                         )


#drop columns that are not needed
Gsheetnames <- Gsheetnames[,c(1:15,22,23)]
#read data (data range has to be adapted!-> row names are Yale-ID-1)
#earliest Yale-ID
EY = 16965
#latest Yale-ID
LY = 20100
Gsheet =read_sheet("1dUm-OtDOxvbS9OdCfnBrjBF9H5AKphauzzWxY-h4oUQ", 
                   sheet = "Sample metadata", range = cell_rows(EY:LY), 
                   col_names = FALSE,
                   col_types = "ccncDccDccnnnnc------cc-"
                   )

#assign columnnames
colnames(Gsheet) <- colnames(Gsheetnames)

#only use data without "do not report"
Gsheet <- Gsheet %>% filter(!grepl('do not report', Filter))

########Weekly Lineages
Gsheet = read_sheet("data/GLab_SC2_sequencing_data.xlsx",
                   sheet = "Sample metadata")
Gsheet$`Collection-date` <- openxlsx::convertToDate(Gsheet$`Collection-date`)

######Filter data###### do only use data without the Connecticut tag
Lineages <- Gsheet %>% dplyr::filter(`Filter`%in% c("connecticut"))
######Assign Week######
Lineages <- Lineages %>% mutate(CopyDate = `Collection-date`)

Lineages <- Lineages %>% group_by(`week` = cut(`CopyDate`, "week"))
Lineages$week <- as.Date(Lineages$week)
#delete data without week assigned
Lineages <- Lineages[!is.na(Lineages$`week`),]
####plot sample numbers for sequencing
WeeklyTotal_Seq <- Lineages %>% group_by(week) %>% count()
samples_seq <- ggplot(WeeklyTotal_Seq, aes(x = as.Date(`week`), y=`n`, fill = `n`)) +
  geom_bar(width = 3, stat = 'identity') +
  theme_classic() + scale_x_date(breaks = "1 week") + xlab("Week") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
  ggtitle("Sequenced Samples") + ylab("Samples") +
  scale_fill_gradient2(low='red', mid = 'yellow', high='darkgreen', midpoint = 90)+
  coord_fixed(ratio = 0.2) 
samples_seq
#select only samples with lineage assignments
Lineages <- Lineages[!is.na(Lineages$`Lineage`),]
#Count Lineages per week
Lineages_Sum <- Lineages %>% group_by(week) %>% count(Lineage)
Lineages_count  <- aggregate(Lineages_Sum$n, by=list(Lineages_Sum$week), FUN = sum)
colnames(Lineages_count) <- c("week","Total")
Lineages_Sum  <- right_join(Lineages_Sum, Lineages_count, by = "week")
Lineages_Sum <- Lineages_Sum %>% mutate(`freq` = ((`n`/`Total`)*100))


#Plot
#colors for plot
#if new lineages are assigned, palette needs to be updated
#check with 
unique(Lineages_Sum$Lineage)
OraBlu <- c("AY.103" = "#A233FF",
            "AY.25.1" = "#FF33AE",
            "AY.103"  = "#FF3345",
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
            "BA.2.3" = "#cc1719", 
            "BA.2.6" = "#ffb370",
            "BA.2.7" = "#df1e32", 
            "BA.2.9" = "#7a0310",
            "BA.2.12" = "#FD8D3C", 
            "BA.2.13" = "#E31A1C",
            "BA.2.14" = "#ea4648",
            "BA.2.22" = "#5d1c1c",
            "BA.2.23" = "#f29091", 
            "BA.2.26" = "#ea6d46",
            "BA.3" = "#663399")

LineagePlot <- ggplot(Lineages_Sum, aes(x = as.Date(`week`), y = `freq`, alluvium = `Lineage`)) + 
  geom_alluvium(aes(col = `Lineage`, fill = `Lineage`), alpha = 0.2,show.legend=FALSE) +
  geom_col(aes(fill = `Lineage` ), width = 3,show.legend=FALSE) +
  theme_classic() + scale_x_date(breaks = "1 week") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black",face ="bold"))+
  ggtitle("Lineages per week") + ylab("% of all samples") +
  coord_fixed(ratio = 0.2)+
  scale_fill_manual(values=c(OraBlu)) +
  scale_color_manual(values=c(OraBlu))
LineagePlot

#fetch legend for plotting together with bubble plots
legendL <- get_legend(ggplot(Lineages_Sum, aes(x = as.Date(`week`), y = `freq`))+ 
                        geom_col(aes(fill = `Lineage` ), width = 3)+
                        guides(fill=guide_legend(ncol=2))+
                        scale_fill_manual(values=c(OraBlu))+
                        theme(legend.text = element_text(size = 10),legend.title = element_blank())
)

####Visualize per Week######
#####Visualize Top 100 clonotype representation
library(packcircles)
###change data of last 4 weeks here
Week1 = c("2022-04-04")
Week2 = c("2022-04-11")
Week3 = c("2022-04-18")
Week4 = c("2022-04-25")
#1 Week
Week <- Lineages_Sum %>% dplyr::filter(`week` == Week1)
packing <- circleProgressiveLayout(Week$freq)
dat.gg <- circleLayoutVertices(packing) 
Fill <- data.frame(lapply(Week, function(x) rep(Week$Lineage, each = 26)))
Fill <- Fill[,1]
dat.w <- cbind(dat.gg, Fill)
sample_name <- c(rep(c(Week1),length.out = (nrow(dat.w))))
Week1 <- cbind(dat.w, sample_name)
rm(Week, Fill, dat.gg, dat.w, packing)

#2 Week
Week <- Lineages_Sum %>% dplyr::filter(`week` == Week2)
packing <- circleProgressiveLayout(Week$freq)
dat.gg <- circleLayoutVertices(packing) 
Fill <- data.frame(lapply(Week, function(x) rep(Week$Lineage, each = 26)))
Fill <- Fill[,1]
dat.w <- cbind(dat.gg, Fill)
sample_name <- c(rep(c(Week2),length.out = (nrow(dat.w))))
Week2 <- cbind(dat.w, sample_name)
rm(Week, Fill, dat.gg, dat.w, packing)

#3 Week
Week <- Lineages_Sum %>% dplyr::filter(`week` == Week3)
packing <- circleProgressiveLayout(Week$freq)
dat.gg <- circleLayoutVertices(packing) 
Fill <- data.frame(lapply(Week, function(x) rep(Week$Lineage, each = 26)))
Fill <- Fill[,1]
dat.w <- cbind(dat.gg, Fill)
sample_name <- c(rep(c(Week3),length.out = (nrow(dat.w))))
Week3 <- cbind(dat.w, sample_name)
rm(Week, Fill, dat.gg, dat.w, packing)

#4 Week
Week <- Lineages_Sum %>% dplyr::filter(`week` == Week4)
packing <- circleProgressiveLayout(Week$freq)
dat.gg <- circleLayoutVertices(packing) 
Fill <- data.frame(lapply(Week, function(x) rep(Week$Lineage, each = 26)))
Fill <- Fill[,1]
dat.w <- cbind(dat.gg, Fill)
sample_name <- c(rep(c(Week4),length.out = (nrow(dat.w))))
Week4 <- cbind(dat.w, sample_name)
rm(Week, Fill, dat.gg, dat.w, packing)

#Summarize
dat.sum <- rbind(Week1,Week2,Week3,Week4)

#Plot
#2x2
#WeeklyBubbles <- ggplot(dat.sum)+
#  geom_polygon(aes(x,y, group = id, fill = factor(Fill)), colour = "black", alpha = 1)+
# scale_y_reverse()+ theme_void()+
#  scale_fill_manual(values = OraBlu)+
#  coord_fixed(ratio = 1/1)+
#  facet_wrap(.~sample_name, nrow = 2) + labs(title = "Overview last 4 weeks") +
#  guides(fill=guide_legend(ncol= 2, title = NULL))+ coord_fixed(ratio = 1/1)
#WeeklyBubbles

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

######Plot All#####
L <- plot_grid(samples_seq, NULL, LineagePlot, legendL, WeeklyBubbles_1L, NULL, 
               rel_widths = c(1,0.2),rel_heights = c(0.9,1.1,0.8),
               ncol = 2, nrow = 3, axis = "b", align = "v")
L

#Write csv files
write.csv(Lineages_Sum,"data/lineages_frequency.csv")
write.csv(Lineages,"data/lineages.csv")
write.csv(dat.sum,"data/dat_sum.csv")

########Ct Values by Lineages#############
######Assign BA.1 vs BA.2 ######
BA1vs2 <- Gsheet %>% mutate(BA1V2 = Lineage) 
BA1vs2$BA1V2 <- substr(BA1vs2$BA1V2, 1,4)
#include only samples with lineage assignments
BA1vs2 <- BA1vs2 %>% dplyr::filter(`BA1V2` %in% c("BA.1", "BA.2"))

#Count
countBA1 <- count(BA1vs2 %>% dplyr::filter(`BA1V2` == c("BA.1")))
countBA2 <- count(BA1vs2 %>% dplyr::filter(`BA1V2` == c("BA.2")))

#######Plot#########
intercept <- BA1vs2 %>% dplyr::filter(`BA1V2` == c("BA.1"))
intercept <- intercept$`Yale-N1(FAM)`
intercept <- mean(intercept)

BA1_2 <- ggplot(data = BA1vs2, aes(x= `BA1V2`, y = `Yale-N1(FAM)`))+
  geom_hline(yintercept = intercept,colour = "lightgrey")+
  geom_boxplot(width=0.4,outlier.shape = NA,colour = "#666666")+
  geom_jitter(aes(col =`BA1V2`), alpha = 0.2,size = 0.2, stroke = 2,shape = 21, width = 0.15)+
  labs(title="BA.1 vs BA.2 Ct value (N)",x=NULL, y = "Ct (N)")+
  scale_color_brewer(palette="Dark2")+
  coord_fixed(ratio = 0.06)+
  scale_y_reverse(breaks = seq(10, 40, by = 5))+
  theme_classic()+theme(legend.position="none",axis.text = element_text(size = 10, color = "black",face ="bold"), 
                        axis.title=element_text(size=12, face ="bold"),title=element_text(size=14, face ="bold"))
BA1_2 

####by vaccination status#####
###filter out NA for vaccination status####
BA1vs2 <- BA1vs2 %>% dplyr::filter(`Vaccinated` %in% c("Yes", "No"))
BA1_2_vac <- ggplot(data = BA1vs2, aes(x= `BA1V2`, y = `Yale-N1(FAM)`,col=`Vaccinated`, group = interaction(BA1V2,Vaccinated)))+
  geom_hline(yintercept = intercept,colour = "lightgrey")+
  geom_boxplot(width=0.65,outlier.shape = NA, colour = "#666666")+
  geom_point(position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.6,),alpha = 0.4,size = 0.2,stroke =2, shape = 21)+
  labs(title="Ct value (N) by Vaccination Status",x=NULL, y = "Ct (N)")+
  scale_color_brewer(palette="Set1")+
  coord_fixed(ratio = 0.06)+
  scale_y_reverse(breaks = seq(10, 43, by = 5))+
  theme_classic()+theme(legend.position="none",
    axis.text = element_text(size = 10, color = "black",face ="bold"), 
    axis.title=element_text(size=12, face ="bold"),
    title=element_text(size=14, face ="bold"))
BA1_2_vac
#plot_grid(BA1_2,BA1_2_vac,rel_widths = c(1, 1.27))

#####Plot BA.2 sublineages######
filter2 <- BA1vs2 %>% dplyr::filter(`BA1V2` %in% c("BA.2"))
#only choose most abundant BA.2 sublineages 
filter2 <- BA1vs2 %>% dplyr::filter(`Lineage` %in% c("BA.2","BA.2.9","BA.2.3","BA.2.7","BA.2.12.1"))

#assign colors,if more sublineages has to be changed accordingly
BApal <- c("#FD8D3C", "#FC4E2A", "#E31A1C" ,"#BD0026", "#800026")
#Plot
BA2_sub <- ggplot(data = filter2, aes(x= `Lineage`, y = `Yale-N1(FAM)`))+
  geom_hline(yintercept = intercept,colour = "lightgrey")+
  geom_boxplot(width=0.5,outlier.shape = NA,colour = "#666666")+
  geom_jitter(alpha = 0.5,size = 0.5, stroke = 2,shape = 21, width = 0.15, aes(col=`Lineage`))+
  labs(title="BA.2 sublineages Ct value (N)",x=NULL, y = "Ct (N)")+
  scale_color_manual(values =BApal)+
  coord_fixed(ratio = 0.15)+
  scale_y_reverse(breaks = seq(10, 43, by = 5))+
  theme_classic()+theme(legend.position="none",axis.text = element_text(
    size = 10, color = "black",face ="bold"), axis.title=element_text(
      size=12, face ="bold"),title=element_text(size=14, face ="bold"))
BA2_sub

#ggplot(data = filter2, aes(x= `Lineage`, y = `Yale-N1(FAM)`))+
#  geom_hline(yintercept = intercept,colour = "lightgrey")+
#  geom_boxplot(width=0.5,outlier.shape = NA)+
#  geom_jitter(alpha = 0.5,size = 0.5, stroke = 2,shape = 21, 
#              width = 0.1, aes(col=`Vaccinated`))+
#  labs(title="BA.2 sublineages Ct value (N)",x=NULL, y = "Ct (N)")+
#  scale_color_brewer(palette="Spectral",n = 5, direction = -1)+
#  coord_fixed(ylim=c(13,41),ratio = 0.15)+
#  scale_y_reverse(breaks = seq(10, 43, by = 5))+
#  theme_classic()+theme(#legend.position="none",
#                        axis.text = element_text(size = 10, color = "black",face ="bold"), 
#                        axis.title=element_text(size=12, face ="bold"),
#                        title=element_text(size=14, face ="bold"))

####by vaccination status#####
BA2_vac <- ggplot(data = filter2, aes(x= `Lineage`, y = `Yale-N1(FAM)`, col=`Vaccinated`,group = interaction(Lineage,Vaccinated)))+
  geom_hline(yintercept = intercept,colour = "lightgrey")+
  geom_boxplot(width=0.7,outlier.shape = NA,show.legend=FALSE,colour = "#666666", aes(fill = NULL))+
  geom_point(position=position_jitterdodge(
    jitter.width = 0.15,dodge.width = 0.7,),
    alpha = 0.4,size = 0.2,stroke =2, shape = 21,
    show.legend=FALSE)+
  labs(title="Ct value (N) by Vaccination Status",x=NULL, y = "Ct (N)")+
  scale_color_brewer(palette = "Set1")+
  coord_fixed(ratio = 0.15)+
  scale_y_reverse(breaks = seq(10, 43, by = 5))+
  theme_classic() + theme(
    axis.text = element_text(size = 10, color = "black",face ="bold"), 
    axis.title=element_text(size=12, face ="bold"),
    title=element_text(size=14, face ="bold"),
    legend.title=element_text(size=14))
BA2_vac

legend <- get_legend(ggplot(data = filter2, aes(x= `Lineage`, y = `Yale-N1(FAM)`, col=`Vaccinated`))+
                       geom_hline(yintercept = intercept,colour = "lightgrey")+
                       geom_boxplot(width=0.65,outlier.shape = NA)+
                       #geom_violin(scale = "width",adjust = 0.5)+
                       geom_point(position=position_jitterdodge(
                         jitter.width = 0.2,dodge.width = 0.5,),
                         alpha = 0.4,size = 0.2,stroke =2, shape = 21)+
                       #geom_jitter(alpha = 0.5,size = 0.5, stroke = 2,shape = 21)+
                       #stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="pointrange")+
                       labs(title="Ct value (N) by Vaccination Status",x=NULL, y = "Ct (N)")+
                       scale_color_brewer(palette="Set1")+
                       coord_fixed(ratio = 0.15)+
                       scale_y_reverse(breaks = seq(10, 43, by = 5))+
                       theme_classic() + theme(
                         axis.text = element_text(size = 10, color = "black"), 
                         axis.title=element_text(size=12, face ="bold"),
                         title=element_text(size=14, face ="bold"),
                         legend.title=element_text(size=10)
                       ))
C <- plot_grid(BA1_2,BA2_sub, NULL, BA1_2_vac,BA2_vac,legend, rel_widths = c(1, 1,0.45), ncol = 3, align = "h", axis = "l")
C
