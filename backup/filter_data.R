#PACKAGE INSTALL ####
packages = c("readxl", "googlesheets4", "ggalluvial", "readr", "dplyr","tidyverse","ggplot2","openxlsx")

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

#####Load Latest GSheet data############
# Authenticate using token. If no browser opens, the authentication works.
gs4_auth(cache = ".secrets", email = "k.pham@yale.edu")

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

######Filter data###### do only use data without the Connecticut tag
Lineages <- Gsheet %>% dplyr::filter(`Filter`%in% c("connecticut"))
######Assign Week######
Lineages <- Lineages %>% mutate(CopyDate = `Collection-date`)

Lineages <- Lineages %>% group_by(`week` = cut(`CopyDate`, "week"))
Lineages$week <- as.Date(Lineages$week)
#delete data without week assigned
Lineages <- Lineages[!is.na(Lineages$`week`),]

#write csv
write.csv(Lineages,"data/lineages.csv")

#select only samples with lineage assignments
Lineages <- Lineages[!is.na(Lineages$`Lineage`),]
#Count Lineages per week
Lineages_Sum <- Lineages %>% group_by(week) %>% count(Lineage)
Lineages_count  <- aggregate(Lineages_Sum$n, by=list(Lineages_Sum$week), FUN = sum)
colnames(Lineages_count) <- c("week","Total")
Lineages_Sum  <- right_join(Lineages_Sum, Lineages_count, by = "week")
Lineages_Sum <- Lineages_Sum %>% mutate(`freq` = ((`n`/`Total`)*100))

#write csv
write.csv(Lineages_Sum,"data/lineages_freq.csv")

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

write.csv(dat.sum,"data/dat_sum.csv")