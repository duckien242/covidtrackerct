#PACKAGE INSTALL ####
packages = c("outbreakinfo", "ggalluvial", "readr", "dplyr","tidyverse","ggplot2","RColorBrewer","cowplot","packcircles")

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

## Read metadata from gisaid
path <- getwd()
filename <- "gisaid_hcov-19.tsv"
metadata <- read_tsv(file = paste0(path, "/data/", filename))
metadata_recent = metadata[metadata$`Collection date` >= "2022-04-04",] # filter by most recent 6 weeks
#metadata_CT = metadata_recent[metadata_recent$`division` == "Connecticut",] # filter by CT
metadata_CT =  metadata_recent[!is.na(metadata_recent$`Lineage`),]
######Assign Week######
metadata_CT$date <- as.Date(metadata_CT$`Collection date`)
metadata_CT <- metadata_CT %>% mutate(CopyDate = `Collection date`)
metadata_CT <- metadata_CT %>% group_by(`week` = cut(`CopyDate`, "week"))
metadata_CT$week <- as.Date(metadata_CT$week)
#delete data without week assigned
metadata_CT <- metadata_CT[!is.na(metadata_CT$`week`),]

#Count Lineages per week
Lineages_Sum <- metadata_CT %>% group_by(week) %>% count(Lineage)
Lineages_count  <- aggregate(Lineages_Sum$n, by=list(Lineages_Sum$week), FUN = sum)
colnames(Lineages_count) <- c("week","Total")
Lineages_Sum  <- right_join(Lineages_Sum, Lineages_count, by = "week")
Lineages_Sum <- Lineages_Sum %>% mutate(`freq` = ((`n`/`Total`)*100))

Lineages_Sum = Lineages_Sum %>%
  add_column(colors = NA)
Lineages_Sum$colors = ifelse(grepl("(BA.1)",Lineages_Sum$Lineage),"#0e6cc9",Lineages_Sum$colors)
Lineages_Sum$colors = ifelse(grepl("(BA.1.1)",Lineages_Sum$Lineage),"#0097A7",Lineages_Sum$colors)
Lineages_Sum$colors = ifelse(grepl("(BA.2)",Lineages_Sum$Lineage),"#e88b43",Lineages_Sum$colors)
Lineages_Sum$colors = ifelse(Lineages_Sum$Lineage == "BA.2","#e69730",Lineages_Sum$colors)
Lineages_Sum$colors = ifelse(grepl("(AY)",Lineages_Sum$Lineage),"#A233FF",Lineages_Sum$colors)
Lineages_Sum$colors = ifelse(grepl("(BA.2.12.1)",Lineages_Sum$Lineage),"#FF6666",Lineages_Sum$colors)
Lineages_Sum$colors = ifelse(grepl("(XE)",Lineages_Sum$Lineage),"black",Lineages_Sum$colors)
Lineages_Sum$colors = ifelse(grepl("(Unassigned)",Lineages_Sum$Lineage),"grey",Lineages_Sum$colors)
Lineages_Sum$colors = ifelse(grepl("(BA.4)",Lineages_Sum$Lineage),"#b85fc9",Lineages_Sum$colors)
Lineages_Sum$colors = ifelse(grepl("(BA.5)",Lineages_Sum$Lineage),"#4aa17d",Lineages_Sum$colors)

# Attempt at plotting
LineagePlot <- ggplot(Lineages_Sum, aes(x = as.Date(`week`), y = `freq`, alluvium = `Lineage`)) + 
  geom_alluvium(aes(col = `Lineage`, fill = `Lineage`), alpha = 0.2,show.legend=FALSE) +
  geom_col(aes(fill = `Lineage` ), width = 3,show.legend=FALSE) +
  theme_classic() + 
  scale_x_date(breaks = "1 week") + 
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black",face ="bold"))+
  ggtitle("Lineages per week") + ylab("% of all samples") +
  coord_fixed(ratio = 0.2)+
  scale_fill_manual(values= Lineages_Sum$colors) +
  scale_color_manual(values= Lineages_Sum$colors)

LineagePlot

########Ct Values by Lineages#############
######Assign BA.1 vs BA.2 ######
BA1vs2 <- metadata_CT %>% mutate(BA1V2 = Lineage) 
BA1vs2$BA1V2 <- substr(BA1vs2$BA1V2, 1,4)
#include only samples with lineage assignments
BA1vs2 <- BA1vs2 %>% dplyr::filter(`BA1V2` %in% c("BA.1", "BA.5"))

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

