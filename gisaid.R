#PACKAGE INSTALL ####
packages = c("dplyr","tidyr","ggplot2","RColorBrewer","shiny","shinythemes","plotly","DT","lubridate","sf","data.table")
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

# Set recent dates
date_3_months = as.Date("2022-04-13")
date_3_weeks = as.Date("2022-06-15") # CHANGE TO MOST RECENT 3 WEEKS
date_gisaid = as.Date("2022-06-24") # FOR FILTERING GISAID DOWNLOAD

## Read metadata from gisaid
metadata_CT_old = read.table("data/metadata_CT_all.tsv", sep = "\t", header = TRUE)
metadata_CT_old$`Collection.date` <- as.Date(metadata_CT_old$`Collection.date`)
new_genomes_title<-list.files('data',pattern = "gisaid",full.names = TRUE,recursive=TRUE,include.dirs=TRUE) # this will now be a list of all files with full paths.
new_genomes = read.table(new_genomes_title, sep = "\t", header = TRUE)
new_genomes$`Collection.date` <- as.Date(new_genomes$`Collection.date`)
metadata_CT_all = rbind(metadata_CT_old,new_genomes)
metadata_CT_all = subset(metadata_CT_all, !duplicated(subset(metadata_CT_all, select=c(`Accession.ID`)))) # remove duplicate sequences
write_tsv(metadata_CT_all,"data/metadata_CT_all.tsv") # write tsv file for next week
write_tsv(metadata_CT_all,"backup/metadata_CT_all.tsv") # write tsv file for next week

metadata_CT_all$pango_lineage <- metadata_CT_all$Lineage
metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "BA.2.12.1 (marker override: BA.2.12 + Spike_L452Q => BA.2.12.1)"] = "BA.2.12.1"

metadata_CT_all <- metadata_CT_all %>% mutate(CopyDate = `Collection.date`)
metadata_CT_all <- metadata_CT_all %>% group_by(`week` = cut(`CopyDate`, "week"))
metadata_CT_all$week <- as.Date(metadata_CT_all$week)
#delete data without week assigned
metadata_CT_all <- metadata_CT_all[!is.na(metadata_CT_all$`week`),]

metadata_CT_recent = metadata_CT_all[metadata_CT_all$`Collection.date` >= date_3_months,] # filter by last 3 months
metadata_CT_recent =  metadata_CT_recent[!is.na(metadata_CT_recent$`pango_lineage`),]
# fix weird BA.2.12.1 lineage name call

# Lineages daily frequency
length(unique(metadata_CT_recent$Lineage))
metadata_CT_recent$`CopyDate` = as.character(metadata_CT_recent$`CopyDate`)

metadata_CT_recent = metadata_CT_recent[-c(which(metadata_CT_recent$`Collection.date` == max(metadata_CT_recent$`Collection.date`))),] # omit latest 3 days
metadata_CT_recent = metadata_CT_recent[-c(which(metadata_CT_recent$`Collection.date` == max(metadata_CT_recent$`Collection.date`))),] # omit latest 3 days
metadata_CT_recent = metadata_CT_recent[-c(which(metadata_CT_recent$`Collection.date` == max(metadata_CT_recent$`Collection.date`))),] # omit latest 3 days

lineages_daily_draft <- metadata_CT_recent %>% complete(`Collection.date`, nesting(`pango_lineage`), fill = list(CopyDate = 0))
lineages_daily_draft$`CopyDate` = ifelse(lineages_daily_draft$`CopyDate` == 0,0,1)

lineages_sum <- metadata_CT_recent %>% group_by(`pango_lineage`) %>% count(pango_lineage)  %>%
  filter(n > 5)
rownames(lineages_sum) = lineages_sum$pango_lineage

lineages_daily_draft$`pango_lineage` = ifelse(lineages_daily_draft$`pango_lineage`%in% rownames(lineages_sum),lineages_daily_draft$`pango_lineage`,"Other")

lineages_daily <- lineages_daily_draft %>% group_by(`Collection.date`,`pango_lineage`) %>% 
  summarise_at(vars(`CopyDate`),list(n = sum)) %>%
  mutate(freq = round(100*n/sum(n),2)) 
#  filter(`pango_lineage`%in% rownames(lineages_sum))


# plot
lineage_col = c("#000000","#f90808","#f96708","#f9a508",
                "#f1d326","#ffff0d","#daff0d","#ebe7e7",
                "#ffe0d8","#ffe7d8","#fff2d8","#fff9d8","#d3d7fe","#f6fed3","#e3fed3",
                "#fffed8","#f9ffd8","#edffd8","#d8ffdc", "#f5e5fa","#faf5e5","#faede5","#d3fcfe",
                "#d8fff4","#d8fdff","#d8f3ff","#d8e2ff","#d0cdff","#f0cdff","#fce2fc","#eaede9",
                "#f1f1e4","#f1fffd","#f1faff","#f5f1ff","#fff1f1","#fffbf1")
lineages_daily_sorted = lineages_daily[order(-lineages_daily$freq),]
lineages_names = unique(lineages_daily_sorted$`pango_lineage`)
names(lineage_col) <- lineages_names
rm(x,col,lineages_names,lineages_daily_draft)

lineage_plot = ggplot(lineages_daily,aes(x=`Collection.date`,y=freq, group=pango_lineage, color = pango_lineage)) +
  geom_line() +
  labs(
    title = "Daily lineage frequencies in Connecticut (last 3 months)",
    y = "Lineage frequency (%)",
    x = "Time (days)",
    subtitle = "% for lattest week shown subject to change"
  )+
  theme_light() + 
  scale_color_manual(values=lineage_col) +
 # scale_x_discrete(breaks = c("2022-04-01","2022-05-01","2022-06-01")) +
  #scale_color_brewer(type = "seq", palette = "Spectral") +
  theme(
    plot.title = element_text(face = "bold", size = 12) ,
      axis.text.x = element_text(hjust = 1, size = 10, color = "black")
  ) 
ggplotly(lineage_plot)


# mix and match
names = read.csv("data/variant_names.csv")
# metadata_CT_all$who_variants = NA
# metadata_CT_all$who_variants[names$pango_variants == metadata_CT_all$Lineage] <- names$who_variants[names$pango_variants == metadata_CT_all$Lineage]
metadata_CT_who = merge(metadata_CT_all,names[,c("pango_lineage","who_variants")], by = "pango_lineage", all.x = TRUE)
#metadata_CT_all$who_names <- names[match(metadata_CT_all$Lineage,names$pango_variants),2]

cumulative <- metadata_CT_who %>% 
  select(`Collection.date`,`who_variants`,`pango_lineage`) %>%
  count(`who_variants`)
colnames(cumulative) = c("WHO.label","Cumulative.sequenced.cases.")
cumulative[which(is.na(cumulative$WHO.label)),1] <- "Other"

last_3_weeks <- metadata_CT_who %>% 
  select(`Collection.date`,`who_variants`,`pango_lineage`) %>%
  subset(`Collection.date` >= date_3_weeks) %>% # update date every 3 weeks
  count(`who_variants`) %>%
  mutate(freq = round(100*n/sum(n),2)) 
colnames(last_3_weeks) = c("WHO.label","Total.sequenced.from.past.3.weeks..","Percent.sequenced.from.past.3.weeks..")
last_3_weeks[which(is.na(last_3_weeks$WHO.label)),1] <- "Other"

metadata_CT_who[which(is.na(metadata_CT_who$who_variants)),length(colnames(metadata_CT_who))] <- 'Other'
metadata_CT_who[which(metadata_CT_who$who_variants == "Eta"),length(colnames(metadata_CT_who))] <- 'Other'
metadata_CT_who[which(metadata_CT_who$who_variants == "Kappa"),length(colnames(metadata_CT_who))] <- 'Other'
metadata_CT_who[which(metadata_CT_who$who_variants == "Zeta"),length(colnames(metadata_CT_who))] <- 'Other'
metadata_CT_who[which(metadata_CT_who$who_variants == "Lambda"),length(colnames(metadata_CT_who))] <- 'Other'

who_weekly <- metadata_CT_who %>% group_by(`week`) %>% count(who_variants) %>%
  mutate(freq = round(n/sum(n),2)) 
who_weekly = who_weekly[c(which(who_weekly$week > '2020-12-31')),] # omit first few days
who_weekly = who_weekly[-c(which(who_weekly$week == max(who_weekly$week))),] # omit latest week
who_weekly$`Variant names` = who_weekly$who_variants

who_colors = c("Other" = "#c4c1c0",
               "Alpha" = "#1160f6",
               "Delta" = "#11adf6",
               "Beta" = "#c81ac8",
               "Gamma" = "#c81ac8",
               "Epsilon" = "#13d84f",
               "Iota" = "#1ac86e",
               "Lambda" = "#bfe59f",
               "Mu" = "#9e5220",
               "Omicron (BA.1)" = "#f1fa30",
               "Omicron (BA.3)" = "#e58b00",
               "Omicron (BA.4)" = "#f47110", 
               "Omicron (BA.5)" = "#d04805"
               )
who_plot = ggplot(who_weekly,aes(x=`week`,y=freq, group = who_variants, color = who_variants)) +
  geom_line(size = 0.6) +
  labs(
    title = "Weekly variant frequencies in Connecticut (since January 2021)",
    y = "Variant frequency"
  ) + 
 scale_color_manual(values = c(who_colors)) +
  theme(
    plot.title = element_text(face = "bold", size = 12)
  ) +
  theme_light()
ggplotly(who_plot)


# attempts at Table 1
table1_old = read.csv("data/Table1.csv")
table1_old = table1_old[,-c(1)]

tab = merge(table1_old,cumulative[,c("WHO.label","Cumulative.sequenced.cases.")], by = "WHO.label", all.x = TRUE)

table1_new = merge(tab,last_3_weeks[,c("WHO.label","Total.sequenced.from.past.3.weeks..","Percent.sequenced.from.past.3.weeks..")], by = "WHO.label", all.x = TRUE)
table1_new[is.na(table1_new)] <- 0
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


###### %seq attempt ##########
cases_path = list.files(path = "data/",pattern = "time_series*",full.names = TRUE,recursive=TRUE,include.dirs=TRUE)
cases = fread(cases_path)
cases_CT_t = cases %>%
  filter(Province_State == "Connecticut") %>%
  select(-c(UID,iso2,iso3,code3,FIPS,Admin2,Country_Region,Lat,Long_,Combined_Key,Province_State)) 
cases_CT = tibble(apply(cases_CT_t,2,sum),mdy(colnames(cases_CT_t)))
colnames(cases_CT) = c("case_cumulative","date")

cases_CT_filt = cases_CT %>% 
  mutate(CopyDate = `date`) %>%
  mutate(case_count = case_cumulative - lag(case_cumulative)) %>%
  group_by(`week` = cut(`CopyDate`, "week")) %>%
  select(-c(`CopyDate`,`case_cumulative`)) 
cases_CT_filt.dt = data.table(cases_CT_filt)
cases_CT_filt.dt = cases_CT_filt.dt[,list(case_count=sum(case_count)), by='week']
cases_CT_filt.dt = cases_CT_filt.dt[-c(1:which(cases_CT_filt.dt$week == "2020-12-28")),]
cases_CT_filt.dt$week = as.Date(cases_CT_filt.dt$week)

seq_CT = metadata_CT_all %>%
  group_by(`week`,`Collection.date`) %>%
  count()
seq_CT_filt.dt = data.table(seq_CT)
seq_CT_filt.dt = seq_CT_filt.dt[,list(seq_count=sum(n)), by='week']
seq_CT_filt.dt = seq_CT_filt.dt[-c(1:which(seq_CT_filt.dt$week == "2020-12-28")),]

subsampler = merge(cases_CT_filt.dt,seq_CT_filt.dt, by = 'week',all.x = TRUE)
subsampler[is.na(subsampler)] <- 0
subsampler = subsampler %>%
  mutate(percent = round(seq_count * 100/ case_count,2)) %>%
  select(-c(seq_count))

rm(cases,cases_CT,cases_CT_filt,seq_CT,cases_path)

ylim.prim <- c(0, 70000)   
ylim.sec <- c(0, 100)    
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

subsampler_plot <- ggplot(subsampler, aes(x=`week`)) +
  geom_col(aes(y=case_count), color = "darkblue", fill = "#e0f4ff") + 
  geom_line(aes(y=percent*b+a), color = "red", size = 1.5) +
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
subsampler_plot


########Ct Values by Lineages#############
glab = read.csv("data/GLab_SC2_sequencing_data - Sample metadata.csv")
######Assign BA.1 vs BA.2 ######
glab_filt <- glab %>% mutate(Ctvals = Lineage) 
glab_filt$Ctvals <- substr(glab_filt$Ctvals, 1,4)
#include only samples with lineage assignments
Ct_values <- glab_filt %>% dplyr::filter(`Ctvals` %in% c("BA.1", "BA.2","BA.4","BA.5"))

#######Plot#########
intercept_char <- Ct_values %>% dplyr::filter(`Ctvals` == c("BA.1"))
intercept <- as.numeric(intercept_char$`Yale.N1.FAM.`)
intercept <- mean(intercept,na.rm = TRUE)

Ct_graph <- ggplot(data = Ct_values, aes(x= `Ctvals`, y = as.numeric(`Yale.N1.FAM.`)))+
  geom_hline(yintercept = intercept,colour = "lightgrey")+
  geom_boxplot(width=0.4,outlier.shape = NA,colour = "#666666")+
  geom_jitter(aes(col =`Ctvals`), alpha = 0.2,size = 0.2, stroke = 2,shape = 21, width = 0.15)+
  labs(title="CR cycle threshold (CT) values by lineage",x=NULL, y = "Ct (N)")+
  scale_color_brewer(palette="Dark2")+
  coord_fixed(ratio = 0.06)+
  scale_y_reverse(breaks = seq(10, 40, by = 5))+
  theme_classic()+theme(legend.position="none",axis.text = element_text(size = 10, color = "black",face ="bold"), 
                        axis.title=element_text(size=12, face ="bold"),title=element_text(size=14, face ="bold"))
ggplotly(Ct_graph) 


####### Mapping attempt #########
ct_map = sf::st_read("data/tl_2019_09_cousub/tl_2019_09_cousub.shp")
m <- leaflet() %>%
  addPolygons(data=ct_map, col="black", weight = 1, layerId = ~id, label = ct_map$`NAME`, 
              highlight = highlightOptions(color = "blue",weight = 2, bringToFront = F, opacity = 0.7))
