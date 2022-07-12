#PACKAGE INSTALL ####
packages = c("dplyr","ggplot2","RColorBrewer","shiny","shinythemes","plotly","DT")
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

## Read metadata from gisaid
path <- getwd()
filename <- "gisaid_hcov-19.tsv"
#metadata <- read_tsv(file = paste0(path, "/data/", filename))
#metadata <- read_tsv(file = "metadata_nextstrain.tsv")
#metadata_CT_all = metadata[metadata$`division` == "Connecticut",]

metadata_CT_all = read.csv("data/metadata_CT_all.csv")
metadata_CT_all = metadata_CT_all[,-c(1)]

metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "BA.2.12.1 (marker override: BA.2.12 + Spike_L452Q => BA.2.12.1)"] = "BA.2.12.1"
metadata_CT_recent = metadata_CT_all[metadata_CT_all$`date` >= "2022-03-22",] # filter by last 3 months
metadata_CT_recent =  metadata_CT_recent[!is.na(metadata_CT_recent$`pango_lineage`),]

metadata_CT_all$`date` <- as.Date(metadata_CT_all$`date`)
metadata_CT_all <- metadata_CT_all %>% mutate(CopyDate = `date`)
metadata_CT_all <- metadata_CT_all %>% group_by(`week` = cut(`CopyDate`, "week"))
metadata_CT_all$week <- as.Date(metadata_CT_all$week)
#delete data without week assigned
metadata_CT_all <- metadata_CT_all[!is.na(metadata_CT_all$`week`),]

#write.csv(metadata_CT_all,"data/metadata_CT_all.csv")
# fix weird BA.2.12.1 lineage name call

# Lineages daily frequency
length(unique(metadata_CT_recent$pango_lineage))

lineages_daily <- metadata_CT_recent %>% group_by(`date`) %>% count(pango_lineage) %>%
  mutate(freq = round(n/sum(n),2)) %>%
  subset(freq > 0.05)

# set cut off for date
date = aggregate(lineages_daily$n, by = list(lineages_daily$`date`), FUN = sum)

# plot
lineage_plot = ggplot(lineages_daily,aes(x=`date`,y=freq, group=pango_lineage, color = pango_lineage)) +
  geom_line() +
  labs(
    title = "Lineage frequencies in Connecticut in the last 3 months",
    y = "Lineage frequency"
  ) + 
  scale_color_brewer(type = "seq", palette = "Spectral") +
  theme(
    plot.title = element_text(face = "bold", size = 12)
  ) +
  theme_light()
ggplotly(lineage_plot)


# mix and match
names = read.csv("data/variant_names.csv")
# metadata_CT_all$who_variants = NA
# metadata_CT_all$who_variants[names$pango_variants == metadata_CT_all$pango_lineage] <- names$who_variants[names$pango_variants == metadata_CT_all$pango_lineage]
metadata_CT_who = merge(metadata_CT_all,names[,c("pango_lineage","who_variants")], by = "pango_lineage", all.x = TRUE)
#metadata_CT_all$who_names <- names[match(metadata_CT_all$pango_lineage,names$pango_variants),2]
metadata_CT_who[which(is.na(metadata_CT_who$who_variants)),31] <- 'Other'
metadata_CT_who[which(metadata_CT_who$who_variants == "Eta"),31] <- 'Other'
metadata_CT_who[which(metadata_CT_who$who_variants == "Kappa"),31] <- 'Other'
metadata_CT_who[which(metadata_CT_who$who_variants == "Zeta"),31] <- 'Other'
metadata_CT_who[which(metadata_CT_who$who_variants == "Lambda"),31] <- 'Other'

who_weekly <- metadata_CT_who %>% group_by(`week`) %>% count(who_variants) %>%
  mutate(freq = round(n/sum(n),2)) 
who_weekly = who_weekly[c(which(who_weekly$week > '2020-12-31')),] # omit first few days

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
  geom_line() +
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
cumulative <- metadata_CT_who %>% 
  select(`date`,`who_variants`,`pango_lineage`) %>%
  group_by() %>% count(`who_variants`)
colnames(cumulative) = c("WHO.label","Cumulative.sequenced.cases.")

tab = merge(table1_old,cumulative[,c("WHO.label","Cumulative.sequenced.cases.")], by = "WHO.label", all.x = TRUE)

last_3_weeks <- metadata_CT_who %>% 
  select(`date`,`who_variants`,`pango_lineage`) %>%
  subset(`date` >= "2022-06-02") %>% # update date every 3 weeks
  count(`who_variants`) %>%
  mutate(freq = round(100*n/sum(n),2)) 
colnames(last_3_weeks) = c("WHO.label","Total.sequenced.from.past.3.weeks..","Percent.sequenced.from.past.3.weeks..")

table1_new = merge(tab,last_3_weeks[,c("WHO.label","Total.sequenced.from.past.3.weeks..","Percent.sequenced.from.past.3.weeks..")], by = "WHO.label", all.x = TRUE)
table1_new[is.na(table1_new)] <- 0
table1_new$`Percent.change.from.previous.report.y` = table1_new$`Percent.sequenced.from.past.3.weeks...y` - table1_new$`Percent.sequenced.from.past.3.weeks...x`
table1_new = table1_new[,-c(4:7)]
colnames(table1_new) = c("WHO label",
                          "Pango lineage",
                         "CDC classification",
                         "Cumulative sequenced cases*",
                         "Total sequenced from past 3 weeks**",
                         "Percent sequenced from past 3 weeks**",
                         "Percent change from previous report")
write.csv(table1_new,"data/Table1_new.csv")


# attempts at subsampler
cases = read.table("data/matrix_cases_epiweeks.tsv", sep = "\t", header = TRUE)
cases_CT = t(cases[which(cases$code == "CT"),])
cases_CT_filt = cases_CT[-c(1:3)]

sampling = read.table("data/weekly_sampling_proportions.tsv", sep = "\t", header = TRUE)
sampling_CT = t(sampling[which(sampling$code == "CT"),])
sampling_CT_filt = sampling_CT[-c(1:3)]

week_CT = sort(c(unique(metadata_CT_all$`week`),"2022-06-13","2022-06-20","2022-06-27"))

subsampler_draft = data.frame(week_CT,cases_CT_filt,sampling_CT_filt)
subsampler = subsampler_draft[-c(1:41),]
subsampler$cases_CT_filt = as.numeric(subsampler$cases_CT_filt)
subsampler$sampling_CT_filt = 100*as.numeric(subsampler$sampling_CT_filt)

rm(cases_CT,cases_CT_filt,sampling,sampling_CT,sampling_CT_filt,week_CT)

ylim.prim <- c(0, 70000)   
ylim.sec <- c(0, 100)    
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

subsampler_plot <- ggplot(subsampler, aes(x=`week_CT`)) +
  geom_col(aes(y=cases_CT_filt), color = "darkblue") + 
  geom_line(aes(y=sampling_CT_filt*b+a), color = "darkred") +
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

