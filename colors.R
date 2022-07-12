#PACKAGE INSTALL ####
packages = c("ggalluvial", "readr", "dplyr","tidyverse","ggplot2","randomcoloR","cowplot","openxlsx")

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

########Weekly Lineages
Lineages = read.csv("data/lineages.csv")
Lineages_Sum = read.csv("data/lineages_freq.csv")
dat.sum = read.csv("data/dat_sum.csv")

########Assigns random colors
#build a pallete mapping using 'Seq' column's value in all available dataframes
pal <- c(randomColor(count=length(unique(Lineages_Sum$`Lineage`))))
pal_seq_mapping <- data.frame(variant=unique(Lineages_Sum$`Lineage`), color=pal)
BA1 = "BA.1"
BA1.1 = "BA.1.1"
BA2 = "BA.2"
BA2.12.1 = "BA.2.12.1"
BA2.8 = "BA.2.8"
BA3 = "BA.3"
BA4 = "BA.4"
BA5 = "BA.5"
XE = "XE"
if(grepl("BA.1", Lineages_Sum$`Lineage`, fixed = TRUE) == TRUE){
  pal_seq_mapping$color = "#1b55b1"
}else {
  pal_seq_mapping$color = "#851bb1"
}

grepl(BA1.1, Lineages_Sum$`Lineage`, fixed = TRUE)

########Plot
LineagePlot <- ggplot(Lineages_Sum, aes(x = as.Date(`week`), y = `freq`, alluvium = `Lineage`)) + 
  geom_alluvium(aes(col = `Lineage`, fill=`Lineage`), alpha = 0.2,show.legend=FALSE) +
  geom_col(aes(fill=`Lineage`), width = 3,show.legend=FALSE) +
  theme_classic() + 
  scale_x_date(breaks = "1 week") + 
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black",face ="bold"))+
  ggtitle("Lineages per week") + ylab("% of all samples") +
  coord_fixed(ratio = 0.2)+
  scale_fill_manual(values=c(pal_seq_mapping[match(Lineages_Sum$`Lineage`, pal_seq_mapping$variant),"color"])) +
  scale_color_manual(values=c(pal_seq_mapping[match(Lineages_Sum$`Lineage`, pal_seq_mapping$variant),"color"]))
LineagePlot
