###############
# Data selection and extraction
###############



install.packages("BiocManager") #one-time installation
BiocManager::install("HMP16SData") #one-time installation
library(HMP16SData)


V35_stool <-
  V35() %>%
  subset(select = HMP_BODY_SUBSITE == "Stool") 

V35_stool_phyloseq <- as_phyloseq(V35_stool)



knitr::opts_chunk$set(message = FALSE)

library(HMP16SData)
library(phyloseq)
library(magrittr)
library(ggplot2)
library(tibble)
library(dplyr)

library(dendextend)
library(circlize)
library(curatedMetagenomicData)
library(gridExtra)
library(cowplot)
library(readr)
library(haven)


###################
#This should export the subset data into the file types I need

V35_participant_data <-
  V35() %>% 
  colData() %>%
  as.data.frame() %>%
  rownames_to_column(var = "PSN")
  
V35_OTU_counts <-
    V35() %>% 
    assay() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "PSN")
    
V35_data <-
  merge.data.frame(V35_participant_data, V35_OTU_counts, by = "PSN")

V35_dictionary <-
  V35() %>%
  rowData() %>%
  t.data.frame() %>%
  as.data.frame()

colnames(V35_data) <-
  colnames(V35_data) %>%
  gsub(pattern = "\\.", replacement = "_", x = .)

colnames(V35_dictionary) <-
  colnames(V35_dictionary) %>%
  gsub(pattern = "\\.", replacement = "_", x = .)


#Saving as external files
library(foreign)
write.foreign(V35_data, "hmpdata.txt", "hmpdata.sas",   package="SAS")
write.foreign(V35_dictionary, "hmptax.txt", "hmptax.sas",   package="SAS")




##################
##############
#HMP Data Introduction
#######################

library(readxl)
library(data.table)
library(magrittr)
library(plyr)
library(dplyr)


####Read in HMP data
hmp <- read_excel("hmpdata.txt")

setDT(hmp)

##get total numbers of observations. original number of observations was 4743
tmp = hmp[HMP_BODY_SUBSITE=="Stool"]
#219 stool samples total

test <- count(tmp, vars="RSID", wt_var="VISITNO")
m1 <- tmp %>% group_by(RSID) %>% filter(row_number()==1)
#in total, there are now 210 FIRST stool samples (which can be from visit 1 or visit 2)
#16384 varialbes, minus the demographics--so now there were 
# > table(m1$VISITNO)
# 
# 1   2 
# 191  19 

# % DETECT--we want 10%+ detect

x <- colnames(tmp)
m1[m1 == 0] <- NA

#this checks % non-detect in each column. we are going for anything 90% or less non-detect (or 10%+ detect)
m2 <- m1[lapply( m1[1:210,], function(x) sum(is.na(x)) / length(x) ) < 0.9 ]
m2[is.na(m2)]=0



#to double check that the % detect is working, and it is checking %detect at 10%+
test_detect <- m1[lapply( m1[1:210,], function(x) sum(!is.na(x)) / length(x) ) > 0.1 ]



m2 %>%
  select(everything()) %>%  # replace to your needs
  summarise_all(funs(sum(is.na(.))))
##maybe test this, above, with the remaining 10% to see if they actually are >90% nondetect

##find just the character vars to figure out exactly how many OTUs there are

tmp %>% 
  select_if(is.character) %>% 
  head(2)
#5, below, and I think also PSID and RSID so 7 total--so 16384-7=16377
#  SEX RUN_CENTER HMP_BODY_SITE HMP_BODY_SUBSITE SRS_SAMPLE_ID

#export dataset

write.csv(m2, file=, row.names=FALSE)

