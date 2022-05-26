############################################################################
#Using the HHEAR Study #1977 Public Microbiome Data
#for an application to the Human Microbiome Project simulation study
#using WQSRS on Microbiome Data
#By Moira Bixby
#Written 2/4/22
#Updated 05/25/22
############################################################################
library(caret)
library(dplyr)
library(tidyverse)
library(magrittr)
library(readxl)
library(MASS)
library(sjPlot)
library("Hmisc")
library(gt)
library(gtsummary)
library(ggpubr)
library(tidyverse)
library(expss)
library(gWQS)
library(naniar)

##this filepath will need to be updated
setwd("C:\\Users\\bixbym01\\OneDrive - The Mount Sinai Hospital\\HMP Microbiome\\Michels Microbiome\\")

###the asvs file will need to be replaced with the the processed dataset that we create
asvs <- read.csv("Final Processed Data\\micheals_HMP_asv_table.csv") #already 10% detect
epi <- read.csv("1977_EPI_DATA.csv")


x <- colnames(asvs[,2:501]) 

asvs2 <- asvs
#relative abundance
asvs3 <-t(apply(asvs2[2:501],1, function(x) x/sum(x)))
# tmp1<-apply(asvs2[2:501],1, function(x) x/sum(x))
# colSums(tmp1[,1:4])
asvs3 <- asvs3 %>% data.frame()


asvs3[asvs3 == 0] <- NA

#going to make this <= because in some cases, the median may be 1--if so, then all values will be the factor 2
for (i in 1:ncol(asvs3)) {
  m <- median(asvs3[,i], na.rm=T)
  asvs3[,i] <- ifelse(asvs3[,i] > m, 2, 1)
}

asvs4 <- cbind(asvs2[,1], asvs3)

asvs4[is.na(asvs4)] <- 0

tax <- colnames(asvs3)


#merge the epi data with the lab data


asvs <- names(asvs3) #128 

###
asvs4  <- asvs4 %>% rename(PID = `asvs2[, 1]`)
epi <- epi %>% rename(PID = pid)
stool <- merge(epi, asvs4, by="PID")


stool <- stool %>% mutate(
  age_stool =as.numeric(age_stool),
  lactdurdays =as.numeric(lactdurdays),
  # antro_baz_depo=as.factor(antro_baz_depo),
  birth_mode=factor(birth_mode,
                    levels = c(0,1),
                    labels = c("C-section", "Vaginal")), 
  patb_5_depo=factor(patb_5_depo,
                     levels = c(0,1),
                     labels = c("No", "Yes")),
  # bmi=factor(antro_baz_depo,
  #            levels=c(-2, -1, 0, 1, 2, 3),
  #            labels=c("Normal", "Normal", "Normal", "Normal",
  #                     "Overweight/obese", "Overweight/obese")),
  # bmi=relevel(bmi, ref="Normal"),
  bmi = case_when(
    antro_baz_depo <= 1  ~ 0,
    antro_baz_depo >1 ~ 1
  ),
  meducat=na_if(meducat, 2),
  meducat=factor(meducat,
                 levels = c(0,1),
                 labels = c("High school or less", "More than high school")), 
                 meducat= relevel(meducat, ref="High school or less")) 

# patb_5_depo=relevel(patb_5_depo="No"), 
#birth_mode=relevel(birth_mode="Vaginal"))



stool = apply_labels(stool,
                     age_stool = "Age (years)",
                     bmi = "BMI category",
                     birth_mode = "Birth mode (vaginal versus c-section)",
                     lactdurdays = "Number of days breastfed as infant",
                     patb_5_depo="Were antibiotics used in the past 6 months (yes vs no)", 
                     meducat="Maternal education")

stool %<>% filter(!is.na(meducat), !is.na(patb_5_depo))


# 
# #for demographics table
# tmp1 <- stool
# # tmp2 <- tmp1 %>% filter(!is.na(foodins))
# table1 <- tmp1 %>% select(age_stool, lactdurdays, bmi, birth_mode,  patb_5_depo, meducat) %>%  filter(!is.na(patb_5_depo)) %>% tbl_summary(by="bmi",
#                                                                                                                                  type = list(c(age_stool, lactdurdays) ~ "continuous"),
#                                                                                                                                  statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
#                                                                                                                                  digits = list(all_continuous() ~1, all_categorical() ~ 0),  missing= "no")   %>% add_p(list(all_continuous() ~ "t.test"), spvalue_fun = function(x) style_pvalue(x, digits = 3))  %>% add_n() %>% as_hux_table()
# 
# 
# table1 %>% as_gt() %>%             # convert to gt table
#   gt::gtsave(             # save table as image
#     filename = "HMP_Michels_Table1_052422.rtf"
#   )
# cov <- glm((bmi) ~ meducat + birth_mode + age_stool + lactdurdays + patb_5_depo, data=stool, family=binomial)
# tab_model(cov)
##referent group is normal weight for 0




ptm <- proc.time()

wqsrs_pos1 <- gwqs(bmi ~ wqs +  meducat + birth_mode + age_stool + lactdurdays + patb_5_depo, data=stool, mix_name=asvs, b=2000, rs = T, n_vars = 22, b1_pos = TRUE, b1_constr = FALSE, q = NULL, validation = 0.6, family=binomial(), seed = 1234) 
proc.time() - ptm
gwqs_summary_tab(wqsrs_pos1)
sum(wqsrs_pos1$bres$b1 >0 ) #772
sum(wqsrs_pos1$bres$b1 <0 ) #1228

##test with positive constraints


##negative constraints

wqsrs_neg1 <- gwqs(bmi ~ wqs +  meducat + birth_mode + age_stool + lactdurdays + patb_5_depo, data=stool, mix_name=asvs, b=2000, rs = TRUE, n_vars = 22, b1_pos = FALSE, b1_constr = TRUE, q = NULL, validation = 0.6, family = binomial(), seed = 1234,
                   plan_strategy = "multisession") 
summary(wqsrs_neg1) 
gwqs_summary_tab(wqsrs_neg1)

################################################
#####


#create the pre-set holdouts for the repeated holdout analysis; this chooses the subjects to go in the 60% validation dataset (and therefore also the training)--based off of the interaction variable of all fo the categorical variables. This is done so in the repeated holdouts, we avoid the issue of having subjects in only one category of a variable, wherein the analysis does not work/run.  
group_list <- createDataPartition(stool$bmi, p = 0.6, list = T, times = 30) #this is for # of repeated holdouts
##put the group_list into the rh= for gwqsrh


ptm <- proc.time()

rhrs_neg2 <- gwqsrh(bmi ~ wqs + meducat + as.factor(birth_mode) + age_stool + lactdurdays + patb_5_depo, data=stool, mix_name=asvs, b=2000, rs = TRUE, n_vars = 22, b1_pos = FALSE, b1_constr = TRUE, q = NULL, validation = 0.6, rh= group_list, family = binomial(), seed = 564)
proc.time() - ptm

write.csv(rhrs_neg2$wmat, "H:\\Michels_BMI_Negcont_WtsRSRH30_051922.csv")
summary(rhrs_neg2) 
gwqs_summary_tab(rhrs_neg2, sumtype = "norm")


sum(rhrs_neg2$coefmat[,2] >0 ) #2
sum(rhrs_neg2$coefmat[,2] < 0 ) #28