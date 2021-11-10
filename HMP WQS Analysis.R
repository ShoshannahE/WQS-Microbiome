# Analysis of HMP data and binary Y using simulated data - June 2021
#rm(list = ls())
# install the packages needed for the analysis
install.packages("gWQS") 
# import the packages just installed
library(gWQS) 
install.packages("devtools") ## if not already installed
library(devtools)
install_github("atahk/pscl")
install.packages("Lasso")
library(Lasso)
# define the path

directory_path = "/Users/gennic01/Desktop/RESEARCH/Shoshannah Eggers/HMP study/"

# import the dataset
#dataset = read.csv(paste0(directory_path, "HMP_10%Detect_RELATIVEAbundance_10exps_RANKS_0_to_2.csv"))
dataset = read.csv(paste0(directory_path, "HMP_10%Detect_relabun_0to3.csv"))
dim(dataset)
names(dataset)

dataset1 = subset(dataset, sim_exp==1)
dim(dataset1)
#make mixture (defined by the column names of the columns of the OTUs from 'dataset')
mixture <- names(dataset[,27:894])
mixture

strong = c("q425", "q450")
strong
medium = c("q545", "q616", "q845", "q58", "q63", "q133", "q141", "q235")
medium
weak = c("q378", "q584", "q737", "q6", "q127", "q202", "q239", "q262", "q359", "q391")
weak
truth = c(strong, medium, weak)
truth

summary(dataset1$female)
summary(dataset1$sim_smk)
summary(dataset1$rv0)
################################################################
# correlation among the chemicals
install.packages("corrplot")
library("corrplot")

dataset_ = dataset1[complete.cases(dataset1),]
forcorr = data.frame(dataset_[,truth])
head(forcorr)
corrplot(cor(forcorr), method="circle", sig.level = c(.001, .01, .05), type="upper")

corrplot(cor(forcorr, method="spearman"), method="shade", type="upper")


## bootstrap WQS using only the variables in truth
results1bs  = gwqs(sim_smk ~ wqs +female , 
                   mix_name = truth, 
                   data = dataset1, q = NULL, 
                   validation = 0.6,  b = 100,  
                   b1_pos = TRUE, b1_constr = TRUE, family = "binomial", 
                   seed = 123456, plan_strategy="multisession", signal="abst",
                   control = list(trace = FALSE, maxit = 3000, reltol = 1e-14))

gwqs_summary_tab(results1bs)  
gwqs_barplot(results1bs) 
gwqs_scatterplot(results1bs)  
gwqs_weights_tab(results1bs) 

## random subset WQS of all OTUs in mixture, with positive constraint, calling 1000 
## random subsets of size 30
results1 = gwqs(sim_smk ~ wqs +female , 
                mix_name = mixture, 
                data = dataset1, q = NULL, 
                validation = 0.6,  b = 1000, rs=TRUE, n_var=30,
                b1_pos = TRUE, b1_constr = TRUE, family = "binomial", 
                seed = 123456, plan_strategy="multisession",
                control = list(trace = FALSE, maxit = 3000, reltol = 1e-14))

gwqs_summary_tab(results1)  
gwqs_barplot(results1) 
gwqs_scatterplot(results1)  
gwqs_weights_tab(results1) 

########### repeated holdout ####################################
# 30 repeated holdouts of the random subset WQS using default signal function (i.e., t2)

results1rh = gwqsrh(sim_smk ~ wqs +female , 
                    mix_name = mixture, 
                    data = dataset1, q = NULL, 
                    validation = 0.6,  b = 1000, rs=TRUE, n_var=30,
                    b1_pos = TRUE, b1_constr = TRUE, family = "binomial", 
                    seed = 123456, plan_strategy="multisession",
                    control = list(trace = FALSE, maxit = 3000, reltol = 1e-14),
                    rh=30)

gwqs_summary_tab(results1rh)
gwqsrh_boxplot(results1rh)
gwqs_weights_tab(results1rh)
summary(results1rh$fit)
results1rh$final_weights

### checking the negative direction
results1rhneg = gwqsrh(sim_smk ~ wqs +female , 
                       mix_name = mixture, 
                       data = dataset1, q = NULL, 
                       validation = 0.6,  b = 1000, rs=TRUE, n_var=30,
                       b1_pos = FALSE, b1_constr = TRUE, family = "binomial", 
                       seed = 123456, plan_strategy="multisession",
                       control = list(trace = FALSE, maxit = 3000, reltol = 1e-14),
                       rh=30)

gwqs_summary_tab(results1rhneg)
gwqsrh_boxplot(results1rhneg)
gwqs_weights_tab(results1rhneg)
summary(results1rhneg$fit)


## compare to using exp(t) as the signal function
results1rh_sigexp = gwqsrh(sim_smk ~ wqs +female , 
                           mix_name = mixture, 
                           data = dataset1, q = NULL, 
                           validation = 0.6,  b = 1000, rs=TRUE, n_var=30,
                           b1_pos = TRUE, b1_constr = TRUE, family = "binomial", 
                           seed = 123456, plan_strategy="multisession",
                           control = list(trace = FALSE, maxit = 3000, reltol = 1e-14),
                           rh=30, signal="expt")

gwqs_summary_tab(results1rh_sigexp)
gwqsrh_boxplot(results1rh_sigexp)
gwqs_weights_tab(results1rh_sigexp)
summary(results1rh_sigexp$fit)
results1rh_sigexp$final_weights

write.csv(results1rh_sigexp$final_weights, paste0(directory_path, "WQSwtsrh_sigexp.csv"))
write.csv(results1rh_sigexp$coefmat, paste0(directory_path, "coefmatrh_sigexp.csv"))
write.csv(results1rh_sigexp$wmat, paste0(directory_path, "wmatrh_sigexp.csv"))

############################################################
## compare to null case
############################################################
results1rh_RV0sigexp = gwqsrh(rv0 ~ wqs +female , 
                              mix_name = mixture, 
                              data = dataset1, q = NULL, 
                              validation = 0.6,  b = 1000, rs=TRUE, n_var=30,
                              b1_pos = TRUE, b1_constr = TRUE, family = "binomial", 
                              seed = 123456, plan_strategy="multisession",
                              control = list(trace = FALSE, maxit = 3000, reltol = 1e-14),
                              rh=30, signal="expt")

gwqs_summary_tab(results1rh_RV0sigexp)
gwqsrh_boxplot(results1rh_RV0sigexp)
gwqs_weights_tab(results1rh_RV0sigexp)
summary(results1rh_RV0sigexp$fit) 

write.csv(results1rh_RV0sigexp$coefmat, paste0(directory_path, "coefmatRV0rh_sigexp.csv"))

##############################################################
## compare to using abs(t) for signal function
results1rh_sigabst = gwqsrh(sim_smk ~ wqs +female , 
                            mix_name = mixture, 
                            data = dataset1, q = NULL, 
                            validation = 0.6,  b = 1000, rs=TRUE, n_var=30,
                            b1_pos = TRUE, b1_constr = TRUE, family = "binomial", 
                            seed = 123456, plan_strategy="multisession",
                            control = list(trace = FALSE, maxit = 3000, reltol = 1e-14),
                            rh=30, signal="abst")

gwqs_summary_tab(results1rh_sigabst)
gwqsrh_boxplot(results1rh_sigabst)
gwqs_weights_tab(results1rh_sigabst)
summary(results1rh_sigabst$fit)
results1rh_sigabst$final_weights

write.csv(results1rh_sigabst$final_weights, paste0(directory_path, "WQSwtsrh_sigabst.csv"))



dim(results1rh$wmat)  

#install.packages("dplyr")
#library("dplyr")
rhval=30
## ncol should be the number of rh times 2 (sens and spec) + 3 (cutpoint and two avgs)
# nrow should be the number of cutpoints considered
minval=0.0005
maxval=.002
incval = (maxval-minval)/39  # there are 40 cutpoints

dataforsensspec = results1rh_sigexp$wmat
#dataforsensspec = results1rh_sigabst$wmat

results = matrix(nrow = 40, ncol=63)
#results[,1]=seq(.0001, .0020,by=.00005)
results[,1]=seq(minval, maxval, by=incval)
results[,1]

# truth, strong, medium, weak
for(i in 1:40){
  cutpoint = results[i,1]
  dataforcalc =dataforsensspec[,truth]>cutpoint
  rowSums(dataforcalc)
  Sensitivity = rowSums(dataforcalc)/20 # compute row sums from truth columns
  results[i,2:(rhval+1)]= Sensitivity
  
  dataforcalc2 = dataforsensspec[, ! mixture %in% truth]<cutpoint    
  Specificity = rowSums(dataforcalc2)/848 # compute row sums without truth columns
  results[i,(rhval+2):(2*rhval+1)] = Specificity
}

SensRes = rowMeans(results[,2:(rhval+1)]) 
results[,(rhval*2+2)] = SensRes
SpecRes = rowMeans(results[, (rhval+2):(2*rhval+1)])
results[,(rhval*2+3)] = SpecRes

resultsdf = as.data.frame(results)
resultsdf$cutpoint= results[,1]
resultsdf$AvgSens = results[,62]
resultsdf$AvgSpec = results[,63]
head(resultsdf)

resultsdf[,c("cutpoint", "AvgSens", "AvgSpec")]

write.csv(resultsdf, paste0(directory_path, "cutpoint_sensspec_sigexp.csv"))
#write.csv(resultsdf, paste0(directory_path, "cutpoint_sensspec_sigabst.csv"))


######################################################

library(broom)
library(tidyr)
library(dplyr)
library(plyr)
# Break up dataset by sim_exp, then fit the specified model to each piece and
# return a list
models <- dlply(dataset, "sim_exp", function(df) 
  results = gwqs(sim_smk ~ wqs +female , 
                 mix_name = mixture, 
                 data = df, q = NULL, 
                 validation = 0.6,  b = 10, rs=TRUE, n_var=30,
                 b1_pos = TRUE, b1_constr = TRUE, family = "binomial", 
                 seed = 123456, plan_strategy="multisession",
                 control = list(trace = FALSE, maxit = 3000, reltol = 1e-14)))




# Apply coef to each model and return a data frame
ldply(models, coef, confint.lm=TRUE)

# Print the summary of each model
l_ply(models, summary, .print = TRUE)
l_ply(models, gwqs_summary_tab, .print=TRUE)
l_ply(models, final_weights)

summary(models[1])


simul_wqs = for(i in 1:2){
  dataseti = subset(dataset, sim_exp==i)
  #positive constraints -  random subset WQS
  results1[i]  = gwqs(sim_smk ~ wqs +female , 
                      mix_name = mixture, 
                      data = dataseti, q = NULL, 
                      validation = 0.6,  b = 20, rs=TRUE, n_var=30,
                      b1_pos = TRUE, b1_constr = TRUE, family = "binomial", 
                      seed = 123456, plan_strategy="multisession",
                      control = list(trace = FALSE, maxit = 3000, reltol = 1e-14))
}
results1[2]
simul_wqs



summary(results1$fit) 

for(i in 1:2){gwqs_summary_tab(results1[1])}  

gwqs_barplot(results1) 
gwqs_scatterplot(results1)  
gwqs_weights_tab(results1) 

results1$final_weights[truth,]
sum(results1$final_weights[truth,2]>(1/868))
sensitivity =  sum(results1$final_weights[truth,2]>(1/868))/20
sensitivity

sum(results1$final_weights[strong,2]>(1/868))
sum(results1$final_weights[medium,2]>(1/868))
sum(results1$final_weights[weak,2]>(1/868))


specificity = (sum(results1$final_weights[mixture,2]<(1/868))-5)/848
specificity
}

summary(simul_wqs)




