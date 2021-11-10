#HMP paper - Comparison Methods


#####Installing and loading packages#######
#Don't actually need all of these packages for this code, but I always install them anyway.

BiocManager::install(c("phyloseq", "DESeq2"))
install.packages('pacman')
pacman::p_load("devtools","here", "BiocManager", "MASS", "lattice", "lmtest", "plyr" , "pscl", "assertr", "tidyverse", "reshape2", "janitor", "dada2", "phyloseq", "qiime2R", "ggplot2", "DESeq2", "GUniFrac", "vegan", "LDM", "scales", "grid", "ggsci", "ggpubr", "viridis", "picante", "knitr", "gridExtra")
library(pacman)
pacman::p_loaded()

install.packages('rgr')
library(rgr)
install.packages(ape)
library(ape)
install.packages(dplyr)
library(dplyr)
install.packages(ggplot2)
library(ggplot2)
install.packages(gplots)
library(gplots)
install.packages(lme4)
library(lme4)
install.packages(phangorn)
library(phangorn)
install.packages(plotly)
library(plotly)
install.packages(tidyr)
library(tidyr)
install.packages(magrittr)
library(magrittr)
install.packages(reshape2)
library(reshape2)
install.packages(randomForest)
library(randomForest)




######Data Processing######

#loading the OTU file
OTU = read.csv("HMP_RELAbun.csv")

#loading the taxonomu file
tax = read.csv('HMP_Tax.csv')

#loading the simulated vars file
sim = read.csv('forShannah070821.csv')

#OTU file also has a bunch of meta vars, so need to make those into their own data set
# selecting meta variables
myvars <- c("PSN", "RSID", "VISITNO", "SEX", "RUN_CENTER", "HMP_BODY_SITE", "HMP_BODY_SUBSITE", "SRS_SAMPLE_ID")
meta <- OTU[myvars]

#need to filter the sim file to only include sim_exp = 1
sim_1 <- sim %>%
  
  # keep obs if sim_exp = 1
  filter(sim_exp == 1)

#now need to merge the meta file with the sim file
meta.sim <- merge(meta,sim_1,by="RSID")

#setting the row names to match in meta and otu files
row.names(OTU) = OTU$RSID
row.names(meta.sim) = meta.sim$RSID

#making a sim factor var
meta.sim$sim_fac <- factor(meta.sim$sim_smk, labels = c("No", "Yes"))
#making a random factor var
meta.sim$rv_fac <- factor(meta.sim$rv0, labels = c("No", "Yes"))

#now removing the non-otu vars from the otu file
myvars <- names(OTU) %in% c("PSN", "RSID", "VISITNO", "SEX", "RUN_CENTER", "HMP_BODY_SITE", "HMP_BODY_SUBSITE", "SRS_SAMPLE_ID")
OTU.clean <- OTU[!myvars]

#setting the row names in the tax file by the OTU number
row.names(tax) = tax$OTU

#removing all the OTUs that aren't in our OTU file
tax.clean = tax[row.names(tax) %in% colnames(OTU.clean),]

#separate the taxonomy into seperate columns
tax.clean = separate(tax.clean, taxonomy, into = c("Root","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep=";")

#removing unnecessary columns
tax.clean = tax.clean[,-which(names(tax.clean) %in% c("OTU", "Root", "Species", "Strain"))]

#setting a random number seed in case we need it later
set.seed(2020)


#
OTU.physeq = otu_table(as.matrix(OTU.clean), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(tax.clean))
meta.physeq = sample_data(meta.sim)

physeq.alpha = phyloseq(OTU.physeq, tax.physeq, meta.physeq)




########ALPHA DIVERSITY###############
#generate alpha diversity columns in phyloseq object
sample_data(physeq.alpha)$shannon.physeq <- estimate_richness(physeq.alpha, measures="Shannon")
sample_data(physeq.alpha)$InvSimpson.physeq <- estimate_richness(physeq.alpha, measures="InvSimpson")

plot_richness(physeq.alpha, measures="Shannon")
plot_richness(physeq.alpha, measures="InvSimpson")

#plotting shannon by sim_fac
plot_richness(physeq.alpha, "sim_fac", measures="Shannon")
#plotting shannon by rv_fac
plot_richness(physeq.alpha, "rv_fac", measures="Shannon")



# #calculating with vegan adds the column to meta file
meta.sim$shannon.vegan <- diversity(OTU.clean, "shannon")
meta.sim$invsimpson.vegan <- diversity(OTU.clean, index="invsimpson")

meta.sim <- as.data.frame(sample_data(physeq.alpha))

#finding the max, min, mean of alpha diversity
mean(meta.sim$shannon) 
min(meta.sim$shannon)  
max(meta.sim$shannon)  


#running a linear regression with the alpha diversity calculated in VEGAN
lmshannon = lm(shannon.vegan ~ sim_smk + SEX, data = meta.sim)
summary(lmshannon)


#running a linear regression comparison with the random variable
lmshannon = lm(shannon.vegan ~ rv0 + SEX, data = meta.sim)
summary(lmshannon)





######Beta Diversity#######

#Here we're looking at all the obs, colored by simulation variable
bray_dist <- phyloseq::distance(physeq.alpha, method="bray")
bray_ord <- ordinate(physeq.alpha, method="MDS", distance=bray_dist)
p1_bray <- plot_ordination(physeq.alpha, bray_ord, color="sim_fac") + theme(aspect.ratio=1) + 
  theme_bw()+
  scale_color_viridis(discrete=TRUE)+ 
  labs(color="Simulated Variable")+
  stat_ellipse(size=1)
p1_bray

  
#Plotting the random var
p2_bray <- plot_ordination(physeq.alpha, bray_ord, color="rv_fac") + theme(aspect.ratio=1) + 
  theme_bw()+
  scale_color_viridis(discrete=TRUE)+ 
  labs(color="Random Variable")+
  stat_ellipse(size=1)
p2_bray




#PERMANOVA
#This tells us if the beta diversity is significantly different between groups
#Calculate bray curtis distance and save as a matrix
BC.dist=vegdist(OTU.clean, distance="bray")

#Run PERMANOVA on distances by the simulated variable.
adonis2(BC.dist ~ sim_smk + SEX , data = meta.sim, by='margin', permutations = 9999)

#Run PERMANOVA on distances by the random variable.
adonis2(BC.dist ~ rv0 + SEX , data = meta.sim, by='margin', permutations = 9999)




#Running the same with Aitchison distance.
#Looks like I have to CLR transform the OTU table first, and then run the Euclidian distance.

OTU.CLR <- clr(OTU.clean)
A.dist=vegdist(OTU.clean, method="euclidean")

#Run PERMANOVA on distances by the simulated variable.
adonis2(A.dist ~ sim_smk + SEX , data = meta.sim, by='margin', permutations = 9999)

#Run PERMANOVA on distances by the random variable.
adonis2(A.dist ~ rv0 + SEX , data = meta.sim, by='margin', permutations = 9999)




#####SIMPER#####
#This tells us which OTUs contributed the most variance to the differences in beta diversity between the two comparison groups.
#running simper by the simulated var
simper(OTU.clean, meta.sim$sim_smk, permutations=999)



#running simper by the random var
simper(OTU.clean, meta.sim$rv0, permutations=999)




#####RANDOM FOREST - Sim#####

# Make a dataframe of training data with OTUs as column and samples as rows
predictors <- data.frame(otu_table(physeq.alpha))
dim(predictors)
predictors$SEX <- as.factor(meta.sim$SEX)

# Make one column for our outcome/response variable 
response <- as.factor(sample_data(physeq.alpha)$sim_smk)

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors)

# Random Forest Model
rf.model <- randomForest(response~., data = rf.data, ntree = 100)
print(rf.model)


# What variables are stored in the output?
names(rf.model)

# Make a data frame with predictor names and their importance
imp <- importance(rf.model)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 20 predictors
imp.20 <- imp.sort[1:20, ]
# here's the whole list
imp.869 <- imp.sort[1:869, ]
# I want to export this data file;
write.csv(imp.869, file = 'RandomForestOTUImportance.csv')


# ggplot
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying Simulated Variable Groups")

ggplot(imp.869, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying Simulated Variable Groups")


# What are those OTUs?
otunames <- imp.20$predictors
r <- rownames(tax_table(physeq.alpha)) %in% otunames
kable(tax_table(physeq.alpha)[r, ])




#####RANDOM FOREST - rand#####

# Make a dataframe of training data with OTUs as column and samples as rows
predictors <- data.frame(otu_table(physeq.alpha))
dim(predictors)
predictors$SEX <- as.factor(meta.sim$SEX)

# Make one column for our outcome/response variable 
response <- as.factor(sample_data(physeq.alpha)$rv0)

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors)

# Random Forest Model
rf.model <- randomForest(response~., data = rf.data, ntree = 100)
print(rf.model)

# What variables are stored in the output?
names(rf.model)


# Make a data frame with predictor names and their importance
imp <- importance(rf.model)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 20 predictors
imp.r.20 <- imp.sort[1:20, ]
# here's the whole list
imp.r.869 <- imp.sort[1:869, ]
# I want to export this data file;
write.csv(imp.r.869, file = 'RandomForestOTUImportance_random.csv')


# ggplot
ggplot(imp.r.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying Simulated Variable Groups")

ggplot(imp.r.869, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying Simulated Variable Groups")


# What are those OTUs?
otunames <- imp.r.20$predictors
r <- rownames(tax_table(physeq.alpha)) %in% otunames
kable(tax_table(physeq.alpha)[r, ])

