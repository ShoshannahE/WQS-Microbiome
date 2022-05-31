###############
##  PACKAGES  
###############
#If you haven't already installed these, run line 5
#install.packages("pacman")

#These loads the specific subpackages we need, and checks to make sure the right packages are loaded.
library(pacman)
pacman::p_load("devtools","here", "BiocManager", "assertr", "tidyverse", "reshape2", "janitor", "dada2", "phyloseq", "qiime2R", "ggplot2", "DESeq2", "GUniFrac", "vegan", "LDM", "scales", "grid", "ggsci", "ggpubr", "viridis", "picante", "knitr", "gridExtra")
pacman::p_loaded()



###############
##  FILE NAMES & DIRECTORY
###############
# files for import (pulled directly from EIGC-provided files)
inputfiles <- list(
  feat = here::here("table.qza"), 
  taxo = here::here("gg_97_taxonomy.qza"), 
  treeroot= here::here("rooted_tree_masked_alignment.qza"),
  meta = here::here("metadata.txt")
)


###############
##  IMPORT DATA
###############
getwd()
# convert imported CHEAR files to phyloseq object 
phyobj<-qza_to_phyloseq(features=inputfiles$feat,
                        taxonomy=inputfiles$taxo,
                        tree=inputfiles$treeroot,
                        metadata=inputfiles$meta)

#I will export the OTU and taxonomy files to CSVs so that I can work with them outside of phyloseq.

# define individual objects
taxtable<-as.data.frame(tax_table(phyobj)) #get taxonomy
phytree<-phy_tree(phyobj) #get phytree
metadata<-sample_data(phyobj) #get sample metadata
asvs<-otu_table(phyobj) #get asvs



# review imported data
phyobj
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8186 taxa and 311 samples ]
# sample_data() Sample Data:       [ 311 samples by 16 sample variables ]
# tax_table()   Taxonomy Table:    [ 8186 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 8186 tips and 8024 internal nodes ]



##########################
##  SUBSET EXPERIMENTAL SAMPLES
##########################
# Subset sample to exclude controls 
sub.phyobj <-subset_samples(phyobj, SubjectGroup %in% (c("A", "B", "C", "Missing")))
sub.phyobj
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8186 taxa and 279 samples ]
# sample_data() Sample Data:       [ 279 samples by 16 sample variables ]
# tax_table()   Taxonomy Table:    [ 8186 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 8186 tips and 8024 internal nodes ]


###########
# convert phyloseq objects to data frames
OTU = as.data.frame(otu_table(sub.phyobj))
TAX = as.data.frame(tax_table(sub.phyobj))

##########

#the OTU table needs to be transposed
library(data.table)
otu_transpose <- transpose(OTU)
rownames(otu_transpose) <- colnames(OTU)
colnames(otu_transpose) <- rownames(OTU)


# export the data frames as CSV files
write.csv(otu_transpose,'micheals_asv_rows.csv')
write.csv(TAX,'micheals_tax_table.csv')


#I also gave the ASVs new numbers in the OTU and tax files because some downstream steps were having issues with the super long ASV names.

##### Moira's Processing Code #####
library(readxl)
library(dplyr)
library(tidyverse)
library(magrittr)
library(stringr)
library(reshape2)
library(lubridate)

#this code will merge the specimen IDs (SID)  to participant IDs (PIDs)
#Written by MB 3/8/22

##upload the microbioem dataset
asvs <- read.csv("micheals_asv_rows.csv")
##in this instance, the "index" should be the name of your ID colun--you didn't have an ID from what I saw so you should name the column:
names(asvs)[1] <- "index"
asvs2 <- asvs %>% mutate(
  SID= str_extract(index, "C[:punct:][:alnum:][:alnum:][:alnum:][:alnum:][:alnum:][:punct:]ST[:punct:][:alnum:][:alnum:]"))

#upload mapping file
map <- read.csv("1977_StoolMap_030822.csv")
#for year
all <- merge(asvs2, map, by="SID") ## you'll be able to merge with the epi data via the PID
##for duplicates--each participant only has 1 sample unless they have a duplicate, so this code works
ASVs <- names(all[,3:8188]) #the column #s of the asvs --this will be based off your dataset
all_small <- all[c("PID", ASVs)]

#make the data long
all2 <- melt( all_small, id.vars="PID")
#average based on PID & microbiome
all3 <- all2 %>%
  group_by(PID, variable) %>%
  summarise(avg=mean(value))
#make datset wide again, with the newly averaged values of the duplicate pairs
micro <- dcast(all3, PID ~ variable, value.var="avg")


#read in the epi file and change the pid lable
epi <- read.csv('1977_EPI_DATA.csv')
epi <- rename(epi, PID=pid)

#####


#####
#subsetting the epi data to only the variables I need.
subset <- epi[, c('PID',"age_stool", "antro_baz_depo", "birth_mode", "lactdurdays", "patb_5_depo", "meducat")]

#adding labels
subset = apply_labels(subset,
                     age_stool = "Age (years)",
                     antro_baz_depo = "BMI z-score",
                     birth_mode = "Birth mode (vaginal versus c-section)",
                     lactdurdays = "Number of days breastfed as infant",
                     patb_5_depo="Were antibiotics used in the past 6 months (yes vs no)", 
                     meducat="Maternal education")

#removing observations with NA
subset.clean <- subset[complete.cases(subset), ]
#169 obs with full data
#####



#now I need to read in the taxonomy table with ASV numbers.
tax <- read.csv('micheals_tax_table.csv')

row.names(subset.clean) = subset.clean$PID
row.names(micro) = micro$PID
row.names(tax) = tax$ASV


OTU.clean = micro[row.names(micro) %in% row.names(subset.clean),]
epi.clean = subset.clean[row.names(subset.clean) %in% row.names(OTU.clean),]
tax.clean = tax[row.names(tax) %in% colnames(OTU.clean),]
tax.clean = tax.clean[,-which(names(tax.clean) %in% c("ASV"))]
OTU.clean = OTU.clean[,-which(names(OTU.clean) %in% c("PID"))]



#making sure the data are all in the same order
OTU.clean = OTU.clean[order(row.names(OTU.clean)),]
epi.clean = epi.clean[order(row.names(epi.clean)),]




#Making these files into a phyloseq object
#have to reformat he files first
OTU.UF = otu_table(as.matrix(OTU.clean), taxa_are_rows=FALSE)
tax.UF = tax_table(as.matrix(tax.clean))
meta.UF = sample_data(epi.clean)

physeq = phyloseq(OTU.UF, tax.UF, meta.UF)
physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8186 taxa and 169 samples ]
# sample_data() Sample Data:       [ 169 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 8186 taxa by 7 taxonomic ranks ]



##########################
#########  FILTERING
##########################

### TAXANOMIC FILTERING ###
# prune ASVs that are not present in at least 2 samples 
po <- prune_taxa(taxa_sums(physeq) >= 2, physeq) # taxa_sums returns the total number of individuals observed from each species/taxa/OTU.
po 
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5804 taxa and 169 samples ]
# sample_data() Sample Data:       [ 169 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 5804 taxa by 7 taxonomic ranks ]


# What phyla are represented? 
get_taxa_unique(po, "Phylum")
# [1] NA                "Spirochaetes"    "Firmicutes"      "Cyanobacteria"   "Elusimicrobia"
# [6] "Actinobacteria"  "TM7"             "Bacteroidetes"   "Tenericutes"     "Euryarchaeota"
# [11] "Fusobacteria"    "Proteobacteria"  "Lentisphaerae"   "Chlamydiae"      "Verrucomicrobia"
# [16] "Synergistetes"   

# READ COUNTS
# Create table, number of features for each phyla 
table.features_phyla <- as.data.frame(table(tax_table(po)[,"Phylum"], exclude=NULL))
View(table.features_phyla) #Majority are Firmicutes & Bacteroidetes; there are 561 with missing phylum


# Remove NA and ambiguous phylum annotation 
subpo1 <- subset_taxa(po, !is.na(Phylum) & !Phylum %in% c("", "NA"))

table.features_phyla2 <- as.data.frame(table(tax_table(subpo1)[,"Phylum"], exclude=NULL))
View(table.features_phyla2)
subpo1 
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5243 taxa and 169 samples ]
# sample_data() Sample Data:       [ 169 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 5243 taxa by 7 taxonomic ranks ]


# drop samples below a read count of 10,000 
subpo2 <- prune_samples(sample_sums(subpo1) >= 10000, subpo1)
subpo2 
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5243 taxa and 165 samples ]
# sample_data() Sample Data:       [ 165 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 5243 taxa by 7 taxonomic ranks ]

#now setting a prevalence threshold of 10%
# convert phyloseq objects to data frames
OTU = as.data.frame(otu_table(subpo2))
TAX = as.data.frame(tax_table(subpo2))
META = as.data.frame(sample_data(subpo2))

x <- colnames(OTU)
OTU[OTU == 0] <- NA

#this checks % non-detect in each column. we are going for anything 90% or less non-detect (or 10%+ detect)
otus_10det <- OTU[lapply( OTU[1:165,], function(x) sum(is.na(x)) / length(x) ) < 0.9 ]

#this takes it down to 500 otus
otus_10det[is.na(otus_10det)] = 0


OTU2 = otus_10det[row.names(otus_10det) %in% row.names(META),]
META2 = META[row.names(META) %in% row.names(OTU2),]
TAX2 = TAX[row.names(TAX) %in% colnames(OTU2),]



##########
# export the data frames as CSV files
write.csv(OTU2,'micheals_HMP_asv_table.csv')
write.csv(TAX2,'micheals_HMP_tax_table.csv')
write.csv(META2,'micheals_HMP_epi_table.csv')


