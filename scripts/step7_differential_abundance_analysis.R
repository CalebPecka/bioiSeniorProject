# Set your working directory to the current file directory
# This command can only be used while working in an RStudio environment.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in the level 7 (species level) csv file
level7 <- read.csv("../results/level-7.csv")

# Read in the metadata file
fullMetadata <- read.delim("../data/10532_20180602-081142.txt", sep = "\t")
mappingFiles <- read.delim("../data/merged_mapping_file.txt", sep = "\t")

library(dplyr)
mappingNew <- left_join(mappingFiles, fullMetadata, by = "sample_name")

mappingNewShort <- mappingNew[,which(colnames(mappingNew) %in% 
                                       c("sample_name", 
                                         "barcode", 
                                         "run", 
                                         "bmi", 
                                         "bbt_donor_id",
                                         "description", 
                                         "experimental_group", 
                                         "height", 
                                         "responder", 
                                         "sex", 
                                         "time_point", 
                                         "week", 
                                         "weight",
                                         "host_subject_id"))]

# Read in the qiime2 alpha diversity
alphaDiversity <- read.delim("../results/alpha-diversity-shannon.tsv", 
                             sep = "\t")

# Join the metadata and alpha diversity
mappingNew <- left_join(mappingNewShort, alphaDiversity, 
                        by = c("sample_name" = "X"))

# Read in the clinical data
clinical <- read.csv("../data/clinical_metadata.csv")

# Merge the clinical data with the patient metadata
clinicalNew <- left_join(mappingNew, clinical, 
                         by = c("host_subject_id" = "ï..Code"))

################################################################################
# Condition 0:
# Is the patient neurotypical or ASD symptomatic?
################################################################################
neurotypical <- clinicalNew[which(clinicalNew$experimental_group == 
                                    "neurotypical"),]$sample_name
ASD_symptomatic <- clinicalNew[which(clinicalNew$experimental_group == 
                                       "autism"),]$sample_name

################################################################################
# Condition 1:
# Did the patient's GI symptoms respond to FMT?
################################################################################
GI_responders <- clinicalNew[which(clinicalNew$responder == 
                                     "Responder"),]$sample_name
GI_non_responders <- clinicalNew[which(clinicalNew$responder == 
                                         "Non-responder"),]$sample_name

################################################################################
# Condition 2:
# Did the patient's microbiome recover from FMT?
################################################################################
# We will inspect two possible approaches to answering this question.

# Option 1: Did alpha diversity increase after week18 compared to initial?
#      Rational: Increased diversity is a common metric of gut health
#               Qin J, Li R, Raes J, Arumugam M, Burgdorf KS, Manichanh C, 
#               Nielsen T, Pons N, Levenez F, Yamada T, Mende DR, Li J, Xu J, 
#               Li S, Li D, Cao J, Wang B, Liang H, Zheng H, Xie Y, Tap J, 
#               Lepage P, Bertalan M, Batto J-M, Hansen T, Le Paslier D, 
#               Linneberg A, Nielsen HB, Pelletier E, Renault P, 
#               Sicheritz-Ponten T, Turner K, Zhu H, Yu C, Li S, Jian M, Zhou Y,
#               Li Y, Zhang X, Li S, Qin N, Yang H, Wang J, Brunak S, DorÃ© J, 
#               Guarner F, Kristiansen K, Pedersen O, Parkhill J, Weissenbach J,
#               MetaHIT Consortium, et al.. 2010. A human gut microbial gene 
#               catalogue established by metagenomic sequencing. Nature 
#               464:59â€“65. doi:10.1038/nature08821.

# This segment of code determines which patients saw in increase in alpha
# diversity after week18 
increasedDiversityPatients <- c()
decreasedDiversityPatients <- c()

# First, we collect a list of the patient ID's with autism diagnosis
for (i in unique(clinicalNew[which(clinicalNew$experimental_group == 
                                   "autism"),]$host_subject_id))
{
  tmp <- clinicalNew[which(clinicalNew$host_subject_id == i),]
  tmp <- tmp[which(!is.na(tmp$shannon_entropy)),]
  maximumWeek <- max(as.numeric(tmp$week)[which(!is.na(as.numeric(tmp$week)))])
  minimumWeek <- min(as.numeric(tmp$week)[which(!is.na(as.numeric(tmp$week)))])
  # Extract week 18 diversity
  week18 <- tmp[which(tmp$week == maximumWeek),]$shannon_entropy 
  # Extract week 0 diversity
  week0 <- tmp[which(tmp$week == minimumWeek),]$shannon_entropy 
  
  # Each sample had multiple runs. Here, we average the runs.
  week18_mean <- mean(week18, na.rm = T)
  week0_mean <- mean(week0, na.rm = T)
  
  # Check for NA values
  if (!is.na(week0_mean) && !is.na(week18_mean)){
    # If the patient's alpha diversity increased, add them to the list
    if (week18_mean > week0_mean){ 
      increasedDiversityPatients <- c(increasedDiversityPatients, i)
    }
    else{
      decreasedDiversityPatients <- c(decreasedDiversityPatients, i)
    }
  }
}

# Now we collect the list of all the responding sample IDs based on the barcode
diversity_option1_responders <- clinicalNew[which(clinicalNew$host_subject_id 
                                  %in% increasedDiversityPatients),]$sample_name
diversity_option1_non_responders <-clinicalNew[which(clinicalNew$host_subject_id
                                  %in% decreasedDiversityPatients),]$sample_name

# Option 2: Did the patient's microbiome become more similar to the donor ID
# that supplied the FMT?
#
# This is calculated using Jaccard distance between the sample ID's, where a
# score of zero is an identical microbial profile.
jaccard <- read.delim("../results/distance-matrix.tsv", sep = "\t", 
                      row.names = 1)
colnames(jaccard) <- gsub("X", "", colnames(jaccard))

increasedDistancePatients <- c()
decreasedDistancePatients <- c()

# First, we collect a list of the patient ID's with autism diagnosis
for (i in unique(clinicalNew[which(clinicalNew$experimental_group == 
                                   "autism"),]$host_subject_id)){
  tmp <- clinicalNew[which(clinicalNew$host_subject_id == i),]
  donor <- names(sort(table(tmp$bbt_donor_id), decreasing = TRUE)[1])
  
  # Collect the donor IDs - multiple runs are possible
  donor_IDs <- clinicalNew[which(clinicalNew$host_subject_id == 
                                   donor),]$sample_name
  
  # Collect the autism patient IDs
  week0_IDs <- tmp[which(tmp$week == 0),]$sample_name
  week18_IDs <- tmp[which(tmp$week == 18),]$sample_name
  
  week0_jaccard_indicies <- c()
  week18_jaccard_indicies <- c()
  for (d_ID in donor_IDs){
    for (ID_0 in week0_IDs){
      addNew <- unlist(jaccard[which(row.names(jaccard) %in% c(d_ID)),
                        which(colnames(jaccard) %in% c(ID_0))])
      week0_jaccard_indicies <- c(week0_jaccard_indicies, addNew)
    }
    for (ID_18 in week18_IDs){
      addNew <- unlist(jaccard[which(row.names(jaccard) %in% c(d_ID)),
                        which(colnames(jaccard) %in% c(ID_18))])
      week18_jaccard_indicies <- c(week18_jaccard_indicies, addNew)
    }
  }
  
  # Each sample had multiple runs. Here, we average the runs.
  week0_mean <- mean(as.numeric(week0_jaccard_indicies), na.rm = T)
  week18_mean <- mean(as.numeric(week18_jaccard_indicies), na.rm = T)
  
  # Check for NA values
  if (!is.na(week0_mean) && !is.na(week18_mean)){
    if (week18_mean > week0_mean){ # If the patient's alpha diversity increased,
      # add them to the list
      increasedDistancePatients <- c(increasedDiversityPatients, i)
    }
    else{
      decreasedDistancePatients <- c(decreasedDiversityPatients, i)
    }
  }
}

# Now we collect the list of all the responding sample IDs based on the barcode
diversity_option2_responders <- clinicalNew[which(clinicalNew$host_subject_id 
                                  %in% decreasedDistancePatients),]$sample_name
diversity_option2_non_responders <- clinicalNew[which(clinicalNew$
                  host_subject_id %in% increasedDistancePatients),]$sample_name

################################################################################
# Condition 3:
# Did the patient's clinical/behavioral symptoms respond to FMT?
################################################################################

# SRS official scoring criteria
# total; 74+ is severe, 45-73 is mild ASD, <45 is normal
clinicalNew$SRS.week.18 <- as.numeric(clinicalNew$SRS.week.18)
clinical_responders <- clinicalNew[which(clinicalNew$SRS.week.18 <= 73),
                                   ]$sample_name
clinical_non_responders <- clinicalNew[which(clinicalNew$SRS.week.18 > 73),
                                       ]$sample_name
    # Note, only 1 patient started week 0 with SRS less than 73. This patient 
    # still had a  13% reduction in SRS
    # symptom scores, so I opted to include them in the list of clinical 
    # responders.

################################################################################
# Figure creation                                                              #
################################################################################
library(ggplot2)
library(ggpubr)

set.seed(0)

# Clinical symptoms depending on whether or not the patient's alpha diversity 
# increased after week18
experimental_group <- c()
week <- c()
C1 <- as.numeric(clinicalNew[which(clinicalNew$week == "18" & 
      clinicalNew$sample_name %in% diversity_option1_responders),]$SRS.week.18)
experimental_group <- c(experimental_group, rep("R Week 18", length(C1)))
week <- c(week, rep("18", length(C1)))
C2 <- as.numeric(clinicalNew[which(clinicalNew$week == "0" & 
      clinicalNew$sample_name %in% diversity_option1_responders),]$SRS.initial)
experimental_group <- c(experimental_group, rep("R Week 0", length(C2)))
week <- c(week, rep("0", length(C2)))
C3 <- as.numeric(clinicalNew[which(clinicalNew$week == "18" & 
  clinicalNew$sample_name %in% diversity_option1_non_responders),]$SRS.week.18)
experimental_group <- c(experimental_group, rep("NR Week 18", length(C3)))
week <- c(week, rep("18", length(C3)))
C4 <- as.numeric(clinicalNew[which(clinicalNew$week == "0" & 
  clinicalNew$sample_name %in% diversity_option1_non_responders),]$SRS.initial)
experimental_group <- c(experimental_group, rep("NR Week 0", length(C4)))
week <- c(week, rep("0", length(C4)))

SRS.score <- c(C1, C2, C3, C4)

dat <- as.data.frame(SRS.score)
dat$Group <- experimental_group
dat$Week <- week

my_comparisons <- list( c("NR Week 0", "NR Week 18"), 
                        c("R Week 0", "R Week 18") )

# This plot showcases how patient's with increased alpha diversity after FMT see
# better clinical improvements than ASD patient's whose alpha diversity does NOT
# increase.
ggplot(dat, aes(x = Group, y = SRS.score)) +
  geom_boxplot(scale = "width", adjust = 1, width = 0.5) + 
  aes(fill = Group) + 
  stat_compare_means(comparisons = my_comparisons, label.y = c(150, 150))+
  stat_compare_means(label.y = 170, label.x = 2.15) +
  geom_hline(yintercept=73, linetype="dashed", color = "red", size = 1)

# We're also interested in determining whether the patient's diet is 
# contributing to the lack of improved alpha diversity after FMT. If the patient
# does not consume food that supports the bacteria, the bacteria will not 
# survive in their new environment.
food <- read.csv("../../ASD_diet.csv", row.names = 1)

R <- unique(clinicalNew[which(clinicalNew$week == "18" & clinicalNew$sample_name
                          %in% diversity_option1_responders),]$host_subject_id)
NR <- unique(clinicalNew[which(clinicalNew$week == "18" & 
  clinicalNew$sample_name %in% diversity_option1_non_responders),
  ]$host_subject_id)

both_R_and_NR <- c(R, NR)

food_dat <- food[which(row.names(food) %in% both_R_and_NR),]

food_dat$group <- rep("NR", nrow(food_dat))
food_dat[which(row.names(food_dat) %in% R),]$group <- "R"

clinicalNewSmall <-clinicalNew[which(!duplicated(clinicalNew$host_subject_id)),]
food_dat$name <- row.names(food_dat)
food_dat <- left_join(food_dat, clinicalNewSmall, by = 
                        c("name" = "host_subject_id"))

# These boxplots compare the diets of responders versus non-responders to FMT
# micrbiome recovery
boxplot(food_dat$DT_FIBE ~ food_dat$group, main = "Fiber Diet")
boxplot(food_dat$DT_CARB ~ food_dat$group, main = "Carb Diet")
boxplot(food_dat$DT_TFAT ~ food_dat$group, main = "Fat Diet")
boxplot(food_dat$DT_PROT ~ food_dat$group, main = "Protein Diet")
boxplot(food_dat$DT_KCAL ~ food_dat$group, main = "Kilocalories Consumed")

# Now we look at metadata related to their age and weight
food_dat$weight <- as.numeric(food_dat$weight)
boxplot(food_dat$age ~ food_dat$group, main = "Age")
boxplot(food_dat$weight ~ food_dat$group, main = "Weight")
boxplot(food_dat$BMI ~ food_dat$group, main = "BMI")

# Now that we have a rough idea what the data looks like, we can make a neater
# plot that ties all this information together, while simultaneously normalizing
# food intake based on kilocalories consumed.

# The following patients recovered microbiome AND improved clinical symptoms
alpha.Recovery.clinical.improvement <- unique(clinicalNew[which(
  clinicalNew$sample_name %in% diversity_option1_responders & 
    clinicalNew$SRS.week.18 <= 73),]$host_subject_id)

# The following patients recovered microbiome BUT DID NOT improve clinical 
# symptoms above the 74 SRS threshold
alpha.Recovery.clinical.stagnant <- unique(clinicalNew[which(
  clinicalNew$sample_name %in% diversity_option1_responders & 
    clinicalNew$SRS.week.18 >= 74),]$host_subject_id)

# And the final group is patients whose microbiomes DID NOT recovery
alpha.noRecovery <- decreasedDiversityPatients

# The data is restructured for ggpplot
all.groups <- c(alpha.Recovery.clinical.improvement, 
                alpha.Recovery.clinical.stagnant, alpha.noRecovery)
food_dat <- food[which(row.names(food) %in% all.groups),]
food_dat$group <- rep("SRS.improve", nrow(food_dat))
food_dat[which(row.names(food_dat) %in% alpha.Recovery.clinical.stagnant),
         ]$group <- "SRS.stagnant"
food_dat[which(row.names(food_dat) %in% alpha.noRecovery),]$group <- 
  "noRecovery"
clinicalNewSmall <- clinicalNew[which(
  !duplicated(clinicalNew$host_subject_id)),]
food_dat$name <- row.names(food_dat)
food_dat <- left_join(food_dat, clinicalNewSmall, by = 
                        c("name" = "host_subject_id"))

# Each food group is normalized by the number of kcalories consumed by the
# patient, multiplied by 1000
food_subset <- food_dat[,c(1:6)]
food_subset$DT_FIBE <- food_subset$DT_FIBE/food_subset$DT_KCAL*1000/
  mean(food_subset$DT_FIBE)
food_subset$DT_CARB <- food_subset$DT_CARB/food_subset$DT_KCAL*1000/
  mean(food_subset$DT_CARB)
food_subset$DT_TFAT <- food_subset$DT_TFAT/food_subset$DT_KCAL*1000/
  mean(food_subset$DT_TFAT)
food_subset$DT_PROT <- food_subset$DT_PROT/food_subset$DT_KCAL*1000/
  mean(food_subset$DT_PROT)
food_subset$DT_KCAL <- NULL

Food.Group <- c()
Recovery.Subset <- c()
Value <- c()

# Create a data frame
for (i in 1:4){
  for (j in 1:nrow(food_subset)){
    Food.Group <- c(Food.Group, colnames(food_subset)[i])
    Recovery.Subset <- c(Recovery.Subset, food_subset$group[j])
    Value <- c(Value, food_subset[j, i])
  }
}

food_ggplot <- as.data.frame(Food.Group)
food_ggplot$Recovery.Subset <- Recovery.Subset
food_ggplot$Value <- Value

# Grouped boxplot
ggplot(food_ggplot, aes(x = Food.Group, y = Value, fill = Recovery.Subset)) + 
  geom_boxplot() + stat_compare_means(label.y = c(1.4, 1.5, 1.4, 1.5))

################################################################################
# The same plots are made, this time removing the differentiation between      #
# clinical recoveries in patients whose microbiome recovered.                  #
################################################################################

# Patients whose microbiome DID recover
alpha.Recovery <- increasedDiversityPatients

# The data is restructured for ggpplot
all.groups <- c(alpha.Recovery, alpha.noRecovery)
food_dat <- food[which(row.names(food) %in% all.groups),]
food_dat$group <- rep("Recovery", nrow(food_dat))
food_dat[which(row.names(food_dat) %in% alpha.noRecovery),]$group <- 
  "noRecovery"
clinicalNewSmall <- clinicalNew[which(
  !duplicated(clinicalNew$host_subject_id)),]
food_dat$name <- row.names(food_dat)
food_dat <- left_join(food_dat, clinicalNewSmall, by = 
                        c("name" = "host_subject_id"))

# Each food group is normalized by the number of kcalories consumed by the
# patient, multiplied by 1000
food_subset <- food_dat[,c(1:6)]
food_subset$DT_FIBE <- food_subset$DT_FIBE/food_subset$DT_KCAL*1000/
  mean(food_subset$DT_FIBE)
food_subset$DT_CARB <- food_subset$DT_CARB/food_subset$DT_KCAL*1000/
  mean(food_subset$DT_CARB)
food_subset$DT_TFAT <- food_subset$DT_TFAT/food_subset$DT_KCAL*1000/
  mean(food_subset$DT_TFAT)
food_subset$DT_PROT <- food_subset$DT_PROT/food_subset$DT_KCAL*1000/
  mean(food_subset$DT_PROT)
food_subset$DT_KCAL <- NULL

Food.Group <- c()
Recovery.Subset <- c()
Value <- c()

# Create a data frame
for (i in 1:4){
  for (j in 1:nrow(food_subset)){
    Food.Group <- c(Food.Group, colnames(food_subset)[i])
    Recovery.Subset <- c(Recovery.Subset, food_subset$group[j])
    Value <- c(Value, food_subset[j, i])
  }
}

food_ggplot <- as.data.frame(Food.Group)
food_ggplot$Recovery.Subset <- Recovery.Subset
food_ggplot$Value <- Value

# Grouped boxplot
ggplot(food_ggplot, aes(x = Food.Group, y = Value, fill = Recovery.Subset)) + 
  geom_boxplot() + stat_compare_means(label.y = c(1.4, 1.5, 1.4, 1.5), 
                                      method = "wilcox")

################################################################################
# No we're going to investigate which bacteria are different between these     #
# recovery subgroups                                                           #
################################################################################

library("ggpubr")

R <- unique(clinicalNew[which(clinicalNew$week >= "16" & 
                                clinicalNew$sample_name %in% 
                                diversity_option1_responders),]$sample_name)
NR <- unique(clinicalNew[which(clinicalNew$week >= "16" & 
                                 clinicalNew$sample_name %in% 
                                diversity_option1_non_responders),]$sample_name)

bacteria <- c()
pVal <- c()
meanR <- c()
meanNR <- c()

for (column in 4:580){
  Rgroup <- level7[which(level7$index %in% R),column]
  NRgroup <- level7[which(level7$index %in% NR),column]
  results <- wilcox.test(Rgroup, NRgroup)
  
  if (!is.na(results[["p.value"]])){
      bacteria <- c(bacteria, colnames(level7)[column])
      pVal <- c(pVal, results[["p.value"]])
      meanR <- c(meanR, mean(Rgroup))
      meanNR <- c(meanNR, mean(NRgroup))
  }
}

sigBacteria <- as.data.frame(bacteria)
sigBacteria$pVal <- pVal
sigBacteria$meanR <- meanR
sigBacteria$meanNR <- meanNR
sigBacteria <- sigBacteria[grep("g__", sigBacteria$bacteria),]
#sigBacteria <- sigBacteria[which(sigBacteria$meanR + sigBacteria$meanNR > 10),]
#sigBacteria$pVal.adjust <- p.adjust(sigBacteria$pVal, method = "BH", n = 
# nrow(sigBacteria))

sigBacteria$name <- gsub('.__', '', gsub('.*g__', '', sigBacteria$bacteria))
sigBacteria$baseMean <- sigBacteria$meanR
sigBacteria$log2FoldChange <- log2(sigBacteria$meanNR/(sigBacteria$meanR))
sigBacteria$padj <- sigBacteria$pVal
sigBacteria$detection_call <- rep(1, nrow(sigBacteria))
sigBacteria <- sigBacteria[,c(5:ncol(sigBacteria))]

ggmaplot(sigBacteria, main = expression("Bacterial Abundance in FMT Microbiome 
                                        Recovery"), 
         fdr = 0.05, fc = 2, size = 2,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(sigBacteria$name),
         legend = "top", top = 10,
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())

# No statistically significant bacterial abundance differences are observed if
# we separate patients based on clinical recovery
R <- unique(clinicalNew[which(clinicalNew$week >= "16" & clinicalNew$
        host_subject_id %in% alpha.Recovery.clinical.improvement),]$sample_name)
NR <- unique(clinicalNew[which(clinicalNew$week >= "16" & 
              clinicalNew$host_subject_id %in% 
                alpha.Recovery.clinical.stagnant),]$sample_name)

bacteria <- c()
pVal <- c()
meanR <- c()
meanNR <- c()

for (column in 4:580){
  Rgroup <- level7[which(level7$index %in% R),column]
  NRgroup <- level7[which(level7$index %in% NR),column]
  results <- wilcox.test(Rgroup, NRgroup)
  
  if (!is.na(results[["p.value"]])){
    bacteria <- c(bacteria, colnames(level7)[column])
    pVal <- c(pVal, results[["p.value"]])
    meanR <- c(meanR, mean(Rgroup))
    meanNR <- c(meanNR, mean(NRgroup))
  }
}

sigBacteria <- as.data.frame(bacteria)
sigBacteria$pVal <- pVal
sigBacteria$meanR <- meanR
sigBacteria$meanNR <- meanNR
sigBacteria <- sigBacteria[grep("g__", sigBacteria$bacteria),]
#sigBacteria <- sigBacteria[which(sigBacteria$meanR + sigBacteria$meanNR > 10),]
#sigBacteria$pVal.adjust <- p.adjust(sigBacteria$pVal, method = "BH", n = 
# nrow(sigBacteria))

sigBacteria$name <- gsub('.__', '', gsub('.*g__', '', sigBacteria$bacteria))
sigBacteria$baseMean <- sigBacteria$meanR
sigBacteria$log2FoldChange <- log2(sigBacteria$meanNR/(sigBacteria$meanR))
sigBacteria$padj <- sigBacteria$pVal
sigBacteria$detection_call <- rep(1, nrow(sigBacteria))
sigBacteria <- sigBacteria[,c(5:ncol(sigBacteria))]

ggmaplot(sigBacteria, main = expression("Bacterial Abundance in FMT Clinical 
                                        Recovery"), 
         fdr = 0.05, fc = 2, size = 2,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(sigBacteria$name),
         legend = "top", top = 10,
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())

# Use this file if you want to check plasma metabolites, or uncomment the line
# beneath if you want to examine stool metabolites.
sample_metabolite_data <- read.csv("../data/plasma_metabolites.csv")
sample_metabolite_data <- read.csv("../data/stool_metabolites.csv")

# First, we will compare the sample_metabolite_data metabolites of those who did
# and did not increase microbiome alpha diversity
R <- alpha.Recovery
NR <- alpha.noRecovery

R <- alpha.Recovery.clinical.improvement
NR <- alpha.Recovery.clinical.stagnant

metabolite <- c()
pVal <- c()
meanR <- c()
meanNR <- c()

for (i in 1:nrow(sample_metabolite_data)){
  Rgroup <- unlist(sample_metabolite_data[i, 
                              which(colnames(sample_metabolite_data) %in% R)])
  NRgroup <- unlist(sample_metabolite_data[i, 
                              which(colnames(sample_metabolite_data) %in% NR)])
  
  if (!is.na(sum(Rgroup)) && !is.na(sum(NRgroup))){
    results <- wilcox.test(Rgroup, NRgroup)
    
    if (!is.na(results[["p.value"]])){
      metabolite <- c(metabolite, sample_metabolite_data$ï..BIOCHEMICAL[i])
      pVal <- c(pVal, results[["p.value"]])
      meanR <- c(meanR, mean(Rgroup))
      meanNR <- c(meanNR, mean(NRgroup))
    }
  }
  
}

sigMetabolite <- as.data.frame(metabolite)
sigMetabolite$pVal <- pVal
sigMetabolite$meanR <- meanR
sigMetabolite$meanNR <- meanNR
sigMetabolite$foldChange <- sigMetabolite$meanR / sigMetabolite$meanNR

# These values can be uncommented for the table values of the paper
#sigMetabolite <- sigMetabolite[order(sigMetabolite$pVal),]
#sigMetabolite <- sigMetabolite[c(1:100),]
#sigMetabolite$meanR <- NULL
#sigMetabolite$meanNR <- NULL
#sigMetabolite$pVal.adjust <- p.adjust(sigMetabolite$pVal, 
#                                      method = "BH", n = nrow(sigMetabolite))

sigMetabolite$name <- gsub('.__', '', gsub('.*g__', '', 
                                           sigMetabolite$metabolite))
sigMetabolite$baseMean <- sigMetabolite$meanR
sigMetabolite$log2FoldChange <- log2(sigMetabolite$meanNR/(sigMetabolite$meanR))
sigMetabolite$padj <- sigMetabolite$pVal
sigMetabolite$detection_call <- rep(1, nrow(sigMetabolite))
sigMetabolite <- sigMetabolite[,c(6:ncol(sigMetabolite))]

ggmaplot(sigMetabolite, main = expression("Stool Metabolite Abundance in 
                                          FMT Clinical Recovery"), 
         fdr = 0.1, fc = 2, size = 2,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(sigMetabolite$name),
         legend = "top", top = 10,
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())

################################################################################
# We can also check the same data with sub-pathways
R <- alpha.Recovery
NR <- alpha.noRecovery

# R <- alpha.Recovery.clinical.improvement
# NR <- alpha.Recovery.clinical.stagnant

metabolite <- c()
pVal <- c()
meanR <- c()
meanNR <- c()

for (i in unique(sample_metabolite_data$SUB_PATHWAY)){
  Rgroup <- unlist(sample_metabolite_data[which(sample_metabolite_data$
            SUB_PATHWAY == i), which(colnames(sample_metabolite_data) %in% R)])
  NRgroup <- unlist(sample_metabolite_data[which(sample_metabolite_data$
            SUB_PATHWAY == i), which(colnames(sample_metabolite_data) %in% NR)])
  
  if (!is.na(sum(Rgroup)) && !is.na(sum(NRgroup))){
    results <- wilcox.test(Rgroup, NRgroup)
    
    if (!is.na(results[["p.value"]])){
      metabolite <- c(metabolite, i)
      pVal <- c(pVal, results[["p.value"]])
      meanR <- c(meanR, mean(Rgroup))
      meanNR <- c(meanNR, mean(NRgroup))
    }
  }
}  

sigMetabolite <- as.data.frame(metabolite)
sigMetabolite$pVal <- pVal
sigMetabolite$meanR <- meanR
sigMetabolite$meanNR <- meanNR
sigMetabolite$pVal.adjust <- p.adjust(sigMetabolite$pVal, 
                                      method = "BH", n = nrow(sigMetabolite)) 




