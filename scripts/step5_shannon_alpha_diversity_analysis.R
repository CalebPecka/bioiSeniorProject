# Set your working directory to the current file directory
# This command can only be used while working in an RStudio environment.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
library(dplyr)
mappingNew <- left_join(mappingNewShort, alphaDiversity, 
                        by = c("sample_name" = "X"))

# Read in the clinical data
clinical <- read.csv("../data/clinical_metadata.csv")

# Merge the clinical data with the patient metadata
clinicalNew <- left_join(mappingNew, clinical, 
                         by = c("host_subject_id" = "ï..Code"))

################################################################################
# Perform linear model analysis based on alpha diversity                       #
################################################################################
library(GGally)
library(ggfortify)
library(car)

# First we need to reorganize our data to only include potentially relevant
# information in our model
data <- clinicalNew[,c(14,4,6,7,10,12,13,15,16,17,21,22,23,25,30,33,38,40)]
data <- data[which(is.na(data$shannon_entropy) == F),]

# Make a consistent NA value
data[data == ""] <- NA
data[data == "na"] <- NA

# Convert stuff to numeric
data[,1] <- as.numeric(data[,1])
data[,2] <- as.numeric(data[,2])
data[,4] <- as.numeric(data[,4])
data[,6] <- as.numeric(data[,6])
data[,7] <- as.numeric(data[,7])
data[,8] <- as.numeric(data[,8])
data[,10] <- as.numeric(data[,10])
data[,14] <- as.numeric(data[,14])
data[,15] <- as.numeric(sub("%", "", data[,15]))
data[,16] <- as.numeric(sub("%", "", data[,16]))
data[,17] <- as.numeric(sub("%", "", data[,17]))
data[,18] <- as.numeric(sub("%", "", data[,18]))


# Check for multicollinearity within our model.
ggpairs(data)

###############################################
# Violin plots
# Lets make some violin plots to check some
# important relations between these variables
###############################################
library(ggplot2)

set.seed(0)

# Violin plots show that alpha diversity does not change depending on autism
# versus neurotypical
C1 <- data[which(data$experimental_group == "neurotypical"), c(1,3)]
C2 <- data[which(data$experimental_group == "autism"), c(1,3)]

sig <- wilcox.test(C1$shannon_entropy, C2$shannon_entropy, paired = F)

dat <- rbind(C1, C2)
ggplot(dat, aes(x = experimental_group, y = shannon_entropy)) +
  geom_violin(scale = "width", adjust = 1, width = 0.5) + 
  geom_boxplot(width = 0.1, fill="white", outlier.shape = NA) +
  labs(title = paste0("Alpha Diversity Comparison in Autism\n", 
                      "Wilcox p-value: ", sig[["p.value"]]), 
       x = "Experimental Group", y = "Shannon Diversity") +
  theme(plot.title = element_text(hjust = 0.5)) +
  aes(fill = experimental_group)

# This data compares the change in alpha diversity over the course over time
# In general, autistic patients start with similar alpha diversity distributions
# that become larger than the neurotypical distribution with time.
C1 <- data[which(data$experimental_group == "neurotypical" & data$week), c(1,3)]
C2 <- data[which(data$experimental_group == "autism" & data$week == 0), c(1,3)]

sig <- wilcox.test(C1$shannon_entropy, C2$shannon_entropy, paired = F)

dat <- rbind(C1, C2)
ggplot(dat, aes(x = experimental_group, y = shannon_entropy)) +
  geom_violin(scale = "width", adjust = 1, width = 0.5) + 
  geom_boxplot(width = 0.1, fill="white", outlier.shape = NA) +
  labs(title = paste0("Alpha Diversity Comparison in Week 0 Autism\n", 
                      "Wilcox p-value: ", sig[["p.value"]]), 
       x = "Experimental Group", y = "Shannon Diversity") +
  theme(plot.title = element_text(hjust = 0.5)) +
  aes(fill = experimental_group)

C1 <- data[which(data$experimental_group == "neurotypical" & data$week), c(1,3)]
C2 <- data[which(data$experimental_group == "autism" & data$week == 18), c(1,3)]

sig <- wilcox.test(C1$shannon_entropy, C2$shannon_entropy, paired = F)

dat <- rbind(C1, C2)
ggplot(dat, aes(x = experimental_group, y = shannon_entropy)) +
  geom_violin(scale = "width", adjust = 1, width = 0.5) + 
  geom_boxplot(width = 0.1, fill="white", outlier.shape = NA) +
  labs(title = paste0("Alpha Diversity Comparison in Week 18 Autism\n", 
                      "Wilcox p-value: ", sig[["p.value"]]), 
       x = "Experimental Group", y = "Shannon Diversity") +
  theme(plot.title = element_text(hjust = 0.5)) +
  aes(fill = experimental_group)

# This observation is especially noticeable when comparing neurotypical
# individuals at week 0 only.
C1 <- data[which(data$experimental_group == "neurotypical" & data$week == 0), 
           c(1,3)]
C2 <- data[which(data$experimental_group == "autism" & data$week == 0), c(1,3)]

sig <- wilcox.test(C1$shannon_entropy, C2$shannon_entropy, paired = F)

dat <- rbind(C1, C2)
ggplot(dat, aes(x = experimental_group, y = shannon_entropy)) +
  geom_violin(scale = "width", adjust = 1, width = 0.5) + 
  geom_boxplot(width = 0.1, fill="white", outlier.shape = NA) +
  labs(title = paste0("Alpha Diversity Comparison in Week 0 Autism, 
                      Week 0 Neurotypical\n", "Wilcox p-value: ", 
                      sig[["p.value"]]), 
       x = "Experimental Group", y = "Shannon Diversity") +
  theme(plot.title = element_text(hjust = 0.5)) +
  aes(fill = experimental_group)

C1 <- data[which(data$experimental_group == "neurotypical" & data$week == 0), 
           c(1,3)]
C2 <- data[which(data$experimental_group == "autism" & data$week == 18), c(1,3)]

sig <- wilcox.test(C1$shannon_entropy, C2$shannon_entropy, paired = F)

dat <- rbind(C1, C2)
ggplot(dat, aes(x = experimental_group, y = shannon_entropy)) +
  geom_violin(scale = "width", adjust = 1, width = 0.5) + 
  geom_boxplot(width = 0.1, fill="white", outlier.shape = NA) +
  labs(title = paste0("Alpha Diversity Comparison in Week 18 Autism, 
                      Week 0 Neurotypical\n", "Wilcox p-value: ", 
                      sig[["p.value"]]), 
       x = "Experimental Group", y = "Shannon Diversity") +
  theme(plot.title = element_text(hjust = 0.5)) +
  aes(fill = experimental_group)




########################################
# Based on these observations, we need to create
# a more sophisticated model for microbial diversity
# in neurotypical individuals
########################################
data_neurotypical <- data[which(data$experimental_group == "neurotypical" & 
                                  data$week == 18),]

lm_alpha_with_age <- lm(shannon_entropy ~ age + gender + weight..pounds., 
                        data_neurotypical)
#lm_alpha_with_age <- lm(shannon_entropy ~ age, data_neurotypical)
summary(lm_alpha_with_age)

library(ggplot2)
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], 
                               y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
ggplotRegression(lm_alpha_with_age)

# Now we check if this relationship is affected by outcome
data_asd <- data[which(data$experimental_group == "autism"),]

lm_alpha_with_age <- lm(shannon_entropy ~ age + weight..pounds., data_asd)
#lm_alpha_with_age <- lm(shannon_entropy ~ age, data_asd)

summary(lm_alpha_with_age)

library(ggplot2)
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], 
                               y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
ggplotRegression(lm_alpha_with_age)

