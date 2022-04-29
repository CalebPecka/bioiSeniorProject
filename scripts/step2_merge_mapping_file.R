# Set your working directory to the current file directory
# This command can only be used while working in an RStudio environment.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
                                         "weight"))]

write.table(mappingNewShort, "../data/mapping_file_full_metadata.tsv", 
            quote = F, sep = "\t", row.names = F)
