#Data pre-processing to transaction table

library(tidyverse)
library(data.table)
library(curatedMetagenomicData)
library(phyloseq)
library(mia)
library(vegan)

###Healthy
sampleMetadata <- curatedMetagenomicData::sampleMetadata %>%
  filter(age >= 18 &
           age <= 65 &
           body_site == "stool" &
           antibiotics_current_use == 'no' &
           sequencing_platform == 'IlluminaHiSeq' &
           study_condition == 'control')

sampleMetadata <- sampleMetadata %>%
  group_by(study_name, subject_id) %>% 
  ungroup()


###Disease
unique(curatedMetagenomicData::sampleMetadata$study_condition)

#CRC - 368
#IBD - 768
#IGT -199
#T2D - 164

Disease <- c("CRC", "IBD", "IGT", "T2D")
sampleMetadata_disease <- curatedMetagenomicData::sampleMetadata %>%
  filter(age >= 18 &
           age <= 65 &
           body_site == "stool" &
           #antibiotics_current_use == 'no' &
           sequencing_platform == 'IlluminaHiSeq' &
           study_condition %in% Disease)

sampleMetadata_disease <- sampleMetadata_disease %>%
  group_by(study_name, subject_id) %>% 
  ungroup()

#================taxonomic profile==========================
stool_profile <- curatedMetagenomicData::returnSamples(rbind(sampleMetadata,sampleMetadata_disease), 
                                                       dataType = "relative_abundance")
species_stool <- splitByRanks(stool_profile, rank = "species")[[1]]
stool_phyloseq <- makePhyloseqFromTreeSummarizedExperiment(species_stool, abund_values = "relative_abundance")

present_proportion <- apply(stool_phyloseq@otu_table,1,function(i) length(which(i != 0))/length(i))


OTU <- as.data.frame(stool_phyloseq@otu_table)

OTU_H <- as.data.frame(stool_phyloseq@otu_table) |> select(sampleMetadata$sample_id)
transaction_H <- apply(OTU_H,2,function(i) ifelse(i>0,1,0))
transaction_H_filter_10 <- transaction_H[names(which(present_proportion >= 0.1)),]
dim(transaction_H_filter_10)
fwrite(t(transaction_H_filter_10),"transaction_H_filter_10.csv")

transaction_H_filter_10 <- fread("transaction_H_filter_10.csv")
transaction_H_filter_10_false <- transaction_H_filter_10
transaction_H_filter_10_false <- as.data.frame(ifelse(transaction_H_filter_10_false == 0,1,0))
colnames(transaction_H_filter_10_false) <- paste0("-",colnames(transaction_H_filter_10_false))
transaction_H_filter_10_pos_neg <- cbind(transaction_H_filter_10,transaction_H_filter_10_false)
fwrite(transaction_H_filter_10_pos_neg,"transaction_H_filter_10+-.csv")


#IBD can be replaced by any disease names in c("CRC", "IBD", "IGT", "T2D")

OTU_IBD <- as.data.frame(stool_phyloseq@otu_table) |> 
  select(sampleMetadata_disease[which(sampleMetadata_disease$study_condition == "IBD"),]$sample_id)
transaction_IBD <- apply(OTU_IBD,2,function(i) ifelse(i>0,1,0))
transaction_IBD_filter_10 <- transaction_IBD[names(which(present_proportion >= 0.1)),]
dim(transaction_IBD_filter_10)
fwrite(t(transaction_IBD_filter_10),"transaction_IBD_filter_10.csv")


transaction_IBD_filter_10 <- fread("transaction_IBD_filter_10.csv")
transaction_IBD_filter_10_false <- transaction_IBD_filter_10
transaction_IBD_filter_10_false <- as.data.frame(ifelse(transaction_IBD_filter_10_false == 0,1,0))
colnames(transaction_IBD_filter_10_false) <- paste0("-",colnames(transaction_IBD_filter_10_false))
transaction_IBD_filter_10_pos_neg <- cbind(transaction_IBD_filter_10,transaction_IBD_filter_10_false)
fwrite(transaction_IBD_filter_10_pos_neg,"transaction_IBD_filter_10+-.csv")


