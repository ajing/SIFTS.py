# dual (or triple) mutational distribution from COSMIC

# The dual mutataion distribution on UniProt disease nsSNP and polymorphism nsSNP

## load data from humsavar.txt
humsavar <- read.table("/users/ajing/SNPDist/Data/humsavar_tru.txt.out_test", sep = "\t", header = F, quote = "", na.string = "-")
humsavar <- read.table("humsavar_tru.txt.out_test", sep = "\t", header = F, quote = "", na.string = "-")


colnames(humsavar) <- c("genename", "uniprot_id", "var_id", "aa_before", "resnum", "aa_after", "var_type", "snp_id", "disease")

dist_each_uniprot <- merge(humsavar, humsavar, by = "uniprot_id")
dist_each_uniprot$dist <- dist_each_uniprot$resnum.x - dist_each_uniprot$resnum.y
dist_each_uniprot <- subset(dist_each_uniprot, dist > 0)
dist_each_uniprot$var_type_c <- apply(dist_each_uniprot, 1, function(x){paste(sort(c(x[["var_type.x"]],x[["var_type.y"]])), collapse = " to ")})

# The mutation distribution between disease nsSNP and polymorphism nsSNP
require(ggplot2)
ggplot(subset(dist_each_uniprot, var_type_c %in% c("Disease to Disease", "Polymorphism to Polymorphism", "Unclassified to Unclassified")), aes(x = dist, color = var_type_c)) + scale_x_log10() + geom_density()
ggsave("snp_dist.jpg")

# The mutation distribution between cancer disease nsSNP and polymorphism nsSNP
ggplot(rbind(subset(dist_each_uniprot, var_type_c %in% c("Disease to Disease", "Polymorphism to Polymorphism", "Unclassified to Unclassified") & grepl("cancer", disease.x) & grepl("cancer", disease.y)), subset(dist_each_uniprot, var_type_c %in% c("Polymorphism to Polymorphism"))), aes(x = dist, color = var_type_c) ) + scale_x_log10() + geom_density()
ggsave("snp_dist_cancer.jpg")




#################################################
#### For cosmic
#################################################
cosmic <- read.table("/users/ajing/SNPDist/Data/CosmicCompleteExport.tsv", sep = "\t", header = T, quote = "", na.string = "")

# Mutation in coding region
library("stringr")
cosmic_missense <- subset(cosmic, Mutation.Description == "Substitution - Missense")
cosmic_missense$Mutation.AA.Before <- apply(cosmic_missense, 1, function(x){str_sub(x[["Mutation.AA"]], 3, 3)}) 
cosmic_missense$Mutation.AA.After <- apply(cosmic_missense, 1, function(x){str_sub(x[["Mutation.AA"]], -1, -1)}) 
cosmic_missense$Mutation.AA.resnum <- apply(cosmic_missense, 1, function(x){as.numeric(str_sub(x[["Mutation.AA"]], 4, -2))}) 

aa_list <- c("R", "Q", "K", "V", "C", "Y", "D", "P", "I", "L", "A", "N", "H", "M", "T", "F", "W", "S", "E", "G", "Z")
cosmic_missense <- subset(cosmic_missense, !is.na(Mutation.AA.resnum) & Mutation.AA.Before %in% aa_list & Mutation.AA.After %in% aa_list)

cosmic_missense$Mutation.AA.Before <- factor(cosmic_missense$Mutation.AA.Before)
cosmic_missense$Mutation.AA.After <- factor(cosmic_missense$Mutation.AA.After)

# only for adenocarcinoma
dist_each_gene <- merge(subset(cosmic_missense, Histology.subtype == "neoplasm"), subset(cosmic_missense, Histology.subtype == "neoplasm"), by = "Gene.name")

getdist <- function(dist_each_uniprot) {
  dist_each_uniprot$dist <- dist_each_uniprot$Mutation.AA.resnum.x - dist_each_uniprot$Mutation.AA.resnum.y
  dist_each_uniprot <- subset(dist_each_uniprot, dist > 0)
  dist_each_uniprot$var_type_c <- apply(dist_each_uniprot, 1, function(x){paste(sort(c(x[["FATHMM.prediction.x"]],x[["FATHMM.prediction.y"]])), collapse = " to ")})
  dist_each_uniprot$var_type_c <- factor(dist_each_uniprot$var_type_c)
  dist_each_uniprot
}
dist_each_gene <- getdist(dist_each_gene)

# The mutation distribution between disease nsSNP and polymorphism nsSNP
require(ggplot2)
ggplot(subset(dist_each_gene, var_type_c %in% c("CANCER to CANCER", "PASSENGER/OTHER to PASSENGER/OTHER")), aes(x = dist, color = var_type_c)) + scale_x_log10() + geom_density()
ggplot(subset(dist_each_gene, var_type_c %in% c("CANCER to CANCER","CANCER to PASSENGER/OTHER", "PASSENGER/OTHER to PASSENGER/OTHER")), aes(x = dist, color = var_type_c)) + scale_x_log10() + geom_density() + theme(legend.position = "top", legend.text = element_text(size = 8))
ggsave("snp_dist_cosmic.pdf")



###################################################################
###########         dual mutation frequency distribution ###########
###################################################################

# 1. I would suggest p53, Ras (HRas, NRas, and KRas), and EGFR to start with as they were involved in lots of distinct cancers with known structural oncogenic mechanism. While Ras and EGFR give gain-of-function (GOF) mutations, p53 has both the GOF and loss-of-function (LOF) mutations.
# 2. by looking for some stand-out frequency of dual mutations from individual cancer patient mutation data in COSMIC (or other database) it might give us some hint of latent mutations in p53, Ras, and EGFR.


# My understanding: 1. the dual mutations is for individual cancer patient. The cancer can be any cancer. But, I could probably start with some specific cancer
#                   2. start with specific gene like p53, Ras (HRas, NRas, and KRas), and EGFR

cosmic <- read.table("/users/ajing/SNPDist/Data/CosmicCompleteExport.tsv", sep = "\t", header = T, quote = "", na.string = "")

# begin with p53
cosmic_p53 <- subset(cosmic, Gene.name %in% c("TP53", "TP53_ENST00000269305", "TP53_ENST00000413465", "TP53_ENST00000414315", "TP53_ENST00000420246", "TP53_ENST00000455263", "TP53_ENST00000545858"))

save(cosmic_p53, file = "cosmic_p53.RData")

# the mutation could be occuring in gene level or amino acid level
# In this case, for different trainscripts, I will count it as different mutations.
p53_join <- merge(subset(cosmic_p53, !is.na(Mutation.CDS) | !is.na(Mutation.AA)), subset(cosmic_p53, !is.na(Mutation.CDS) | !is.na(Mutation.AA)), by = c("ID_sample", "ID_STUDY", "Gene.name"))

library(dplyr)
get_summary_data <- function(p53_join) {
p53_join <- subset(p53_join, Mutation.ID.x != Mutation.ID.y)
p53_join$Mutation.CDS.C <- apply(p53_join, 1, function(x){paste(sort(c(x[["Mutation.CDS.x"]],x[["Mutation.CDS.y"]])), collapse = " and ")})
p53_join$Mutation.AA.C <- apply(p53_join, 1, function(x){paste(sort(c(x[["Mutation.AA.x"]],x[["Mutation.AA.y"]])), collapse = " and ")})
p53_join$Mutation.CDS.C <- factor(p53_join$Mutation.CDS.C)
p53_join$Mutation.AA.C <- factor(p53_join$Mutation.AA.C)

p53_aa_sum <- subset(p53_join, Mutation.AA.x != Mutation.AA.y, select = Mutation.AA.C) %>% 
              group_by(Mutation.AA.C) %>%
              summarise(total_num = length(Mutation.AA.C))
p53_aa_sum$Mutation.AA.C = reorder(p53_aa_sum$Mutation.AA.C, -p53_aa_sum$total_num)

p53_cds_sum <- subset(p53_join, Mutation.CDS.x != Mutation.CDS.y, select = Mutation.CDS.C) %>% 
  group_by(Mutation.CDS.C) %>%
  summarise(total_num = length(Mutation.CDS.C))
p53_cds_sum$Mutation.CDS.C = reorder(p53_cds_sum$Mutation.CDS.C, -p53_cds_sum$total_num)

list("p53_aa_sum" = p53_aa_sum, "p53_cds_sum" = p53_cds_sum)
}
s_result <- get_summary_data(p53_join)

p53_aa_sum <- s_result$p53_aa_sum
ggplot(data=subset(p53_aa_sum, total_num >= tail(sort(p53_aa_sum$total_num), 21)[1]), aes(x=Mutation.AA.C,y=total_num)) + geom_bar(stat="identity") + theme(axis.text.x = element_text( size=8, angle=30))
ggsave("sample_gene_aa.jpg")
p53_aa_sum <- s_result$p53_cds_sum
ggplot(data=subset(p53_cds_sum, total_num >= tail(sort(p53_cds_sum$total_num), 21)[1]), aes(x=Mutation.CDS.C,y=total_num)) + geom_bar(stat="identity") + theme(axis.text.x = element_text( size=8, angle=30))
ggsave("sample_gene_cds.jpg")

# If we limit to occurance on single transcript
# the mutation could be occuring in gene level or amino acid level
# In this case, for different trainscripts, I will count it as different mutations.
p53_join_t <- merge(subset(cosmic_p53, !is.na(Mutation.CDS) | !is.na(Mutation.AA)), subset(cosmic_p53, !is.na(Mutation.CDS) | !is.na(Mutation.AA)), by = c("ID_sample", "ID_STUDY", "Accession.Number"))

s_result <- get_summary_data(p53_join_t)

p53_aa_sum <- s_result$p53_aa_sum
ggplot(data=subset(p53_aa_sum, total_num >= tail(sort(p53_aa_sum$total_num), 21)[1]), aes(x=Mutation.AA.C,y=total_num)) + geom_bar(stat="identity") + theme(axis.text.x = element_text( size=8, angle=30))
ggsave("sample_transcript_aa.jpg")
p53_aa_sum <- s_result$p53_cds_sum
ggplot(data=subset(p53_cds_sum, total_num >= tail(sort(p53_cds_sum$total_num), 21)[1]), aes(x=Mutation.CDS.C,y=total_num)) + geom_bar(stat="identity") + theme(axis.text.x = element_text( size=8, angle=30))
ggsave("sample_transcript_cds.jpg")


# This is for each gene, because they have different name of Gene.name for transcript also. So, this part is more correct for gene mutation distribution
p53_join_one <- merge(subset(cosmic_p53, !is.na(Mutation.CDS) | !is.na(Mutation.AA)), subset(cosmic_p53, !is.na(Mutation.CDS) | !is.na(Mutation.AA)), by = c("ID_sample", "ID_STUDY"))

s_result <- get_summary_data(p53_join_one)

p53_aa_sum <- s_result$p53_aa_sum
ggplot(data=subset(p53_aa_sum, total_num >= tail(sort(p53_aa_sum$total_num), 21)[1]), aes(x=Mutation.AA.C,y=total_num)) + geom_bar(stat="identity") + theme(axis.text.x = element_text( size=8, angle=30))
ggsave("sample_onegene_aa.jpg")
p53_cds_sum <- s_result$p53_cds_sum
ggplot(data=subset(p53_cds_sum, total_num >= tail(sort(p53_cds_sum$total_num), 21)[1]), aes(x=Mutation.CDS.C,y=total_num)) + geom_bar(stat="identity") + theme(axis.text.x = element_text( size=8, angle=30))
ggsave("sample_onegene_cds.jpg")


####################################################################################
##########   Why no consistency between transcript level and gene level############
####################################################################################
result <- with(subset(cosmic_p53, Mutation.AA %in% c("p.R175H", "p.R43H")), table(factor(Accession.Number), factor(Mutation.AA)))
result <- with(subset(cosmic_p53, Mutation.AA %in% c("p.R175H", "p.R43H", "p.R82H")), table(factor(Accession.Number), factor(Mutation.AA)))

save.image("p53.RData")

write.csv(result, file = "tmp.csv")


##############################################################################################
###########    dual mutation frequency distribution for only one transcript and    ###########
##############################################################################################
# If yes, we should look for dual mutations within each transcript instead of mixed transcripts in searching for latent mutations. Also, statistically, the significance of dual mutations should be normalized by the occurrence of individual mutation.

# The number of co-occured mutations at each transcript for each patient
transcript_co_mut_num <- cosmic_missense %>% group_by(ID_sample, ID_STUDY, Accession.Number) %>% summarise(total_num = length(Mutation.AA))

ggplot(transcript_co_mut_num, aes(x = total_num)) + geom_histogram()
ggsave("mut_dist_co_trans.jpg")

# distribution of occurrence for individual mutation
transcript_mut_num <- cosmic_missense %>% group_by(ID_sample, ID_STUDY, Accession.Number, Mutation.AA) %>% summarise(total_num = n())

ggplot(transcript_mut_num, aes(x = total_num)) + geom_histogram()
ggsave("mut_dist_trans.jpg")

## for just p53
mut_count_p53 <- subset(cosmic_p53, !is.na(Mutation.CDS) | !is.na(Mutation.AA)) %>% group_by(Accession.Number, Mutation.AA) %>% summarise(mut_count = n())

ggplot(mut_count_p53, aes(x = mut_count)) + geom_histogram()  + scale_x_log10()
ggsave("mut_dist_p53.jpg")

# The dual mutation for each transcript divided by min of occurrence of individual mutation

# recalculaet the data, because some redundancy involved for HGNC column
cosmic_p53_nr <- subset(cosmic_p53, !is.na(Mutation.CDS) | !is.na(Mutation.AA)) %>% distinct(ID_sample, ID_STUDY, Accession.Number, Mutation.AA)

p53_join_t <- merge(cosmic_p53_nr, cosmic_p53_nr, by = c("ID_sample", "ID_STUDY", "Accession.Number"))

## For normalization
p53_join_norm <- merge(p53_join_t, mut_count_p53, by.x = c("Accession.Number", "Mutation.AA.x"), by.y = c("Accession.Number", "Mutation.AA"))
p53_join_norm <- merge(p53_join_norm, mut_count_p53, by.x = c("Accession.Number", "Mutation.AA.y"), by.y = c("Accession.Number", "Mutation.AA"))
p53_join_norm$mut_count.min <- do.call(pmin, subset(p53_join_norm, select = c(mut_count.x, mut_count.y)))

get_summary_data_norm <- function(p53_join) {
  p53_join <- subset(p53_join, Mutation.ID.x != Mutation.ID.y)
  p53_join$Mutation.CDS.C <- apply(p53_join, 1, function(x){paste(sort(c(x[["Mutation.CDS.x"]],x[["Mutation.CDS.y"]])), collapse = " and ")})
  p53_join$Mutation.AA.C <- apply(p53_join, 1, function(x){paste(sort(c(x[["Mutation.AA.x"]],x[["Mutation.AA.y"]])), collapse = " and ")})
  p53_join$Mutation.CDS.C <- factor(p53_join$Mutation.CDS.C)
  p53_join$Mutation.AA.C <- factor(p53_join$Mutation.AA.C)
  
  p53_aa_sum <- subset(p53_join, Mutation.AA.x != Mutation.AA.y, select = c(Accession.Number, Mutation.AA.C, mut_count.min)) %>% 
    group_by(Accession.Number, Mutation.AA.C) %>%
    select(mut_count.min) %>%
    summarise(total_score =  n() / (2 * first(mut_count.min)), mut_count_min = first(mut_count.min))
  p53_aa_sum$Mutation.AA.C = reorder(p53_aa_sum$Mutation.AA.C, -p53_aa_sum$total_score)
  
#  p53_cds_sum <- subset(p53_join, Mutation.CDS.x != Mutation.CDS.y, select = Mutation.CDS.C) %>% 
#    group_by(Mutation.CDS.C) %>%
#    summarise(total_num = length(Mutation.CDS.C) / )
#  p53_cds_sum$Mutation.CDS.C = reorder(p53_cds_sum$Mutation.CDS.C, -p53_cds_sum$total_num)
  
#  list("p53_aa_sum" = p53_aa_sum, "p53_cds_sum" = p53_cds_sum)

  list("p53_aa_sum" = p53_aa_sum)
}
s_result <- get_summary_data_norm(p53_join_norm)

p53_aa_sum <- merge(s_result$p53_aa_sum, mut_count_p53)

ggplot(data=subset(s_result$p53_aa_sum, mut_count_min > 1 & total_score > 0.5), aes(x=Mutation.AA.C,y=total_score)) + geom_bar(stat="identity") + theme(axis.text.x = element_text( size=8, angle=30)) + facet_grid(Accession.Number ~.)
ggsave("transcript_aa_norm.jpg")


##############################################################################################
###########    dual mutation frequency impact on protein structures/functions      ###########
##############################################################################################
subset(p53_cds_sum, total_num >= tail(sort(p53_cds_sum$total_num), 21)[1])
