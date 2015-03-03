# Script for normalized SNP data set

protein_annotate<- read.table("./Data/SNPOnStruct_final.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
# with aabefore and aaafter
protein_annotate<- read.table("./Data/SNPOnStruct_final_3.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
# remove redundancy in gene level.
protein_annotate<- read.table("./Data/SNPOnStruct_final_7.txt", sep = "\t", header = F, quote = "", na.string = "\\N")


#
colnames(protein_annotate)<-c("pdbid", "biounit", "ModelID", "chainid", "resnam", "resnum", "structcode", "ssa", "rsa", "UniProtID", "uniprot_resnam", "uniprot_resnum", "bs_biounit", "bs_p_chainid", "bs_p_resnum", "ligandName", "BindingSiteComment", "distance", "genename", "SwissProt_AC", "FTID", "AABefore", "prot_resnum", "AAAfter", "VarType", "dbSNPID", "DiseaseName", "AABeforeProp", "AAAfterProp", "proteinname", "reviewed", "gene_name_acc", "interproname_acc")

protein_annotate$location = "Core"

#MAX_ACC["GLY"] = 84.0
protein_annotate$location[protein_annotate$ssa > 5] = "Surface"

protein_annotate[ which(protein_annotate$distance <= 4.0 | !is.na(protein_annotate$BindingSiteComment)),'location'] = "Binding Site"

table(protein_annotate$VarType, protein_annotate$location)


# only keep UniProt with at least one nsSNP
snp_num = ddply(protein_annotate, .(UniProtID), summarise, N = length(unique(uniprot_resnum[!is.na(VarType)])))
protein_annotate_withsnp = subset(protein_annotate, UniProtID %in% factor(unique(subset(protein_annotate, !is.na(VarType))$UniProtID)))
length(unique(protein_annotate_withsnp$UniProtID))

> table(protein_annotate_withsnp$location)

Binding Site         Core      Surface 
17010        49292       154993 

> length(unique(protein_annotate[protein_annotate$VarType == "Unclassified",]$UniProtID))
[1] 150
> length(unique(protein_annotate[protein_annotate$VarType == "Disease",]$UniProtID))
[1] 246
> length(unique(protein_annotate[protein_annotate$VarType == "Polymorphism",]$UniProtID))
[1] 581


> dim(unique(protein_annotate[protein_annotate$VarType == "Polymorphism",][,c("UniProtID", "uniprot_r
esnum")]))
[1] 1656    2
> dim(unique(protein_annotate[protein_annotate$VarType == "Disease",][,c("UniProtID", "uniprot_resnum
")]))
[1] 3006    2
> dim(unique(protein_annotate[protein_annotate$VarType == "Unclassified",][,c("UniProtID", "uniprot_r
esnum")]))
[1] 813   2


> dim(unique(protein_annotate[protein_annotate$VarType == "Unclassified",]))
[1] 1766   34
> dim(unique(protein_annotate[protein_annotate$VarType == "Polymorphism",]))
[1] 1717   34
> dim(unique(protein_annotate[protein_annotate$VarType == "Disease",]))
[1] 3709   34



# The first table
get_stat_eachtype <- function(p_annotate, snp_type){
  allres_table <- table(p_annotate$location)
  loc_table = table(subset(p_annotate, VarType == snp_type)$location)
  print("Observed")
  print(loc_table)
  print(sum(loc_table))
  print("Expected:")
  exp_table = allres_table / sum(allres_table) * sum(loc_table)
  print(exp_table)
  print("O/E")
  print(loc_table/exp_table)
}

ggplot(protein_annotate, aes(x=ssa)) +
  geom_histogram() + facet_grid(.~location)

ggplot(subset(protein_annotate, location %in% c("Surface", "Binding Site")), aes(x=ssa, color = location)) +
  geom_density() + facet_grid(location~.)

ggsave(filename = "tmp.pdf")

summary(subset(protein_annotate, location == "Binding Site" & VarType =="Polymorphism"))

# odds ratio is to quantify how strongly the presence or absence of property A is associated with the presence or absense of property B in a given population.
# If the OR is greater than 1, then having "A" is considered to be "associated" with having "B" in the sense that the having of "B" raises (relative to not-having "B") the odds of having "A".

library("fmsb")
result = oddsratio(, conf.level = 0.95)
result$p.value
result$estimate
result$conf.int

install.packages(epitools)
library("epitools")

deal_with_na <- function(x) {
  x[is.na(x)] = F
  factor(x)
}

odds_ratio_stat <- function(p_annotate, vartype){
table_res <- table(deal_with_na(p_annotate$VarType == vartype), deal_with_na(p_annotate$location == "Core"))
table_res <- apply(table_res, 1:2, as.numeric)
rownames(table_res) <- c(paste("Not", vartype), vartype)
colnames(table_res) <- c("Not Core", "Core")
print(oddsratio.wald(table_res))

table_res <- table(deal_with_na(subset(p_annotate, location %in% c("Surface", "Binding Site"))$VarType == vartype),
        deal_with_na(subset(p_annotate, location %in% c("Surface", "Binding Site"))$location == "Surface"))
table_res <- apply(table_res, 1:2, as.numeric)
rownames(table_res) <- c(paste("Not", vartype), vartype)
colnames(table_res) <- c("Binding Site", "Surface")
print(oddsratio.wald(table_res))

table_res <- table(deal_with_na(subset(p_annotate, location %in% c("Core", "Binding Site"))$VarType == vartype),
        deal_with_na(subset(p_annotate, location %in% c("Core", "Binding Site"))$location == "Core"))
table_res <- apply(table_res, 1:2, as.numeric)
rownames(table_res) <- c(paste("Not", vartype), vartype)
colnames(table_res) <- c("Binding Site", "Core")
print(oddsratio.wald(table_res))
}

# Disease
odds_ratio_stat(protein_annotate, "Disease")

# Polymorphism
odds_ratio_stat(protein_annotate, "Polymorphism")

# Polymorphism without antigen
odds_ratio_stat(subset(protein_annotate, !grepl("*histocompatibility antigen*", x = proteinname)), "Polymorphism")

# Unclassified
odds_ratio_stat(protein_annotate, "Unclassified")

# Fisher's Exact Test
# examine the significance of the association (contingency) between the two kinds of classification.
# null hypothesis that disease and non-disease are equally likely to be in the structure core,
fish_result = fisher.test(table(DF$Admit, DF$Gender))
fish_result$p.value
fish_result$estimate
fish_result$conf.int

# For amino acid mutated in binding site
# the subset I am interested in is subset(protein_annotate, 'location' == "Binding Site")
p_annotate_bs <- subset(protein_annotate, location == "Binding Site")
p_annotate_bs <- subset(protein_annotate, location %in% c("Surface", "Binding Site"))
amino_names <- factor(unique(p_annotate_bs$resnam)[1:20])
each = as.character(amino_names[1])
vartype = "Disease"
vartypes   = c("Disease", "Polymorphism", "Unclassified")
ratio_amino <- data.frame(merge(amino_names, vartypes, all = TRUE))
colnames(ratio_amino) <- c("amino_names", "vartypes")
for (vartype in vartypes) {
  for (each in amino_names){
    print(each)
    print(vartype)
    each = as.character(each)
    table_res <- table(deal_with_na(p_annotate_bs$VarType == vartype),
                       deal_with_na(as.character(p_annotate_bs$resnam) == each))
    table_res <- apply(table_res, 1:2, as.numeric)
    rownames(table_res) <- c(paste("Not", vartype), vartype)
    colnames(table_res) <- c(paste("Not", each), each)
    
    fish_result = fisher.test(table(data.frame(residueName = (p_annotate_bs$resnam == each), Disease = p_annotate_bs$VarType == vartype)))
    print(fish_result)
    print(ratio_amino$amino_names == each & ratio_amino$vartypes == vartype)
    ratio_amino[ratio_amino$amino_names == each & ratio_amino$vartypes == vartype,"mean"] = fish_result$estimate
    ratio_amino[ratio_amino$amino_names == each & ratio_amino$vartypes == vartype,"ci_low"] = fish_result$conf.int[1]
    ratio_amino[ratio_amino$amino_names == each & ratio_amino$vartypes == vartype,"ci_up"] = fish_result$conf.int[2]
    ratio_amino[ratio_amino$amino_names == each & ratio_amino$vartypes == vartype,"pvalue"] = fish_result$p.value
  }
}


# plots odds ratio with error bar
ggplot(dfwc.between, aes(x=condition, y=value, group=1)) +
  geom_bar() +
  geom_errorbar(width=.1, aes(ymin=value-ci, ymax=value+ci), colour="red")


library(ggplot2)

# http://stackoverflow.com/questions/13386177/how-to-create-odds-ratio-and-95-ci-plot-in-r
# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_%28ggplot2%29/
# http://stackoverflow.com/questions/14069629/plotting-confidence-intervals
ggplot(subset(ratio_amino, vartypes %in% c("Disease", "Polymorphism")), aes(x = amino_names, y = mean, ymin = ci_low, ymax = ci_up)) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.30))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10()
ggsave(filename = "tmp.pdf", height=9, width=12) 

# Also odds ratio for each type of amino acids
# for disease and polymorphism seperately
ggplot(dat, aes(x = pollut, y = or, ymin = lcl, ymax = ucl)) + geom_pointrange(aes(col = factor(lag)), position=position_dodge(width=0.30)) + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("")

# odds ratio comparison between disease and polymorphism
ggplot(dat, aes(x = pollut, y = or, ymin = lcl, ymax = ucl)) + geom_pointrange(aes(col = factor(lag)), position=position_dodge(width=0.30)) + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("")


# Distribution of BLOSSUM Matrix
# 3x3 chis-sq test



