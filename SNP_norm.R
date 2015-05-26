# Script for normalized SNP data set

protein_annotate<- read.table("./Data/SNPOnStruct_final.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
protein_annotate<- read.table("./Data/SNPOnStruct_final_3.txt", sep = "\t", header = F, quote = "", na.string = "\\N")

#
colnames(protein_annotate)<-c("pdbid", "biounit", "ModelID", "chainid", "resnam", "resnum", "structcode", "ssa", "rsa", "resolution", "UniProtID", "uniprot_resnum", "BindingSiteComment", "distance", "ligandName", "genename", "FTID", "AABefore", "AAAfter", "VarType", "dbSNPID", "DiseaseName", "AABeforeProp", "AAAfterProp", "proteinname", "reviewed", "interproname")

protein_annotate$location = "Core"

#MAX_ACC["GLY"] = 84.0
protein_annotate$location[protein_annotate$ssa > 5] = "Surface"

protein_annotate[ which(protein_annotate$distance <= 4.0 | !is.na(protein_annotate$BindingSiteComment)),'location'] = "Binding Site"

table(protein_annotate$VarType, protein_annotate$location)

# The first table
get_stat_eachtype <- function(protein_annotate, snp_type){
  allres_table <- table(protein_annotate$location)
  loc_table = table(subset(protein_annotate, VarType == snp_type)$location)
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

odds_ratio_stat <- function(protein_annotate, vartype){
table_res <- table(deal_with_na(protein_annotate$VarType == vartype), deal_with_na(protein_annotate$location == "Core"))
table_res <- apply(table_res, 1:2, as.numeric)
rownames(table_res) <- c(paste("Not", vartype), vartype)
colnames(table_res) <- c("Not Core", "Core")
print(oddsratio.wald(table_res))

table_res <- table(deal_with_na(subset(protein_annotate, location %in% c("Surface", "Binding Site"))$VarType == vartype), deal_with_na(subset(protein_annotate, location %in% c("Surface", "Binding Site"))$location == "Surface"))
table_res <- apply(table_res, 1:2, as.numeric)
rownames(table_res) <- c(paste("Not", vartype), vartype)
colnames(table_res) <- c("Binding Site", "Surface")
print(oddsratio.wald(table_res))

table_res <- table(deal_with_na(subset(protein_annotate, location %in% c("Core", "Binding Site"))$VarType == vartype), deal_with_na(subset(protein_annotate, location %in% c("Core", "Binding Site"))$location == "Core"))
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
odds_ratio_stat(subset(protein_annotate, !grepl("*antigen*", x = proteinname)), "Polymorphism")

# Unclassified
odds_ratio_stat(protein_annotate, "Unclassified")

# Fisher's Exact Test
# examine the significance of the association (contingency) between the two kinds of classification.
# null hypothesis that disease and non-disease are equally likely to be in the structure core,
fish_result = fisher.test(table(DF$Admit, DF$Gender))
fish_result$p.value
fish_result$estimate
fish_result$conf.int

# For each amino acid
amino_names <- unique(DF$Gender)
ratio_amino <- data.frame(amino_names)
for (each in amino_names){
  fish_result = fisher.test(table(DF$Admit, DF$Gender))

}

# plots odds ratio with error bar
ggplot(dfwc.between, aes(x=condition, y=value, group=1)) +
  geom_bar() +
  geom_errorbar(width=.1, aes(ymin=value-ci, ymax=value+ci), colour="red")


library(ggplot2)

# http://stackoverflow.com/questions/13386177/how-to-create-odds-ratio-and-95-ci-plot-in-r
# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_%28ggplot2%29/
# http://stackoverflow.com/questions/14069629/plotting-confidence-intervals
ggplot(dat, aes(x = pollut, y = or, ymin = lcl, ymax = ucl)) + geom_pointrange(aes(col = factor(lag)), position=position_dodge(width=0.30)) + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("")

# Also odds ratio for each type of amino acids
# for disease and polymorphism seperately
ggplot(dat, aes(x = pollut, y = or, ymin = lcl, ymax = ucl)) + geom_pointrange(aes(col = factor(lag)), position=position_dodge(width=0.30)) + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("")

# odds ratio comparison between disease and polymorphism
ggplot(dat, aes(x = pollut, y = or, ymin = lcl, ymax = ucl)) + geom_pointrange(aes(col = factor(lag)), position=position_dodge(width=0.30)) + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("")


# Distribution of BLOSSUM Matrix
# 3x3 chis-sq test



