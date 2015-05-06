
############## data preparing for Brandon
install.packages("sqldf")
library(sqldf)
s01 <- sqldf("select pdbid, chainid from protein_annotate_onlysnp group by UniProtID")
write.table(s01, file = "unique_pdb.txt", quote = F, row.names = F)

# get uniprot with disease snp in binding site
uniprot_bs_ds <- sqldf('select UniProtID from protein_annotate_onlysnp where VarType = "Disease" and location = "Binding Site" group by UniProtID')

with(subset(protein_annotate_onlysnp, VarType == "Disease"), table(factor(UniProtID), location))

# odds ratio for each uniprot
odds_ratio_stat <- function(p_annotate, vartype){
  table_res <- table(data.frame(location = deal_with_na(subset(p_annotate, location %in% c("Surface", "Binding Site"))$VarType == vartype),
                     Disease = deal_with_na(subset(p_annotate, location %in% c("Surface", "Binding Site"))$location == "Surface")))
  #table_res <- apply(table_res, 1:2, as.numeric)
 # rownames(table_res) <- c(paste("Not", vartype), vartype)
  #colnames(table_res) <- c("Binding Site", "Surface")
  #print(table_res)
  fisher.test(table_res)
}

snp_uniprot = ddply(subset(protein_annotate_withsnp, UniProtID %in% uniprot_bs_ds[,1]), .(UniProtID), summarise, pvalue = odds_ratio_stat(data.frame(location = location, VarType = VarType), "Disease")$p.value, oddsratio = 1/ odds_ratio_stat(data.frame(location = location, VarType = VarType), "Disease")$estimate)

snp_uniprot[order(snp_uniprot$pvalue),]
# test whether this is right
with(subset(protein_annotate_withsnp, UniProtID == "P01112"), table(location, VarType))

VarType
location       Disease Polymorphism Unclassified
Binding Site      16            0            2
Core               0            0            0
Surface            1            0            0


with(subset(protein_annotate_withsnp, UniProtID == "P38398"), table(location, VarType))

VarType
location       Disease Polymorphism Unclassified
Binding Site       4            0            2
Core               2            3           15
Surface            0            1           22

# top 40
snp_uniprot_t40 <- snp_uniprot[order(snp_uniprot$pvalue),][1:40,]
# 
s01 <- sqldf("select pdbid, chainid, snp_uniprot_t40.UniProtID, snp_uniprot_t40.pvalue from protein_annotate_withsnp inner join snp_uniprot_t40 on protein_annotate_withsnp.UniProtID = snp_uniprot_t40.UniProtID group by protein_annotate_withsnp.UniProtID")
write.table(s01, file = "top_40_pdb.txt", quote = F, row.names = F)


