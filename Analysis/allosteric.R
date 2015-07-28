# allosteric
allo_annotate<- read.table("./Data/SNPOnStruct_biolip.txt", sep = "\t", header = F, quote = "", na.string = "\\N")

colnames(allo_annotate)<-c("pdbid", "biounit", "ModelID", "chainid", "resnam", "resnum", "structcode", "ssa", "rsa", "UniProtID", "uniprot_resnam", "uniprot_resnum", "bs_biounit", "bs_p_chainid", "bs_p_resnum", "ligandName", "allo_pubmedId", "allo_pubmedTitle", "allo_organism", "allo_bsid", "genename", "SwissProt_AC", "FTID", "AABefore", "prot_resnum", "AAAfter", "VarType", "dbSNPID", "DiseaseName", "uniprot_resnam_3d", "AABeforeProp", "AAAfterProp", "proteinname", "reviewed", "gene_name_acc", "interproname_acc")

allo_annotate$location = "Core"

#MAX_ACC["GLY"] = 84.0
allo_annotate$location[allo_annotate$ssa > 5] = "Surface"

allo_annotate[ which(allo_annotate$distance <= 4.0), 'location'] = "Binding Site"

allo_annotate[ which(!is.na(allo_annotate$allo_pubmedId)), 'location'] = "Allo"
# For allosteric site
allo_annotate[]

table(allo_annotate$VarType, allo_annotate$location)

# only keep UniProt with at least one nsSNP
snp_num = ddply(allo_annotate, .(UniProtID), summarise, N = length(unique(uniprot_resnum[!is.na(VarType)])))
allo_annotate_withsnp = subset(allo_annotate, UniProtID %in% factor(unique(subset(allo_annotate, !is.na(VarType))$UniProtID)))
length(unique(allo_annotate_withsnp$UniProtID))

> table(allo_annotate_withsnp$location)


length(unique(allo_annotate[allo_annotate$VarType == "Unclassified",]$UniProtID))

length(unique(allo_annotate[allo_annotate$VarType == "Disease",]$UniProtID))

length(unique(allo_annotate[allo_annotate$VarType == "Polymorphism",]$UniProtID))

# allosteric sites or other sites
fish_bs <- function(p_annotate_bs, vartype, loc){
  fish_result = fisher.test(table(data.frame(Disease = p_annotate_bs$VarType == vartype, Allosite = p_annotate_bs$allosite == loc)))
  fish_result
}
fish_bs(subset(allo_site, allosite %in% c("Binding Site", "Allo")), "Disease", "Allo")
fish_bs(subset(allo_site, allosite %in% c("Binding Site", "Allo")), "Polymorphism", "Allo")
fish_bs(subset(allo_site, allosite %in% c("Binding Site", "Allo")), "Unclassified", "Allo")


# comparing allosteric sites and other ligand binding sites for each protein
each_uniprot = by(subset(allo_site, allosite %in% c("Binding Site", "Allo")), subset(allo_site, allosite %in% c("Binding Site", "Allo"))[, "UniProtID"], function(x) {fish_bs(x, "Disease", "Allo")})


snp_num = ddply(allo_site, .(UniProtID), summarise, fish_odds = fisher.test(table(data.frame(Disease = VarType, Allosite = allosite)))$estimate, fish_pvalue = fisher.test(table(data.frame(Disease = VarType, Allosite = allosite)))$p.value)
