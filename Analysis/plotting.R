#
#
#

polybind <- read.table("./Data/polybind.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
colnames(polybind) = c("PDBID", "chainID", "resnum", "UniProtID", "uniprot_resnum", "VarType", "DiseaseName", "ligName", "distance")

paste(subset(polybind, UniProtID == "P00918" & PDBID == "4MO8" & distance < 4)$resnum, collapse = "+")


length(unique(subset(protein_annotate_withsnp, VarType %in% c("Polymorphism", "Disease"))$UniProtID))
length(unique(subset(protein_annotate_withsnp, VarType %in% c("Polymorphism", "Unclassified"))$UniProtID))
length(unique(subset(protein_annotate_withsnp, VarType %in% c("Unclassified", "Disease"))$UniProtID))

library(dplyr)
# get the proteins with all 3 var type

disease =  245
polymor =  580
unclass =  149
dis_pol =  disease + polymor - length(unique(subset(protein_annotate_withsnp, VarType %in% c("Polymorphism", "Disease"))$UniProtID))
pol_unc =  polymor + unclass - length(unique(subset(protein_annotate_withsnp, VarType %in% c("Polymorphism", "Unclassified"))$UniProtID))
dis_unc =  disease + unclass - length(unique(subset(protein_annotate_withsnp, VarType %in% c("Unclassified", "Disease"))$UniProtID))
dis_unc_poly = 41

install.packages("VennDiagram")
library(VennDiagram)

pdf("vennplot.pdf")
venn.plot <- draw.triple.venn(disease, polymor, unclass, dis_pol, pol_unc, dis_unc, dis_unc_poly, category = c("Disease", "Polymorphism", "Unclassified"), fill = c("red", "blue", "orange"), cex = 2, cat.cex=2, cat.pos = c(-40, 32, 180), cat.col=c("red", "blue", "orange"))
dev.off()


114 / length(unique(subset(protein_annotate_withsnp, VarType %in% c("Disease"))$UniProtID))