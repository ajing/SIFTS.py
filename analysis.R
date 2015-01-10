surface <- read.table("./Result/avedist2surface.txt_small", sep = "\t", header = F)

colnames(surface)<-c("pdbid", "chainid", "resnum", "resnam", "surf_dist", "surf_dist_alpha")

surface <- melt(surface, id = c("pdbid","chainid","resnam","surf_dist","sur_dist_alpha"))
library("reshape2")
surface <- melt(surface, id = c("pdbid","chainid","resnum", "resnam"))
qplot(value, data = surface, color = variable, geom="density", log="x")
newsurf<-subset(surface, value != "None")
subset(newsurf, is.na(as.numeric(value)))
surface <- newsurf


qplot(value, data = surface, geom="density", log="x") + facet_wrap(resnam~.)
aanam <- c('ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL')
qplot(value, data = subset(surface,resnam %in% aanam[1:11]), color = resnam, geom="density") + coord_cartesian(xlim = c(1,10))
qplot(value, data = subset(surface,resnam %in% aanam[12:22]), color = resnam, geom="density")  + coord_trans(x="log2")+ coord_cartesian(xlim = c(1,10))


library("scales")
qplot(value, data = subset(surface,resnam %in% aanam[12:22]), color = resnam, geom="density")  + scale_x_continuous(trans=log2_trans())+ coord_cartesian(xlim = c(1,10))
snp_dist <-read.csv("./Result/snp_dist.txt",sep="\t", header=F)
colnames(snp_dist)<-c("pdbid","snpid1","chainid1","resnam1","resnum1","snpid2","chainid2","resnam2","resnum2", "spa_dist","seq_dist")

snp_dist_nofalse <- subset(snp_dist, seq_dist != 'False')
snp_dist_nofalse[,"seq_dist"] <- as.numeric(levels(snp_dist_nofalse$seq_dist)[snp_dist_nofalse$seq_dist])
ggplot(snp_dist_nofalse, aes(x=spa_dist, y=seq_dist)) + geom_bin2d()

humsavar <- read.table("./Data/humsavar_tru.txt.out", sep="\t",header=F, quote="")
humsavar <- subset(humsavar, V6 !="-")
humsavar <- humsavar[,c("V5","V6")]
colnames(humsavar)<-c("type", "snpid")

snp_dist_nof_var<-merge(x = snp_dist_nofalse, y = humsavar, by.x = "snpid1", by.y="snpid", all.x=TRUE)
colnames(snp_dist_nof_var)[12]<- "type1"
snp_dist_nof_var<-merge(x = snp_dist_nof_var, y = humsavar, by.x = "snpid2", by.y="snpid", all.x=TRUE)
colnames(snp_dist_nof_var)[13]<- "type2"
ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "unclassified"), aes(x=spa_dist, y=seq_dist)) + geom_bin2d() + facet_grid(type1 ~.)
ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified"), aes(x=spa_dist, y=seq_dist)) + stat_density2d(geom="tile",aes(fill=..density..), contour=FALSE) + facet_grid(type1 + type2 ~.)

snp_dist_nof_var$type = paste(snp_dist_nof_var$type1,snp_dist_nof_var$type2)
snp_dist_nof_var[snp_dist_nof_var$type1 != snp_dist_nof_var$type2 & snp_dist_nof_var$type1 != "Unclassified" & snp_dist_nof_var$type2 != "Unclassified","type"] <- "Polymorephism Disease"
snp_dist_nof_var[snp_dist_nof_var$type1 != snp_dist_nof_var$type2 & snp_dist_nof_var$type1 != "Polymorphism" & snp_dist_nof_var$type2 != "Polymorphism","type"] <- "Unclassified Disease"
snp_dist_nof_var[snp_dist_nof_var$type1 != snp_dist_nof_var$type2 & snp_dist_nof_var$type1 != "Disease" & snp_dist_nof_var$type2 != "Disease","type"] <- "Polymorphism Unclassified"
snp_dist_nof_var$type = as.factor(snp_dist_nof_var$type)

ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified"), aes(x=spa_dist, y=seq_dist)) + stat_density2d(geom="tile",aes(fill=..density..), contour=FALSE) +  scale_fill_gradient(low = "white", high = "blue") + facet_grid(type ~.) + coord_cartesian(xlim = c(0, 50))
ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified"), aes(x=spa_dist, y=seq_dist))  + facet_grid(type ~.) + stat_density2d(geom="tile",aes(fill=..density..), contour=FALSE) +  scale_fill_gradient(low = "white", high = "blue")+ scale_x_continuous(limits = c(0, 50), expand=c(0,0)) + scale_y_continuous(limits = c(0,200), expand=c(0,0))

ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified" & type %in% c("Polymorphism Polymorphism")), aes(x=spa_dist, y=seq_dist))  + facet_grid(type ~.) + stat_density2d(geom="tile",aes(fill=..density..), contour=FALSE) +  scale_fill_gradient(low = "white", high = "blue")+ scale_x_continuous(limits = c(0, 50), expand=c(0,0)) + scale_y_continuous(limits=c(0,200), expand=c(0,0))

dist_bs <- read.table("./Result/distaa_biolip.txt_filtered", sep = "\t", header =F)
colnames(dist_bs)<-c("pdbid","bs_chainid","bsid","resnam","res_chainid","resnum", "distance", "snpid")
dist_bs$pdbid<-toupper(dist_bs$pdbid)

disease_spot_melt_dist<-merge(x = disease_spot_melt, y = dist_bs[,c("pdbid","snpid","res_chainid", "resnam","resnum","distance")], by.x = c("pdbid","snpid","chainid", "resnam","resnum"), by.y = c("pdbid","snpid","res_chainid", "resnam","resnum"), all.x=TRUE)

c_disease_spot <- ddply(subset(disease_spot_melt_dist, !is.na(distance)), "cond", spot, rating.mean=mean(distance))
qplot(distance, data=subset(disease_spot_melt_dist, !is.na(distance)), geom="density", facets=spot~., alpha=I(0.3)) + geom_vline(data=cdf, aes(xintercept=rating.mean,  colour=cond),linetype="dashed", size=1)+ xlim(c(0,100)) 

pdf("eachPDB.pdf")
ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified"), aes(x=spa_dist, y=seq_dist))  + facet_grid(type ~ pdbid) + stat_density2d(geom="tile",aes(fill=..density..), contour=FALSE) +  scale_fill_gradient(low = "white", high = "blue")+ scale_x_continuous(limits = c(0, 50), expand=c(0,0)) + scale_y_continuous(limits=c(0,200), expand=c(0,0))
dev.off()

for (pdb in snp_dist_nof_var$pdbid) {
  filename = paste(pdb,".pdb", sep="")
  pdf(filename)
  ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified" & pdbid == pdb), aes(x=spa_dist, y=seq_dist))  + facet_grid(type ~.) + stat_density2d(geom="tile",aes(fill=..density..), contour=FALSE) +  scale_fill_gradient(low = "white", high = "blue")+ scale_x_continuous(limits = c(0, 50), expand=c(0,0)) + scale_y_continuous(limits=c(0,200), expand=c(0,0))
  dev.off()
}
