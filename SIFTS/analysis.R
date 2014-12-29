filename = "../pair.txt"
#filename = "pair.txt"
xy = read.table(filename, sep = "\t")
x  = xy[,1]
y  = xy[,2]

install.packages("hexbin")

library(hexbin)
hbo <- hexbin(x,y, xlab = "spatial distance", ylab = "sequence distance")
plot(hbo)
jpeg("hbo.jpg")
plot(hbo)
dev.off()

hist(x)
histpercent("pairdist.jpg", x, "Pair Spatial Distance", "Distance in Angstrom", "Probability")

content <- data.frame
content <- read.table("order_join.txt", sep="\t", quote = "")
pie(content[,2], labels = content[,1])

distfile <- "../onlydist.txt_noinf"
distance <- read.table(distfile, header=F)

distance_dat <- data.matrix(distance)
distance_num <- as.numeric(distance[,"V1"])


histpercent <- function(filename, datavector, gtitle, xlabel, ylabel){
jpeg(filename)
disthist = hist(datavector)
disthist$counts = disthist$counts / sum(disthist$count)
plot(disthist, main = gtitle, ylab = ylabel, xlab = xlabel)
dev.off()
}

histpercent("bsdist_hist.jpg", distance[distance[,1]<250,1], "Distance to Binding Site", "Distance in Angstrom", "Probability")
histpercent("bsdist_hist.jpg", distance[,1], "Distance to Binding Site",  "Distance in Angstrom", "Probability")

# For distance between residues
filename = "../distaa.txt_filtered"
distdata = read.table(filename)

histpercent("bsdist_hist_allchains_snp.jpg", as.numeric(distdata[,"V7"]), "Distance to Binding Site", "Distance in Angstrom", "Probability")
histpercent("bsdist_hist_allchains_less10_snp.jpg", as.numeric(distdata[distdata[,"V7"]< 10,"V7"]), "Distance to Binding Site", "Distance in Angstrom", "Probability")
histpercent("bsdist_hist_allchains_less20_snp.jpg", as.numeric(distdata[distdata[,"V7"]< 20,"V7"]), "Distance to Binding Site", "Distance in Angstrom", "Probability")

filename = "../distaa.txt"
distdata = read.table(filename)
histpercent("bsdist_hist_allchains.jpg", as.numeric(distdata[,"V7"]), "Distance to Binding Site", "Distance in Angstrom", "Probability")
histpercent("bsdist_hist_allchains_less10.jpg", as.numeric(distdata[distdata[,"V7"]< 10,"V7"]), "Distance to Binding Site", "Distance in Angstrom", "Probability")
histpercent("bsdist_hist_allchains_less20.jpg", as.numeric(distdata[distdata[,"V7"]< 20,"V7"]), "Distance to Binding Site", "Distance in Angstrom", "Probability")


disease_snp <- read.table("../disease_snp_dist.txt", header=F)
disease_snp <- data.matrix(disease_snp)

datavector <- disease_snp
plot(disthist,  xlim=c(0,200))

histpercent("disease_dist.jpg", disease_snp, "Distance from Disease SNP to BindingSite", "Probability", "Distance in Angstrom")
histpercent("disease_dist_less10.jpg", disease_snp[disease_snp < 10], "Distance from Disease SNP to BindingSite", "Probability", "Distance in Angstrom")

library(ggplot2)
filename = "../distaa.txt"
distdata = read.table(filename)
all_dist = data.frame(distance = as.numeric(distdata[,"V7"]))

disease_snp <- read.table("../disease_snp_dist.txt", header=F)
disease_snp <- as.numeric(unlist(disease_snp))

polymore_snp <- read.table("/tmp/polymorephism_snp_dist.txt", header=F)
polymore_snp <- as.numeric(unlist(polymore_snp))

filename = "../distaa.txt_filtered"
distdata = read.table(filename)
snp_dist = data.frame(distance = as.numeric(distdata[,"V7"]))
disease_dist = data.frame(distance = disease_snp)
snp_dist$type = "all_snp"
disease_dist$type = "disease_snp"


all_dist$type = "all_pair"
jpeg("all_hist.jpg")
ggplot(rbind(all_dist, snp_dist,disease_dist), aes(distance , fill = type)) + geom_density(alpha = 0.2)
dev.off()


polymore_dist = data.frame(distance = polymore_snp)
polymore_dist$type = "polymorephism_snp"
jpeg("polymore_hist.jpg")
ggplot(rbind(snp_dist,polymore_dist), aes(distance , fill = type)) + geom_density(alpha = 0.2)
dev.off()

unclass_snp <- read.table("/tmp/unclassified_snp_dist.txt", header=F)
unclass_snp <- as.numeric(unlist(unclass_snp))
unclass_dist = data.frame(distance = unclass_snp)
unclass_dist$type = "unclassified_snp"
jpeg("unclass_hist.jpg")
ggplot(rbind(snp_dist,unclass_dist), aes(distance , fill = type)) + geom_density(alpha = 0.2)
dev.off()

# The first column is family size (here is the number of pairs.. so very big number) second is max sequence length third column is standard deviation of distance
dist_var <- read.table("../Data/dist_variance.txt", sep="\t", quote= "", header=F)
plot(dist_var[,1], dist_var[,3], xlab = "family size", ylab = "standard deviation of distance", col = dist_var[,2])

install.packages("plot3D")
library(plot3D)
scatter3D(dist_var[,1], dist_var[,2], dist_var[,3])
pdf("dist_sd_family_relation.pdf")
png("dist_sd_family_relation.png")
scatter2D(dist_var[,1], dist_var[,3], colvar=dist_var[,2], log = "x", xlab = "The number of distances for a family (to binding sites)", ylab = "standard deviation of distances", clab = "Sequence Length")
index = dist_var[,3] > 60 & dist_var[,1] / dist_var[,2] ** 2 > 0.1
text(dist_var[index,1], dist_var[index,3], dist_var[index,4], cex=0.7, pos=1, col="black")


# a figure for biolip for all SNP, disease SNP and polymorephism SNP
library(ggplot2)
filename = "../distaa_biolip.txt_filtered"
distdata = read.table(filename)
snp_dist = data.frame(distance = as.numeric(distdata[,"V7"]))
snp_dist$type = "all_snp"

disease_snp <- read.table("/tmp/disease_snp_dist_2.txt", header=F)
disease_snp <- as.numeric(unlist(disease_snp))
disease_dist = data.frame(distance = disease_snp)
disease_dist$type = "disease_snp"

polymore_snp <- read.table("/tmp/polymorephism_snp_dist_2.txt", header=F)
polymore_snp <- as.numeric(unlist(polymore_snp))
polymore_dist = data.frame(distance = polymore_snp)
polymore_dist$type = "polymorephism_snp"

jpeg("biolip_hist.jpg")
ggplot(rbind(snp_dist, disease_dist, polymore_dist), aes(distance , fill = type)) + geom_density(alpha = 0.2) + xlim(0,50)
dev.off()


# distance variance summary
install.packages("plot3D")
library(plot3D)
filename = "/tmp/dist_summary_2.txt"
dist_var = read.table(filename, sep="\t", header=F)
jpeg("dist_var_all.jpg")
scatter2D(dist_var[,1], dist_var[,2], colvar=dist_var[,3], log = "x", xlab = "Average Distance", ylab = "Standard Deviation of Distances", clab = "Family Size(PDBs)", pch=19)
dev.off()


filename = "/tmp/disease_dist_summary.txt"
disease_dist_var = read.table(filename, sep="\t", header=F)
disease_dist_var = cbind(disease_dist_var, type = 0)
filename = "/tmp/polym_dist_summary.txt"
poly_dist_var = read.table(filename, sep="\t", header=F)
poly_dist_var = cbind(poly_dist_var, type=1)
dist_var = rbind(disease_dist_var, poly_dist_var)
jpeg("disease_var_2.jpg")
scatter2D(dist_var[,1], dist_var[,2], col=as.numeric(lapply(dist_var[,5],  function(x){ifelse(x==0, 100, 2)})), log = "x", xlab = "Average Distance (Disease)", ylab = "Standard Deviation of Distances", clab = "Family Size(PDBs)", pch = as.numeric(lapply(dist_var[,5],  function(x){ifelse(x==0, 0, 18)})))
legend("topleft", legend=c("Disease", "Polymorephism"), col = c(100,2), pch = c(0,18))
dev.off()


filename = "/tmp/polym_dist_summary.txt"
poly_dist_var = read.table(filename, sep="\t", header=F)
jpeg("dist_var_polymore.jpg")
scatter2D(dist_var[,1], dist_var[,2], colvar=dist_var[,3], log = "x", xlab = "Average Distance (Polymorephism)", ylab = "Standard Deviation of Distances", clab = "Family Size(PDBs)" ,pch = 19)
dev.off()


filename = "/tmp/std_pair_4.txt"
dist_var = read.table(filename, sep="\t", header=F)
jpeg("dist_var_pair.jpg")
scatter2D(dist_var[,1], dist_var[,2], colvar=dist_var[,3], xlab = "Standard Deviation of Distance (Diseases)", ylab = "Standard Deviation of Distances (Polymorephism)", clab = "Family Size(PDBs)" ,pch = 19)
index = dist_var[,1] > 20 & dist_var[,2] < 7
text(dist_var[index,1], dist_var[index,2], dist_var[index,4], cex=0.7, pos=1, col="black")
dev.off()



## IMPORTANT
filename = "/tmp/std_pair_5.txt"
dist_var = read.table(filename, sep="\t", header=F)
jpeg("dist_mean_pair.jpg")
scatter2D(dist_var[,1], dist_var[,2], colvar=dist_var[,5], xlab = "Diseases", ylab = "Polymorephism", clab = "Family Size(PDBs)" ,pch = 19, main = "Average Distance")
abline(0,1)
index = dist_var[,1] < 20 & dist_var[,2] >50
text(dist_var[index,1], dist_var[index,2], dist_var[index,6], cex=0.7, pos=1, col="black")
dev.off()


filename = "/tmp/std_pair_5.txt"
dist_var = read.table(filename, sep="\t", header=F)
jpeg("dist_var_pair.jpg")
scatter2D(dist_var[,3], dist_var[,4], colvar=dist_var[,5], xlab = "Diseases", ylab = "Polymorephism", clab = "Family Size(PDBs)" ,pch = 19, main = "Standard Deviation of Distance")
abline(0,1)
index = dist_var[,3] < 20 & dist_var[,4] >50
text(dist_var[index,3], dist_var[index,4], dist_var[index,6], cex=0.7, pos=1, col="black")
dev.off()




library(ggplot2)

filename = "/tmp/std_pair_5.txt"
dist_var = read.table(filename, sep="\t", header=F)
index = dist_var[,1] < 20 & dist_var[,2] >50

interestingfig <- function(interest_idx){
filename = "/tmp/std_pair_5.txt"
dist_var = read.table(filename, sep="\t", header=F)
index = dist_var[,1] < 20 & dist_var[,2] >50
interested = dist_var[index,6][interest_idx]
filename = "/tmp/snp_dist_desc_2.txt"
dist_var = read.table(filename, sep="\t", quote = "", header=F)
index    = (as.character(dist_var[,5]) == as.character(interested))
dist_var = dist_var[index, c(1,2)]
all_dist = dist_var
all_dist$distance = dist_var[,1]
all_dist$type = dist_var[,2]
print(dim(all_dist))
print(colnames(all_dist))
atitle = paste("/users/ajing/SNPDist/", "interest_fam_", interest_idx, ".jpg", sep="")
jpeg(atitle)
ggplot(all_dist, aes(distance , fill = type)) + geom_density(alpha = 0.2) + labs(title=interested) + xlim(0,50)
dev.off()
}
interestingfig(1)
interestingfig(2)
interestingfig(3)
interestingfig(4)
interestingfig(5)
interestingfig(6)






library(MASS)  # in case it is not already loaded
filename = "/tmp/disease_dist_summary.txt"
disease_dist_var = read.table(filename, sep="\t", header=F)
X <- disease_dist_var[,c(1,2)]

filename = "/tmp/polym_dist_summary.txt"
poly_dist_var = read.table(filename, sep="\t", header=F)
X <- poly_dist_var[,c(1,2)]

## some pretty colors
library(RColorBrewer)
k <- 5
my.cols <- rev(brewer.pal(k, "RdYlBu"))

## compute 2D kernel density, see MASS book, pp. 130-131
z <- kde2d(X[,1], X[,2], n=50)

#levels <- sapply(prob, function(x) {approx(c1, sz, xout = 1 - x)$y })

jpeg("contour.jpg")
plot(X, xlab = "Average Distance", ylab = "Standard Deviation", xlim = c(0,185), ylim=c(0,80),pch=19, cex=.4)
#contour(z, levels= round (levels,7), add=TRUE, col = "red")
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=4)
abline(h=mean(X[,2]), v=mean(X[,1]), lwd=2)
legend("topleft", paste("R=", round(cor(X)[1,2],2)), bty="n")
dev.off()

library(hexbin)
hbo <- hexbin(X[,1],X[,2], xlab = "Average Distance", ylab = "Standard Deviation")
plot(hbo)
jpeg("hbo.jpg")
plot(hbo)
dev.off()



histpercent <- function(filename, datavector, gtitle, xlabel, ylabel){
jpeg(filename)
disthist = hist(datavector, breaks = 10)
disthist$counts = disthist$counts / sum(disthist$count)
plot(disthist, main = gtitle, ylab = ylabel, xlab = xlabel)
dev.off()
}

X <- disease_dist_var[,1]
histpercent("bsdist_hist.jpg", X[X<50], "Average Distance Distribution", "Average Distance in Angstrom", "Probability")


filename = "/tmp/snpinfamily.txt"
content <- read.table(filename, sep="\t", quote = "")
jpeg("snp_family_hist.jpg")
hist(content[,1], xlab = "The number of SNP in a family", ylab = "Frequency", main = "Histogram of SNP number in each family")
dev.off()
