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
