filename = "../pair_small.txt"
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
histpercent("pairdist.jpg", x, "Pair Spatial Distance", "Probability", "Distance in Angstrom")

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

histpercent("bsdist_hist.jpg", distance[,1], "Distance to Binding Site", "Probability", "Distance in Angstrom")
