
xy = read.table(filename, sep = "\t")
x  = xy[,1]
y  = xy[,2]

library(hexbin)
hbo <- hexbin(x,y)
plot(hbo)
