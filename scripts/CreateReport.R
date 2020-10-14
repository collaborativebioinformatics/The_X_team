t <- read.table("MaternalPaternal.tab")

library(stringr)

head(t)
hom <- which(str_detect(t$V1, "Hom"))
het <- which(str_detect(t$V1, "Het"))
low <- which(str_detect(t$V1, "Lowfreq"))
maxLen <- 3e4
shortSV <- which(t$V9 < 3e4 & t$V9 >= 50)
matConf <- which(t$V7 < 10) | t$V7 < 0.05*t$V9)
patConf <- which(t$V15 < 10 | t$V15 < 0.05*t$V17)

matConfStrict <- which(t$V7 < 5) | t$V7 < 0.05*t$V9)
patConfStrict <- which(t$V15 < 5 | t$V15 < 0.05*t$V17)

plot(sort(t$V9*0.05), log="y")
o <- order(t$V9*.005)
points(t$V7[o], col='red')
both <- union(matConf,patConf)
notConf <- setdiff(seq(1,length(t$V15)), both)
br <- seq(0,maxLen,100)
hist(t$V9[intersect(shortSV,intersect(hom,both))],breaks=br)
hist(t$V9[intersect(shortSV,intersect(het,both))],breaks=br)
hist(t$V9[intersect(shortSV,intersect(low,both))],breaks=br)
hist(t$V9[intersect(shortSV,intersect(hom,notConf))],breaks=100)


