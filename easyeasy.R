library(dplyr)
library(moiR)
library(data.table)
library(devtools)
#load_all("~/napML")
# install("~/napML")
library(napML)


bedfile<-"easy.bed"
mapfile<-"easy.map"
famfile<-"easy.fam"

map<-fread(mapfile)
fam<-fread(famfile)
bed<-readbed(bedfile, N = nrow(fam), p = nrow(map), myrows = 1:nrow(fam),mycols=1:nrow(map))

y<-fam[,6]
N<-nrow(fam)
p=nrow(map)
nsnps=500 # define first X SNs to analyse


r2<-nap(bedfile,famfile,mapfile,
        myrows=1:N,
        mycols=1:nsnps,
        s = NULL,
        mod=1,epi=1,iter=100,k=2)

plot(r2$w,r2$y)
