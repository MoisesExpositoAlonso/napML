require(BB)
library(moiR)
library(data.table)
library(devtools)
library(bigmemory)
load_all(".")

bedfile<-"easy.bed"
mapfile<-"easy.map"
famfile<-"easy.fam"

map<-fread(mapfile)
fam<-fread(famfile)


y<-fam[,6]
N<-nrow(fam)
p=nrow(map)

library(microbenchmark)
nsnps=1000
sstart<-rnorm(nsnps,0,0.01)

# X<-readbed(bedfile,N,p , 1:N,1:nsnps)
# G<-attach.big.matrix("example.desc")


napSPG_R(bedfile,N,p,1:N,1:nsnps,fn(y),sstart,mod=1,epi=1,iter=100)
# r<-napSPG(fn(y),A = G,h_ = 1:nrow(G),m_=1:nsnps,s = sstart,iter = 3)

d<-read.table("solution.txt")
sinf<-iseltr(d[1:nsnps,1])
iptr(d[(nsnps+1):(nsnps+3),1])
iptr05(d[nsnps+4,1])

strue=ssimC(1:nsnps,0.01)
plot(sinf,sstart)
plot(sinf,strue)


microbenchmark(times = 3,
  memory=napSPG_(fn(y),A = readbed(bedfile,N,p , 1:N,1:nsnps),s = sstart,iter=50) ,
  filebacked=napSPG(fn(y),A = G,h_ = 1:nrow(G),m_=1:nsnps,s = sstart,iter=50),
  onlyC=napSPG_R(bedfile,N,p,1:N,1:nsnps,fn(y),sstart,mod=1,epi=1,iter=50)
)
