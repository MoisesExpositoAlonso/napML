

fieldcodes=function(sites=c('m','t'),water=c('h','l'),rep=c('i','p')){
  tmp=expand.grid(sites,water,rep)
  paste(tmp[,1],tmp[,2],tmp[,3],sep="")

}
