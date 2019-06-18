require(devtools)
require(dplyr)
require(cowplot)
require(latex2exp)
require(data.table)
require(moiR)
require(Rcpp)
require(argparser)
devtools::load_all("~/napML")
setwd("~/nap")

cat("\n")
p <- arg_parser("Run Non-Additive Polygenic ML model:")
# Input files
p <- add_argument(p, "--param", help="path to .param.txt",default="output/b0.01_a0.01_p0_svar0.01_epi1_mod2_h21.param.txt")
p <- add_argument(p, "--fam", help="path to .fam", default="databig/b0.01_a0.01_p0_svar0.01_epi1_mod2_h21/example.fam")
p <- add_argument(p, "--map", help="path to .map", default="databig/b0.01_a0.01_p0_svar0.01_epi1_mod2_h21/example.map")
p <- add_argument(p, "--bed", help="path to .bed without extension", default = "databig/example.bed")
# dros
# p <- add_argument(p, "--param", help="path to .param.txt",default="output/drosophila.param.txt")
# p <- add_argument(p, "--fam", help="path to .fam", default="databig/dgrp2.fam")
# p <- add_argument(p, "--map", help="path to .map", default="databig/dgrp2.map")
# p <- add_argument(p, "--bed", help="path to .bed without extension", default = "databig/dgrp2.bed")
# ara
# p <- add_argument(p, "--param", help="path to .param.txt",default="output/rFitness_mhp.param.txt")
# p <- add_argument(p, "--fam", help="path to .fam", default="databig/rFitness_mhp/515g.fam")
# p <- add_argument(p, "--map", help="path to .map", default="databig/515g.map")
# p <- add_argument(p, "--bed", help="path to .bed without extension", default = "databig/515g.bed")
# yeast
# p <- add_argument(p, "--param", help="path to .param.txt",default="output/YPDHU.param.txt")
# p <- add_argument(p, "--fam", help="path to .fam", default="databig/YPDHU/yeast.fam")
# p <- add_argument(p, "--map", help="path to .map", default="databig/YPDHU/yeast.map")
# p <- add_argument(p, "--bed", help="path to .bed without extension", default = "databig/YPDHU/yeast.bed")
# p <- add_argument(p, "--param", help="path to .param.txt",default="output/YPD6AU.param.txt")
# p <- add_argument(p, "--fam", help="path to .fam", default="databig/YPD6AU/yeast.fam")
# p <- add_argument(p, "--map", help="path to .map", default="databig/YPD6AU/yeast.map")
# p <- add_argument(p, "--bed", help="path to .bed without extension", default = "databig/YPD6AU/yeast.bed")
# Analysis conditions
p <- add_argument(p, "--l", help="number of SNPs to analyze", default=500)
p <- add_argument(p, "--n", help="fraction of individuals to analyze", default=0.9)
p <- add_argument(p, "--e", help="global epistasis power ('free' to be inference)", default=1)
p <- add_argument(p, "--m", help="fitness mode: 1 for additive, 2 for multiplicative", default=2)
p <- add_argument(p, "--i", help="maximum optimization iterations", default=1000)
p <- add_argument(p, "--k", help="cross-validations of the dataset, k=10 default", default=10)
p <- add_argument(p, "--r", help="random start of parameters?", default=F)
# Output behaviour
p <- add_argument(p, "--d", help="dry run?", default=FALSE)
p <- add_argument(p, "--o", help="overwrite?", default=FALSE)
argv<-parse_args(p)

#### Create filenames

dname<-dirname(argv$p)
dataname<-dirname(argv$f)
bname<-gsub(x=basename(argv$p),pattern=".param.txt",replacement="")

finalbase=paste0(dname,"/",bname,".results.",
            ifelse(argv$m==1,"a_","m_"),
            paste0("e",argv$e,"_"),
            paste0("l",argv$l)
          )
finalfile<-paste0(finalbase,".tsv")
finallog<-paste0(finalbase,".log")
finalrda<-paste0(finalbase,".rda")
finaltsv<-paste0(finalbase,".tsv")
finaliplot<-paste0(finalbase,"_indplot.pdf")
finalsplot<-paste0(finalbase,"_splot.pdf")
successfile<-paste0(finalbase,".success")

#### Smart stop
if(all(file.exists(c(successfile)) ) ){
  prevAIC<-fn(head(read.table(successfile),1))
  if(abs(prevAIC) < 1e+200){
    if(argv$o==T){
      message("Overwriting mode!")
    }else{
      stop("Result files already exist and run AIC does not indicate failed optimization! Stopping now")
    }
  }else{
      message("Result files already exist but AIC indicates a failed optimization run")
  }
}



#### BSLMM starting point
cat("Reading BSLMM GWA to define starting conditions \n")
cat("(second column of .parm.txt and .map need to match) \n")

snps<-read_and_top_gwa(argv$param,argv$map,argv$l)
s<-snps$s
m<-snps$m
# if(argv$r)
s<-rnorm(length(m),0,0.1)

#### OPtimization
cat("Optimization  \n")
iter=argv$i
mod=argv$m
epi=argv$e
k=argv$k

bedfile<-argv$bed
famfile<-argv$fam
mapfile<-argv$map

# production
res<-nap(bedfile,famfile,mapfile,
           myrows = NULL,
           mycols=m,
           s=s,
           mod=mod,
           epi=epi,
           iter=iter,
           k=k
           )

#### write ####

d<-data.frame(phenotype=bname,
            method=c("nap","nap_nozero"),
            mod=mod,
            epi=epi,
            r=c(res$rpearson,res$rnonzero)
            )
print(d)
write.table(row.names = F, quote = F, sep = "\\t",
            x=d,
            file=finaltsv
            )

#### save ####
saveRDS(file = finalrda,object = res)

#### plots ####
# plotresults<-plot_grid(labels=c("ML","BSLMM"),ncol=1,
#                       indplot(res$ycv,res$winf)$pind  ,
#                       indplot(res$y,res$wgwa)$pind
#                       )
# save_plot(finaliplot,
#           plotresults,base_height = 2*5,base_width = 5)
plotresults<-indplot(res$y,res$w)$pind
plotresults

save_plot(finaliplot,
          plotresults,base_height = 5,base_width = 5)

#### save success file
sink(successfile)
cat(res$AIC,"\n")
sink()
