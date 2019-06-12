####************************************************************************####
#### WRAPS ####
# leadingSNP<-function(dat, significancecol='significance',window=5e4){
#
# colnames(dat) [which(colnames(dat) == significancecol)] <- 'significance'
#
# dat   <- dat %>%
#             mutate(posround=round(pos / window)*window) %>%
#             mutate(posround=paste(chr,posround,sep='_'))  %>%
#             group_by(posround) %>%
#             summarize(SNPtop = head(SNP[significance == min(significance)],n=1)  ) ->
#             selected
#
# return(selected$SNPtop)
# }

# top_gwa<-function(as){
#
#   res<-lapply(as, FUN = function(a){
#                 message("Reading ", a)
#                 tmp<-.read_assoc(a)
#                 tmp<-.selecttop(tmp)
#                 return(tmp) }
#               ) %>%
#     do.call(rbind,.)
#   return(res)
# }
#
# sub_gwa<-function(as,SNPs){
#   res<-lapply(as, FUN = function(a){
#                 message("Reading ", a)
#                 tmp<-.read_assoc(a)
#                 tmp$SNP<-paste(tmp$chr,tmp$pos,sep='_')
#                 tmp<-dplyr::filter(tmp,SNP %in% SNPs)
#                 }
#               ) %>%
#     do.call(rbind,.)
#   return(res)
# }
#
# all_gwa<-function(as){
#   res<-lapply(as, FUN = function(a){
#                 message("Reading ", a)
#                 d<-.read_assoc(a)
#                   return(d)
#                 }
#               ) %>%
#     do.call(rbind,.)
#   return(res)
# }
#
#
# all_gwa_onlyeff<-function(as,mycol='effect'){
#   data(map)
#   res<-lapply(as, FUN = function(a){
#                 message("Reading ", a)
#                 d<-.read_assoc(a)
#                 myname<-d$env[1]
#                 tmp <- data.frame(d[,'SNP'],d[,mycol])
#                 colnames(tmp) <- c('SNP',myname)
#                 tmp2<-merge(map[,'SNP'],tmp, by='SNP',all.x=T)
#                 return(tmp2)
#                 }
#               ) %>%
#     do.call(cbind,.)
#   ressnp<-res [,'SNP']
#   res<-dplyr::select(data.frame(res),-contains('SNP'))
#   res<-data.frame(ressnp,res)
#
#   fin<-merge(map,res, by='SNP',all.x=T)
#
#   return(fin)
# }


# wrap_gwa<-function(a,...){
#   gwares=.read_assoc(a)
#   aname=.cleanname(a)
#
#   p<-gws::ggmanhattan(gwares,stat.col = 'p_lrt', stat.type = 'p',subset = 5000,...)
#   es<-gws::ggmanhattan(gwares,stat.col = 'beta', stat.type = 'es',subset = 5000,...)
#
#   comb=plot_grid(p, es,nrow=2)
#
#   save_plot(file.path('figs/',paste0(aname,".pdf") ),plot = comb,base_width = 12,base_height = 4 ,useDingbats = F)
#
#   return(list(p,es))
# }




####************************************************************************####
#### READING ####
bedto012<-function(outplink){
  command1<-sprintf("plink --noweb  --bfile %s --recode --tab --out %s", outplink,outplink)
  system(command1)
  # prettymap $outplink".map"
  # prettybim $outplink".bim"

  command2<-sprintf("plink --noweb --file %s  --recodeA --out %s", outplink,outplink)
  system(command2)

  command3<-sprintf("cat %s.raw | cut -d ' '  -f 7-  > %s.012", outplink,outplink)
  system(command3)
}

#' Get the base name
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
.cleanname<-function(x){
  sapply(x, function(i){
    strsplit(i, split =  '/', fixed=TRUE) %>% unlist %>% tail(1) %>%
    strsplit(split =  '.', fixed=TRUE) %>% unlist %>%
    head(1)
  })
}

#' Find gemma files in a folder
#'
#' @param folder
#' @param type
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
.find_gemmas<-function(folder='output/', type='heritability',...){
  if(type=='heritability'){
    filelist=list.files(folder,pattern = '.hyp.',...)
  }else if(type=='association_bslmm'){
    filelist=list.files(folder,pattern = '.param.',...)
  }else if(type=='association_lm'){
    filelist=list.files(folder,pattern = '.assoc.',...)
  }else if(type=='association_lmm'){
    filelist=list.files(folder,pattern = '.assoc.',...)
  }else{stop("argument not recognized")}
  return(filelist)
}


#' Read different files from GEMMA
#'
#' @param folder
#' @param name
#' @param what
#'
#' @return
#' @export
#'
#' @examples
.read_gemma<-function(folder="output", name, what='heritability'){
  ## Define functions
    .read_hyperparameters<-function(outfile, hyperparameter=1, MCMCchain=1000,quantiles= c(0.5,0.025,0.975) ){
    hyp<-read.table(outfile,header=TRUE, stringsAsFactors = FALSE) %>% tail(n=MCMCchain)
    # hist(hyp[,hyperparameter], main = outfile)
    quantile(hyp[,hyperparameter],p=quantiles)
  }

  .read_assoc<-function(outfile="output/fecundity_mhi.assoc.txt"){
    d <- data.table::fread(outfile,header=TRUE,stringsAsFactors = FALSE)
    # d <- read.table(outfile,header=TRUE,stringsAsFactors = FALSE)
    d<-dplyr::rename(d,pos=ps)
    d<-mutate(d, SNP= paste(chr,pos,sep="_"),
              env= .cleanname(outfile),
              cumpos=1:nrow(d))

    if (grepl('.param.', outfile)){
      d<-mutate(d,effect= alpha + (gamma*beta) )
      d<-mutate(d,sparseeffect= (gamma*beta) )
      d$BAY<- d$gamma !=0

    }else if(grepl('.assoc.', outfile)){
      d<-mutate(d,effect= beta)
      d$FDR<- d$p_lrt < fdr_level(d$p_lrt)
      d$BONF<- d$p_lrt < 1/nrow(d)
    }

    return(d)
  }
  ## Run
  if(what=='heritability'){
    res=.read_hyperparameters(file.path(folder,paste0(name,".hyp.txt")))
  }else if(what=='lm'){
    res=.read_assoc(file.path(folder,paste0(name,".assoc.txt")))
  }else if(what=="bslmm"){
    res=.read_assoc(file.path(folder,paste0(name,".param.txt")))
  }else if(what=="bv"){
    res=read.table(file.path(folder,paste0(name,".bv.txt")))
  }
    return(res)
}




####************************************************************************####

#' Create a directory with all files needed for GEMMA
#' Generates the plink necessary fam file to run gemma
#'
#' @param data
#' @param id
#' @param phenotype
#' @param naval
#' @export
genplink<-function(data,
                   id='id',
                   phenotype,
                   naval='-9',
                   out='test',
                   path='../gemma/',
                   famfile='515g.fam',
                   bimfile='515g.bim',
                   bedfile='515g.bed',
                   cleardir=F,
                   makerelative=F
                   ){

# load the fam
data(fam)
fam<-data.frame(fam)

idname<-ifelse('id' %in% colnames(fam),'id','sample.ID')

# merge with the dataset
myplink<-merge(fam, data[,c('id',phenotype)] , by.x=idname,by.y='id' ,all.x = T)

myplink<-data.frame(
                    sample.ID= myplink[,idname],
                    family.ID= myplink[,'family.ID'],
                    paternal.ID= myplink[,'paternal.ID'],
                    maternal.ID= myplink[,'maternal.ID'],
                    sex= myplink[,'sex'],
                    affection= myplink[,phenotype]
                    )

# make relative if needed
if(makerelative==TRUE){
    myplink[,'affection'] = myplink[,'affection']/mean(myplink[,'affection'],na.rm=T)
}

# check that the NA are as naval
myplink[,'affection'][is.na(myplink[,'affection'])]<-naval # sixth column is the phenotype

# output name
fampath=file.path(path, out)
famfile=file.path(fampath,'/515g.fam')

# create directory
catchexit=system(paste('mkdir',fampath))

if(catchexit !=0){
message("the directory already exists!")
  if(cleardir==TRUE){
    message('flag cleandir=T, so removing directory')
    clear_gemma_dir(out)
    system(paste('mkdir',fampath))
  }else{
    stop('provide permission to remove directory with cleardir=TRUE')
  }
}
# hard link the genome matrix files for gemma
system(paste('ln', file.path(path, bedfile),fampath))
system(paste('ln', file.path(path, bimfile),fampath))


# out
write.table(myplink,file=famfile,col.names=F,row.names=F,quote=F,sep=" ")
print(head(myplink))

}


clear_gemma_dir<-function(out,gemmapath='../gemma',force=F){
  fampath=file.path(gemmapath,out)
  if(force) system(paste('rm -rf',fampath))
  else message("Run with force=T to run and erase: ", paste('rm -rf',fampath))
}


####************************************************************************####
#### RUNNING ####
#' Call GEMMA from within R
#'
#' @param out
#' @param plinkbase
#' @param plinkfolder
#' @param type
#' @param maf
#' @param background
#' @param dryrun
#'
#' @export
run_gemma<-function(
                    plinkbase,
                    plinkfolder='.',
                    out=".",
                    type='bslmm',
                    maf=0.0,
                    background=TRUE,
                    dryrun=FALSE
                    ){

background=ifelse(background==F, " ", " &")

if(type=='bslmm'){
  command= paste0('nice ~/bin/gemma -bfile ', file.path(plinkfolder, plinkbase) ,'  -bslmm ')
}else if(type=='lm'){
  command= paste0('nice  ~/bin/gemma -bfile ', file.path(plinkfolder, plinkbase) ,'  -lm 4  ')
}else if(type=='lmm'){
  command0= paste0('nice  ~/bin/gemma -bfile ', file.path(plinkfolder, plinkbase) ,'  -gk 1  -o ' ,out)
  command= paste0('nice  ~/bin/gemma -bfile ', file.path(plinkfolder, plinkbase) ,'  -k ', paste0(out,'.cXX.txt') , '  -lmm 4 ')
}
if( !type %in% c('bslmm','lm','lmm')) stop('type of GWA not recognized')

# Add maf
command<-paste(command, "-maf", maf)
# Add out
command<-paste(command, ' -o ' ,out)
# Add background
command<-paste(command, background)

# Run gemma
if(dryrun==TRUE){
  message("Dry run, only printing command")
  if(type=='lmm') print(command0)
  print(command)
}else{
message(paste('running GEMMA command: ',command))
  if(type=='lmm') system(command0)
  system(command)
}

return(TRUE)
}

#' Title
#'
#' @param bfile
#' @param emp
#'
#' @return
#' @export
#'
#' @examples
pred_gemma<-function(
                     predict=1,
                     bfile="../gemma/rFitness_mhi/515g",
                     epm="output/rFitness_mhi.param.txt",
                     emu="output/rFitness_mhi.log.txt",
                     o="rFitness_mhi-predict",
                     dryrun=F
                     ){
  command<-paste("nice  ~/bin/gemma -bfile",bfile,
                 "-predict",predict,
                 "-epm",epm,
                 "-emu",emu,
                 "-o",o
                )
  # gemma -bfile ../gemma/rFitness_mhi/515g -epm output/rFitness_mhi.param.txt -emu output/rFitness_mhi.log.txt  -o rFitness_mhi -predict 1
  if(dryrun) system(command)
  print(command)
  return(TRUE)
}


####************************************************************************####
#' Get the FDR threshold from a set of p-values
#'
#' @param pval
#'
#' @return
#' @details
#' The function will do a p-value adjustment using p.adjust(method='fdr') and
#' then will search for the p-value in the vector provided that corresponds to
#' that element that is below the alpha level. If there is none, it will report
#' the minimum p-value after the FDR transformation, which will be above the
#' alpha level and then indicates that none of the provided p-values would be
#' significant after the FDR correction.
#'
#' @export
#'
#' @examples
fdr_level<-function(pval,al=0.05){
  stopifnot(class(pval) == 'numeric')

  pad= p.adjust(pval,method = 'fdr')

  tmp=data.frame(pad=pad,pval=pval)
  tmp=dplyr::arrange(tmp,pad)
  head(tmp)
  tail(tmp)


  signis=pval[pad<al]

  whatmin=tmp$pad[tmp$pad == min(tmp$pad) ]

  res=ifelse(length(signis) !=0,
                max(signis),
                min(pad))    # if there is none that is significative, just the smallest value

  return(res)
}
