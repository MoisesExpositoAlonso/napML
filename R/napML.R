

accuracies<-function(y,x){
  lmo<-lm(y~x)
  a<-fn(coefficients(lmo)[1])
  b<-fn(coefficients(lmo)[2])
  lmos<-summary(lmo)
  r2<-fn(lmos$adj.r.squared)
  rho<-fn(cor(y,x,method="spearman"))
  r<-fn(cor(y,x,method="pearson"))
  return(unlist(list("a"=a,"b"=b,"R2"=r2,"rho"=rho,"r"=r)))
}

read_and_top_gwa<-function(parfile,mapfile,nloci=500){
    gammas<-.read_gemma(folder=dirname(parfile),
                        name=gsub(pattern = ".param.txt", "",basename(parfile)),
                        what = "bslmm")
    if(file.exists(mapfile)){
      map<-data.table::fread(mapfile,verbose = F)
    }else{
      map<-data.table::fread(gsub(x=mapfile, pattern = ".map",replacement = ".bim"))
    }

    try(mapping<-match(moiR::fc(gammas[,2]), moiR::fc(map[,2]) ),silent = T)
    if(!exists("mapping")) stop(".map/.bim file do not share SNP names with .param.txt")

    bslmm<-rep(0,nrow(map))
    bslmm[mapping]<-gammas$effect

    # Define SNPs to analyse
    m<-which(rank(max(abs(bslmm))-abs(bslmm)) <= nloci )

    # Propose starting point of selection coefficients based on BSLMM
    s<-rep(0,length(m))
    s<-bslmm[m]
    return(list(s=s,m=m,bslmm=bslmm))
}

ptr<-function(p) sapply(p,function(i) 1+(-log((1-i)/(i))))
iptr<-function(h) sapply(h,function(i) (1/(1+exp(1-i))))
ptr05<-function(p) sapply(p,function(i) 1+(-log((0.5-i)/(i)))) # this is for the p parameter
iptr05<-function(h) sapply(h,function(i) (0.5/(1+exp(1-i))))
ptr2<-function(p) sapply(p,function(i) 1+(-log((2-i)/(i)))) # this is for the p parameter
iptr2<-function(h) sapply(h,function(i) (2/(1+exp(1-i))))
seltr<-function(s) sapply(s,function(i) 1+ log(1+i))
iseltr<-function(x) sapply(x,function(i) exp(i-1) -1 )


####************************************************************************####
#### Rcpp BM implementation ####
# likelihoodR<-function(y, w, b,a, p){
#   mymin<-.Machine$double.min.exp
#
#   tmp<-ifelse(y==0,
#               p  + (1-p) *  pnorm(0,w,a+w*b,TRUE,FALSE),
#               (1-p) * dnorm(y,w,a+w*b,FALSE)
#              )
#
#   tmp<-log(tmp)
#   tmp[is.na(tmp)]<-mymin
#   tmp[is.infinite(tmp) & tmp<0 ] <- mymin
#   tmp[is.infinite(tmp) & tmp>0 ] <- -mymin # seems to work best
#   LL<-sum(tmp)
#
#   if(is.infinite(LL) & LL<0) LL <- mymin
#   if(is.na(LL)) LL<-mymin
#   return(LL)
# }



# lik.nap_<-function(y,A,par,mod,epi,mu){
#    s=iseltr(par[1:ncol(A)])
#    e=ifelse(epi=="free",par["e"],epi)
#    mu=ifelse(mu=="free",par["mu"],mu)
#    b=iptr(par["b"])
#    a=iptr(par["a"])
#    p=iptr(par["p"])
#
#    w = wC(X=A, s=s,
#             mode=mod,epi=e)
#
#   -likelihoodC(y=y/mu,w=w,b=b,a=a,p=p)
# }

# napSPG_<-function(y,A,s,mod=1,epi=1,mu=1,iter=20){
#   # starting parameters
#   parstart<-par<-list(
#                       "s"=seltr(s),
#                       "b"=ptr(0.1),
#                       "a"=ptr(0.1),
#                       "p"=ptr(0.1)
#                       )
#   if(epi=="free") parstart$e=1
#   if(mu=="free") parstart$mu=1
#   parstart<-unlist(parstart)
#
#   # Test starting parameters
#   w1<-wC(X = A,
#                   s=s,
#                   mode=mod,
#                   epi=epi)
#   # if NAs make zero
#   if(any(is.na(y) || y == -9)){
#     cat("All NAs to zero \n")
#     y[is.na(y)]<-0
#   }
#
#   cat("Initial prediction of fitness: ",cor.test(w1,y,na.rm=T)$estimate,"\n")
#   lik1<-lik.nap_(y,A,parstart,mod,epi,mu)
#   cat("Initial likelihood of fitness: ",lik1,"\n")
#
#   # Optimization
#   start_time <- Sys.time()
#   cat("\n === Spectral Projected Gradient SPG === \n")
#   cat("Starting",iter,"iterations ... \n")
#   r<-spg(
#           fn= lik.nap_,
#           y=y,
#           A=A,
#           mod=mod,
#           epi=epi,
#           mu=mu,
#           par = parstart,
#           control = list(trace=1,maxit=iter)
#           )
#   diftime<-Sys.time() - start_time
#   cat("Finished after",as.numeric(diftime, units = "mins"), "mins \n")
# return(r)
# }

napGR <- function(par, func=lik.nap,...) {
            f<-func(par=par,...)
            df <- rep(NA, length(par))
            for (i in 1:length(par)) {
                dx <- par
                dx[i] <- dx[i] + 1.0e-07
                df[i] <- (func(par=dx, ...) - f)/1.0e-07
            }
            df
}

lik.nap<-function(y,h_,m_,A,par,mod,epi,mu){
   s=iseltr(par[1:length(m_)])
   e=ifelse(epi=="free",par["e"],epi)
   mu=ifelse(mu=="free",par["mu"],mu)
   b=iptr(par["b"])
   a=iptr(par["a"])
   p=iptr(par["p"])

   w = wCBM(A=A@address, s=s,
            mycols =m_,myrows = h_,
            mode=mod,epi=e)

  -likelihoodC(y=y/mu,w=w,b=b,a=a,p=p)
}
napR<-function(bedfile,famfile,mapfile,
                  myrows=NULL,mycols=NULL,
                  s=NULL,mod=1,epi=1,iter=20,k=1,
                  log=NULL,
                  verbose=1){
  require(BB)
  cat("Input files\n")
  cat(mapfile,"\n")
  cat(famfile,"\n")
  cat(bedfile,"\n")
  descfile<-gsub(x = bedfile,pattern = ".bed",replacement = ".desc")
  if(file.exists(descfile)){
    A<-bigmemory::attach.big.matrix(descfile)
  }else{
    stop("transform .bed file into .desc file! not found!")
  }
  ###=======================================================================###
  # Read other files
  cat(" Reading .map file\n")
  if(file.exists(mapfile)){
      map<-data.table::fread(mapfile,verbose = F)
    }else{
      map<-data.table::fread(gsub(x=mapfile, pattern = ".map",replacement = ".bim"))
    }
  cat(" Reading .fam file\n")
  fam<-data.table::fread(famfile)

  y<-moiR::fn(fam[,6])
  N<-nrow(fam)
  p=nrow(map)

  if(is.null(mycols)) mycols<-1:p
  if(is.null(myrows)) myrows<-1:N
  nsnps<-length(mycols)

  if(any(is.na(y)) | any(y == -9)){
    cat("All NAs to zero \n")
    y[is.na(y)]<-0
    y[ y == -9 ]<-0
  }
  ###=======================================================================###
  # starting parameters
  if(is.null(s)) s<-rnorm(nsnps,0,0.01)
  parstart<-par<-list(
                      "s"=seltr(s),
                      "b"=ptr(0.1),
                      "a"=ptr(0.1),
                      "p"=ptr05(0.1),
                      "mu"=ptr2(1)
                      )
  parstart<-unlist(parstart)
  ###=======================================================================###
  # Optimization
  start_time <- Sys.time()
  cat("\n === Spectral Projected Gradient SPG === \n")
  cat("Starting",iter,"iterations ... \n")
  r<-spg(
          fn= lik.nap,
          gr=napGR, ## CAREFUL!
          y=y,
          A=A,
          h_=myrows,
          m_=mycols,
          mod=mod,
          epi=epi,
          mu="free",
          par = parstart,
          control = list(trace=1,maxit=iter)
          )
  diftime<-Sys.time() - start_time
  cat("Finished after",as.numeric(diftime, units = "mins"), "mins \n")
  rr<-parseoptim(r)
  rr$y<-y[myrows]
  rr$w<- rr$mu*wC(A[myrows,mycols],s=rr$s,mode = mod,epi = epi)

  return(r)
}



napSPG<-function(y,h_,m_,A,s,mod=1,epi=1,mu="free",iter=20){
  # starting parameters
  parstart<-par<-list(
                      "s"=seltr(s),
                      "b"=ptr(0.1),
                      "a"=ptr(0.1),
                      "p"=ptr(0.1)
                      )
  if(epi=="free") parstart$e=1
  if(mu=="free") parstart$mu=1
  parstart<-unlist(parstart)

  # Test starting parameters
  w1<-wCBM(A = A@address,
                  s=s,
                  mycols=m_,
                  myrows= h_,
                  mode=mod,
                  epi=epi)
  # if NAs make zero
  if(any(is.na(y))){
    cat("All NAs to zero \n")
    y[is.na(y)]<-0
  }

  cat("Initial prediction of fitness: ",cor.test(w1,y,na.rm=T)$estimate,"\n")
  lik1<-lik.nap(y,h_,m_,A,parstart,mod,epi,mu)
  cat("Initial likelihood of fitness: ",lik1,"\n")

  # Optimization
  start_time <- Sys.time()
  cat("\n === Spectral Projected Gradient SPG === \n")
  cat("Starting",iter,"iterations ... \n")
  r<-spg(
          fn= lik.nap,
          y=y,
          A=A,
          h_=h_,
          m_=m_,
          mod=mod,
          epi=epi,
          mu=mu,
          par = parstart,
          control = list(trace=1,maxit=iter)
          )
  diftime<-Sys.time() - start_time
  cat("Finished after",as.numeric(diftime, units = "mins"), "mins \n")
return(r)
}


# napBFGS<-function(y,h_,m_,A,s,mod=1,epi=1,mu=1,iter=20){
#   # starting parameters
#   parstart<-par<-list(
#                       "s"=seltr(s),
#                       "b"=ptr(0.1),
#                       "a"=ptr(0.1),
#                       "p"=ptr05(0.1)
#                       )
#   if(epi=="free") parstart$e=1
#   if(mu=="free") parstart$mu=1
#   parstart<-unlist(parstart)
#
#   # Optimization
#   start_time <- Sys.time()
#   cat("BFGS \n")
#   cat("Starting",iter,"iterations ... \n")
#   r<-optim(
#           fn= lik.nap,
#           y=y,
#           A=A,
#           h_=h_,
#           m_=m_,
#           mod=mod,
#           epi=epi,
#           mu=mu,
#           par = parstart,
#           method="BFGS",
#           control = list(trace=1,maxit=iter)
#           )
#   diftime<-Sys.time() - start_time
#   cat("Finished after",as.numeric(diftime, units = "mins"), "mins \n")
# return(r)
# }



parseoptim<-function(r){
  res<-list()
  res$s<-iseltr(unlist(r$par[grepl("s",names(r$par))]))
  if("mu" %in% names(r$par)) res$mu<-iptr(r$par["mu"])
  if("b" %in% names(r$par)) res$b<-iptr(r$par["b"])
  if("a" %in% names(r$par)) res$a<-iptr(r$par["a"])
  if("p" %in% names(r$par)) res$p<-iptr(r$par["p"])
  if("e" %in% names(r$par)) res$e<-iptr(r$par["e"])
  L= -r$value
  res$LL= L
  res$AIC= 2* (length(r$par)) -2 * L
  return( res)
}

summaryseoptim<-function(r,y,G){
    winf<-sapply(1:k,function(i){wCBM(A = G@address,
                                  s=r$s,
                                  mycols=m,
                                  myrows= h[which(cv==i)],
                                  mode=mod,
                                  epi=epi)
                                }
              ) %>% unlist
    wgwa=wCBM(A = G@address,
                      s=s,
                      mycols=m,
                      myrows= h,
                      mode=1,
                      epi=1
          )
    ytest<-sapply(1:k, function(i) y[cv==i]) %>% unlist
    rnap<-cor(winf,ytest)
    rgwa<-cor(wgwa,y)
    rnap_no<-cor(winf[ytest!=0],ytest[ytest!=0])
    rgwa_no<-cor(wgwa[y!=0],y[y!=0])
}


writeaccuracy<-function(pheno,mod,epi,rgwa,rnap,rgwano,rnapno,finaltsv){
  d<-data.frame(phenotype=pheno,
              method=c("gwa","nap","gwa_nozero","nap_nozero"),
              mod=mod,
              epi=epi,
              r=c(rgwa,rnap,rgwano,rnapno)
              )
  write.table(row.names = F, quote = F, sep = "\t",
              x=d,
              file=finaltsv
              )
  return(d)
}

accuracyresults<-function(pheno,mod,epi,wgwa,winf,wtest,y,ytest,finaltsv){
  yno<-which(y!=0)
  ytestno<-which(ytest!=0)
  d<-data.frame(phenotype=pheno,
              method=c("gwa","nap","nap_t",
                       "gwa_nozero","nap_nozero","nap_t_nozero"),
              mod=mod,
              epi=epi
              )
  d_<-rbind(
    accuracies(wgwa,y),
    accuracies(winf,y),
    accuracies(wtest,ytest),
    accuracies(wgwa[yno],y[yno]),
    accuracies(winf[yno],y[yno]),
    accuracies(wtest[ytestno],ytest[ytestno])
  )
  d<-data.frame(cbind(d,d_))
  write.table(row.names = F, quote = F, sep = "\t",
              x=sapply(d,as.character),
              file=finaltsv
              )
  return(d)
}




startcluster<-function(no_cores=NULL){
require("parallel")
message("starting cluster for parallel operations")

# Calculate the number of cores
if(is.null(no_cores)){
  no_cores <- detectCores() - 1
  if(is.na(no_cores)) no_cores <- 2
  if(no_cores>4){no_cores=10}
}

 # Initiate cluster
# cl <- makeCluster(no_cores,type="FORK")
# Now we just call the parallel version of lapply, parLapply:
cl <- makeCluster(no_cores)

return(cl)
}

stopcluster<-function(cl){
  stopCluster(cl)
  message("cluster closed")

}
kcv<-function(y, h, k){
  mapp<-sample(1:k,size = length(h),replace = T)
}
napCV<-function(y, m , h, G,bmfile=NULL, s,mod, epi, iter,k=10, parallel=TRUE,cores=NULL,is_sim_data=FALSE){
    # stopping conditions
    stopifnot(k>1)
    # Get NAs out
    # h<-h[which(!is.na(y))]
    # y<-y[which(!is.na(y))]
    cat("All NAs to zero \n")
    y[is.na(y)]<-0
    ## Cross-Validation design
    cv<-kcv(y,h,k)
    ## Parallel optimization
    # function to evaluate
    ex<-expression(parseoptim(napSPG(y = y[which(cv!=i)],
                            m_ =m,
                            h_ = h[which(cv!=i)],
                            A = G,
                            s=s,
                            mod=mod,
                            epi=epi,
                            iter=iter
                            ) )
          )
    # CV optimization
    if(parallel==TRUE ){

    cat("Starting analyses in parallel \n")
    require(parallel)
      # set up cluster and source the file on each core
      # if you don't do this you will likely get a strange error about a NULL value passed as symbol address, because reasons
      cl <- makeSOCKcluster(4)
      registerDoSNOW(cl)
      clusterExport(cl=cl, list("ex","y","h","m","G","A","s","cv","mod","epi","iter","bmfile",
                                "parseoptim","napSPG","lik.nap","spg",
                                "iseltr","seltr","ptr","iptr","iptr","ptr"),
                    envir=environment())
      cat("Compile C++ in each node \n")
      clusterEvalQ(cl,{
        Rcpp::sourceCpp(here::here("napML/src/ML.cpp"))
        G<-A<-bigmemory::attach.big.matrix(bmfile)
      })
      cat("Running optimizations \n")
      res<-foreach(i=1:k,
                   .noexport = c("wCBM","likelihoodC")) %dopar% {
           r<-eval(ex)
      }
      stopCluster(cl);rm(cl);gc()

    }else{#
      res<-lapply(1:k,function(i){eval(ex)})
    }
    #### analyze results
    G<-A<-bigmemory::attach.big.matrix(bmfile)
    winf<-sapply(1:k,function(i){wCBM(A = G@address,
                                  s=res[[i]]$s,
                                  mycols=m,
                                  myrows= h[which(cv==i)],
                                  mode=mod,
                                  epi=epi)
                                }
              ) %>% unlist
    wgwa=wCBM(A = G@address,
                      s=s,
                      mycols=m,
                      myrows= h,
                      mode=1,
                      epi=1
          )
    ytest<-sapply(1:k, function(i) y[cv==i]) %>% unlist
    rnap<-cor(winf,ytest)
    rgwa<-cor(wgwa,y)
    rnap_no<-cor(winf[ytest!=0],ytest[ytest!=0])
    rgwa_no<-cor(wgwa[y!=0],y[y!=0])

    sinf<-lapply(1:k,function(i){res[[i]]$s}) %>%
      do.call(what = rbind,.) %>%
      apply(.,2,mean)

    #### AIC
    AICall<-min(sapply(1:k,function(i) res[[i]]$AIC))

    return(list(run=res,s0=s,sinf=sinf,
                winf=winf, ycv=ytest,
                wgwa=wgwa,y=y,
                r=rnap,rgwa=rgwa,
                rno=rnap_no,rgwano=rgwa_no,
                AIC=AICall))
}

####************************************************************************####
#### C++ implementation ####
# napCV<-function(y, m , h, G,bmfile=NULL, s,mod, epi, iter,k=10, parallel=TRUE,cores=NULL,is_sim_data=FALSE){

nap<-function(bedfile,famfile,mapfile,
                  myrows=NULL,mycols=NULL,
                  s=NULL,mod=1,epi=1,iter=20,k=1,
                  log=NULL,
                  verbose=1){
  cat("Input files\n")
  cat(mapfile,"\n")
  cat(famfile,"\n")
  cat(bedfile,"\n")
  ###=======================================================================###
  # Read other files
  cat(" Reading .map file\n")
  if(file.exists(mapfile)){
      map<-data.table::fread(mapfile,verbose = F)
    }else{
      map<-data.table::fread(gsub(x=mapfile, pattern = ".map",replacement = ".bim"))
    }
  cat(" Reading .fam file\n")
  fam<-data.table::fread(famfile)

  y<-moiR::fn(fam[,6])
  N<-nrow(fam)
  p=nrow(map)

  if(is.null(mycols)) mycols<-1:p
  if(is.null(myrows)) myrows<-1:N
  nsnps<-length(mycols)

  if(any(is.na(y)) | any(y == -9)){
    cat("All NAs to zero \n")
    y[is.na(y)]<-0
    y[ y == -9 ]<-0
  }

  ###=======================================================================###
  # starting parameters
  if(is.null(s)) s<-rnorm(nsnps,0,0.01)
  parstart<-par<-list(
                      "s"=seltr(s),
                      "b"=ptr(0.1),
                      "a"=ptr(0.1),
                      "p"=ptr05(0.1),
                      "mu"=ptr2(1)
                      )
  parstart<-unlist(parstart)
  ###=======================================================================###
  # Optimization
  if(!is.null(log)) sink(log)
  start_time <- Sys.time()
    ### Cross-Validation design ###
    if(k>1){
      cat("Dividing data for",k, "fold crossvalidation \n" )
      cv<-sample(1:k,size = length(myrows),replace = T)
      rcv<-lapply(1:k,function(i){
                                  cat("=>Runing group",i," \n" )
                                  napspgC(
                                  bedfile,
                                  N, p,myrows[which(cv!=i)],mycols,
                                  y[which(cv!=i)],
                                  parstart,
                                  epi,
                                  mod,
                                  iter,
                                  verbose=F)
                                    })
      r_<-list()
      r_$w<-lapply(1:k,function(i){rcv[[i]]$w[which(cv==i)]}) %>% unlist
      r_$y<-sapply(1:k, function(i) y[which(cv==i)]) %>% unlist()
    }
    ### Run on 100% data
    r<-napspgC(
              bedfile,
              N, p,myrows,mycols,
              y[myrows],
              parstart,
              epi,
              mod,
              iter,
              verbose)
    r$AIC <- 2*length(parstart) -2 *(-r$f) # function value given negative likelihood, so need double neg
    r$y<-y[myrows]
    if(k>1){
      r$y_in<-r$y
      r$w_in<-r$w # keep also in-sample prediction
      r$y<-r_$y
      r$w<-r_$w
    }
    r$myrows<-myrows
    r$mycols<-mycols

  diftime<-Sys.time() - start_time
  cat("Finished after",as.numeric(diftime, units = "mins"), "mins \n")
  if(!is.null(log)) sink(log)
  ###=======================================================================###
  # Accuracy
  r$w[r$w<0] <-0
  r$rpearson<-cor.test((r$w),(r$y),na.rm=T)$estimate
  r$rnonzero<-cor.test(r$w[r$y!=0],r$y[r$y!=0],na.rm=T)$estimate
  ###=======================================================================###
  return(r)
}


napcall<-function(bedfiles,famfiles,mapfiles,parfiles,
                   es=c(0.9,1,1.1),mods=c(1,2),loci=100,
                   iter=100,
                   k=10,
                   dryrun=T,
                   override=F,
                   background=T,
                   maxbackground=10
                  ){
  counter<-1
  for(i in 1:length(bedfiles)){
    for(e in es){
      for(mod in mods){
        for(l in loci){
        background<-ifelse(system("ps | grep R | wc -l")<20,
                           " &"," ")
        if((counter %% maxbackground)==0) background<-" "
        command<-paste("nice Rscript ~/napML/napcall.R",
                       "--param", parfiles[i],
                       "--fam", famfiles[i],
                       "--map", mapfiles[i],
                       "--bed", bedfiles[i],
                       "--e", e,
                       "--m", mod,
                       "--l", l,
                       "--o", override,
                       "--i", iter ,
                       "--k", k ,
                       " ", background
                       )
        message(command)
        if(!dryrun){system(command)}
        counter=counter+1
        }
      }
    }
  }
}

