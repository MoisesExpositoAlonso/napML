

allelesim <- function(
mu=0,
nu=0,
m=0,
wAA=1,
wAa=0.75,
waa=0.5,
Fi=0,
p0=0.5,
psource=0.1,
tmax=100,
d=0,
N=1000,
rep=50
  ) {

  sapply( 1:rep , FUN=function(rep){
    p <- c()
    p[1] <- p0

    for (t in 1:(tmax-1)) {
      # mutation first
      pp <- (1-mu)*p[t] + (1-p[t])*nu

      # next, migration
      ppp <- (1-m)*pp + m*psource

      # then selection
      w_hat=( wAA*ppp^2 + wAa*2*ppp*(1-ppp) + waa*(1-ppp)^2 )
      print(w_hat)
      if ( w_hat > 0 & !is.na(w_hat) ) {
        p[t+1] <- ( wAA*ppp^2 + wAa*ppp*(1-ppp) ) / w_hat
      } else {
        p[t+1] <- NA
      }
      # then imbreeding (this equation is general)
      fAA <- (p[t+1]^2 * (1-Fi)) + p[t+1]*Fi
      fAa <- 2*p[t+1]*(1-p[t+1]) * (1-Fi)
      faa <- (1-p[t+1])^2 * (1-Fi) + (1-p[t+1])* Fi
      # no imbreeding
      # fAA <- p[t+1]^2
      # fAa <- 2*p[t+1]*(1-p[t+1])
      # faa <- (1-p[t+1])^2

      # then drift
      NAA <- round(N*fAA)
      NAa <- round(N*fAa)
      Naa <- round(N*faa)

      # if (NAA <= 0) {
      #   NAA <- 0
      #   NAA.prime <- 0
      # } else {
      #   NAA.prime <- sum(rbinom(n=NAA, size=1, prob=d))
      # }
      # if (NAa <= 0) {
      #   NAa <- 0
      #   NAa.prime <- 0
      # } else {
      #   NAa.prime <- sum(rbinom(n=NAa, size=1, prob=d))
      # }
      # if (Naa <= 0) {
      #   Naa <- 0
      #   Naa.prime <- 0
      # } else {
      #   Naa.prime <- sum(rbinom(n=Naa, size=1, prob=d))
      # }
      # N.prime <- NAA.prime + NAa.prime + Naa.prime
      #
      # if (N.prime <= 0) {
      #   p[t+1] <- NA
      # } else {
      # p[t+1] <- (2*NAA.prime + NAa.prime) / (2*N.prime)
      # }

      p[t+1] <- (2*NAA + NAa) / (2*Naa)

    } #end t loop


    return(p)
}) #end sapply
}

############ funcitons
mafsim<-function(p, type='uniform',rate=1, cutoff=1){
  stopifnot(type %in% c('uniform','exponential'))

  if(type=='uniform'){
    runif(p,0,0.5)
  }else if(type=='exponential'){
    es<-rexp(p,rate = 1)
    # es/max(es)
    0.5* es/ (cutoff*max(es))
  }
}

# Xsim<-function(n=100,p=1000, maf){
#   X<-cbind(sapply(1:p,function(i) sample(c(-1,+1),size=n,replace=T,prob = c(1-maf[i],maf[i] ) )))
# return(X)
# }
Xsim<-function(n=100,p=1000, maf){
  X<-cbind(sapply(1:p,function(i) sample(c(0,1),size=n,replace=T,prob = c(1-maf[i],maf[i] ) )))
return(X)
}

XsimLD<-function(n=100,p=1000, maf,r2=0.4){
  require(MASS)
  R<-diag(p)
  R[R==1]<-1
  R[R==0]<-r2
  X_<- ((mvrnorm(n = n,Sigma = R,mu=maf,empirical = T)) )
  X_<- lapply(1:p,function(i) {X_[,i]>maf[i]} ) %>% do.call(what = cbind,.)
  # X_<- ((mvrnorm(n = n,Sigma = R,mu=rep(0.5,p),empirical = T)) > 0.5)*1
  # X_<- lapply(1:p,function(i) {X_[,i]>0.5} ) %>% do.call(what = cbind,.)
  X_=X_*1
  cor(X_)
return(X_)
}


#
# Specify the parameters.
#


PropoS<-function(n, m=1,svar=0.1){
  s = exp( rnorm(n,0,svar) ) - 1
 return(s)
}

ssim<-function(nsnp,svar=0.1){
  exp( rnorm(nsnp,0,svar) ) - 1
}


wsim<-function(X,s,mode=1, epi=1,mu=1){
  wC(X,s,mode,epi,mu)
}


sampleW<-function(Eys,a,b,p,rep=1){
  Yobs<-c()
  for(i in 1:length(Eys)){
    Yobs<-c(Yobs,rnorm(rep,Eys[i], abs(a+(Eys[i]*b)) ))
  }
  Yobs[Yobs<0] <-0
  if(p!=0){
    Yobs[sample(1:length(Yobs),size = ceiling(p*length(Yobs)) ) ]<-0
  }
return(Yobs)
}

sampleW2<-function(Eys,a,b,p,rep=1){
  Yobs<-c()
  for(i in 1:length(Eys)){
    Yobs<-c(Yobs, Eys[i]+rnorm(rep,0, abs(a+(Eys[i]*b)) ))
  }
  Yobs[Yobs<0] <-0
  if(p!=0){
    Yobs[sample(1:length(Yobs),size = ceiling(p*length(Yobs)) ) ]<-0
  }
return(Yobs)
}




#####**********************************************************************#####
multievo<-function(p,
                   n,
                   rate,
                   svar,
                   a,
                   b,
                   mu,
                   mode,
                   No,
                   Nmax,
                   d,
                   tmax,
                   Replicates=2){

    pop<-lapply(1:Replicates,function(i){
                maf=mafsim(p,'exponential',1)
                X=Xsim(n,p,maf)
                s=PropoS(p,svar)
                w=  wC(X,s,mode,1,mu)

                po<-multigenpopsimC(fitness=w,
                                     No=No,
                                     Nmax=Nmax,
                                     a=a,
                                     b=b,
                                     t=tmax,
                                     d=d
                                     )

                return(po)
      })

  return(pop)
}




#####**********************************************************************#####

popsizeplot<-function(obj){

  dat_<-sapply(obj,function(x) apply(as.matrix(x),2,sum) )
  dat<-data.frame(popsize=as.numeric(dat_), generations=1:nrow(dat_), replicate=sort(rep(1:ncol(dat_), nrow(dat_)))  )

spline_int <- as.data.frame(spline(dat$popsize~ dat$generations))

p<-ggplot(dat) +
  geom_line(aes(y=popsize,x=generations,group=replicate), color='grey')+
  geom_line(data = spline_int, aes(y=y,x=x),lwd=2)+
      ylab('N(t)') +
      xlab("Generations")
  return(p)
}


####  Allele frequencies
allelefreqplot<-function(obj){

  dat_<-sapply(obj,function(x) apply(as.matrix(x),2,sum) )
  dat<-data.frame(popsize=as.numeric(dat_), generations=1:nrow(dat_), replicate=sort(rep(1:ncol(dat_), nrow(dat_)))  )

spline_int <- as.data.frame(spline(dat$popsize~ dat$generations))

p<-ggplot(dat) +
  geom_line(aes(y=popsize,x=generations,group=replicate), color='grey')+
  geom_line(data = spline_int, aes(y=y,x=x),lwd=2)
  # stat_summary(fun.y = mean, geom="line", aes(y=popsize,x=generations))
  # stat_summary(aes(y=popsize,x=generations))
  # stat_smooth( method = lm, formula = y ~ poly(x, nrow(dat_)-1 ), aes(y=popsize,x=generations))
  # stat_smooth(aes(y=popsize,x=generations), se=F,col='black',span=10)

  # p<-ggplot() +
  #     ylim(c(0,max(as.numeric(as.matrix(dat))))) +
  #     xlim(c(1, input$tmax))+
  #     ylab('N(t)') +
  #     xlab("Generations")
  # for(rep in 1:ncol(dat)){
  #   newdat<-data.frame(y=dat[,rep],
  #                      x=1:input$tmax)
  #   p <- p+ geom_line(data=newdat,aes(y=y,x=x))
  # }
  return(p)
}
