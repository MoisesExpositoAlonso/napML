saccuracy <-function(s,sinf){
  ## Bias and accuracy at the individual level
  lmobj<-summary(lm(s~sinf))
  accuracy<-lmobj$r.squared %>% round(.,digits=3)

  if(dim(lmobj$coefficients)[1] ==1){
    bias2<-"na"
  }else{
    bias1<-coefficients(lmobj)[2,1]  %>% round(.,digits=3)
    bias2<-coefficients(lmobj)[2,2]  %>% round(.,digits=3)
  }
  return(c("R2"=accuracy,
              "b"=bias1,
              "a"=bias2)
         )
}

iaccuracy<-function(y,h,inf){
    ## Prepare data
    yinf=inf[h]

  ## Bias and accuracy at the selsection level
    lmobj<-summary(lm(y~yinf))
    accuracy<-lmobj$r.squared %>% round(.,digits=3)
    if(dim(lmobj$coefficients)[1] ==1){
      bias<-"na"
    }else{
      bias1<-coefficients(lmobj)[2,1]  %>% round(.,digits=3)
      bias2<-coefficients(lmobj)[2,2]  %>% round(.,digits=3)
    }

  return(c("R2"=accuracy,
            "b"=bias1,
            "a"=bias2)
       )
}




scorplot <-function(s,sinf,sinf_range=NULL){

    ## Bias and accuracy at the selsection level
    lmobj<-summary(lm(s~sinf))
    # lmobj<-summary(lm(sinf ~s))
    accuracy<-lmobj$r.squared %>% round(.,digits=3)
    if(dim(lmobj$coefficients)[1] ==1){
      bias<-"na"
    }else{
      bias<-coefficients(lmobj)[2,1]  %>% round(.,digits=3)
      bias2<-coefficients(lmobj)[2,2]  %>% round(.,digits=3)
    }

    if( is.null(sinf_range)) sinf_range =cbind(sinf,sinf)

    ## PLot selection coefficients
  psel<-qplot(y=s,x=sinf,
              xlim=c(range(c(s,sinf_range))),
              ylim=c(range(c(s,sinf_range)))
              ) +
  # geom_segment(aes(x=s,xend=s,y=sinf_range[,1],yend=sinf_range[,2]) )+
      geom_abline(intercept = 0,slope = 1,lty="dotted")+
      xlab("Inferred selection coefficients")+
      ylab("True selection coefficients")+
      ggtitle(TeX(paste("$R^2$ = ",accuracy, ", $\\beta = $",bias,", $\\a = $",bias2 )))
    # print(psel)
    if(!is.null(sinf_range)) {
      psel<-psel+geom_segment(aes(x=s,xend=s,y=sinf_range[,1],yend=sinf_range[,2]) )
    }
  return(list(psel=psel,accuracy=accuracy,bias=bias))
}

ssummary<-function(parchain){
  require(coda)
  r<-as.mcmc(parchain)
  colnames(r) <- paste0("SNP",1:dim(parchain)[2])
  summary(r)
}


shist <- function(parchain){

  if(ncol(parchain)>5) parchain=parchain[,1:5 ]

  parchain<-data.frame(parchain)
  colnames(parchain)<- paste0("SNP",1:dim(parchain)[2])

  pnames<-paste0("SNP",1:dim(parchain)[2])

  plotlist<-list()

  for(i in pnames){
    panel<-plot_grid(
      qplot(parchain[,i],
            x=1:nrow(parchain),
            # ylim=c(0,1),
            geom='line',xlab='iterations',ylab=i), #+geom_hline(yintercept = sdat[[i]], col='grey',lty='dashed')   ,
      qplot(parchain[,i],geom="density",
            xlab=i,
            # xlim=c(0,1),
            fill=I(transparent("black"))) #+geom_vline(xintercept = sdat[[i]], lty="dashed",col='grey')

    )
    plotlist[[i]]<-panel
  }

  bigpanel<-plot_grid(plotlist=plotlist,ncol=1)
  return(bigpanel)
}


indplot<-function(y,inf,h=NULL){
  if(is.null(h)) h=1:length(y)
  ## Prepare data
  toplot<-data.frame(x=inf[h], y=y)
  # ## Bias and accuracy at the individual level
  # lmobj2<-summary(lm(toplot$y~toplot$x))
  # R2<-lmobj2$r.squared %>% round(.,digits=3)
  #
  # if(dim(lmobj2$coefficients)[1] ==1){
  #   a<-"na"
  # }else{
  #   a<-coefficients(lmobj2)[1,1]  %>% round(.,digits=3)
  #   b<-coefficients(lmobj2)[2,1]  %>% round(.,digits=3)
  # }
  ac<-accuracies(y = toplot$y, x=toplot$x) %>% round(.,digits=3)
  pind<-ggplot(data=toplot,aes(y=y,x=x)) +
    geom_point(col="grey") +
    stat_smooth(aes(y=y , x=x), se=F,
                method="glm",lty="dashed",col="black")+
          ylim(range(c(toplot$x,toplot$y))) +
          xlim(range(c(toplot$x,toplot$y)))+
          xlab("Inferred individual fitness")+
          ylab("True individual fitness")+
      geom_abline(intercept = 0,slope = 1,lty="dotted")+
      geom_hline(yintercept = 0,lty="dotted")+
      geom_vline(xintercept = 0,lty="dotted")+
      ggtitle(TeX(paste("$R^2$ = ",ac["R2"],
                        ", $\\beta = $",ac["b"],
                        ", $\\a = $",ac["a"] )))
    return(list(pind=pind,accuracy=ac))
}

# iplotreal<-function(true, inf){
#
#   toplot<-data.frame(x=true, y=inf)
#   toplot$xnozero <- toplot$x
#   toplot$xnozero[toplot$xnozero==0 ] <- NA
#
#   return(iplot_(toplot))
# }
#
# iplot_<-function(toplot){
#     ## Bias and accuracy at the individual level
#     lmobj2<-summary(lm(toplot$y~toplot$xnozero))
#     # lmobj2<-summary(lm(toplot$xnozero ~toplot$y))
#     accuracy2<-lmobj2$r.squared %>% round(.,digits=3)
#     if(dim(lmobj2$coefficients)[1] ==1){
#     bias2<-"na"
#     }else{
#     bias2<-coefficients(lmobj2)[2,1]  %>% round(.,digits=3)
#     }
#
#
#     pind<-ggplot(data=toplot,aes(y=y,x=x)) +
#       geom_point(col="grey") +
#       stat_smooth(aes(y=y , x=xnozero), se=F,
#                   method="glm",lty="dashed",col="black")+
#             ylim(range(c(toplot$x,toplot$y))) +
#             xlim(range(c(toplot$x,toplot$y)))+
#             ylab("Inferred individual fitness")+
#             xlab("True individual fitness")+
#         geom_abline(intercept = 0,slope = 1,lty="dotted")+
#         geom_hline(yintercept = 0,lty="dotted")+
#         geom_vline(xintercept = 0,lty="dotted")+
#         ggtitle(TeX(paste("$R^2$ = ",accuracy2, ", $\\beta = $",bias2)))
#   return(list(pind=pind,accuracy2=accuracy2,bias2=bias2))
# }
parnames<-function() return(c("b","a","p","mu","epi","svar","ss"))

parametersummary<-function(parchain,pnames=c("b","a","p","mu","epi","svar","ss")){
  require(coda)
  r<-as.mcmc(parchain)
  colnames(r) <- pnames
  summary(r)
}

paramplot<-function(parchain,pnames=c("b","a","p","mu","epi","svar","ss"),
                    truevalues=NULL){

  parchain<-data.frame(parchain)
  colnames(parchain)<-pnames

  plotlist<-list()

  for(i in pnames){

    p1<-qplot(parchain[,i],
            x=1:nrow(parchain),
            # ylim=c(0,1),
            geom='line',xlab='iterations',ylab=i) #+geom_hline(yintercept = sdat[[i]], col='grey',lty='dashed')   ,
    p2<-qplot(parchain[,i],geom="density",
            xlab=i,
            # xlim=c(0,1),
            fill=I(transparent("black"))) #+geom_vline(xintercept = sdat[[i]], lty="dashed",col='grey')

    if(!is.null(truevalues)){
      p2+geom_vline(xintercept = truevalues[i],color="red")
    }

    plotlist[[i]]<-plot_grid(p1,p2)
  }

  bigpanel<-plot_grid(plotlist=plotlist,ncol=1)

  return(bigpanel)
}

plotNAP<-function(x){
  p1<-paramplot(x$parchain)
  p0<-posteriorplot
  if(ncol(x$chain)>5) toplot=x$chain[,1:5]
  else toplot=x$chain
  p2<-shist(x$chain)

  plot_grid(
            plot_grid(p0,p2,ncol=1),
            p1,
            ncol=2
  )
}

posteriorplot<-function(post,mylab="Posterior"){
  qplot(post,
            x=1:nrow(post),
            geom='line',xlab='iterations',ylab=mylab) #+
      # geom_vline(xintercept=(1:nrow(x$accept))[x$accept==1], col=transparent("darkgreen",0.2))
}
