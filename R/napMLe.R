###### Optimization
napoptime<-function(s,lik.nap,y,hall,m,G,mod,maxit = 1000){
  parstart<-list(
                "s"=s,
                "b"=0.5,
                "a"=0.5,
                "p"=0.1,
                "e"=1
                ) %>% unlist
  parlow<- list(
                "s"=s-sd(s)*2,
                "b"=0.0,
                "a"=0.0,
                "p"=0,
                "e"=0.5
                ) %>% unlist
  parhigh<- list(
                "s"=s+sd(s)*2,
                "b"=1,
                "a"=1,
                "p"=1,
                "e"=1.5
                )%>% unlist

  cat("Likelihood with starting parameters ")
  cat(lik.nap(y=y,
          h_=hall,
          m_=m,
          A=G,
          par=c(s,parstart["b"],parstart["a"],parstart["p"]),
          mod=mod,
          e=parstart["e"]
          )
  )
  cat("\n")
  cat("Spectral Projected Gradient method SPG ... \n")
  r1 <- tryCatch(
                  # spg( # wrappter for spg
                  BBoptim( # wrappter for spg
                      par = parstart ,
                      fn=lik.nap.e,
                      y=y,
                      h_=hall,
                      m_=m,
                      A=G,
                      mod=mod,
                      lower = parlow,
                      upper=parhigh,
                      control=list(maxit=maxit)
                      )
                   ,
                   error = function(e){list(value = NA)}
                   )
  cat("Low-storage Broyden-Fletcher-Goldfarb-Shanno L-BFGS ... \n")
  r2 <- tryCatch(
                  lbfgs(
                      x0 = parstart ,
                  # optimx(
                  #     par = parstart ,
                      fn=lik.nap.e,
                      y=y,
                      h_=hall,
                      m_=m,
                      A=G,
                      mod=mod,
                      lower = parlow,
                      upper=parhigh,
                      control=list(maxit=maxit),
                      method="L-BFGS-B"
                      )
                   ,
                   error = function(e){list(value = NA)}
                   )
  results<-c(r1$value,r2$value)
  if(!all(is.na(results))){
    decidemode<-which(results == min(results,na.rm = T))
    cat(paste("Best model is ", c("spg","lbfgsb")[decidemode],"\n"))
    r<-list(r1,r2)[head(decidemode)]
  }else{
    # stop("Neither optimization routines worked!")
    r<-list(value = NA)
  }
  return(r)
}

####************************************************************************####
lik.nap.e<-function(y,h_,m_,A,par,mod,e,debug=F){
  h=h_
  m=m_
  w=wCBM(A = A@address,
            s=par[1:length(m)],
            mycols =m,
            myrows =h,
            mode=mod,
            epi=par[length(m)+5])
  LL<-likelihoodR(y = y,
            w = w,
            b=par[length(m)+1],
            a=par[length(m)+2],
            p=par[length(m)+3]
            )
  if(debug) print(LL)
  if(debug) print(par[seq(length(m)+1,length(par))])
  if(debug) print(par[1:5])

  return( - LL) # DO NOT FORGET MINUS
}

