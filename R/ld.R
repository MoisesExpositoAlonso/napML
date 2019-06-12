getD<-function(xs){
  apply(xs,1,function(i)(i[1]*i[4]) - (i[2]*i[3]))
}
getR<-function(xs){
  apply(xs,1,function(i)(i[1]*i[4]) / (i[2]*i[3]))
}
getRchange<-function(w){
  apply(w,1,function(i)(i[1]*i[4]) / (i[2]*i[3]))
}
getr2<-function(){
  p=(x1+x2)
  q=(x3+x4)
  
  r2=(x1*x4 - x2*x3)/sqrt(p*q*(1-p)*(1-q))
  r2
}
getR2change<-function(w){
  apply(w,1,function(i) {
    w_hat<- i[1]*x1+i[2]*x2+i[3]*x3+i[4]*x4
    p=(x1+x2)
    q=(x3+x4)
    R0=(x1*x4 - x2*x3)/sqrt(p*q*(1-p)*(1-q))
    x1p<-(i[1]*x1) /w_hat
    x2p<-(i[2]*x2) /w_hat
    x3p<-(i[3]*x3) /w_hat
    x4p<-(i[4]*x4) /w_hat
    p_prime=(x1p+x2p)
    q_prime=(x3p+x4p)
    R1=(x1p*x4p - x2p*x3p)/sqrt(p_prime*q_prime*(1-p_prime)*(1-q_prime))
    Rchange=R1-R0
  })
}

getDchange<-function(w){
  apply(w,1,function(i) {
    
    D1=(x1*x4) - (x2*x3)
    
    w_hat<- i[1]*x1+i[2]*x2+i[3]*x3+i[4]*x4
    x1p<-(i[1]*x1) /w_hat
    x2p<-(i[2]*x2) /w_hat
    x3p<-(i[3]*x3) /w_hat
    x4p<-(i[4]*x4) /w_hat
    D2=(x1p*x4p) - (x2p*x3p)
   
    D2-D1
  })
}


################################################################################

cleanld<-function(LD){

LD[LD>10]<-10
LD[LD< (-10)]<- (-10)
LD[LD==Inf]<-10
LD[LD == (-Inf)]<- (-10)
LD[is.na(LD)]<-1
LD<-LD[lower.tri(LD)]

return(LD)
}



ld.lik<-function(theta,y){

e=theta[1]
sel=theta[-1]

ldexpected<- ldCexpect(sel, e) %>% cleanld()

mu<-mean(ldexpected)

sigma2<-var(ldexpected)

n<-length(y)
logl<- -.5*n*log(2*pi) -.5*n*log(sigma2) -(1/(2*sigma2))*sum((y-mu)^2)
return(-logl)
}
