

#install.packages(c("sampling", rockchalk", "matrixcalc", "lme4","numDeriv", "RColorBrewer","xtable"))
packages<-c("sampling", "rockchalk", "matrixcalc", "lme4", "numDeriv", "RColorBrewer","xtable" )
lapply(packages, library,character.only = TRUE)
dyn.load("RIFourOrdPiTwostagePoisson.so")


#setting: notation
N1=100 #  number of sampling cluster in the first stage (population level)
N2=100 ##number of elements in each sampling cluster (population level)
latitude<-1:N2
longitude<-1:N1
population<-expand.grid(lat=latitude,long=longitude)
population$PSU<-population$long
overlap=ceiling(N2*3/4)

model_cluster<-function(population, overlap){
   population$cluster<-numeric(nrow(population))
   
   id<-ifelse(population$lat<=overlap, 
              population$long, 
              ((population$long+population$lat-overlap) %% N1)+1
   )
   population$cluster<-id
   population	
}

population<-model_cluster(population, overlap)
T=length(unique(population$cluster))


##check
population<-model_cluster(population, overlap)
table(population$cluster==population$PSU)
table(table(population$cluster))
table(table(population$PSU))

#Model: parameter from random intercept model
truebeta1=1
truebeta2=3
truesigma2=1
truetau2=0.5
truevalue<-c(truebeta1,truebeta2, truesigma2, truetau2)
names(truevalue)<-c("beta1", "beta2", "sigma2", "tau2")


##Population data
population$u<-rnorm(T,s=sqrt(truetau2))[population$cluster]
population$x<-rnorm(N1*N2)+rnorm(T)[population$cluster]
population$y<-with(population, truebeta1+truebeta2*x+u+rnorm(N1*N2,s=sqrt(truesigma2)))
population$r=with(population, x*(y-truebeta1-truebeta2*x))
population$ID_unit=with(population, 1:(N1*N2))

range(population$r)

#uninformative two-stage sampling design (first-stage: Poisson, Second-stage:SRSWOR)
##uninformative Poisson first-stage 
pi1=runif(N1) ## sampling inclusion probability for sampling cluster
FirststagePoisson=UPpoisson(pi1)
n1=sum(FirststagePoisson) # number of PSU in the first-stage sample

FirststagePoissonSample=subset(population, population$PSU%in% which(FirststagePoisson==1))

#uninformative SRSWOR second-stage  
n2=ceiling(N2/10)##number of elements in each sampling cluster ( sample level)
SecondstageSRSWOR=unlist(lapply(rep(n2,n1), function(v) return(srswor(v, N2))))
TwostagePoissonSample<-FirststagePoissonSample[c(which(SecondstageSRSWOR==1)), ]


#informative two-stage sampling design (first-stage: Poisson, Second-stage:SRSWOR)
#informative Poisson first-stage sample
param1=c(0.05, 0.45)
pi1informative= function(r, sc,  param, N1){
   a=rep(NA, N1)
   b=rep(NA, N1)
   for (i in 1: N1){
      a[i]=mean(r[sc==i])
      b[i]=(param[1]*exp(-param[2]*a[i]))/(1 +param[1]*exp(-param[2]*a[i]))
   }
   b}
pi1is=pi1informative(population$r,population$PSU, param1 ,N1)
FirststagePoissonis=UPpoisson(pi1is)
n1is=sum(FirststagePoissonis) # number of PSU in the first-stage sample

FirststagePoissonSampleis=subset(population, population$PSU%in% which(FirststagePoissonis==1))


#informative second-stage sample(SRSWOR)
#param2=c(0.05, 0.45)
#n2informative= function(r, sc, param, N2){
#   a=rep(NA, length=length(unique(population$sc)))
#   b=rep(NA, length=length(unique(population$sc)))
#   for (i in unique(sc)){
#      a[i]=mean(r[sc==i])
#      b[i]=2*ceiling((param[1]*exp(-param[2]*a[i]))/(1 +param[1]*exp(-param[2]*a[i]))*N2/2)
#   }
#   b
#}

#n2pop=n2informative(population$r,population$PSU, param2 ,N2)
#n2is=n2pop*FirststagePoissonis
#SecondstageSRSWORis=unlist(lapply(n2is[c(which(n2is!=0))], function(v) return(srswor(v, N2))))
#TwostagePoissonSampleis=FirststagePoissonSampleis[c(which(SecondstageSRSWORis==1)), ]

n2is=n2*FirststagePoissonis
SecondstageSRSWORis=unlist(lapply(n2is[c(which(n2is!=0))], function(v) return(srswor(v, N2))))
TwostagePoissonSampleis=FirststagePoissonSampleis[c(which(SecondstageSRSWORis==1)), ]



#check
summary(population$y)
summary(TwostagePoissonSample$y)
summary(TwostagePoissonSampleis$y)

#Estimation: full-likelihood
##Census estimator 
lmer(y~(1|cluster)+x,data=population)

##uninformative two-stage sampling design (Poisson)
lmer(y~(1|cluster)+x,data=TwostagePoissonSample)

##informative two-stage sampling design (Poisson)
lmer(y~(1|cluster)+x,data=TwostagePoissonSampleis)

# Estimation: pairwise likelihood (without weight)
l2<-function(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2){
   st<-sigma2+tau2
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   iftau<-ifelse(g1==g2,tau2,0)
   det<-st^2-iftau^2
   (-(r1*r1*st-2*r1*r2*iftau+r2*r2*st)/2/det-log(det)/2)
}	


dalpha<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2,tau2){
   st<-sigma2+tau2
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   iftau<-ifelse(g1==g2,tau2,0)
   det<-st^2-iftau^2
   ((r1*st-r2*iftau-r1*iftau+r2*st)/det)
}	

dbeta<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2,tau2){
   st<-sigma2+tau2
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   iftau<-ifelse(g1==g2,tau2,0)
   det<-st^2-iftau^2
   ((r1*st*x1-r2*iftau*x1-r1*iftau*x2+r2*st*x2)/det)
}	


dsigma2<-function(y1,y2, g1,g2, x1,x2,alpha, beta, sigma2,tau2){
   st<-sigma2+tau2
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   iftau<-ifelse(g1==g2,tau2,0)
   det<-st^2-iftau^2
   (-( (r1*r1+r2*r2)/2/det - (r1*r1*st-2*r1*r2*iftau+r2*r2*st)*(2*st)/2/det/det)- st/det)
}	


dtau2<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2){
   st<-sigma2+tau2
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   iftau<-ifelse(g1==g2,tau2,0)
   ifistau<-1*(g1==g2)
   det<-st^2-iftau^2
   ddet<-2*(st-iftau)
   (-( (r1*r1-2*r1*r2*ifistau+r2*r2)/2/det - (r1*r1*st-2*r1*r2*iftau+r2*r2*st)*(ddet)/2/det/det) -ddet/det/2)
}	



##Define the pairwise score function and checking the pairwise score at PML (without weight)
pairscore_PL<-function(y,g,x, theta){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   incrementda=dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                      sigma2=theta[3],tau2=theta[4])
   incrementdb=dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                     sigma2=theta[3],tau2=theta[4])
   incrementds=dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                       sigma2=theta[3],tau2=theta[4])
   incrementdt=dtau2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                     sigma2=theta[3],tau2=theta[4])
   
   c(sum(incrementda), sum(incrementdb), sum(incrementds),sum(incrementdt))/T
}

cenfit_PL<-function(y,g,x, pars){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   func1<-function(theta){
      increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                   sigma2=theta[3],tau2=theta[4])
      sum(increment)/T
   }
   gr<-function(theta){
      incrementda=dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                         sigma2=theta[3],tau2=theta[4])
      incrementdb=dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                        sigma2=theta[3],tau2=theta[4])
      incrementds=dsigma2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                          sigma2=theta[3],tau2=theta[4])
      incrementdt=dtau2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                        sigma2=theta[3],tau2=theta[4])
      c(sum(incrementda), sum(incrementdb), sum(incrementds), sum(incrementdt))/T
   }
   optim(pars,func1, gr, method="L-BFGS-B", lower=c(-Inf, -Inf, 0.1, 0),control=list(fnscale=-1,parscale=c(1/n^0.2,1/n^0.2,1/n^0.2, 1/n^0.2)))
}


###census PL estimator 
cen_PL<- cenfit_PL(y=population$y, g=population$cluster, x=population$x, pars=truevalue)
cen_PL[[1]]-truevalue
pairscore_PL(y=population$y,g=population$cluster, x=population$x, cen_PL[[1]])

#optimization problem for PL (without weight)
fit_PL<-function(y,g,x, pars){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   func1<-function(theta){
      increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                   sigma2=theta[3],tau2=theta[4])
      sum(increment)/T
   }
   gr<-function(theta){
      incrementda=dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                         sigma2=theta[3],tau2=theta[4])
      incrementdb=dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                        sigma2=theta[3],tau2=theta[4])
      incrementds=dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                         sigma2=theta[3],tau2=theta[4])
      incrementdt=dtau2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                        sigma2=theta[3],tau2=theta[4])
      c(sum(incrementda), sum(incrementdb), sum(incrementds), sum(incrementdt))/T
   }
   optim(pars,func1, gr, method="L-BFGS-B", lower=c(-Inf, -Inf, 0.1, 0) ,control=list(fnscale=-1, parscale=c(1/n^(1/2),1/n^(1/2),1/n^(1/2),1/n^(1/2))))
}


###uninformative PL
estimator_PL<-fit_PL(TwostagePoissonSample$y, TwostagePoissonSample$cluster, TwostagePoissonSample$x, pars=truevalue)
estimator_PL[[1]]-truevalue

###uninformative at the estimated value
pairscore_PL(y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x, estimator_PL[[1]])
###uninformative at the true value
pairscore_PL(y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x, truevalue)

fitis_PL<-function(y,g,x, pars){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   func1<-function(theta){
      increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                   sigma2=theta[3],tau2=theta[4])
      sum(increment)/T
   }
   gr<-function(theta){
      incrementda=dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                         sigma2=theta[3],tau2=theta[4])
      incrementdb=dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                        sigma2=theta[3],tau2=theta[4])
      incrementds=dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                          sigma2=theta[3],tau2=theta[4])
      incrementdt=dtau2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                        sigma2=theta[3],tau2=theta[4])
      c(sum(incrementda), sum(incrementdb), sum(incrementds), sum(incrementdt))/T
   }
   optim(pars,func1, gr, method="L-BFGS-B", lower=c(-Inf, -Inf, 0.1, 0),control=list(fnscale=-1))
}

###informative sampling (without weight)
estimatoris_PL<- fitis_PL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster, TwostagePoissonSampleis$x, pars=truevalue)
estimatoris_PL[[1]]-truevalue

###informative sampling at the estimated value 
pairscore_PL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster, TwostagePoissonSampleis$x, estimatoris_PL[[1]])
###informative sampling at the true value 
pairscore_PL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster, TwostagePoissonSampleis$x, truevalue)


# second-order inclusion probability
C2<-function(pos1, pos2,sc1, sc2,fss, n2infor,N2){
   .C("SecOrdPi",as.integer(pos1), as.integer(pos2),as.integer(sc1), as.integer(sc2), as.double(fss), as.double(n2infor),as.double(N2),length(pos1),rval=numeric(length(pos1)))$rval
}

# fourth-order inclusion probability
C4<-function(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,fss, n2infor,N2){
   .C("FourOrdPi",as.integer(pos1), as.integer(pos2),as.integer(pos3), as.integer(pos4),as.integer(sc1), as.integer(sc2),
      as.integer(sc3), as.integer(sc4), as.double(fss), as.double(n2infor),as.double(N2),length(pos1),rval=numeric(length(pos1)))$rval	
   
}

##Define the second-order inclusion probability
SecOrdPi<-function(pos1, pos2,sc1, sc2,fss, n2infor,N2){
   #pi<-SecOrdPiInternal(pos1, pos2,sc1, sc2,n1, N1, n2infor,N2)	
   Cpi<-C2(pos1, pos2,sc1, sc2,fss, n2infor,N2)
   #if ((pi-Cpi)/(pi+Cpi)>1e-10) stop(paste(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,":",pi,Cpi,sep=","))
   Cpi
}



##Define the  fourth-order inclusion probability
FouOrdPi<-function(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,fss, n2infor,N2){
   #pi<-FouOrdPiInternal(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,n2infor,N2)	
   Cpi<-C4(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,fss, n2infor,N2)
   #if ((pi-Cpi)/(pi+Cpi)>1e-10) stop(paste(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,":",pi,Cpi,sep=","))
   Cpi
}

#Define the fourth-order Delta
FouOrdDel=function(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,fss, n2infor, N2){
   FouOrdPi(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4, fss, n2infor,N2)-
      SecOrdPi(pos1, pos2,sc1, sc2, fss, n2infor,N2)*
      SecOrdPi(pos3, pos4,sc3, sc4, fss, n2infor,N2)
}

# Estimation: weighted pairwise likeliood 
wl2<-function(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2, pos1, pos2,sc1, sc2, fss,  n2infor,N2){
   1/SecOrdPi(pos1, pos2,sc1, sc2, fss,  n2infor, N2)*(l2(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2))
}	

pairscore_WPL<-function(y,g,x, theta, pos, sc,fss,  n2infor, N2){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   wij<-1/SecOrdPi(pos[i], pos[j],sc[i], sc[j], fss,  n2infor,N2)
   incrementda=wij*dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                          sigma2=theta[3],tau2=theta[4])
   incrementdb=wij*dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                         sigma2=theta[3],tau2=theta[4])
   incrementds=wij*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                           sigma2=theta[3],tau2=theta[4])
   incrementdt=wij*dtau2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                         sigma2=theta[3],tau2=theta[4])
   
   
   c(sum(incrementda), sum(incrementdb), sum(incrementds),sum(incrementdt))/T
}

#optimization (WPL)
fit_WPL<-function(y,g,x, pos, sc,fss,  n2infor, N2,  pars){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   func1<-function(theta){
      wincrement=wl2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                     sigma2=theta[3],tau2=theta[4], pos[i], pos[j], sc[i], sc[j], fss,  n2infor,N2)
      sum(wincrement)/T
   }
   gr<-function(theta){
      wij<-1/SecOrdPi(pos[i], pos[j],sc[i], sc[j], fss,  n2infor,N2)
      wincrementda=wij*dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                              sigma2=theta[3],tau2=theta[4])
      wincrementdb=wij*dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                             sigma2=theta[3],tau2=theta[4])
      wincrementds=wij*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                              sigma2=theta[3],tau2=theta[4])
      wincrementdt=wij*dtau2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                              sigma2=theta[3],tau2=theta[4])
      c(sum(wincrementda), sum(wincrementdb), sum(wincrementds), sum(wincrementdt))/T
   }
   optim(pars,func1,gr,  method="L-BFGS-B", lower=c(-Inf, -Inf, 0.1, 0), control=list(fnscale=-1,parscale=c(1/n^(2/3),1/n^(2/3),1/n^(2/3),1/n^(2/3))))
}

##uninformative sampling WPL (with weight)
estimator_WPL=fit_WPL(y=TwostagePoissonSample$y, g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
                      pos=TwostagePoissonSample$ID_unit, sc=TwostagePoissonSample$PSU, fss=pi1, n2infor=FirststagePoisson*n2, N2, 
                      pars=truevalue)
estimator_WPL[[1]]-truevalue

## WPL uninformative sampling at the estimated value 
pairscore_WPL(TwostagePoissonSample$y, TwostagePoissonSample$cluster, TwostagePoissonSample$x,estimator_WPL[[1]],pos=TwostagePoissonSample$ID_unit,
              sc=TwostagePoissonSample$PSU, fss=pi1, n2infor=FirststagePoisson*n2,N2)
## WPL uninformative sampling at the true value
pairscore_WPL(TwostagePoissonSample$y, TwostagePoissonSample$cluster, TwostagePoissonSample$x,truevalue,pos=TwostagePoissonSample$ID_unit,
              sc=TwostagePoissonSample$PSU, fss=pi1, n2infor=FirststagePoisson*n2,N2)


fitis_WPL<-function(y,g,x, pos, sc,fss,  n2infor, N2,  pars){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   func1<-function(theta){
      wincrement=wl2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                     sigma2=theta[3],tau2=theta[4], pos[i], pos[j], sc[i], sc[j], fss,  n2infor,N2)
      sum(wincrement)/T
   }
   gr<-function(theta){
      wij<-1/SecOrdPi(pos[i], pos[j],sc[i], sc[j], fss,  n2infor,N2)
      wincrementda=wij*dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                              sigma2=theta[3],tau2=theta[4])
      wincrementdb=wij*dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                             sigma2=theta[3],tau2=theta[4])
      wincrementds=wij*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                               sigma2=theta[3],tau2=theta[4])
      wincrementdt=wij*dtau2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                             sigma2=theta[3],tau2=theta[4])
      c(sum(wincrementda), sum(wincrementdb), sum(wincrementds), sum(wincrementdt))/T
   }
   optim(pars,func1,gr,  method="L-BFGS-B", lower=c(-Inf, -Inf, 0.1, 0),
         control=list(fnscale=-1,parscale=c(1/n^2,1/n^2,1/n^0.5,1/n^0.5)))
}


##informative sampling WPL (with weight)
estimatoris_WPL=fitis_WPL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster,TwostagePoissonSampleis$x,
                        TwostagePoissonSampleis$ID_unit, TwostagePoissonSampleis$PSU, fss=pi1is,  n2infor=n2is , N2,  
                        pars=truevalue)
estimatoris_WPL[[1]]-truevalue

## WPL informative sampling at the estimated value
pairscore_WPL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster,TwostagePoissonSampleis$x,estimatoris_WPL[[1]],  pos=TwostagePoissonSampleis$ID_unit, sc=TwostagePoissonSampleis$PSU, 
              fss=pi1is,  n2infor=n2is,N2)
## WPL informative sampling at the true value 
pairscore_WPL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster,TwostagePoissonSampleis$x,truevalue,  pos=TwostagePoissonSampleis$ID_unit, sc=TwostagePoissonSampleis$PSU, 
              fss=pi1is,  n2infor=n2is,N2)



#variance estimation for PL under stratified SRSWORS
#uninformative 
#Calculate Hessian matrix H for PL (bread for uninformative sampling design)
estH_PL=-jacobian(function(theta){with(TwostagePoissonSample,
                                       pairscore_PL(y,cluster,x,theta))}, x=estimator_PL[[1]],method="complex")



#Calculate  variance matrix J  for PL (meat for uninformative sampling design)
fast_J_PL<-function(y,g,x,pos, sc,fss, n2infor,N2, theta){
   n<-length(y)
   sum=0
   
   kl<-expand.grid(1:n,1:n)
   kl<-kl[kl[,1]<kl[,2],]
   kl<-kl[g[kl[,1]]==g[kl[,2,]],]
   k<-kl[,1]
   l<-kl[,2]
   
   for (i in 1:(n-1)){
      cat(i)
      js <- (i+1):n
      js <- js[g[js] %in% g[i]]
      for(j in js){
         cat(".")
         incrementdaij=dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=theta[3],
                              tau2=theta[4])
         incrementdbij=dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],sigma2=theta[3],
                             tau2=theta[4])
         incrementdsij=dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],sigma2=theta[3],
                              tau2=theta[4])
         incrementdtij=dtau2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=theta[3],
                              tau2=theta[4])
         psij=c(incrementdaij, incrementdbij, incrementdsij, incrementdtij)
         
         ## k,l vectorised: probably can't afford memory to do that for ijkl 
         ii <-rep(i, length(k))
         jj<-rep(j,length(k))
         incrementdakl=dalpha(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                              sigma2=theta[3],tau2=theta[4])
         incrementdbkl=dbeta(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                             sigma2=theta[3],tau2=theta[4])
         incrementdskl=dsigma2(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                              sigma2=theta[3],tau2=theta[4])
         incrementdtkl=dtau2(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                              sigma2=theta[3],tau2=theta[4])
         pskl=cbind(incrementdakl, incrementdbkl, incrementdskl, incrementdtkl)
         sumpskl<-colSums( FouOrdDel(pos[ii], pos[jj], pos[k], pos[l], sc[ii], sc[jj], sc[k], sc[l],fss,n2infor,N2)*(FouOrdPi(pos[ii], pos[jj], pos[k], pos[l], sc[ii], sc[jj], sc[k], sc[l],fss,n2infor,N2))^(-1)* pskl)
         psijkl<-tcrossprod(psij,sumpskl)
         sum=sum+psijkl
      }
   }
   rval<-sum/(T^2)
   ##attr(rval, "pairs")<-keep ##debug
   rval
}

estJ_PL=fast_J_PL(y=TwostagePoissonSample$y, g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
                  pos=TwostagePoissonSample$ID_unit, sc=TwostagePoissonSample$PSU,fss=pi1,
                  n2infor=FirststagePoisson*n2, N2, theta=estimator_PL[[1]] )

#sanwich estimator (uninformative sampling )
sanestimator_PL= solve(estH_PL)%*% estJ_PL%*% solve(t(estH_PL))

#Informative
#Calculate Hessian matrix H for PL (bread for informative sampling design)
estHis_PL=-jacobian(function(theta){with(TwostagePoissonSampleis,
                                         pairscore_PL(y,cluster,x,theta))}, x=estimatoris_PL[[1]],method="complex")


#Calculate  variance matrix J  for PL (meat for informative sampling design)
estJis_PL=fast_J_PL(y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster,x=TwostagePoissonSampleis$x,pos=TwostagePoissonSampleis$ID_unit,  
                    sc=TwostagePoissonSampleis$PSU, fss=pi1is, n2infor=n2is, N2, theta=estimatoris_PL[[1]] )
#sanwich estimator (informative sampling )
sanestimatoris_PL = solve(estHis_PL)%*% estJis_PL%*% t(solve(estHis_PL))


#variance estimation for WPL under two-stage SRSWORS
##define H as in page 96 of my thesis as \hat{H}(\est)
#define weighted pairwise likelihood WPL 

##Hessian for WPL uninformative sampling
estH_WPL=-jacobian(function(theta){with(TwostagePoissonSample,
                                        pairscore_WPL(y,cluster,x,theta,ID_unit,PSU, pi1, FirststagePoisson*n2, N2))}, x=estimator_WPL[[1]],method="complex")

##define \hat{J}(\theta) as in page 97 of my thesis and  evaluate at the WPLE
fast_J_WPL<-function(y,g,x,  pos,  sc, fss, n2infor,N2, theta){
   n<-length(y)
   sum=0
   
   kl<-expand.grid(1:n,1:n)
   kl<-kl[kl[,1]<kl[,2],]
   kl<-kl[g[kl[,1]]==g[kl[,2,]],]
   k<-kl[,1]
   l<-kl[,2]
   
   for (i in 1:(n-1)){
      cat(i)
      js <- (i+1):n
      js <- js[g[js] %in% g[i]]
      for(j in js){
         cat(".")
         wij<-1/SecOrdPi(pos[i], pos[j],sc[i], sc[j], fss, n2infor,N2)
         incrementdaij=wij*dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=theta[3],
                                  tau2=theta[4])
         incrementdbij=wij*dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],sigma2=theta[3],
                                 tau2=theta[4])
         incrementdsij=wij*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],sigma2=theta[3],
                                  tau2=theta[4])
         incrementdtij=wij*dtau2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=theta[3],
                                  tau2=theta[4])
         wpsij=c(incrementdaij, incrementdbij, incrementdsij, incrementdtij)
         
         ## k,l vectorised: probably can't afford memory to do that for ijkl 
         ii <-rep(i, length(k))
         jj<-rep(j,length(k))
         wkl<-1/SecOrdPi(pos[k], pos[l],sc[k], sc[l], fss, n2infor,N2)
         incrementdakl=wkl*dalpha(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                  sigma2=theta[3],tau2=theta[4])
         incrementdbkl=wkl*dbeta(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                                 sigma2=theta[3],tau2=theta[4])
         incrementdskl=wkl*dsigma2(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                                 sigma2=theta[3],tau2=theta[4])
         incrementdtkl=wkl*dtau2(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                 sigma2=theta[3],tau2=theta[4])
         wpskl=cbind(incrementdakl, incrementdbkl, incrementdskl, incrementdtkl)
         sumwpskl<-colSums( (1/FouOrdPi( pos[ii], pos[jj], pos[k], pos[l], sc[ii], sc[jj], sc[k], sc[l], fss,  n2infor,N2))*FouOrdDel(pos[ii], pos[jj], pos[k], pos[l], sc[ii], sc[jj], sc[k], sc[l],fss,  n2infor,N2)* wpskl)
         wpsijkl<-tcrossprod(wpsij,sumwpskl)
         sum=sum+wpsijkl
      }
   }
   rval<-sum/(T^2)
   # attr(rval, "pairs")<-keep ##debug
   rval
}

estJ_WPL=fast_J_WPL(y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster,
                    x=TwostagePoissonSample$x, pos=TwostagePoissonSample$ID_unit,  sc=TwostagePoissonSample$PSU, fss=pi1, 
                    n2infor= FirststagePoisson*n2, N2, theta=estimator_WPL[[1]] )

# sanwich estimator H^{-1}J (H^{-1})^\T
##uninformaitve
sanestimator_WPL= solve(estH_WPL)%*% estJ_WPL%*% t(solve(estH_WPL))

##informative sampling
estHis_WPL=-jacobian(function(theta){with(TwostagePoissonSampleis,
                                          pairscore_WPL(y,cluster,x,theta,ID_unit,PSU, pi1is, n2is, N2))}, x=estimatoris_WPL[[1]],method="complex")
estJis_WPL=fast_J_WPL(y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster,
                      x=TwostagePoissonSampleis$x, pos=TwostagePoissonSampleis$ID_unit,  sc=TwostagePoissonSampleis$PSU, fss=pi1is, 
                      n2infor= n2is, N2, theta=estimatoris_WPL[[1]] )
# sanwich estimator H^{-1}J (H^{-1})^\T
##informaitve
sanestimatoris_WPL= solve(estHis_WPL)%*% estJis_WPL%*% t(solve(estHis_WPL))


#simulation

LOTS=150
#Fit from NML,PL, WPL for uninformative sampling
Fit_NML<-matrix(0,nrow=LOTS,ncol=4)
Fit_PL<-matrix(0,nrow=LOTS,ncol=4)
Fit_WPL<-matrix(0,nrow=LOTS, ncol=4)

#Fit from NML, PL, WPL for informative sampling
Fitis_NML<-matrix(0,nrow=LOTS,ncol=4)
Fitis_PL<-matrix(0,nrow=LOTS,ncol=4)
Fitis_WPL<-matrix(0,nrow=LOTS, ncol=4)

#Hessian matrix for PL(without weight) for uninformative sampling
H_PL<-array(0, c(4,4, LOTS))

#Hessian matrix for PL(without weight) for informative sampling
His_PL<-array(0, c(4,4, LOTS))

#Hessian matrix for WPL for uninformative sampling
H_WPL<-array(0, c(4,4, LOTS))

#Hessian matrix for WPL for informative sampling
His_WPL<-array(0, c(4,4, LOTS))

#Variance matrix J for PL for uninformative sampling
J_PL<-array(0, c(4,4, LOTS))

#Variance matrix J for PL for informative sampling
Jis_PL<-array(0, c(4,4, LOTS))

#Variance matrix J for WPL for uninformative sampling
J_WPL<-array(0, c(4,4, LOTS))

#Variance matrix J for WPL for informative sampling
Jis_WPL<-array(0, c(4,4, LOTS))

#Sanwich variance estimator for PL for uninformative sampling
G_PL<-array(0, c(4,4, LOTS))

#Sanwich variance estimator for PL for informative sampling
Gis_PL<-array(0, c(4,4, LOTS))

#Sanwich variance estimator for WPL for uninformative sampling
G_WPL<-array(0, c(4,4, LOTS))

#Sanwich variance estimator for WPL for informative sampling
Gis_WPL<-array(0, c(4,4, LOTS))

#Pairwise score function for PL for informative sampling at the estimated value 
PS_PL<-matrix(0,nrow=LOTS,ncol=4)

#Pairwise score function for PL for informative sampling at the true value 
PS_PL_true<-matrix(0,nrow=LOTS,ncol=4)

#Pairwise score function for PL for informative sampling at the census value 
PS_PL_cen<-matrix(0,nrow=LOTS,ncol=4)

#Pairwise score function for PL for  informative sampling at the estimated value
PSis_PL<-matrix(0,nrow=LOTS,ncol=4)

#Pairwise score function for PL for  informative sampling at the true value 
PSis_PL_true<-matrix(0,nrow=LOTS,ncol=4)

#Pairwise score function for PL for  informative sampling at the census value 
PSis_PL_cen<-matrix(0,nrow=LOTS,ncol=4)

#Pairwise score function for WPL for informative sampling at the estimated value 
PS_WPL<-matrix(0,nrow=LOTS,ncol=4)

#Pairwise score function for WPL for informative sampling at the true value 
PS_WPL_true<-matrix(0,nrow=LOTS,ncol=4)

#Pairwise score function for WPL for informative sampling at the census value 
PS_WPL_cen<-matrix(0,nrow=LOTS,ncol=4)

#Pairwise score function for WPL for  informative sampling at the estimated value 
PSis_WPL<-matrix(0,nrow=LOTS,ncol=4)

#Pairwise score function for WPL for  informative sampling at the true value 
PSis_WPL_true<-matrix(0,nrow=LOTS,ncol=4)

#Pairwise score function for WPL for  informative sampling at the census value 
PSis_WPL_cen<-matrix(0,nrow=LOTS,ncol=4)



#define the squre root of J
sqrt_diagJ_PL=matrix(0, nrow=LOTS, ncol=4)
sqrt_diagJis_PL=matrix(0, nrow=LOTS, ncol=4)
sqrt_diagJ_WPL=matrix(0, nrow=LOTS, ncol=4)
sqrt_diagJis_WPL=matrix(0, nrow=LOTS, ncol=4)


#define the squre root of G
sqrt_diagG_PL=matrix(0, nrow=LOTS, ncol=4)
sqrt_diagGis_PL=matrix(0, nrow=LOTS, ncol=4)
sqrt_diagG_WPL=matrix(0, nrow=LOTS, ncol=4)
sqrt_diagGis_WPL=matrix(0, nrow=LOTS, ncol=4)

##Estimation: NML, PL and WPL 
for(i in 1:LOTS){
   
   cat(i)

   ##uninformative Poisson first-stage 
   pi1=runif(N1) ## sampling inclusion probability for sampling cluster
   FirststagePoisson=UPpoisson(pi1)
   n1=sum(FirststagePoisson) # number of PSU in the first-stage sample
   
   FirststagePoissonSample=subset(population, population$PSU%in% which(FirststagePoisson==1))
   
   #uninformative SRSWOR second-stage  
   n2=ceiling(N2/10)##number of elements in each sampling cluster ( sample level)
   SecondstageSRSWOR=unlist(lapply(rep(n2,n1), function(v) return(srswor(v, N2))))
   TwostagePoissonSample<-FirststagePoissonSample[c(which(SecondstageSRSWOR==1)), ]
   

   #informative two-stage sampling design (first-stage: Poisson, Second-stage:SRSWOR)
   ##informative Poisson first-stage sample
   pi1is=pi1informative(population$r,population$PSU, param1 ,N1)
   FirststagePoissonis=UPpoisson(pi1is)
   n1is=sum(FirststagePoissonis) # number of PSU in the first-stage sample
   
   FirststagePoissonSampleis=subset(population, population$PSU%in% which(FirststagePoissonis==1))
   

   #n2pop=n2informative(population$r,population$PSU, param2 ,N2)
   n2is=n2*FirststagePoissonis
   SecondstageSRSWORis=unlist(lapply(n2is[c(which(n2is!=0))], function(v) return(srswor(v, N2))))
   TwostagePoissonSampleis=FirststagePoissonSampleis[c(which(SecondstageSRSWORis==1)), ]
   
   
   #NML, PL and WPL (uninformative sampling)
   ra<-lmer(y~(1|cluster)+x,data=TwostagePoissonSample)
   rb<-fit_PL(TwostagePoissonSample$y, TwostagePoissonSample$cluster, TwostagePoissonSample$x, pars=truevalue)
   rc<-fit_WPL(y=TwostagePoissonSample$y, g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
               pos=TwostagePoissonSample$ID_unit, sc=TwostagePoissonSample$PSU, fss=pi1, n2infor=FirststagePoisson*n2, N2, 
               pars=truevalue)
   
   #NML, PL and WPL (informative sampling)
   rais<-lmer(y~(1|cluster)+x,data=TwostagePoissonSampleis)
   rbis<-fitis_PL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster, TwostagePoissonSampleis$x,  pars=truevalue)
   rcis<-fitis_WPL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster,TwostagePoissonSampleis$x,
                 TwostagePoissonSampleis$ID_unit, TwostagePoissonSampleis$PSU, fss=pi1is,  n2infor=n2is , N2,  
                 pars=truevalue)
   
   #NML (uniformative sampling)
   Fit_NML[i,1:2]<-fixef(ra)-truevalue[1:2]
   Fit_NML[i,3]<-sigma(ra)^2-truevalue[3]
   Fit_NML[i,4]<-as.numeric(VarCorr(ra))-truevalue[4]
   
   #PL (uninformative sampling)
   Fit_PL[i,1:4]<-rb$par[1:4]-truevalue[1:4]
    
   
   #WPL (uniformative sampling)
   Fit_WPL[i,1:4]<-rc$par[1:4]-truevalue[1:4]
   
   #NML (informative sampling)
   Fitis_NML[i,1:2]<-fixef(rais)-truevalue[1:2]
   Fitis_NML[i,3]<-sigma(rais)^2-truevalue[3]
   Fitis_NML[i,4]<-as.numeric(VarCorr(rais))-truevalue[4]
   
   #PL (informative sampling)
   Fitis_PL[i,1:4]<-rbis$par[1:4]-truevalue[1:4]
   
   #WPLE (informative sampling)
   Fitis_WPL[i,1:4]<-rcis$par[1:4]-truevalue[1:4]
   
   
   #Calculate Hessian matrix H for PL (bread for uninformative sampling design)
   H_PL[,,i]=-jacobian(function(theta){with(TwostagePoissonSample, pairscore_PL(y,cluster,x,theta)) }, x=rb[[1]],method="complex")
      
      #Calculate  variance matrix J  for PL (meat for uniformative sampling design)
      J_PL[, , i]=fast_J_PL(y=TwostagePoissonSample$y, g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
                            pos=TwostagePoissonSample$ID_unit, sc=TwostagePoissonSample$PSU,fss=pi1,
                            n2infor=FirststagePoisson*n2, N2, theta=rb[[1]] )
      
      #sanwich estimator (uninformative sampling )
      G_PL[, ,i] = tryCatch(solve(H_PL[,,i])%*% J_PL[, , i]%*% t(solve(H_PL[,,i])), error=function(e) matrix(NaN, 4,4)) 
      
      #Pairwise score function PL (uninformative sampling)
      PS_PL[i, ]<- pairscore_PL(y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
                                     theta=rb[[1]])
      
      #Pairwise score function PL (uninformative sampling)
      PS_PL_true[i, ]<- pairscore_PL(y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
                                theta=truevalue)
      
      #Pairwise score function PL (uninformative sampling)
      PS_PL_cen[i, ]<- pairscore_PL(y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
                                     theta=cen_PL[[1]])
      
      
      #Calculate Hessian matrix H for PL (bread for informative sampling design)
       His_PL[,,i]=-jacobian(function(theta){with(TwostagePoissonSampleis, pairscore_PL(y,cluster,x,theta))},  x=rbis[[1]], method="complex")
         
         #Calculate  variance matrix J  for PL (meat for informative sampling design)
         Jis_PL[, , i]=fast_J_PL(y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster,x=TwostagePoissonSampleis$x,pos=TwostagePoissonSampleis$ID_unit,  
                                 sc=TwostagePoissonSampleis$PSU, fss=pi1is,    n2infor=n2is, N2, theta=rbis[[1]] )
         
         #sanwich estimator (informative sampling )
         Gis_PL[, ,i] =  tryCatch(solve(His_PL[,,i])%*% Jis_PL[, , i]%*% t(solve(His_PL[,,i])),  error=function(e) matrix(NaN, 4,4)) 
         
         #Pairwise score function PL (informative sampling)
         PSis_PL[i, ]<- pairscore_PL(y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster, x=TwostagePoissonSampleis$x,
                                          theta=rbis[[1]])
         
         #Pairwise score function PL (informative sampling)
         PSis_PL_true[i, ]<- pairscore_PL(y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster, x=TwostagePoissonSampleis$x,
                                     theta=truevalue)
         #Pairwise score function PL (informative sampling)
         PSis_PL_cen[i, ]<- pairscore_PL(y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster, x=TwostagePoissonSampleis$x,
                                          theta=cen_PL[[1]])
         
         #Calculate Hessian matrix H for WPL (bread for uninformative sampling design)
   
         H_WPL[,,i]=-jacobian(function(theta){with(TwostagePoissonSample,
                                                      pairscore_WPL(y,cluster,x,theta, ID_unit, PSU, pi1, FirststagePoisson*n2, N2 ))}, x=rc[[1]],method="complex")
            

            
            #Calculate  variance matrix J  for WPL (meat for uniformative sampling design)
            J_WPL[, , i]=fast_J_WPL(y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster,
                                    x=TwostagePoissonSample$x, pos=TwostagePoissonSample$ID_unit,  sc=TwostagePoissonSample$PSU, fss=pi1, 
                                    n2infor= FirststagePoisson*n2, N2, theta=rc[[1]])
            
            #sanwich estimator (uninformative sampling )
            G_WPL[, ,i] =tryCatch(solve(H_WPL[,,i])%*% J_WPL[, , i]%*% t(solve(H_WPL[,,i])), error=function(e) matrix(NaN, 4,4)) 
            
            #Pairwise score function WPL (uninformative sampling)
            PS_WPL[i, ]<- pairscore_WPL(TwostagePoissonSample$y, TwostagePoissonSample$cluster, TwostagePoissonSample$x,
                                             theta=rc[[1]],   TwostagePoissonSample$ID_unit, TwostagePoissonSample$PSU, fss=pi1, n2infor=FirststagePoisson*n2, N2)
            
            #Pairwise score function WPL (uninformative sampling)
            PS_WPL_true[i, ]<- pairscore_WPL(TwostagePoissonSample$y, TwostagePoissonSample$cluster, TwostagePoissonSample$x,
                                        theta=truevalue,   TwostagePoissonSample$ID_unit, TwostagePoissonSample$PSU, fss=pi1, n2infor=FirststagePoisson*n2, N2)
            
            #Pairwise score function WPL (uninformative sampling)
            PS_WPL_cen[i, ]<- pairscore_WPL(TwostagePoissonSample$y, TwostagePoissonSample$cluster, TwostagePoissonSample$x,
                                             theta=cen_PL[[1]],   TwostagePoissonSample$ID_unit, TwostagePoissonSample$PSU, fss=pi1, n2infor=FirststagePoisson*n2, N2)
            
            
            #Calculate Hessian matrix H  for WPL (bread for informative sampling design)
            ##informative sampling
   
            
            His_WPL[,,i]=-jacobian(function(theta){with(TwostagePoissonSampleis,
                                                           pairscore_WPL(y,cluster,x,theta, ID_unit, PSU, pi1is, n2is, N2))}, x=rcis[[1]],method="complex")
               
               #Calculate Variance matrix J  for WPL (meat for  informative sampling design)
            Jis_WPL[, , i]=fast_J_WPL(y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster,
                                         x=TwostagePoissonSampleis$x, pos=TwostagePoissonSampleis$ID_unit,  sc=TwostagePoissonSampleis$PSU, fss=pi1is, 
                                         n2infor= n2is, N2,  theta=rcis[[1]] )
               
               #sanwich estimator for WPL (informative sampling )
               Gis_WPL[,,i]= tryCatch(solve(His_WPL[, , i])%*% Jis_WPL[, , i]%*% t(solve(His_WPL[, , i])), error=function(e) matrix(NaN, 4,4))
               
               
               #Pairwise score function WPL (informative sampling)
               PSis_WPL[i, ]<- pairscore_WPL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster,TwostagePoissonSampleis$x,
                                                  theta=rcis[[1]],TwostagePoissonSampleis$ID_unit, TwostagePoissonSampleis$PSU, fss=pi1is,  n2infor=n2is , N2)
               
               #Pairwise score function WPL (informative sampling)
               PSis_WPL_true[i, ]<- pairscore_WPL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster,TwostagePoissonSampleis$x,
                                             theta=truevalue,TwostagePoissonSampleis$ID_unit, TwostagePoissonSampleis$PSU, fss=pi1is,  n2infor=n2is , N2)
               #Pairwise score function WPL (informative sampling)
               PSis_WPL_cen[i, ]<- pairscore_WPL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster,TwostagePoissonSampleis$x,
                                                  theta=cen_PL[[1]],TwostagePoissonSampleis$ID_unit, TwostagePoissonSampleis$PSU, fss=pi1is,  n2infor=n2is , N2)
               
               sqrt_diagJ_PL[i,]= sqrt(diag(J_PL[,,i]))
               sqrt_diagJis_PL[i,]= sqrt(diag(Jis_PL[,,i]))
               sqrt_diagJ_WPL[i,]= sqrt(diag(J_WPL[,,i]))
               sqrt_diagJis_WPL[i,]= sqrt(diag(Jis_WPL[,,i]))
               
               sqrt_diagG_PL[i,]= sqrt(diag(G_PL[,,i]))
               sqrt_diagGis_PL[i,]= sqrt(diag(Gis_PL[,,i]))
               sqrt_diagG_WPL[i,]= sqrt(diag(G_WPL[,,i]))
               sqrt_diagGis_WPL[i,]= sqrt(diag(Gis_WPL[,,i]))
               
               
}


#boxplot for uninformative sampling (NML, PL and WPL)
color<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")

#boxplot for uninformative sampling (NML, PL and WPL)
color=c( rep(color, 3))
name=c("alpha_NML", "beta_NML", "sigma^2_NML", "tau^2_NML", "alpha_PL", "beta_PL", "sigma^2_PL", "tau^2_PL", "alpha_WPL", "beta_WPL", "sigma^2_WPL", "tau^2_WPL" )
boxplot(cbind(Fit_NML[,c(1:4)],Fit_PL[,c(1:4)], Fit_WPL[,c(1:4)]) ,   col=color)
abline(h=0)


#boxplot for informative sampling (NML,PL and WPL)
boxplot(cbind(Fitis_NML[,c(1:4)],Fitis_PL[,c(1:4)], Fitis_WPL[,c(1:4)]) ,   col=color)
abline(h=0)

construct_header <- function(df, grp_names, span, align = "c", draw_line = T) {
   if (length(align) == 1) align <- rep(align, length(grp_names))
   if (!all.equal(length(grp_names), length(span), length(align)))
      stop("grp_names and span have to have the same length!")
   
   if (ncol(df) < sum(span)) stop("Span has to be less or equal to the number of columns of df")
   
   header <- mapply(function(s, a, grp) sprintf("\\multicolumn{%i}{%s}{%s}", s, a, grp),
                    span, align, grp_names)
   header <- paste(header, collapse = " & ")
   header <- paste0(header, " \\\\")
   
   if (draw_line) {
      # where do we span the lines:
      min_vals <- c(1, 1 + cumsum(span)[1:(length(span) - 1)])
      max_vals <- cumsum(span)
      line <- ifelse(grp_names == "", "",
                     sprintf("\\cmidrule(lr){%i-%i}", min_vals, max_vals))
      line <- paste(line[line != ""], collapse = " ")
      
      header <- paste0(header, "  ", line, "\n  ")
   }
   
   addtorow <- list(pos = list(-1, -1, nrow(df)),
                    command = c("\\hline\n  ", header, "\\hline\n  "))
   return(addtorow)
}


#bias,  sd (empirical standard deviation)  and estimated sd (G) for uninformative sampling (NML, PL, WPL)
df<- matrix(c(apply(Fit_NML, 2,  median), apply(Fit_NML, 2, mad), apply(Fit_PL, 2, median), apply(Fit_PL, 2, mad),
              apply(sqrt_diagG_PL, 2, function(x) median(x, na.rm=TRUE)),  apply(Fit_WPL, 2, median), apply(Fit_WPL, 2,mad),
              apply(sqrt_diagG_WPL, 2, function(x) median(x,na.rm=TRUE))),ncol=8)
df<-round(df, 2)
df<-cbind(c("beta_0", "beta_1", "sigma^2", "tau^2"), df)
colnames(df)<-c("",c("media bias", "mad"), rep(c("median bias", "mad","esd"), 2))
df
df_header <- construct_header(
   # the data.frame or matrix that should be plotted
   df,
   # the labels of the groups that we want to insert
   grp_names = c("uninformtive", "NML", "PL", "WPL"),
   # the number of columns each group spans
   span = c(1, 2, 3, 3),
   # the alignment of each group, can be a single character (lcr) or a vector
   align = "c"
)
print(xtable(df), add.to.row = df_header, include.rownames = F,  hline.after = F, latex.environments=NULL,booktabs=TRUE)   


#variance estimator (sd of PS and  J) for uninformative sampling (PL, WPL)
vardf<- matrix(c( apply(PS_PL_true, 2, median),apply(PS_PL_true, 2, mad),  apply(sqrt_diagJ_PL, 2, function(x) median(x,na.rm=TRUE)) ,
                  apply(PS_WPL_true, 2, median),apply(PS_WPL_true, 2, mad),  apply(sqrt_diagJ_WPL, 2, function(x) median(x,na.rm=TRUE))),ncol=6)
vardf<-round(vardf,2)
vardf<-cbind(c("beta_0", "beta_1", "sigma^2", "tau^2"), vardf)
colnames(vardf)<-c("parameter", c("median", "mad", "esd","median", "mad", "esd"))
vardf
vardf_header <- construct_header(
   # the data.frame or matrix that should be plotted
   vardf,
   # the labels of the groups that we want to insert
   grp_names = c("",  "Pairwise score", "Weighted pairwise score"),
   # the number of columns each group spans
   span = c(1, 3, 3),
   # the alignment of each group, can be a single character (lcr) or a vector
   align = "c"
)
print(xtable(vardf), add.to.row = vardf_header, include.rownames = F, hline.after = F, latex.environments=NULL,booktabs=TRUE)   


#bias and sd for informative sampling (NML, PL, WPL)
dfis<- matrix(c(apply(Fitis_NML, 2,  median), apply(Fitis_NML, 2, mad), apply(Fitis_PL, 2, median), apply(Fitis_PL, 2, mad), apply(sqrt_diagGis_PL, 2, function(x) median(x, na.rm=TRUE)),  apply(Fitis_WPL, 2, median), apply(Fitis_WPL, 2, mad),
                apply(sqrt_diagGis_WPL, 2, function(x) median(x, na.rm=TRUE))),ncol=8)

dfis=round(dfis, 2)
dfis<-cbind(c("beta_0", "beta_1", "sigma^2", "tau^2"), dfis)
colnames(dfis)<-c("parameter", c("median bias", "mad"), c("median bias", "mad", "esd",
"median bias", "mad", "esd"))
dfis
dfis_header <- construct_header(
   # the data.frame or matrix that should be plotted
   dfis
   # the labels of the groups that we want to insert
   grp_names = c("", "NML", "PL", "WPL"),
   # the number of columns each group spans
   span = c(1, 2, 3, 3),
   # the alignment of each group, can be a single character (lcr) or a vector
   align = "c"
)
print(xtable(dfis), add.to.row = dfis_header, include.rownames = F,  hline.after = F, latex.environments=NULL,booktabs=TRUE) 

#variance estimator for informative sampling (PL, WPL)
vardfis<-matrix(c( apply(PSis_PL_true, 2, median),apply(PSis_PL_true, 2, mad),  apply(sqrt_diagJis_PL, 2, function(x) median(x, na.rm=TRUE)) , apply(PSis_WPL_true, 2, median),apply(PSis_WPL_true, 2, mad),  apply(sqrt_diagJis_WPL, 2, function(x) median(x, na.rm=TRUE))),ncol=6)
vardfis<-round(vardfis,2)
vardfis<-cbind(c("beta_0", "beta_1", "sigma^2", "tau^2"), vardfis)
colnames(vardfis)<-c("parameter", c( "median bias", "mad", "esd","median bias", "mad", "esd" ))
vardfis
vardfis_header <- construct_header(
   # the data.frame or matrix that should be plotted
   vardfis,
   # the labels of the groups that we want to insert
   grp_names = c("",  "Pairwise score", "Weighted pairwise score"),
   # the number of columns each group spans
   span = c(1, 3, 3),
   # the alignment of each group, can be a single character (lcr) or a vector
   align = "c"
)

print(xtable(vardfis), add.to.row = vardfis_header, include.rownames = F, hline.after = F, latex.environments=NULL,booktabs=TRUE)   




