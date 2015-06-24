################################################################################################################################################################
################################################################################################################################################################
#########       A Retrospective Likelihood Approach for Efficient Integration of Multiple Omics Factors in Case-Control Association Studies  ##################
#########                                                     Brunilda-Sophia Balliu                                                          ##################
#########                                                     Leiden November 2013                                                            ##################
################################################################################################################################################################
################################################################################################################################################################

# libraries -----------------------------------------------------------------
library(mvtnorm)  # Multivariate normal density
library(Matrix)   # Nearest positive definite matrix 
library(boot)     # logit 
library(ltm)      # GHQ
 
# Data Simulation Function ------------------------------------------------
SimulateDataSet=function(N=NULL, DesiredNrOfCases=NULL,DesiredNrOfControls=NULL, a, b, v1, v2, v12, pG=pG, pA=pA, Ascertainment, Covariates){
  #N: number of individuals, 
  #a: parameters in the logistic regression
  #b :fixed effect for effect of covariates on joint E and M distribution
  #v1, v2, v12 : variance and covariance of errors on joint E and M distribution
  #p: allele frequency 
  #pA: P(Age>=50)
    if(!Ascertainment) {  
      print('Simulating population data')
      if(any(Covariates=='G')) G=rbinom(N,2,pG)  else G=NULL        
      if(any(Covariates=='S')) S=rbinom(N,1,.5)  else S=NULL
      if(any(Covariates=='A')) A=rbinom(N,1,pA)  else A=NULL
      Sigma=matrix(c(v1,v12,v12,v2),ncol=2,byrow=T) 
      e=rmvnorm(N,mean=c(0,0),sigma=Sigma)
      p=length(Covariates)
      E=rowSums(cbind(t(b[1:p]*rbind(G,S,A)), e[,1]))
      M=rowSums(cbind(t(b[(p+1):(2*p)]*rbind(G,S,A)) , e[,2]));      
      Y=rbinom(N,1,plogis(rowSums(t(a*rbind(1,E,M,G,S,A)))))  
      Data=cbind(Y,E,M,G,S,A)
      Data=data.frame(Data)
      names(Data)=c('Y','E', 'M', Covariates)
    } else{
      print('Simulating case-control data')
      ObsNrCases=0
      while(ObsNrCases<DesiredNrOfCases){
        if(any(Covariates=='G')) G=rbinom(N,2,pG)  else G=NULL        
        if(any(Covariates=='S')) S=rbinom(N,1,.5)  else S=NULL
        if(any(Covariates=='A')) A=rbinom(N,1,pA)  else A=NULL
        Sigma=matrix(c(v1,v12,v12,v2),ncol=2,byrow=T) 
        e=rmvnorm(N,mean=c(0,0),sigma=Sigma)
        p=length(Covariates)
        E=rowSums(cbind(t(b[1:p]*rbind(G,S,A)), e[,1]))
        M=rowSums(cbind(t(b[(p+1):(2*p)]*rbind(G,S,A)) , e[,2]));      
        Y=rbinom(N,1,plogis(rowSums(t(a*rbind(1,E,M,G,S,A)))))  
        Data=cbind(Y,E,M,G,S,A)
        Data=data.frame(Data)
        names(Data)=c('Y','E', 'M', Covariates)
        ObsNrCases=sum(Data[,'Y'])
      }
      DataCases=Data[Data$Y==1,]
      DataControls=Data[Data$Y==0,]
      CaseID=sample(1:nrow(DataCases),DesiredNrOfCases) 
      ControlsID=sample(1:nrow(DataControls),DesiredNrOfControls) 
      Data=rbind(DataCases[CaseID,], DataControls[ControlsID,])        
    } 
  return(Data)
} 

# Data Simulation Function ------------------------------------------------
SimulateDataSetTdist=function(N=NULL, DesiredNrOfCases=NULL,DesiredNrOfControls=NULL, a, b, v1, v2, v12, pG=pG, pA=pA, Ascertainment, Covariates,DF=100){
  #N: number of individuals, 
  #a: parameters in the logistic regression
  #b :fixed effect for effect of covariates on joint E and M distribution
  #v1, v2, v12 : variance and covariance of errors on joint E and M distribution
  #p: allele frequency 
  #pA: P(Age>=50)
    if(!Ascertainment) {  
      print('Simulating population data')
      if(any(Covariates=='G')) G=rbinom(N,2,pG)  else G=NULL        
      if(any(Covariates=='S')) S=rbinom(N,1,.5)  else S=NULL
      if(any(Covariates=='A')) A=rbinom(N,1,pA)  else A=NULL
      Sigma=matrix(c(v1,v12,v12,v2),ncol=2,byrow=T) 
      e=rmvt(N,df = DF, delta=c(0,0),sigma=Sigma)
      p=length(Covariates)
      E=rowSums(cbind(t(b[1:p]*rbind(G,S,A)), e[,1]))
      M=rowSums(cbind(t(b[(p+1):(2*p)]*rbind(G,S,A)) , e[,2]));      
      Y=rbinom(N,1,plogis(rowSums(t(a*rbind(1,E,M,G,S,A)))))  
      Data=cbind(Y,E,M,G,S,A)
      Data=data.frame(Data)
      names(Data)=c('Y','E', 'M', Covariates)
    } else{
      print('Simulating case-control data')
      ObsNrCases=0
      while(ObsNrCases<DesiredNrOfCases){
        if(any(Covariates=='G')) G=rbinom(N,2,pG)  else G=NULL        
        if(any(Covariates=='S')) S=rbinom(N,1,.5)  else S=NULL
        if(any(Covariates=='A')) A=rbinom(N,1,pA)  else A=NULL
        Sigma=matrix(c(v1,v12,v12,v2),ncol=2,byrow=T) 
        e=rmvt(N,df = DF, delta=c(0,0),sigma=Sigma)
        p=length(Covariates)
        E=rowSums(cbind(t(b[1:p]*rbind(G,S,A)), e[,1]))
        M=rowSums(cbind(t(b[(p+1):(2*p)]*rbind(G,S,A)) , e[,2]));      
        Y=rbinom(N,1,plogis(rowSums(t(a*rbind(1,E,M,G,S,A)))))  
        Data=cbind(Y,E,M,G,S,A)
        Data=data.frame(Data)
        names(Data)=c('Y','E', 'M', Covariates)
        ObsNrCases=sum(Data[,'Y'])
      }
      DataCases=Data[Data$Y==1,]
      DataControls=Data[Data$Y==0,]
      CaseID=sample(1:nrow(DataCases),DesiredNrOfCases) 
      ControlsID=sample(1:nrow(DataControls),DesiredNrOfControls) 
      Data=rbind(DataCases[CaseID,], DataControls[ControlsID,])        
    } 
  return(Data)
} 

# Data Simulation Function ------------------------------------------------
SimulateDataSetSN=function(N=NULL, DesiredNrOfCases=NULL,DesiredNrOfControls=NULL, a, b, v1, v2, v12, pG=pG, pA=pA, Ascertainment, Covariates){
  #N: number of individuals, 
  #a: parameters in the logistic regression
  #b :fixed effect for effect of covariates on joint E and M distribution
  #v1, v2, v12 : variance and covariance of errors on joint E and M distribution
  #p: allele frequency 
  #pA: P(Age>=50)
    if(!Ascertainment) {  
      print('Simulating population data')
      if(any(Covariates=='G')) G=rbinom(N,2,pG)  else G=NULL        
      if(any(Covariates=='S')) S=rbinom(N,1,.5)  else S=NULL
      if(any(Covariates=='A')) A=rbinom(N,1,pA)  else A=NULL
      Sigma=matrix(c(v1,v12,v12,v2),ncol=2,byrow=T) 
      e=rmsn(n=N, xi=c(0,0), Omega=Sigma, alpha=c(0,0), tau=0, dp=NULL)
      p=length(Covariates)
      E=rowSums(cbind(t(b[1:p]*rbind(G,S,A)), e[,1]))
      M=rowSums(cbind(t(b[(p+1):(2*p)]*rbind(G,S,A)) , e[,2]));      
      Y=rbinom(N,1,plogis(rowSums(t(a*rbind(1,E,M,G,S,A)))))  
      Data=cbind(Y,E,M,G,S,A)
      Data=data.frame(Data)
      names(Data)=c('Y','E', 'M', Covariates)
    } else{
      print('Simulating case-control data')
      ObsNrCases=0
      while(ObsNrCases<DesiredNrOfCases){
        if(any(Covariates=='G')) G=rbinom(N,2,pG)  else G=NULL        
        if(any(Covariates=='S')) S=rbinom(N,1,.5)  else S=NULL
        if(any(Covariates=='A')) A=rbinom(N,1,pA)  else A=NULL
        Sigma=matrix(c(v1,v12,v12,v2),ncol=2,byrow=T) 
        e=rmsn(n=N, xi=c(0,0), Omega=Sigma, alpha=c(2,2), tau=0, dp=NULL)
        e=e-apply(e,2,mean)
        p=length(Covariates)
        E=rowSums(cbind(t(b[1:p]*rbind(G,S,A)), e[,1]))
        M=rowSums(cbind(t(b[(p+1):(2*p)]*rbind(G,S,A)) , e[,2]));      
        Y=rbinom(N,1,plogis(rowSums(t(a*rbind(1,E,M,G,S,A)))))  
        Data=cbind(Y,E,M,G,S,A)
        Data=data.frame(Data)
        names(Data)=c('Y','E', 'M', Covariates)
        ObsNrCases=sum(Data[,'Y'])
      }
      DataCases=Data[Data$Y==1,]
      DataControls=Data[Data$Y==0,]
      CaseID=sample(1:nrow(DataCases),DesiredNrOfCases) 
      ControlsID=sample(1:nrow(DataControls),DesiredNrOfControls) 
      Data=rbind(DataCases[CaseID,], DataControls[ControlsID,])        
    } 
  return(Data)
} 


# The Penetrance Function ------------------------------------------
Penetrance=function(a, Data, loglik=FALSE, Covariates) {
  #a: parameters penetrance function
  #Data: Data Data.frame as generated from 'SimulateDataSet'
  #loglik: TRUE / FALSE should the loglikelihood be computed?
  #Covariates=FALSE : TRUE / FALSE : are there any covariates, besides E, M and G included?
  p=length(Covariates)
  if(length(a)!=p+3) print('Warning, the number of parameters in P(Y|X) is incorrect!')
  prob=plogis(apply(t(t(cbind(1,Data[c('E','M', Covariates)])) *a),1,sum)) 
  probYigivenXi=ifelse(Data$Y,prob,1-prob)
  if(loglik) return(-sum(log(probYigivenXi))) else return(probYigivenXi)
}

# Joint Distribution of E and M -------------------------------------------
JointEandM=function(b, Sigma, Data, loglik=FALSE, Covariates) {
  #b: parameters for the multivariate normal
  #loglik: TRUE / FALSE should the loglikelihood be computed?
  p=length(Covariates)
  if(length(b)!=2*p) print('Warning, the number of parameters in P(E,M|X) is incorrect!')
  Y=Data[c('E','M')]
  mu=cbind(rowSums(t(b[1:p]*t(Data[,Covariates]))),rowSums(t(b[(p+1):(2*p)]*t(Data[,Covariates]))))
  JointProbs=sapply(1:nrow(Data), function(i) dmvnorm(Y[i,], mean=mu[i,], sigma=Sigma)) 
  if(loglik) return(-sum(log(JointProbs))) else return(JointProbs)
}
JointEandMFixedSigma=function(b, Data, loglik=FALSE, Covariates) {
  #b: parameters for the multivariate normal
  #loglik: TRUE / FALSE should the loglikelihood be computed?
  p=length(Covariates)
  if(length(b)!=2*p) print('Warning, the number of parameters in P(E,M|X) is incorrect!')
  Y=Data[c('E','M')]
  mu=cbind(rowSums(t(b[1:p]*t(Data[,Covariates]))),rowSums(t(b[(p+1):(2*p)]*t(Data[,Covariates]))))
  Sigma=matrix(as.vector(nearPD(var(Y - mu))$mat),ncol=2)
  JointProbs=sapply(1:nrow(Data), function(i) dmvnorm(Y[i,], mean=mu[i,], sigma=Sigma)) 
  if(loglik) return(-sum(log(JointProbs))) else return(JointProbs)
}


# Genotype Distribution ---------------------------------------------------
GenDist<-function(pG, Data, loglik=FALSE){ 
  #p: allele frequency, 
  #loglik: TRUE / FALSE should the loglikelihood be computed?
  if(is.null(pG)) { 
    lik=rep(1,nrow(Data))
  } else {
    I1<-Data$G==2 ; I2<-Data$G==1 ; I3<-Data$G==0
    q1<-pG^2 ; q2<-2*pG*(1-pG) ; q3<-(1-pG)^2
    lik<-(q1^I1)*(q2^I2)*(q3^I3)
    if(loglik) return(lik=-sum(log(lik)))
  }
  return(lik)
}

# Sex Distribution  ---------------------------------------------------
SexDist<-function(pS,Data, loglik=FALSE){ 
  #loglik: TRUE / FALSE should the loglikelihood be computed?
  if(is.null(pS)) {
    lik=rep(1,nrow(Data))   
  } else {
    lik=ifelse(Data$S,pS,1-pS)
    if(loglik) lik=-sum(log(lik)) 
  }
  return(lik)
}

# Age Distribution ---------------------------------------------------
AgeDist<-function(pA,Data, loglik=FALSE){ 
  #loglik: TRUE / FALSE should the loglikelihood be computed?
  if(is.null(pA)) {
    lik=rep(1,nrow(Data))   
  } else {
  lik=ifelse(Data$A,pA,1-pA)
  if(loglik) lik=-sum(log(lik)) 
  }
  return(lik)
}

# P(Y) via  GHQ -----------------------------------------------
GHQ=function(a, b, Sigma, pG, pS, pA, Data, Covariates, npoints=100) {
  p=length(Covariates)
  H=solve(Sigma)
  B=chol(H)
  GH <- ltm:::gauher(npoints)
  w<-as.matrix(expand.grid(GH$w,GH$w)[,c('Var1','Var2')])
  x=as.matrix(expand.grid(GH$x,GH$x)[,c('Var1','Var2')])
  exp_x=exp(x[,1]^2+x[,2]^2)
  if(any(Covariates=='G')) G=0:2 else G=NULL
  if(any(Covariates=='S')) S=0:1 else S=NULL
  if(any(Covariates=='A')) A=0:1 else A=NULL
  X=list(G=G,S=S,A=A)  
  GSA=expand.grid(X[Covariates])
  mu=cbind(rowSums(t(b[1:p]*t(GSA))) , rowSums(t(b[(p+1):(2*p)]*t(GSA))))
  Return=sum(sapply(1:nrow(mu), function(i) {
    theta=mu[i,]
    x_star=t(sqrt(2)*solve(B)%*%t(x) + theta)
    Data.GSA=data.frame(cbind(Y=1,E=x_star[,1],M=x_star[,2],G=GSA[i,'G'],S=GSA[i,'S'],A=GSA[i,'A']))
    names(Data.GSA)=c('Y', 'E', 'M', Covariates)
    sum(2*(det(B)^(-1))*apply(w,1,prod)*exp_x * 
          Penetrance(a=a,Data=Data.GSA, Covariates=Covariates)*
          dmvnorm(x=x_star,mean=theta,sigma=Sigma)* 
          GenDist(pG=pG,Data=Data.GSA) * 
          SexDist(pS=pS,Data=Data.GSA)*
          AgeDist(pA=pA,Data=Data.GSA))
  }))
  return(Return)
}

GHQWithoutSigma=function(a, b, pG, pS, pA, Data, Covariates, npoints=100) {
  Y=Data[,c("E","M")]
  p=length(Covariates)
  mu=cbind(rowSums(t(b[1:p]*t(Data[,Covariates]))),rowSums(t(b[(p+1):(2*p)]*t(Data[,Covariates]))))
  Sigma=matrix(as.vector(nearPD(var(Y-mu))$mat),ncol=2) #matrix(c(v1,v12,v12,v2),ncol=2)
  H=solve(Sigma)
  B=chol(H)
  GH <- ltm:::gauher(npoints)
  w<-as.matrix(expand.grid(GH$w,GH$w)[,c('Var1','Var2')])
  x=as.matrix(expand.grid(GH$x,GH$x)[,c('Var1','Var2')])
  exp_x=exp(x[,1]^2+x[,2]^2)
  if(any(Covariates=='G')) G=0:2 else G=NULL
  if(any(Covariates=='S')) S=0:1 else S=NULL
  if(any(Covariates=='A')) A=0:1 else A=NULL
  X=list(G=G,S=S,A=A)  
  GSA=expand.grid(X[Covariates])
  mu=cbind(rowSums(t(b[1:p]*t(GSA))) , rowSums(t(b[(p+1):(2*p)]*t(GSA))))
  Return=sum(sapply(1:nrow(mu), function(i) {
    theta=mu[i,]
    x_star=t(sqrt(2)*solve(B)%*%t(x) + theta)
    Data.GSA=data.frame(cbind(Y=1,E=x_star[,1],M=x_star[,2],G=GSA[i,'G'],S=GSA[i,'S'],A=GSA[i,'A']))
    names(Data.GSA)=c('Y', 'E', 'M', Covariates)
    sum(2*(det(B)^(-1))*apply(w,1,prod)*exp_x * 
          Penetrance(a=a,Data=Data.GSA, Covariates=Covariates)*
          dmvnorm(x=x_star,mean=theta,sigma=Sigma)* 
          GenDist(pG=pG,Data=Data.GSA) * 
          SexDist(pS=pS,Data=Data.GSA)*
          AgeDist(pA=pA,Data=Data.GSA))
  }))
  return(Return)
}

# P(Y) via Monte Carlo integration -----------------------------------------------
MonteCarlo=function(a, b, Sigma, pG, pS, pA, Data, Covariates, npoints=2000) {
  Y=Data[,c("E","M")]
  p=length(Covariates)
  if(any(Covariates=='G')) G=0:2 else G=NULL
  if(any(Covariates=='S')) S=0:1 else S=NULL
  if(any(Covariates=='A')) A=0:1 else A=NULL
  X=list(G=G,S=S,A=A)  
  GSA=expand.grid(X[Covariates])
  mu=cbind(rowSums(t(b[1:p]*t(GSA))) , rowSums(t(b[(p+1):(2*p)]*t(GSA))))  
  Return=sum(sapply(1:nrow(mu), function(i) {
    x_star=rmvnorm(n=npoints, mean=mu[i,],sigma=Sigma)
    Data.GSA=data.frame(cbind(Y=1,E=x_star[,1],M=x_star[,2],G=GSA[i,'G'],S=GSA[i,'S'],A=GSA[i,'A']))
    names(Data.GSA)=c('Y', 'E', 'M', Covariates)
    sum(Penetrance(a=a,Data=Data.GSA, Covariates=Covariates) * GenDist(pG=pG,Data=Data.GSA) * SexDist(pS=pS,Data=Data.GSA) * AgeDist(pA=pA,Data=Data.GSA))/npoints
  })) #IndS=any(Covariates=='S')
  return(Return)
}

MonteCarloWithoutSigma=function(a, b, pG, pS, pA, Data, Covariates, npoints=2000) {
  Y=Data[,c("E","M")]
  p=length(Covariates)
  mu=cbind(rowSums(t(b[1:p]*t(Data[,Covariates]))),rowSums(t(b[(p+1):(2*p)]*t(Data[,Covariates]))))
  Sigma=matrix(as.vector(nearPD(var(Y-mu))$mat),ncol=2) #matrix(c(v1,v12,v12,v2),ncol=2)
  if(any(Covariates=='G')) G=0:2 else G=NULL
  if(any(Covariates=='S')) S=0:1 else S=NULL
  if(any(Covariates=='A')) A=0:1 else A=NULL
  X=list(G=G,S=S,A=A)  
  GSA=expand.grid(X[Covariates])
  mu=cbind(rowSums(t(b[1:p]*t(GSA))) , rowSums(t(b[(p+1):(2*p)]*t(GSA))))  
  Return=sum(sapply(1:nrow(mu), function(i) {
    x_star=rmvnorm(n=npoints, mean=mu[i,],sigma=Sigma)
    Data.GSA=data.frame(cbind(Y=1,E=x_star[,1],M=x_star[,2],G=GSA[i,'G'],S=GSA[i,'S'],A=GSA[i,'A']))
    names(Data.GSA)=c('Y', 'E', 'M', Covariates)
    sum(Penetrance(a=a,Data=Data.GSA, Covariates=Covariates) * GenDist(pG=pG,Data=Data.GSA) * SexDist(pS=pS,Data=Data.GSA) * AgeDist(pA=pA,Data=Data.GSA))/npoints
  })) #IndS=any(Covariates=='S')
  return(Return)
}

# Minus log likelihood -----------------------------------------------------------
LogLikelihood = function(theta, names.theta, WithSigma=FALSE, Data, Model, Type, Covariates){ 
  # theta : parameters of the likelihood
  # Model = c('NULL','ALTERNATIVE')   
  # Type = c('PROSPECTIVE', 'RETROSPECTIVE', 'JOINT')   
  p=length(Covariates)
  names(theta)=names.theta
#   write.table(x=t(theta[c('VarE','CovEM','VarM')]),file = '~/Desktop/Values.txt',append = T,quote = F,sep = '&',row.names = F, col.names = F)
  if(Type=='PROSPECTIVE') { 
    if( Model=='NULL') { 
      a=rep(NA,times=3+p) ; names(a)=c('a0', 'aE', 'aM',paste('a', Covariates,sep=''))
      a[1:3]=c(theta[1],0,0)
      if(any(Covariates=='G')) a['aG']=0 ;
      if(any(Covariates=='S')) a['aS']=theta['aS'];
      if(any(Covariates=='A')) a['aA']=theta['aA']; 
    }
    if(Model=='ALTERNATIVE'){ 
      a=rep(NA,times=3+p) ; names(a)=c('a0', 'aE', 'aM',paste('a', Covariates,sep=''))
      a=theta[1:(3+p)]
    }  
    b=pG=pA=pS=NULL
  } 
  
  if(Type=='RETROSPECTIVE') {
    names.a=c('a0', 'aE', 'aM',paste('a', Covariates,sep=''))
    a=rep(NA,times=length(names.a)) ; names(a)=names.a
    names.b=c(paste(paste('b', Covariates,sep=''), 'oE' ,sep=''),paste(paste('b', Covariates,sep=''), 'oM' ,sep=''))
    b=rep(NA,times=length(names.b)) ; names(b)=names.b

    pG=pA=pS=NULL
    
    if(Model=='NULL') {  
      a[1:3]=c(theta[1],0,0) 
      if(any(Covariates=='G')) { a['aG']=0 ; b['bGoE']=theta['bGoE'] ; b['bGoM']=theta['bGoM'] ; pG=inv.logit(theta['pG'])}      
    }   
    if(Model=='ALTERNATIVE'){ 
      a=theta[1:length(a)]  
      if(any(Covariates=='G')) { b['bGoE']=theta['bGoE'] ; b['bGoM']=theta['bGoM'] ; pG=inv.logit(theta['pG'])} 
    }
    if(any(Covariates=='S')) {a['aS']=theta['aS']; b['bSoE']=theta['bSoE'] ; b['bSoM']=theta['bSoM']; pS=inv.logit(theta['pS'])}
    if(any(Covariates=='A')) {a['aA']=theta['aA']; b['bAoE']=theta['bAoE'] ; b['bAoM']=theta['bAoM']; pA=inv.logit(theta['pA'])}
  
      if(WithSigma){ 
      sigma_E=theta[which(names.theta=='VarE')]
      sigma_M=theta[which(names.theta=='VarM')]
      cov_EM=theta[which(names.theta=='CovEM')]
      Theta=c(sigma_E, cov_EM, sigma_M)
      Theta[c(1,3)]=exp(Theta[c(1,3)])
      Chol.Sigma=matrix(c(Theta[1], 0,Theta[2], Theta[3]), byrow=T, ncol=2)
      Sigma=Chol.Sigma%*%t(Chol.Sigma)
      Sigma=matrix(as.vector(nearPD(Sigma)$mat),ncol=2) 
    } else {
      Y=Data[,c("E","M")]
      mu=cbind(rowSums(t(b[1:p]*t(Data[,Covariates]))),rowSums(t(b[(p+1):(2*p)]*t(Data[,Covariates]))))
      Sigma=matrix(as.vector(nearPD(var(Y - mu))$mat),ncol=2)  
    }
  } 
  
  if(Type=='PROSPECTIVE') {logLik=sum(log(Penetrance(a=a,Data=Data,Covariates=Covariates)))} 
  
  if(Type=='RETROSPECTIVE'){
      Numerator=sum(
        log(Penetrance(a=a,Data=Data,Covariates=Covariates, loglik=FALSE)) + 
        log(JointEandM(b=b,Sigma=Sigma,Data=Data,Covariates=Covariates, loglik=FALSE)) + 
        log(GenDist(pG=pG,Data=Data, loglik=FALSE)) + 
        log(SexDist(pS=pS,Data=Data, loglik=FALSE)) + 
        log(AgeDist(pA=pA,Data=Data, loglik=FALSE))
        )          
      ProbY1=GHQ(a=a,b=b,Sigma=Sigma,pG=pG,pS=pS, pA=pA,Data=Data, Covariates=Covariates,npoints=100) 
      ProbY=sum(log(ifelse(Data$Y,ProbY1, 1-ProbY1)))
      logLik=Numerator - ProbY
  }
  print(sprintf("%.3f", c(inv.logit(a)[1],exp(a[-1]),b,c(Sigma),inv.logit(pG),logLik)))
  #,inv.logit(pS),inv.logit(pA)
  return(logLik)
}
