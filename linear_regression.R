library(tidyverse)
N=100
k=2
X=matrix(rnorm(N*k),ncol=k)
y = -0.5+0.2*X[,1]+0.1*X[,2]+rnorm(N,sd=0.5)
df=data.frame(X,y)
df
#maximum likelihood
mle = function(par,X,y){
  beta=par[-1]
  sigma2=par[1]
  sigma=sqrt(sigma2)
  N=nrow(X)
  LP=X%*%beta
  mu = LP
  L = dnorm(y,mu,sigma,log=T)
  return(-sum(L,na.rm=T))
}
lse = function(par,X,y){
  beta=par
  LP = X%*%beta
  mu = LP
  L = crossprod(y-mu)
}

X = cbind(1,X)
init = c(1,rep(0,ncol(X)))
names(init) = c("sigma2","intercept","b1","b2")
fitML = optim(par=init,
              fn=mle,
              X=X,
              y=y,
              control=list(reltol=1e-8))
fitLS = optim(par=init[-1],
              fn=lse,
              X=X,
              y=y,
              control=list(reltol=1e-8))
parsML=fitML$par
parsLS=c(sigma2=fitLS$value/(N-k-1),fitLS$par)
fitLM=lm(y~.,df)
fitGLM=glm(y~.,df,family="gaussian")
parsML
parsLS
fitLM
fitGLM
