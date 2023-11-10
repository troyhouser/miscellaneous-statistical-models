library(tidyverse)
N=10000
k=2
X = matrix(rnorm(N*k),ncol=k)
lp = -0.5+0.2*X[,1]+0.5*X[,2]
y=rbinom(N,size=1,prob=plogis(lp))
df=data.frame(X,y)
mle = function(par,X,y){
  beta=par
  N=nrow(X)
  LP=X%*%beta
  mu=plogis(LP)
  L=dbinom(y,1,prob=mu,log=T)
  return(-sum(L,rm.na=T))
}
logregexp = function(par,X,y){
  beta=par
  LP = X%*%beta
  L=sum(exp(-ifelse(y,1,-1)*0.5*LP))
}
X = cbind(1, X)

init = rep(0, ncol(X))
names(init) = c('intercept', 'b1', 'b2')

fit_ml = optim(
  par = init,
  fn  = mle,
  X   = X,
  y   = y,
  control = list(reltol = 1e-8)
)

fit_exp = optim(
  par = init,
  fn  = logregexp,
  X   = X,
  y   = y, 
  control = list(reltol = 1e-15)
)

fit_ml$par
fit_exp$par
fitGLM=glm(y~.,df,family=binomial)
fitGLM
