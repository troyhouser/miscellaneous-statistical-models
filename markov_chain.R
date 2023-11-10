library(tidyverse)
library(markovchain)
A = matrix(c(0.7,0.3,
             0.3,0.7),nrow=2,byrow=T)
dtmcA = new("markovchain",
            transitionMatrix=A,
            states=c("a","b"),
            name = "MC_A")
dtmcA
plot(dtmcA)
initialState = c(0, 1)
steps = 4
finalState = initialState * dtmcA^steps 
finalState
steadyStates(dtmcA)
observed_states = sample(c('a', 'b'), 50, c(.7, .3), replace = TRUE)

createSequenceMatrix(observed_states)
markovchainFit(observed_states)

# take matrix power
mat_power = function(M,N){
  if(N==1) return(M)
  M %*% mat_power(M,N-1)
}

# sequence generator
create_sequence = function(states,len,tmat){
  states_numeric = length(unique(states))
  out = numeric(len)
  out[1] = sample(states_numeric,1,prob=colMeans(tmat))
  for(i in 2:len){
    out[i] = sample(states_numeric,1,prob=tmat[out[i-1],])
  }
  states[out]
}
test_matrix = matrix(rep(2,4),nrow=2)
mat_power(test_matrix,2)
A=matrix(c(0.7,0.3,0.3,0.7),nrow=2,byrow=T)
observed_states = create_sequence(c("a","b"),500,tmat=A)
prop.table(createSequenceMatrix(observed_states),1)
fit=markovchainFit(observed_states)
sum(createSequenceMatrix(observed_states)*log(fit$estimate@transitionMatrix))
markov_ll <- function(par, x) {
  # par should be the c(A) of transition probabilities A
  nstates = length(unique(x))
  
  # create transition matrix
  par = matrix(par, ncol = nstates)
  par = t(apply(par, 1, function(x) x / sum(x)))
  
  # create seq matrix
  seq_mat = table(x[-length(x)], x[-1])
  
  # calculate log likelihood
  ll = sum(seq_mat * log(par))
  
  -ll
}
initpar = rep(1, 4)

fit = optim(
  par = initpar,
  fn  = markov_ll,
  x   = observed_states,
  method  = 'BFGS',
  control = list(reltol = 1e-12)
)
