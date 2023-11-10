library(tidyverse)
N_g = 2
N_1 = 50
N_2 = 50
mu_1 = 1
mu_2 = -0.5
sigma_1 = 1
sigma_2 = 1
y_1 = rnorm(N_1,mu_1,sigma_1)
y_2 = rnorm(N_2,mu_2,sigma_2)
y = c(y_1,y_2)
group_id = as.numeric(gl(2,N_1))
d = data.frame(y,group_id)
stan_data = list(N=length(y),
                 N_g=N_g,
                 group_id=group_id,
                 y=y)
library(rstan)
stanmodelcode = '
data {
  int<lower = 1> N;                              // sample size 
  int<lower = 2> N_g;                            // number of groups
  vector[N] y;                                   // response
  int<lower = 1, upper = N_g> group_id[N];       // group ID
}

transformed data{
  real y_mean;                                   // mean of y; see mu prior
  
  y_mean = mean(y); 
}

parameters {
  vector[2] mu;                                  // estimated group means and sd
  vector<lower = 0>[2] sigma;                    // Kruschke puts upper bound as well; ignored here
  real<lower = 0, upper = 100> nu;               // df for t distribution
}

model {
  // priors
  // note that there is a faster implementation of this for stan, 
  // and that the sd here is more informative than in Kruschke paper
  mu    ~ normal(y_mean, 10);                       
  sigma ~ cauchy(0, 5);
  
  // Based on Kruschke; makes average nu 29 
  // might consider upper bound, as if too large then might as well switch to normal
  nu    ~ exponential(1.0/29);                
  
  // likelihood
  for (n in 1:N) {
    y[n] ~ student_t(nu, mu[group_id[n]], sigma[group_id[n]]);
    
    // compare to normal; remove all nu specifications if you do this;
    //y[n] ~ normal(mu[group_id[n]], sigma[group_id[n]]);           
  }
}

generated quantities {
  vector[N] y_rep;                               // posterior predictive distribution
  real mu_diff;                                  // mean difference
  real cohens_d;                                 // effect size; see footnote 1 in Kruschke paper
  real CLES;                                     // common language effect size
  real CLES2;                                    // a more explicit approach; the mean should roughly equal CLES

  for (n in 1:N) {
    y_rep[n] = student_t_rng(nu, mu[group_id[n]], sigma[group_id[n]]);
  }

  mu_diff  = mu[1] - mu[2];
  cohens_d = mu_diff / sqrt(sum(sigma)/2);
  CLES     = normal_cdf(mu_diff / sqrt(sum(sigma)), 0, 1);
  CLES2    = student_t_rng(nu, mu[1], sigma[1]) - student_t_rng(nu, mu[2], sigma[2]) > 0;
}


'
fit = stan(model_code=stanmodelcode, data=stan_data, iter=12000, warmup=2000, cores=4, thin=10)
print(
  fit,
  pars   = c('mu', 'sigma', 'mu_diff', 'cohens_d', 'CLES', 'CLES2', 'nu'),
  probs  = c(.025, .5, .975), 
  digits = 3
)
y_rep   = extract(fit, par = 'y_rep')$y_rep
mu_diff = extract(fit, par = 'mu_diff')$mu_diff

init = d %>% 
  group_by(group_id) %>% 
  summarise(
    mean = mean(y),
    sd = sd(y),
  )

means = init$mean
sds   = init$sd

mu_1 - mu_2     
cohens_d = extract(fit, par = 'cohens_d')$cohens_d
(mu_1 - mu_2) / sqrt((sigma_1 ^ 2 + sigma_2 ^ 2)/2)      # population
CLES = extract(fit, par='CLES')$CLES
pnorm((mu_1 - mu_2) / sqrt(sigma_1^2 + sigma_2^2))       # population
pnorm((means[1] - means[2]) / sqrt(sum(sds^2))) 
mean(CLES)                                      
t.test(y_1, y_2)

