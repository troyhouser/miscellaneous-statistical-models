library(tidyverse)
# 1 factor
y = matrix(c(22.6, 20.5, 20.8,
             22.6, 21.2, 20.5,
             17.3, 16.2, 16.6,
             21.4, 23.7, 23.2,
             20.9, 22.2, 22.6,
             14.5, 10.5, 12.3,
             20.8, 19.1, 21.3,
             17.4, 18.6, 18.6,
             25.1, 24.8, 24.9,
             14.9, 16.3, 16.6), 
           10, 3, byrow = TRUE)

one_factor_re = function(mu,sigma2_mu,sigma2){
  d=nrow(y)
  ni=ncol(y)
  Sigmai=sigma2*diag(ni)+sigma2_mu*matrix(1,ni,ni)
  l=rep(NA,nrow(y))
  for(i in 1:d){
    l[i] = 0.5*t(y[i,]-mu) %*% chol2inv(chol(Sigmai))%*%(y[i,]-mu)
  }
  ll=-(ni*d)/2*log(2*pi) - d/2*log(det(Sigmai)) - sum(i)
  return(-ll)
}
starts = list(
  mu = mean(y),
  sigma2_mu = var(rowMeans(y)),
  sigma2    = mean(apply(y, 1, var))
)
one_factor_re(mu=starts[[1]],
              sigma2_mu = starts[[2]],
              sigma2 = starts[[3]])

library(lme4)
library(tidyverse)

d = data.frame(y) %>% 
  pivot_longer(everything(), names_to = 'x', values_to = 'value') %>% 
  arrange(x) %>% 
  group_by(x) %>% 
  mutate(group = 1:n())

fit_mer = lmer(value ~ 1 | group, data = d, REML = FALSE)

summary(fit_mer)

#2 factor
y = c(1.39,1.29,1.12,1.16,1.52,1.62,1.88,1.87,1.24,1.18,
      .95,.96,.82,.92,1.18,1.20,1.47,1.41,1.57,1.65)
d = expand.grid(sire = rep(1:5, 2), dam = 1:2)
d = data.frame(d[order(d$sire), ], y)
two_factor_re <- function(mu, eta_alpha, eta_gamma, eta) {
  sigma2_alpha = exp(eta_alpha)
  sigma2_gamma = exp(eta_gamma)
  sigma2 = exp(eta)
  n = length(y)
  Sigma = sigma2 * diag(n) + sigma2_alpha * tcrossprod(Xalpha) + 
    sigma2_gamma * tcrossprod(Xgamma)
  ll = -n / 2 * log(2 * pi) - sum(log(diag(chol(Sigma)))) -
    .5 * t(y - mu) %*% chol2inv(chol(Sigma)) %*% (y - mu)
  return(-ll)
}
starts = list(
  mu = mean(y),
  eta_alpha = var(tapply(y, d$sire, mean)),
  eta_gamma = var(y) / 3,
  eta = var(y) / 3
)

Xalpha = diag(5) %x% rep(1, 4)

Xgamma = diag(10) %x% rep(1, 2)

two_factor_re(
  mu  = starts[[1]],
  eta_alpha = starts[[2]],
  eta_gamma = starts[[3]],
  eta = starts[[4]]
)

fit_mer = lmer(y ~ (1 | sire) + (1 | dam:sire), d, REML = FALSE)

summary(fit_mer)
