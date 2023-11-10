shots = c("goal","goal","goal","miss","miss",
          "goal","goal","miss","miss","goal")
shots_01 = as.numeric(shots=="goal")
N=length(shots)
n_goals = sum(shots=="goal")
n_missed = sum(shots=="miss")

x1 = rbinom(1000,size=10,p=0.5)
x2 = rbinom(1000,size=10,p=0.05)

#prior
theta = seq(from=1/(N+1),
            to=N/(N+1),
            length=10)
p_theta = pmin(theta,1-theta)
p_theta = p_theta/sum(p_theta)
#likelihood
p_data_given_theta = choose(N,n_goals)*
  theta^n_goals * (1-theta)^n_missed
#posterior
p_data = p_data_given_theta*p_theta
p_theta_given_data = p_data_given_theta*p_theta/sum(p_data)
posterior_mean = sum(p_theta_given_data*theta)
posterior_mean

plot(p_theta,type="l",ylim=c(0,0.3))
lines(p_data_given_theta,col="blue")
lines(p_theta_given_data,col="red")
legend("topleft",legend=c("prior","likelihood",
                          "posterior"),col=c("black","blue","red"),pch=15)
