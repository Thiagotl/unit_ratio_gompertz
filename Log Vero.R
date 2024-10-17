## URGo
#initial values
sigma=5.2
mu=0.17
tau=.5
x=.025

# density function
dURGo<-function(x, mu, sigma, tau, log = FALSE)
{
  dURGo<- sigma*log((1-tau)^(-1))*exp(sigma*x/(1-x))/
    ((1-x)^(2)*(exp(sigma*mu/(1-mu))-1))*
    (1-tau)^((exp(sigma*x/(1-x))-1)/(exp(sigma*mu/(1-mu))-1))
}

#log-liklihood
ll<-function(x, mu, sigma, tau)
{
  ll<-log(sigma)+log(-log(1-tau))+((sigma*x)/(1-x))-2*log(1-x)-log(exp((sigma*mu)/(1-mu))-1)+
    (exp((sigma*x)/(1-x))-1)/(exp((sigma*mu)/(1-mu))-1) * log(1-tau)
  
}

log(dURGo(x,mu,sigma, tau))
ll(x, mu, sigma, tau)

#checking for more values
z=c(.34,.41,.59)

log(dURGo(z[1],mu, sigma, tau)*dURGo(z[2],mu, sigma, tau)*dURGo(z[3],mu, sigma, tau))

sum(ll(z,mu, sigma, tau))

