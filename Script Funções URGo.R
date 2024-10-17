#------------------------------------------------------------------------------------------
# density function
dURGo<-function(x, mu = 0.7, sigma =5 , tau = 0.5, log = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(x <= 0) | any(x >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  fy1 <- sigma*log((1-tau)^(-1))*exp(sigma*x/(1-x))/
    ((1-x)^(2)*(exp(sigma*mu/(1-mu))-1))*
    (1-tau)^((exp(sigma*x/(1-x))-1)/(exp(sigma*mu/(1-mu))-1))
  
  
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}
# integrate(dURGo,0,1) # checking the pdf
integrate(dURGo,0,0.99)
#------------------------------------------------------------------------------------------
# cumulative distribution function
pURGo<-function(q, mu = 0.7, sigma = 5, tau = 0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  cdf1<- 1-(1-tau)^((exp(sigma*q/(1-q))-1)/(exp(sigma*mu/(1-mu))-1))
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
  
}
pURGo(.3,sigma = 5)
integrate(dUG,0,.5) # checking the cdf with the pdf
integrate(dURGo,0,.5)
#------------------------------------------------------------------------------------------
# quantile function
qURGo<-function(u,mu = 0.7, sigma = 2.1, tau = 0.5)
{
  q<- log(log(1-u)/log(1- tau)*(exp(sigma - mu)/(1 - mu) -1)+1)/
    sigma + log(log(1-u)/log(1- tau)*(exp(sigma - mu)/(1 - mu) -1)+1)
  q
}
u=pURGo(.5)
qURGo(u,mu=.7,sigma=2.1) # checking the qf with the cdf
#------------------------------------------------------------------------------------------
# inversion method for randon generation
rUW<-function(n,mu,sigma,tau=.5)
{
  u<- runif(n)
  y<- qUW(u,mu =mu, sigma =sigma)
  y
}