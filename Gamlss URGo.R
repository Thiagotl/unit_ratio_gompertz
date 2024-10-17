#log-likelihood of the unit ratio Gompertz 
library(gamlss)
sigma=5.2
mu=0.17
tau=.5
y=.025



URGo <- expression(
  log(sigma)+log(-log(1-tau))+((sigma*y)/(1-y))-2*log(1-y)-log(exp((sigma*mu)/(1-mu))-1)+
    (exp((sigma*y)/(1-y))-1)/(exp((sigma*mu)/(1-mu))-1) * log(1-tau)
)
m1URGo<-D(URGo,"mu")
s1URGo<-D(URGo,"sigma")
ms2URGo<-D(m1URGo,"sigma")
URGo<-function (mu.link = "logit", sigma.link = "log")
{
  mstats <- checklink("mu.link", "URGo", substitute(mu.link),
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "URGo", substitute(sigma.link),
                      c("inverse", "log", "identity", "own"))
  structure(list(family = c("URGo", "Unit-Ratio-Gompertz"),
                 parameters = list(mu = TRUE, sigma = TRUE),
                 nopar = 2,
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 dldm = function(y, mu, sigma) {
                   tau=.5
                   dldm <- eval(m1URGo)
                   dldm
                 },
                 d2ldm2 = function(y,mu, sigma) {
                   tau=.5
                   dldm <- eval(m1URGo)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)
                   d2ldm2
                 },
                 dldd = function(y, mu, sigma) {
                   tau=.5
                   dldd <- eval(s1URGo)
                   dldd
                 },
                 d2ldd2 = function(y,mu, sigma) {
                   tau=.5
                   dldd <- eval(s1URGo)
                   d2ldd2 = -dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu, sigma) {
                   tau=.5
                   dldm <- eval(m1URGo)
                   dldd <- eval(s1URGo)
                   d2ldmdd = -(dldm * dldd)
                   d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
                   d2ldmdd
                 },
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * log(dURGo(y=y, mu=mu, sigma=sigma)),
                 rqres = expression(
                   rqres(pfun = "pURGo", type = "Continuous", y = y, mu = mu, sigma = sigma)
                 ),
                 mu.initial = expression(mu <- rep(mean(y),length(y))),
                 sigma.initial = expression(sigma<- rep(4, length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1),
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0 & y < 1)
  ),
  class = c("gamlss.family", "family"))
}
#------------------------------------------------------------------------------------------

mu = 0.7
sigma =5 
tau = 0.5
# y=2
# density function
dURGo<-function(y, mu = 0.7, sigma =5 , tau = 0.5, log = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(y <= 0) | any(y >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  fy1 <- sigma*log((1-tau)^(-1))*exp(sigma*y/(1-y))/
    ((1-y)^(2)*(exp(sigma*mu/(1-mu))-1))*
    (1-tau)^((exp(sigma*y/(1-y))-1)/(exp(sigma*mu/(1-mu))-1))
  
  
  if(log==FALSE) 
  fy<-fy1 else fy<-log(fy1)
  fy
}
integrate(dURGo,0,0.99) # checking the pdf
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
# pURGo(.5)
# integrate(dURGo,0,.5) # checking the cdf with the pdf
#------------------------------------------------------------------------------------------
# quantile function
qURGo<-function(u,mu = 0.7, sigma = 5, tau = 0.5)
{
        # log(log(1-u)/log(1-tau)*(exp(sigma*mu/(1-mu))-1)+1)
    q<- log(log(1-u)/log(1-tau)*(exp(sigma*mu/(1-mu))-1)+1)/
    (sigma + log(log(1-u)/log(1- tau)*(exp(sigma * mu/(1 - mu)) -1)+1))
  q
}
u=pURGo(.183665,mu=.7,sigma=5)
qURGo(u,mu=.7,sigma=5) # checking the qf with the cdf
#------------------------------------------------------------------------------------------
# inversion method for randon generation
rURGo<-function(n,mu,sigma)
{
  u<- runif(n)
  y<- qURGo(u,mu =mu, sigma =sigma)
  y
}



# h<-qURGo(runif(100),mu =mu, sigma =sigma)
# 
# print(h)









# Checking the results
library(gamlss)

set.seed(10)
n<-100
# Case 1: without regressors
mu_true<-.7
sigma_true<-.15
mu_result<-sigma_result<-c()
for (i in 1:100) {
  y<-rURGo(n,mu_true,sigma_true)
  fit1<-gamlss(y~1, family="URGo", trace = F)
  logit_link<-make.link("logit")
  mu_result[i]<-logit_link$linkinv(fit1$mu.coefficients)
  sigma_result[i]<-exp(fit1$sigma.coefficients)
}
result1<- matrix(c(mu_true, mean(mu_result),
                   sigma_true, mean(sigma_result)),2,2)
colnames(result1)<-c("mu","sigma")
rownames(result1)<-c("true value","mean")
print(round(result1,2))
