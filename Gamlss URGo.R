#log-likelihood of the unit ratio Gompertz 

rm(list = ls())

URGo <- expression(log(sigma)+log(-log(1-tau))+((sigma*y)/(1-y))-2*log(1-y)-log(exp((sigma*mu)/(1-mu))-1)+
                     (exp((sigma*y)/(1-y))-1)/(exp((sigma*mu)/(1-mu))-1) * log(1-tau))

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
                 sigma.initial = expression(sigma<- rep(3.14, length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1),
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0 & y < 1)
  ),
  class = c("gamlss.family", "family"))
}
#------------------------------------------------------------------------------------------

# density function
dURGo<-function(y, mu, sigma, tau = 0.5, log = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(y <= 0) | any(y >= 1)) stop(paste("y must be between 0 and 1", "\n", ""))
  fy1 <- sigma*log((1-tau)^(-1))*exp(sigma*y/(1-y))/
    ((1-y)^(2)*(exp(sigma*mu/(1-mu))-1))*
    (1-tau)^((exp(sigma*y/(1-y))-1)/(exp(sigma*mu/(1-mu))-1))
  
  
  if(log==FALSE) 
  fy<-fy1 else fy<-log(fy1)
  fy
}
#integrate(dURGo,0,0.99) # checking the pdf
#------------------------------------------------------------------------------------------
# cumulative distribution function
pURGo<-function(q, mu, sigma, tau = 0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("y must be between 0 and 1", "\n", ""))
  cdf1<- 1-(1-tau)^((exp(sigma*q/(1-q))-1)/(exp(sigma*mu/(1-mu))-1))
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
  
}
# pURGo(.5)
# integrate(dURGo,0,.5) # checking the cdf with the pdf
#------------------------------------------------------------------------------------------
# quantile function
qURGo<-function(u,mu, sigma, tau = 0.5)
{
        # log(log(1-u)/log(1-tau)*(exp(sigma*mu/(1-mu))-1)+1)
    q<- log(log(1-u)/log(1-tau)*(exp(sigma*mu/(1-mu))-1)+1)/
    (sigma + log(log(1-u)/log(1- tau)*(exp(sigma * mu/(1 - mu)) -1)+1))
  q
}
#u=pURGo(.183665,mu=.7,sigma=5)
#qURGo(u,mu=.7,sigma=5) # checking the qf with the cdf
#------------------------------------------------------------------------------------------
# inversion method for randon generation
rURGo<-function(n,mu,sigma)
{
  u<- runif(n)
  y<- qURGo(u,mu =mu, sigma =sigma)
  y
}

# h<-qURGo(runif(100),mu =mu, sigma =sigma)
# print(h)

# Checking the results
library(gamlss)

set.seed(10)
n<-c(70, 150, 300, 500, 1000)
R<-10000

# Case 1: without regressors
mu_true<-.7
sigma_true<-5
mu_result<-sigma_result<-c()
logit_link<-make.link("logit")



for (i in n){
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  for (i in 1:R){
  y<-rURGo(n,mu_true,sigma_true)
  fit1<-gamlss(y~1, family="URGo", c.crit = 0.001, n.cyc = 700,
               mu.step = .1, sigma.step = .1,trace = F)
  
  mu_result[i]<-logit_link$linkinv(fit1$mu.coefficients)
  sigma_result[i]<-exp(fit1$sigma.coefficients)
  setTxtProgressBar(pb, i)
  }
  result1<- matrix(c(mu_true, mean(mu_result),
                     sigma_true, mean(sigma_result)),2,2)
  colnames(result1)<-c("mu","sigma")
  rownames(result1)<-c("true value","mean")
  print(round(result1,2))
  }


# n<-100 #TA FUNCIONANDO NORMAL PARA ESSE n
# 
# for (i in 1:100){
#   pb <- txtProgressBar(min = 0, max = n, style = 3)
#   y<-rURGo(n,mu_true,sigma_true)
#   fit1<-gamlss(y~1, family="URGo", c.crit = 0.001, n.cyc = 700,
#                  mu.step = .1, sigma.step = .1,trace = F)
# 
#   mu_result[i]<-logit_link$linkinv(fit1$mu.coefficients)
#   sigma_result[i]<-exp(fit1$sigma.coefficients)
#   setTxtProgressBar(pb, i)
# }



result1<- matrix(c(mu_true, mean(mu_result),
                   sigma_true, mean(sigma_result)),2,2)
colnames(result1)<-c("mu","sigma")
rownames(result1)<-c("true value","mean")
print(round(result1,2))




# result1<- matrix(c(mu_true, mean(mu_result),
#                    sigma_true, mean(sigma_result)),2,2)
# colnames(result1)<-c("mu","sigma")
# rownames(result1)<-c("true value","mean")
# print(round(result1,2))



# CHEKING WITH REGRESSORS
n<-1000 #amostra
X<-runif(n)
logit_link<-make.link("logit")
log_link<-make.link("log")
b1<-.7 #mu
b2<-3  #mu
mu_true<-logit_link$linkinv(b1+b2*X)
g1<-5 #sigma
g2<-1.5 #sigma
sigma_true<-log_link$linkinv(g1+g2*X)
R<-10000
mu_result<-sigma_result<-matrix(NA,R,2)

for (i in 1:R) {
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  
  y<-rURGo(n,mu_true,sigma_true)
  fit1<-gamlss(y~X,sigma.formula =~ X, family=URGo(sigma.link = "log"), c.crit = 0.001, n.cyc = 700,
               mu.step = .1, sigma.step = .1, trace = F)
  mu_result[i,]<-fit1$mu.coefficients
  sigma_result[i,]<-fit1$sigma.coefficients
  
  setTxtProgressBar(pb, i)
}

true_values<-c(b1,b2, g1,g2)
mean_values<-c(apply(mu_result,2,mean),
               apply(sigma_result,2,mean))
b_values<-(true_values-mean_values)/true_values*100
eqm_values<-c(apply(mu_result,2,var),
              apply(sigma_result,2,var))+(true_values-mean_values)^2
result1<- cbind(true_values,
                mean_values,
                b_values,
                eqm_values
)
colnames(result1)<-c("true value","mean","bias","eqm")
rownames(result1)<-c("b1","b2","g1","g2")
print(round(result1,2))









