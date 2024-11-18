tau=.5;mu=.85;sigma=.2;x=.63

ll<-expression(
  log(sigma)+log(-log(1-tau))+((sigma*x)/(1-x))-2*log(1-x)-log(exp((sigma*mu)/(1-mu))-1)+(exp((sigma*x)/(1-x))-1)/(exp((sigma*mu)/(1-mu))-1) * log(1-tau)
  )

lmu<-D(ll,"mu")
eval(lmu)

# Resultado mu
exp_sigma_mu <- exp(sigma * mu / (1 - mu))
exp_sigma_x <- exp(sigma * x / (1 - x))
term1 <- -sigma / (1 - mu)^2 * exp_sigma_mu / (exp_sigma_mu - 1)
term2 <- -log(1 - tau) * (exp_sigma_x - 1) * sigma / (1 - mu)^2 * exp_sigma_mu / (exp_sigma_mu - 1)^2
print(term1 + term2)



#Resultado Sigma

lsigma<-D(ll,"sigma")
eval(lsigma)

exp_sigma_mu <- exp(sigma * mu / (1 - mu))
exp_sigma_x <- exp(sigma * x / (1 - x))
term1 <- 1 / sigma
term2 <- x / (1 - x)
term3 <- -mu / (1 - mu) * exp_sigma_mu / (exp_sigma_mu - 1)
term4 <- log(1 - tau) * ((x / (1 - x)) * exp_sigma_x * (exp_sigma_mu - 1) -(exp_sigma_x - 1) * (mu / (1 - mu)) * exp_sigma_mu) / (exp_sigma_mu - 1)^2
print(term1 + term2 + term3 + term4)





##### SEGUNDAS DERIVADAS 

# PARA TESTAR DEPOIS
d2ldmu2<-D(lmu,"mu")
eval(d2ldmu2)

d2ldd2<-D(lsigma,"sigma")
eval(d2ldd2)

# SEGUNDA DERIVADA mu 


exp_sigma_mu <- exp(sigma * mu / (1 - mu))  
exp_sigma_x <- exp(sigma * x / (1 - x))    
denom_mu <- (1 - mu)^2
exp_sigma_mu_1 <- exp_sigma_mu - 1
term1 <- sigma * exp_sigma_mu * (exp_sigma_mu + 1) / (denom_mu^2 * exp_sigma_mu_1^2)
log_tau <- log(1 - tau)
num <- (exp_sigma_x - 1) * exp_sigma_mu * sigma * (exp_sigma_mu + 1)
term2 <- -log_tau * (num / (denom_mu^2 * exp_sigma_mu_1^3))

print(term1 + term2)



# SEGUNDA DERIVDA sigma

exp_sigma_mu <- exp(sigma * mu / (1 - mu))  
exp_sigma_x <- exp(sigma * x / (1 - x))    

denom_mu <- (exp_sigma_mu - 1)
denom_x <- (exp_sigma_mu - 1)^2

term1 <- -1 / (sigma^2)

term3 <- -mu / (1 - mu) * (exp_sigma_mu * (exp_sigma_mu + 1)) / (denom_mu^2)

log_tau <- log(1 - tau)
num <- exp_sigma_x * (exp_sigma_mu - 1) * (x / (1 - x)) - (exp_sigma_x - 1) * exp_sigma_mu * (mu / (1 - mu))
denom <- denom_x^2
term4 <- log_tau * (num / denom)


 
print(term1 + term3 + term4)







