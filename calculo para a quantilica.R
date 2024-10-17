
mu= 0.7
tau= 0.5
sigma= 5
q=0.3




u1 = 1-(1-tau)^((exp(sigma*q/(1-q))-1)/(exp(sigma*mu/(1-mu))-1))


log(1-u1)/log(1-tau)
(exp(sigma*q/(1-q))-1)/(exp(sigma*mu/(1-mu))-1)

(exp(sigma*q/(1-q))-1)
log(1-u1)/log(1-tau)*(exp(sigma*mu/(1-mu))-1)

(exp(sigma*q/(1-q)))
log(1-u1)/log(1-tau)*(exp(sigma*mu/(1-mu))-1)+1

sigma*q/(1-q)
log(log(1-u1)/log(1-tau)*(exp(sigma*mu/(1-mu))-1)+1)

q/(1-q)
(log(log(1-u1)/log(1-tau)*(exp(sigma*mu/(1-mu))-1)+1))/sigma

a=(log(log(1-u1)/log(1-tau)*(exp(sigma*mu/(1-mu))-1)+1))

q
((log(log(1-u1)/log(1-tau)*(exp(sigma*mu/(1-mu))-1)+1))/sigma)*(1-q)

q
((a/sigma)-((a*q)/sigma))

(q+((a*q)/sigma))
(a/sigma)

q*(sigma+a)/sigma 
a/sigma

q  
a/(sigma+a)
