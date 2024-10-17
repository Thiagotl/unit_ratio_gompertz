# Exemplo de cálculo para o vetor escore da Kumaraswamy
# parametros q e phi, de quantil e forma respectivamente
# Renata Rojas Guerra (renata.r.guerra@ufsm.br), january/2023

# função densidade
dkum<-function(y,q,varphi,tau=.5)
{
  d<-log(1-tau)/log(1-q^varphi)*
    varphi*y^(varphi-1)*(1-y^varphi)^(log(1-tau)/log(1-q^varphi)-1)
  d
}

#calculando a log-verossimilhança analítica
loglik<-expression(log(
  log(1-tau)/log(1-q^phi)*phi*y^(phi-1)*
    (1-y^phi)^(log(1-tau)/log(1-q^phi)-1))
  )
# definindo valores para testar
q=.7;phi=1.7;tau=.5;y<-.23
# comparando a loglik analitica e numerica
log(dkum(y,q,phi)) # numerica
eval(loglik) # analitica

# calculando a derivada da loglik com respeito a q
D(loglik,"q") # analitica do R
#avaliando nos valores definidos para testar
eval(D(loglik,"q"))
# analitica simplificada
phi*q^(phi-1)/((1-q^phi)*log(1-q^phi))*
  (log(1-tau)/log(1-q^phi)*log(1-y^phi)+1)

# calculando a derivada da loglik com respeito a phi
D(loglik,"phi") # analitica do R
#avaliando nos valores definidos para testar
eval(D(loglik,"phi"))
# analitica simplificada
1/phi+log(y)+((q^(phi-1))/((1-q^phi)*log(1-q^phi))*
                ((log(1-tau)/log(1-q^phi))*log(1-y^phi)+1))*q*log(q)-
  ((log(1-tau)/log(1-q^phi))-1)*y^phi*log(y)/(1-y^phi)
# podemos simplificar mais
# definimos delta e c como
delta=(log(1-tau)/log(1-q^phi))
ct=(q^(phi-1))/((1-q^phi)*log(1-q^phi))*
  (delta*log(1-y^phi)+1)
# entao
1/phi+log(y)+ct*q*log(q)-
  (delta-1)*y^phi*log(y)/(1-y^phi)
