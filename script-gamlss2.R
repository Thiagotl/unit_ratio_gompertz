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
                 sigma.initial = expression(sigma<- rep(5, length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1),
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0 & y < 1)
  ),
  class = c("gamlss.family", "family"))
}

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------
# quantile function
qURGo<-function(u,mu, sigma, tau = 0.5)
{
  # log(log(1-u)/log(1-tau)*(exp(sigma*mu/(1-mu))-1)+1)
  q<- log(log(1-u)/log(1-tau)*(exp(sigma*mu/(1-mu))-1)+1)/
    (sigma + log(log(1-u)/log(1- tau)*(exp(sigma * mu/(1 - mu)) -1)+1))
  q
}



#-------------------------------------------------------------------------------
# inversion method for randon generation

rURGo<-function(n,mu, sigma){
  tau=0.5
  u<- runif(n)
  y<- qURGo(u,mu=mu,sigma=sigma)
  return(y)
}


library(gamlss)


# Definir parâmetros globais
vn <- c(30, 70, 150, 300) # Tamanhos de amostra
logit_link <- make.link("logit")
log_link <- make.link("log")
b1 <- 0.7  # mu
b2 <- 0.3  # mu
g1 <- .7 # sigma
g2 <- .25  # sigma
R <- 50  # Número de repetições

set.seed(10)

# Função auxiliar para calcular métricas
calculate_metrics <- function(mu_result, sigma_result, true_values) {
  mean_values <- c(
    apply(mu_result, 2, mean, na.rm = TRUE),
    apply(sigma_result, 2, mean, na.rm = TRUE)
  )
  bias_values <- (true_values - mean_values) / true_values * 100
  eqm_values <- c(
    apply(mu_result, 2, var, na.rm = TRUE),
    apply(sigma_result, 2, var, na.rm = TRUE)
  ) + (true_values - mean_values)^2
  
  result <- cbind(
    "True Value" = true_values,
    "Mean" = mean_values,
    "Bias (%)" = bias_values,
    "EQM" = eqm_values
  )
  rownames(result) <- c("b1", "b2", "g1", "g2")
  return(result)
}

# Inicialização de resultados finais
#final_results <- list()
bug_counter <- 0  
final_results <- data.frame()

for (n in vn) {
  
  mu_result <- matrix(NA, R, 2)  # Para b1, b2
  sigma_result <- matrix(NA, R, 2)  # Para g1, g2
  X <- runif(n)
  
  mu_true <- logit_link$linkinv(b1+b2*X)
  # mean(mu_true)
  sigma_true <- log_link$linkinv(g1+g2*X)
  # mean(sigma_true)
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  
  i<-0
  while (i < R) {
    y <- rURGo(n, mu_true, sigma_true)  # Geração dos dados
    
    # Ajustar o modelo 
    fit1 <- try(
        gamlss(
          y ~ X, sigma.formula = ~ X,
          family = URGo(sigma.link = "log"),
          # sigma.start = start,
          c.crit = 0.001,
          n.cyc = 700,
          mu.step = 0.1,
          sigma.step = 0.1,
          trace = FALSE
        ),
        silent = TRUE
      )
    
    if (inherits(fit1, "try-error")) {
      bug_counter <- bug_counter + 1
      next
    }
    
    i<-i+1
    # print(i)
    #  coeficientes ajustados
    mu_result[i, ] <- fit1$mu.coefficients
    sigma_result[i, ] <- fit1$sigma.coefficients
    
    setTxtProgressBar(pb, i)
  }
  # print(sigma_result)
  close(pb)
  
  # Calcular métricas
  true_values <- c(b1, b2, g1, g2)
  result <- calculate_metrics(mu_result, sigma_result, true_values)
  result <- as.data.frame(result)
  result$Sample_Size <- n
  
  # result_file <- paste0("simulation_results_n_", n, ".txt")
  # write.table(round(result, 2), file = result_file, sep = "\t", col.names = NA, quote = FALSE)
  # final_results[[as.character(n)]] <- result
  
  
  final_results <- rbind(final_results, result)
  
  
  # Armazenar e exibir resultados
  #final_results[[as.character(n)]] <- result
  cat("\nTamanho da amostra:", n, "\n")
  print(round(result, 2))
  cat("N de erros no ajuste do modelo:", bug_counter, "\n")
}

# Exibir contadores de erros
# cat("\nNúmero total de erros no ajuste do modelo:", bug_counter, "\n")


## SAIDA TXT -------------------------------------------------------------------
final_results_aligned <- as.data.frame(
  lapply(final_results, function(x) {
    if (is.numeric(x)) {
      format(round(x, 2), nsmall = 2, justify = "right")  # Arredondar e alinhar
    } else {
      format(x, justify = "right")  # Apenas alinhar texto
    }
  })
)

final_results_aligned$Sample_Size <- format(as.integer(final_results$Sample_Size), justify = "right")

write.table(
  final_results_aligned,
  file = "all_simulation_results.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# SAIDA LATEX -----------------------------------------------------------------

# Instalar os pacotes
if (!require(kableExtra)) install.packages("kableExtra")

# Formatar os números no data.frame
final_results_aligned <- as.data.frame(
  lapply(final_results, function(x) {
    if (is.numeric(x)) {
      round(x, 2)  # Arredondar para 2 casas decimais
    } else {
      x
    }
  })
)


library(kableExtra)
kable_latex <- kable(
  final_results_aligned,
  format = "latex",
  booktabs = TRUE,           
  col.names = c("True Value", "Mean", "Bias (%)", "EQM", "Sample Size"),
  align = "c",               
  caption = "Simulation Results"
) %>%
  kable_styling(
    latex_options = c("striped", "hold_position"),  
    full_width = FALSE
  )

# Salvar a tabela em um arquivo LaTeX
output_file <- "all_simulation_results_kable.tex"
cat(kable_latex, file = output_file)

cat("\nTabela LaTeX gerada com kableExtra salva no arquivo 'all_simulation_results_kable.tex'.\n")





## SAÍDA TXT PARA LATEX --------------------------------------------------------

# Carregar o pacote xtable
if (!require(xtable)) install.packages("xtable")
library(xtable)

# Formatar os números no data.frame
final_results_aligned <- as.data.frame(
  lapply(final_results, function(x) {
    if (is.numeric(x)) {
      round(x, 2)  
    } else {
      x
    }
  })
)

# Gerar tabela LaTeX com xtable
latex_table <- xtable(final_results_aligned,
                      caption = "Simulation Results",
                      align = c("c", "c", "c", "c", "c", "c"))

# Salvar em um arquivo .txt
output_file <- "all_simulation_results_xtable.txt"
sink(output_file)
print(latex_table, include.rownames = FALSE, booktabs = TRUE)
sink()

## GRÁFICOS -------------------------------------------------------------------


# PARA VALORES ALEATÓRIOS 
pdf.plot(
  obj = fit1,
  obs = 1:8,
  from = 0.001,  
  to = 0.999,    
  no.points = 201,
  y.axis.lim = 5
)


# PARA MU

pdf.plot(
  family = URGo(sigma.link = "log"),
  mu = c(0.3,0.5,0.7,0.9),
  sigma = 0.7,
  from = 0.001,
  to = 0.999,
  no.points = 201,
  y.axis.lim = 5
)

# PARA SIGMA

pdf.plot(
  family = URGo(sigma.link = "log"),
  mu = 0.5,
  sigma = c(0.3, 0.5, 0.7, 0.9),
  from = 0.001,
  to = 0.999,
  no.points = 201,
  y.axis.lim = 5
)



library(ggplot2)

y_values<-seq(0.001, 0.999, length.out = 100)
sigma_values <- c(0.2, 0.3, 0.5, 0.7)
mu <- 0.5


results <- data.frame(
  y = rep(y_values, length(sigma_values)),
  sigma = rep(sigma_values, each = length(y_values)),
  density = NA
)


for (s in sigma_values) {
  results$density[results$sigma == s] <- dURGo(y_values, mu, s)
}


ggplot(results, aes(x = y, y = density, color = as.factor(sigma))) +
  geom_line(size = 1) +
  labs(
    title = "Função de Densidade para Diferentes Valores de Sigma",
    x = "y",
    y = "Densidade",
    color = "Sigma"
  ) +
  theme_minimal()

