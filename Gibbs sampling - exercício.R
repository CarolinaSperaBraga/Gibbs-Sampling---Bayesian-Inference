# Biblioteca utilizada
library(invgamma)

set.seed(6)

M <- 10000
beta0 <- 1
alpha0 <- 2
sigma20 <- 100
a <- 0.0001
b <- 0.0001
n <- 16
k <- 10
alpha <- alpha0
beta <- beta0
sigma2 <- sigma20
amostra.alpha <- c(alpha0)
amostra.beta <- c(beta0)
amostra.sigma2 <- c(sigma20)
x <- runif(M,0,1)
y <- rnorm(M, alpha+beta*x, sigma2^(1/2))

# Obtenção da amostra
for(i in 1:M){
  # Amostra para alpha
  media_alpha = (sum(x*y)/sigma2 - beta*sum(x)/sigma2)/(sum(x^2)/sigma2 + 1/(100^2))
  se_alpha = ((sum(x^2)/sigma2 + 1/(100^2))^(-1/2))
  alpha <- rnorm(1,media_alpha,se_alpha)
  amostra.alpha <- c(amostra.alpha,alpha)
  
  # Amostra para beta
  media_beta <- (sum(y)/sigma2 + alpha*sum(x)/sigma2)/(M/sigma2 + 1/(100^2))
  se_beta <- ((M/sigma2 + 1/(100^2))^(-1/2))
  beta<- rnorm(1,media_beta, se_beta)
  amostra.beta <- c(amostra.beta,beta)
  
  # Amostra para sigma2
  sigma2 <- rinvgamma(1, M/2 + a, b + (1/2)*sum(y - alpha - beta*x)^2)
  amostra.sigma2 <- c(amostra.sigma2, sigma2)
}


# Análise gráfica
par(mfrow=c(2,3))
ts.plot(amostra.alpha)
ts.plot(amostra.beta)
ts.plot(amostra.sigma2)
## alpha: variável aleatória contínua
hist(amostra.alpha,col="lightgray",border="hotpink")
## beta: variável aleatória contínua
hist(amostra.beta,col="lightgray",border="hotpink")
## sigma2: variável aleatória contínua
hist(amostra.sigma2,col="lightgray",border="hotpink")


# Resumo a posteriori baseado na amostra (am)
require(coda)

tab1 <-round(data.frame(media=c(mean(amostra.alpha), mean(amostra.beta), mean(amostra.sigma2)),
                        mediana=c(median(amostra.alpha), median(amostra.beta), median(amostra.sigma2))),4)
row.names(tab1) <- c("alpha","beta", "sigma2")
tab1


tab2 <- data.frame(IC95.alpha = quantile(amostra.alpha,c(0.025,0.975)),
                   IC95.beta = quantile(amostra.beta,c(0.025,0.975)),
                   IC95.sigma2 = quantile(amostra.sigma2,c(0.025,0.975)))

tab2


tab3 <- data.frame(
  HPDInt.alpha = HPDinterval(as.mcmc(amostra.alpha),0.95)[1,],
  HPDInt.beta = HPDinterval(as.mcmc(amostra.beta),0.95)[1,],
  HPDInt.sigma2 = HPDinterval(as.mcmc(amostra.sigma2),0.95)[1,])
tab3



# Com burn-in e thin
M <- 10000
m <- 1000 #burn-in
k <- 10 #thin
length(seq(m+1,M,k))

# Resumo a posteriori baseado na amostra (am) com burn-in e thin
amostra.alphabt <- amostra.alpha[seq(m+1,M-m,k)]
amostra.betabt <- amostra.beta[seq(m+1,M-m,k)]
amostra.sigma2bt <- amostra.sigma2[seq(m+1,M-m,k)]
tab1 <-round(data.frame(media=c(mean(amostra.alphabt), mean(amostra.betabt), mean(amostra.sigma2bt)),
                        mediana=c(median(amostra.alphabt), median(amostra.betabt), median(amostra.sigma2bt))),4)
row.names(tab1) <- c("alpha","beta","sigma2")
tab1


tab2 <- data.frame(IC95.alpha = quantile(amostra.alphabt,c(0.025,0.975)),
                   IC95.beta = quantile(amostra.betabt,c(0.025,0.975)),
                   IC95.sigma2 = quantile(amostra.sigma2bt,c(0.025,0.975)))

tab2


tab3 <- data.frame(
  HPDInt.alpha = HPDinterval(as.mcmc(amostra.alphabt),0.95)[1,],
  HPDInt.beta = HPDinterval(as.mcmc(amostra.betabt),0.95)[1,],
  HPDInt.sigma2 = HPDinterval(as.mcmc(amostra.sigma2bt),0.95)[1,])
tab3



# Veja os gráficos
# alpha
par(mfrow=c(1,2))
acf(amostra.alpha)
acf(amostra.alphabt)

#beta
par(mfrow=c(1,2))
acf(amostra.beta)
acf(amostra.betabt)

#sigma2
par(mfrow=c(1,2))
acf(amostra.sigma2)
acf(amostra.sigma2bt)
