library(MASS)
library(R2BayesX)
gen.beta.data <- function(n = 100, napaka ){
  # generirano kovariate - originalne spremenljivke brez NA
  # definicijo mu in Sigma lahko damo tudi izven funkcije
  mu <- c(3, 3, 4, 1)
  Sigma <- rbind(c(0.5, 0.1, 0.2, 0.3),
                 c(0.1, 1.0, 0.4, 0.1),
                 c(0.2, 0.4, 1.0, 0.2),
                 c(0.3, 0.1, 0.2, 1.0))
  
  podatki <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  podatki <- cbind(podatki, rbinom(n, size =0:1, prob=0.5))
  colnames(podatki) <- paste("X", 1:(length(mu)+1), sep="")
  #dodam še binarno spremenljivko

  # generiramo ciljno spremenljivko ("odvisno spremenljivko")
  b0 <- -4.5
  b1 <- 0.8
  b2 <- 0.6
  b3 <- 1.0
  b4 <- 0
  b5 <- 1.4
  
  y <- b0 + podatki[, "X1"]*b1 + podatki[, "X2"]*b2 + 
    podatki[, "X3"]*b3 + podatki[, "X4"]*b4 + podatki[, "X5"]*b5 + napaka
  # summary(lm(y ~ podatki))
  
  # zdruzimo podatke
  # data.frame popravi konverzijo formata
  podatki <- data.frame(cbind(podatki, y))
  return(podatki)
}

bete <- data.frame()
sd.bete <- data.frame()
for(i in 1:1000){
  data.lm <- gen.beta.data(n=100, napaka = rnorm(100,0,100))
  lm.model <- lm(y~., data = data.lm)
  bete.hat <- summary(lm.model)$coef[,1]
  sd.bet <- summary(lm.model)$coef[,2]
  bete <- rbind(bete,bete.hat)
  sd.bete <- rbind(sd.bet, sd.bete)
  
  bayesx.norm <- bayesx(y ~X1+X2+X3+X4+X5 , data = data.lm, family = "gaussian", method = "MCMC")
  bayesx.mu <- attr(bayesx.norm$fixed.effects, "sample")[,1] #za vsako beto posebaj
  bayesx.sigma2 <- attr(bayesx.norm$variance, "sample")
  par(mfrow = c(2, 2))
  plot(bayesx.mu, type = "l", main = "Povprecje,\nveriga", xlab = "")
  plot(bayesx.sigma2, type = "l", main = "Varianca,\nveriga", xlab = "")
  hist(bayesx.mu, prob = T, main = "Povprecje,\nrobna aposteriorna porazdelitev")
  lines(density(bayesx.mu), col = "red", lwd = 2)
  
}
colnames(bete) <- paste("B", 0:5, sep="")
colMeans(bete)
  
