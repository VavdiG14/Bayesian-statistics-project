library(MASS)

#Generiranje podatkov za simulacije I - Bayes/Freq

gen.beta.data <- function(n = 100,r=0, napaka ){
  mu <- c(3, 3, 4, 1)
  Sigma <- rbind(c(1, r, r, r),
                 c(r, 1.0, r, r),
                 c(r, r, 1.0, r),
                 c(r, r, r, 1.0))
  podatki <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  #dodam še binarno spremenljivko
  podatki <- cbind(podatki, rbinom(n, size =0:1, prob=0.5)) 
  colnames(podatki) <- paste("X", 1:(length(mu)+1), sep="")
  
  # generiramo ciljno spremenljivko ("odvisno spremenljivko")
  b0 <- -4.5
  b1 <- 0.8
  b2 <- 0.6
  b3 <- 3
  b4 <- 0
  b5 <- 1.4
  
  y <- b0 + podatki[, "X1"]*b1 + podatki[, "X2"]*b2 + 
    podatki[, "X3"]*b3 + podatki[, "X4"]*b4 + podatki[, "X5"]*b5 + napaka
  # zdruzimo podatke
  # data.frame popravi konverzijo formata
  podatki <- data.frame(cbind(podatki, y))
  return(podatki)
}

#Generiranje podatkov za simulacije II - Variable selection
gen.beta.data.2 <- function(n = 100,r=0, napaka ){
  # generirano kovariate - originalne spremenljivke brez NA
  # definicijo mu in Sigma lahko damo tudi izven funkcije
  mu <- c(3, 3, 4, 1)
  Sigma <- rbind(c(1, r, r, r),
                 c(r, 1.0, r, r),
                 c(r, r, 1.0, r),
                 c(r, r, r, 1.0))
  mu.brez.vpliva <- c(6, 17, 22, 12, 5)
  Sigma.brez.vpliva <- rbind(c(1, r, r, r, r),
                             c(r, 1.0, r, r, r),
                             c(r, r, 1.0, r, r),
                             c(r, r, r, 1.0, r),
                             c(r, r, r, r, 1.0))
  
  podatki <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  podatki <- cbind(podatki, rbinom(n, size =0:1, prob=0.5))
  podatki.brez.vpliva <-mvrnorm(n = n, mu = mu.brez.vpliva, Sigma = Sigma.brez.vpliva)
  podatki <- cbind(podatki, podatki.brez.vpliva)
  colnames(podatki) <- paste("X", 1:ncol(podatki), sep="")
  
  # generiramo ciljno spremenljivko ("odvisno spremenljivko")
  b0 <- -4.5
  b1 <- 0.8
  b2 <- 0.6
  b3 <- 1.0
  b4 <- 0
  b5 <- 1.4
  b6 <- 0
  b7 <-0
  b8 <- 0
  b9 <-0
  b10 <- 0
  
  y <- b0 + podatki[, "X1"]*b1 + podatki[, "X2"]*b2 + 
    podatki[, "X3"]*b3 + podatki[, "X4"]*b4 + podatki[, "X5"]*b5 + 
    podatki[, "X6"]*b6 + podatki[, "X7"]*b7 + podatki[, "X8"]*b8 + 
    podatki[, "X9"]*b9 + podatki[, "X10"]*b10 + napaka
  
  # data.frame popravi konverzijo formata
  podatki <- data.frame(cbind(podatki,"Y"= y))
  return(podatki)
}
