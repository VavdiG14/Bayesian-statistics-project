library(basicMCMCplots)
library(nimble)
library(coda)
source("generation_data.R")

#Prvi scenarij
data.NIBLE.N1 <- gen.beta.data.2(100, napaka = rnorm(100))
data.NIBLE.N1.selectY <- data.NIBLE.N1$Y
data.NIBLE.N1.selectX <- sweep(data.NIBLE.N1[,-5], 2, colMeans(data.NIBLE.N1[,-5]), FUN="-")[,-10] #centriramo
data.NIBLE.N1.selectX$X5 <- data.NIBLE.N1$X5
data.NIBLE.N1.selectX <- data.NIBLE.N1.selectX[, c(1,2,3,4,10,5,6,7,8,9)]


codeSelect <- nimbleCode({
  sigma ~ dunif(0, 20)
  psi ~ dunif(0,1)
  beta0 ~ dnorm(0, sd=100)
  for(i in 1:p) {
    z[i] ~ dbern(psi) #indikator za vsak koeficient
    beta[i] ~ dnorm(0, sd = 100)
    zbeta[i] <- z[i] * beta[i]
  }
  for(i in 1:N) {
    y[i] ~ dnorm(beta0 + inprod(X[i, 1:p], zbeta[1:p]), sd = sigma)
  }
})

N <- dim(data.NIBLE.N1.selectX)[1]
p <- dim(data.NIBLE.N1.selectX)[2]
constantsSelect <- list(N = N, p = p)
initsSelect <- list(sigma = 1, psi = 0.5, beta0 = 0,
                    beta = rnorm(p),
                    z = sample(c(0, 1), p, replace = TRUE))

dataSelect <- list(y = data.NIBLE.N1.selectY, X = data.NIBLE.N1.selectX)

RmodelRJ <- nimbleModel(code = codeSelect, constants = constantsSelect, #isto
                        inits = initsSelect, data = dataSelect)

confRJ <- configureMCMC(RmodelRJ) #isto

confRJ$addMonitors('z') #isto

configureRJ(confRJ,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2))

RmcmcRJ <- buildMCMC(confRJ)
CmodelRJ <- compileNimble(RmodelRJ)
CmcmcRJ <- compileNimble(RmcmcRJ, project = CmodelRJ)
samplesRJ <- runMCMC(CmcmcRJ, niter = 12000, nburnin = 2000)


grafGG <- function(model, naslov = ""){
  shrani.rez <- as.data.frame(round(samplesSummary(model),2))
  shrani.rez$ime <- rownames(shrani.rez)
  z.izbrani <- shrani.rez %>% 
    filter(grepl("z",ime)) %>% 
    dplyr::select(Mean) %>% 
    mutate(X = paste("X",1:10, sep =""))
  gg <- ggplot(z.izbrani, aes(x = X, y = Mean)) + geom_bar(stat = "identity")+
    ggtitle(paste("Izbira spremenljivk pri", naslov))
  return(gg)
}


graf.nicneBete <- function(model, imeFile="img/nimble_bete69.png"){
  samplesPlot(model, var = c("beta[10]", "beta[9]", "beta[8]", "beta[7]", "beta[6]", "beta[4]"),
            file = imeFile)
}

graf.nenicneBete <- function(model, imeFile="img/nimble_bete69.png"){
  samplesPlot(model, var = c("beta[1]", "beta[2]", "beta[3]", "beta[5]", "beta[0]"),
              file = imeFile)
}

data.NIBLE.N1 <- gen.beta.data.2(100, napaka = rnorm(100))
data.NIBLE.N1.selectY <- data.NIBLE.N1$Y
data.NIBLE.N1.selectX <- sweep(data.NIBLE.N1[,-5], 2, colMeans(data.NIBLE.N1[,-5]), FUN="-")[,-10] #centriramo
data.NIBLE.N1.selectX$X5 <- data.NIBLE.N1$X5
data.NIBLE.N1.selectX <- data.NIBLE.N1.selectX[, c(1,2,3,4,10,5,6,7,8,9)]
mod.N1 <- get.rjMCMC(data.NIBLE.N1.selectX, data.NIBLE.N1.selectY)

grafGG(mod.N1, "modelu z nakljucno napako: N(0,1)")
graf.nenicneBete(mod.N1, "img/N1_neBeta.png")
graf.nicneBete(mod.N1, "img/N1_Beta.png")
effectiveSize(mod.N1)


data.NIBLE.N3 <- gen.beta.data.2(100, napaka = rnorm(100))
data.NIBLE.N3.selectY <- data.NIBLE.N3$Y
data.NIBLE.N3.selectX <- sweep(data.NIBLE.N3[,-5], 2, colMeans(data.NIBLE.N3[,-5]), FUN="-")[,-10] #centriramo
data.NIBLE.N3.selectX$X5 <- data.NIBLE.N3$X5
data.NIBLE.N3.selectX <- data.NIBLE.N3.selectX[, c(1,2,3,4,10,5,6,7,8,9)]
mod.N3 <- get.rjMCMC(data.NIBLE.N3.selectX, data.NIBLE.N3.selectY)

grafGG(mod.N3, "modelu z nakljucno napako: N(0,3)")
graf.nenicneBete(mod.N3, "img/N3_neBeta.png")
graf.nicneBete(mod.N3, "img/N3_Beta.png")
effectiveSize(mod.N3)


data.NIBLE.H1 <- gen.beta.data.2(100, napaka = rchisq(100, 1))
data.NIBLE.H1.selectY <- data.NIBLE.H1$Y
data.NIBLE.H1.selectX <- sweep(data.NIBLE.H1[,-5], 2, colMeans(data.NIBLE.H1[,-5]), FUN="-")[,-10] #centriramo
data.NIBLE.H1.selectX$X5 <- data.NIBLE.H1$X5
data.NIBLE.H1.selectX <- data.NIBLE.H1.selectX[, c(1,2,3,4,10,5,6,7,8,9)]
mod.H1 <- get.rjMCMC(data.NIBLE.H1.selectX, data.NIBLE.H1.selectY)

grafGG(mod.H1, "modelu z nakljucno napako: Hi_1")
graf.nenicneBete(mod.H1, "img/H1_neBeta.png")
graf.nicneBete(mod.H1, "img/H1_Beta.png")

data.NIBLE.H4 <- gen.beta.data.2(100, napaka = rchisq(100, 4))
data.NIBLE.H4.selectY <- data.NIBLE.H4$Y
data.NIBLE.H4.selectX <- sweep(data.NIBLE.H4[,-5], 2, colMeans(data.NIBLE.H4[,-5]), FUN="-")[,-10] #centriramo
data.NIBLE.H4.selectX$X5 <- data.NIBLE.H4$X5
data.NIBLE.H4.selectX <- data.NIBLE.H4.selectX[, c(1,2,3,4,10,5,6,7,8,9)]
mod.H4 <- get.rjMCMC(data.NIBLE.H4.selectX, data.NIBLE.H4.selectY)


grafGG(mod.H4, "modelu z nakljucno napako: Hi_4")
graf.nenicneBete(mod.H4, "img/H4_neBeta.png")
graf.nicneBete(mod.H4, "img/H4_Beta.png")
