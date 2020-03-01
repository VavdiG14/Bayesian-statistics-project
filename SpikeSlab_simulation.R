library(BoomSpikeSlab)
library(gridExtra)
library(magrittr)
source("generation_data.R")

selection.spike.slab <- function(dataSetSelectX, Y, alpha = 0.2){
  dataX <- cbind(rep(1, nrow(dataSetSelectX)), as.matrix(dataSetSelectX))
  prior <- IndependentSpikeSlabPrior(dataX, Y,
                                     expected.model.size = 5, # expect 3 nonzero predictors
                                     prior.df = .0, # weaker prior than the default
  )
  dataX <- dataX[,-1]
  model <- lm.spike(Y ~ dataX,  niter = 1000, prior = prior)
  return(model)
}


pon <- 100
sampleSize <- 100
vrsta.napake <- c("N(0,1)", "N(0,3)", "Hi_1", "Hi_4")
#vrsta.metode <- c("AIC", "BIC", "alpha")
korelacije <- c(0, 0.8)
zasnova <- expand.grid(korelacije, vrsta.napake)
zasnova <- do.call(rbind, replicate(pon, zasnova, simplify=FALSE)) %>%
  `colnames<-`(c("Korelacija","Napaka"))
matrika.rez <- matrix(NA, nrow = nrow(zasnova), ncol = 13)
i = 1

while(i <= nrow(zasnova)){
  vrsta.napake.i <- zasnova[i, "Napaka"]
  if(vrsta.napake.i == "N(0,1)"){
    error <- rnorm(sampleSize, 0, 1)
  }
  else if(vrsta.napake.i == "N(0,3)"){
    error <- rnorm(sampleSize, 0, 3)
  }
  else if(vrsta.napake.i == "Hi_1"){
    error <- rchisq(sampleSize, df = 1)
  }
  else if(vrsta.napake.i == "Hi_4"){
    error <- rchisq(sampleSize, df = 4)
  }
  r.i <- zasnova[i, "Korelacija"]
  
  data.Spike <- gen.beta.data.2(sampleSize, r = r.i, napaka = error)
  data.Spike.selectY <- data.Spike$Y
  data.Spike.selectX <- sweep(data.Spike[,-5], 2, colMeans(data.Spike[,-5]), FUN="-")[,-10] #centriramo
  data.Spike.selectX$X5 <- data.Spike$X5
  data.Spike.selectX <- data.Spike.selectX[, c(1,2,3,4,10,5,6,7,8,9)]
  
  mod.spikeSlab <- selection.spike.slab(dataSetSelectX = data.Spike.selectX, Y =data.Spike.selectY)
  inc.prob <- summary(mod.spikeSlab)$coef[,"inc.prob"]
  inc.prob.urejen <- inc.prob[order(factor(names(inc.prob), levels = paste("dataXX", 0:10, sep="")))]

  matrika.rez[i,] <- c("R" = r.i,"Napaka" = vrsta.napake.i, inc.prob.urejen)
  i <- i + 1
}

spikeSlab.df <- as.data.frame(matrika.rez[,-13])
spikeSlab.df$V2 <- factor(spikeSlab.df$V2,labels = vrsta.napake)
colnames(spikeSlab.df) <- c("Korelacija", "Napaka", paste("X", 1:10, sep =""))
saveRDS(spikeSlab.df, "data/spikeSlab_simulation.RDS")
