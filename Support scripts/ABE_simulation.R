library(gridExtra)
library(magrittr)
library(abe)
source("generation_data.R")

pon <- 100
sampleSize <- 100
vrsta.napake <- c("N(0,1)", "N(0,3)", "Hi_1", "Hi_4")
#vrsta.metode <- c("AIC", "BIC", "alpha")
korelacije <- c(0, 0.8)

zasnova <- expand.grid(korelacije, vrsta.napake)
zasnova <- do.call(rbind, replicate(pon, zasnova, simplify=FALSE)) %>%
  `colnames<-`(c("Korelacija","Napaka"))
matrika.rez <- matrix(NA, nrow = 3*nrow(zasnova), ncol = 14)
i = 1
j = 1
while(i < nrow(zasnova)){
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
    
  data.ABE <- gen.beta.data.2(sampleSize, r = r.i, napaka = error)
  data.ABE.selectY <- data.ABE$Y
  data.ABE.selectX <- sweep(data.ABE, 2, colMeans(data.ABE), FUN="-")[,-11] #centriramo
  
  abe.AIC <-selection.ABE(data.ABE.selectX, Y =data.ABE.selectY, 
                             p = 10, metrika = "AIC",num.boot = 500)
  abe.BIC <-selection.ABE(data.ABE.selectX, Y =data.ABE.selectY, 
                          p = 10, metrika = "BIC",num.boot = 500)
  abe.alpha <-selection.ABE(data.ABE.selectX, Y =data.ABE.selectY, 
                          p = 10, metrika = "alpha",num.boot = 500)
  
  rez.aic <- summary(abe.AIC)$var.rel.frequencies
  rez.bic <- summary(abe.BIC)$var.rel.frequencies
  rez.alpha <- summary(abe.alpha)$var.rel.frequencies
  matrika.rez[j,] <- c("R" = r.i,  "Metoda" = 1, "Napaka" = vrsta.napake.i, rez.aic)
  matrika.rez[j+1,] <- c("R" = r.i,  "Metoda" = 2, "Napaka" = vrsta.napake.i, rez.bic)
  matrika.rez[j+2,] <- c("R" = r.i,  "Metoda" = 3, "Napaka" = vrsta.napake.i, rez.alpha)
  j <- j + 3 
  i <- i + 1
}


rez.df <- na.omit(as.data.frame(matrika.rez))
rez.df.nov <- as.data.frame(apply(rez.df, 2, as.numeric))
rez.df.nov$Metoda <- factor(rez.df.nov$V2, labels = c("AIC", "BIC", "alpha"))
rez.df.nov$Napaka <- factor(rez.df.nov$V3,labels = vrsta.napake)
rez.df.nov$Korelacija <- rez.df.nov$V1
abe.rez <- rez.df.nov[,4:17]
colnames(abe.rez) <- c( paste("X", 0:10, sep=""),"Metoda", "Napaka","Korelacija")

saveRDS(abe.rez,"data/abe_simulation.RDS")


abe.data <- readRDS("data/abe_simulation.RDS")

aba.data.sim <- abe.data %>%
  group_by(Korelacija, Metoda, Napaka)%>%
  summarise_all(mean) %>% 
  gather(key = "X", "Vrednost", -c(Korelacija, Metoda, Napaka))  %>%
  filter(X != "X0")


ggplot(aba.data.sim) + 
  geom_bar( aes(x = X, y = Vrednost, color = Metoda, fill = Metoda ), position = "dodge2", stat = "identity")+
  facet_grid(cols = vars(Napaka), rows = vars(Korelacija))+
  ggtitle("Pristranskost ocenjenih parametrov", subtitle = "Število simulacij: 1000")

