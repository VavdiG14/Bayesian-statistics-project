library(gridExtra)
library(magrittr)
library(MASS)
library(R2BayesX)
source("generation_data.R")

pon <- 1000
sampleSize <- 100
vrsta.napake <- c("N(0,1)", "N(0,3)", "Hi_1", "Hi_4")
korelacije <- c(0, 0.4, 0.8)
metoda <- c("Freq", "Bayes")
zasnova <- expand.grid(metoda, vrsta.napake,korelacije)
zasnova.lm <- do.call(rbind, replicate(pon, zasnova, simplify=FALSE)) %>% `colnames<-`(c("Metoda","Napaka", "Korelacija"))
zasnova.osnova <- do.call(rbind, replicate(pon, expand.grid(vrsta.napake,korelacije), simplify=FALSE)) %>% `colnames<-`(c("Napaka", "Korelacija"))

skupaj.df <- as.data.frame(matrix(NA, nrow = nrow(zasnova.lm), ncol = 9))
table.sd <- as.data.frame(matrix(NA, nrow = nrow(zasnova.lm), ncol = 9))

start.time <- Sys.time()
i <- 1
j <- 1
while(i <= nrow(zasnova.osnova)){
  napaka.i <- as.character(zasnova.osnova[i, "Napaka"])
  #print(napaka.i)
  if(napaka.i == "N(0,1)"){
    error <- rnorm(sampleSize, 0, 1)
  }
  else if(napaka.i == "N(0,3)"){
    error <- rnorm(sampleSize, 0, 3)
  }
  else if(napaka.i == "Hi_1"){
    error <- rchisq(sampleSize, df = 1)
  }
  else if(napaka.i == "Hi_4"){
    error <- rchisq(sampleSize, df = 4)
  }
  r.i <- zasnova.osnova[i, "Korelacija"]
  data.i <- gen.beta.data(n = sampleSize, r= r.i, napaka = error)
  lm.model <- lm(y~., data = data.i)
  bete.hat <- summary(lm.model)$coef[,1]
  sd.bet <- summary(lm.model)$coef[,2]
  index.i <- which(zasnova.lm$Napaka == napaka.i & zasnova.lm$Korelacija==r.i& zasnova.lm$Metoda == "Freq")
  index.i.bay <- which(zasnova.lm$Napaka == napaka.i & zasnova.lm$Korelacija==r.i& zasnova.lm$Metoda == "Bayes")
  
  
  bayesx.norm <- bayesx(y ~ X1+X2+X3+X4+X5 , data = data.i, family = "gaussian", method = "MCMC")
  bayesx.mu <- attr(bayesx.norm$fixed.effects, "sample") #za vsako beto posebaj
  bayes.mu.mean.bet <- apply(bayesx.mu, 2, mean)
  bayes.mu.sd.bet <- apply(bayesx.mu, 2, sd)
  
  skupaj.df[j, ] <- cbind(zasnova.lm[index.i.bay, ], "Beta0" = bayes.mu.mean.bet[1], 
        "Beta1" = bayes.mu.mean.bet[2], 
        "Beta2" = bayes.mu.mean.bet[3], 
        "Beta3" = bayes.mu.mean.bet[4], 
        "Beta4" = bayes.mu.mean.bet[5], 
        "Beta5" = bayes.mu.mean.bet[6])
  
  table.sd[j, ] <- cbind(zasnova.lm[index.i.bay, ], 
                         "Beta0" = bayes.mu.sd.bet[1], 
                          "Beta1" = bayes.mu.sd.bet[2], 
                          "Beta2" = bayes.mu.sd.bet[3], 
                          "Beta3" = bayes.mu.sd.bet[4], 
                          "Beta4" = bayes.mu.sd.bet[5], 
                          "Beta5" = bayes.mu.sd.bet[6])
  
  skupaj.df[(j+1), ] <- cbind(zasnova.lm[index.i, ], 
                           "Beta0" = bete.hat[1], 
                           "Beta1" = bete.hat[2], 
                          "Beta2" = bete.hat[3], 
                          "Beta3" = bete.hat[4], 
                          "Beta4" = bete.hat[5], 
                          "Beta5" = bete.hat[6])
  
  table.sd[(j+1), ] <- cbind(zasnova.lm[index.i, ], 
                             "Beta0" = sd.bet[1], 
                             "Beta1" = sd.bet[2], 
                             "Beta2" = sd.bet[3], 
                             "Beta3" = sd.bet[4], 
                             "Beta4" = sd.bet[5], 
                             "Beta5" = sd.bet[6])
  j <- j + 2
  i <- i + 1
}
skupaj.df$Metoda <- factor(skupaj.df$V1, labels = metoda)
skupaj.df$Napaka <- factor(skupaj.df$V2,labels = vrsta.napake)
bete.mean <- skupaj.df[,-c(1,2)]
colnames(bete.mean) <- c("Korelacija", "Beta0", "Beta1", "Beta2", "Beta3", "Beta4", "Beta5", "Metoda", "Napaka")
table.sd$Metoda <- factor(table.sd$V1, labels = metoda)
table.sd$Napaka <- factor(table.sd$V2,labels = vrsta.napake)
bete.sd <- table.sd[,-c(1,2)]
colnames(bete.sd) <- c("Korelacija", "SD_Beta0", "SD_Beta1", "SD_Beta2", "SD_Beta3", "SD_Beta4", "SD_Beta5", "Metoda", "Napaka")


end.time <- Sys.time()
end.time - start.time

saveRDS(bete.mean, "data/simulation_I_bete_mean.RDS")
saveRDS(bete.sd, "data/simulation_I_bete_sd.RDS")


library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)

bete.mean.est <- readRDS("data/simulation_I_bete_mean.RDS")
bete.sd.est <-readRDS("data/simulation_I_bete_sd.RDS")

mean.square.error <- function(beta.original, beta.est){
  return(sum((beta.original - beta.est)^2) / length(beta.est))
}

grupiraj.podatke <- bete.mean.est %>% group_by(Korelacija, Metoda, Napaka)%>% 
  summarise_all(mean) %>% 
  gather(key = "Bete", "Vrednost", -c(Korelacija, Metoda, Napaka))

original.beta.df <- data.frame("Bete"= c("Beta0", "Beta1" ,"Beta2","Beta3", "Beta4", "Beta5"),
                            "Orig.vred" = c(-4.5, 0.8, 0.6, 3, 0, 1.4))
ggplot(grupiraj.podatke.gather) + 
  geom_bar( aes(x = Bete, y = Vrednost, color = Metoda, fill = Metoda ), position = "dodge2", stat = "identity")+
  geom_point(data = original.beta.df, aes(x = Bete, y = Orig.vred), shape = 95, size = 3)+
  facet_grid(cols = vars(Napaka), rows = vars(Korelacija))

bias.data <- grupiraj.podatke %>% mutate_each(funs(. - original.beta$.),
                                                starts_with("Beta")) %>%
  gather(key = "Bete", "Vrednost", -c(Korelacija, Metoda, Napaka))


original.beta <- list("Beta0" = -4.5, "Beta1"= 0.8 ,"Beta2" = 0.6,
                      "Beta3" = 3, "Beta4" = 0, "Beta5" = 1.4)

ggplot(bias.data) + 
  geom_bar( aes(x = Bete, y = Vrednost, color = Metoda, fill = Metoda ), position = "dodge2", stat = "identity")+
  facet_grid(cols = vars(Napaka), rows = vars(Korelacija))+
  ggtitle("Pristranskost ocenjenih parametrov", subtitle = "Število simulacij: 1000")
  
#Standardna napaka simulacije
grupiraj.podatke.sd <- bete.mean.est %>% 
  group_by(Korelacija, Metoda, Napaka) %>%
  summarise_all(sd) %>%  
  gather(key = "Bete", "SE", -c(Korelacija, Metoda, Napaka))

grupiraj.podatke.mean <- bete.mean.est %>% 
  group_by(Korelacija, Metoda, Napaka) %>%
  summarise_all(mean)%>%  
  gather(key = "Bete", "Mean", -c(Korelacija, Metoda, Napaka))

podatki.se.mean.join <- grupiraj.podatke.mean %>% left_join(grupiraj.podatke.sd,by = c("Korelacija", "Metoda", "Napaka", "Bete"))


ggplot(data = podatki.se.mean.join) +
  geom_point(aes(x = Bete, y = Mean, color = Metoda ), position = position_dodge(width = 0.3), stat = "identity")+
  geom_errorbar(aes(x = Bete, ymin= Mean-1.96*SE, ymax = Mean+ 1.96*SE,color = Metoda), position = position_dodge(width = 0.3), width = 0.3)+
  facet_grid(cols = vars(Napaka), rows = vars(Korelacija))+
  ggtitle("Intervali zaupanja ocenjenih parametrov", subtitle = "Število simulacij: 1000")

bete.srednje <- bete.mean.est %>% group_by(Korelacija, Metoda, Napaka)%>% 
  summarise(Beta0 = sqrt(mean.square.error(original.beta[[1]], Beta0)),
            Beta1 = sqrt(mean.square.error(original.beta[[2]], Beta1)),
            Beta2 = sqrt(mean.square.error(original.beta[[3]], Beta2)),
            Beta3 = sqrt(mean.square.error(original.beta[[4]], Beta3)),
            Beta4 = sqrt(mean.square.error(original.beta[[5]], Beta4)),
            Beta5 = sqrt(mean.square.error(original.beta[[6]], Beta5))) %>%
  gather(key = "Bete", "SrednjeVrednosti", -c(Korelacija, Metoda, Napaka))


ggplot(data = bete.srednje) +
  geom_point(aes(x = Bete, y = SrednjeVrednosti, color = Metoda ), position = position_dodge(width = 0.3), stat = "identity")+
  facet_grid(cols = vars(Napaka), rows = vars(Korelacija))+
  ggtitle("Koren srednje vrednosti ocenjenih parametrov", subtitle = "Število simulacij: 1000")



grupiraj.podatke.sd <- bete.sd.est %>% group_by(Korelacija, Metoda, Napaka)%>% 
  summarise_all(mean) %>% 
  gather(key = "Bete", "Vrednost", -c(Korelacija, Metoda, Napaka))


ggplot(data = grupiraj.podatke.sd) +
  geom_point(aes(x = Bete, y = Vrednost, color = Metoda ), position = position_dodge(width = 0.3), stat = "identity")+
  facet_grid(cols = vars(Napaka), rows = vars(Korelacija))+
  scale_x_discrete(labels = c("Beta0" , "Beta1" ,"Beta2" ,
                              "Beta3" , "Beta4", "Beta5"))+
  ggtitle("Povprečja standardne napake parametrov", subtitle = "Število simulacij: 1000")
