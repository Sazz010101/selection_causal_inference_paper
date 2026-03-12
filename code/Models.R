####Title: Models.R
####Description: code for running and analysing selection on antler size in adult male red deer using a univariate model, multi-response model 
#with genetic information and an animal multi-response model (genetic information available).



rm(list=ls())

###packages needed###
require(tidyverse)
require(MCMCglmm)
require(MasterBayes)
require(ggpubr)
require(viridis)

####running the models#####
###load data
data<- read.csv("../data/empricial_data_anon.csv", sep= ",", header = T, stringsAsFactors = F) %>% mutate(Year = as.character(Year),
                                                                                                          Year2 = as.numeric(Year),
                                                                                                          animal = as.character(animal),
                                                                                                          StagCode = as.character(StagCode))%>% 
  na.omit() 



####formatting pedigree######

pedigree<-read.csv("../data/Pedigree_anon.csv", sep= ",", header = T, stringsAsFactors = F)%>% mutate(animal = as.character(animal),
                                                                                                      Sire = as.character(Sire),
                                                                                                      Dam = as.character(Dam))


pedigree1 <- pedigree%>% mutate_all(na_if,"")
pedigree1<-orderPed(pedigree1)
pedigree1<-prunePed(pedigree1, data$animal)
pedigree1<-insertPed(pedigree1, data$animal)
pedigree1 = data.frame(pedigree1)


##Priors##


priorG <- list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

priorGG<- list(R = list(V = diag(2), nu = 2),
               G = list(G1 = list(V = diag(2), nu = (2)),
                        G2 = list(V = diag(2), nu = (2)),
                        G3 = list(V = diag(2), nu = (2))))

priorGG2<- list(R = list(V = diag(2), nu = 2),
               G = list(G1 = list(V = diag(2), nu = (2)),
                        G2 = list(V = diag(2), nu = (2))))



####univariate model
model1_uni<-MCMCglmm(rel_ARS ~ std_Weight + 
                       Centered_Age + 
                       I(Centered_Age^2)+
                       Year2,
                       random = ~StagCode + Year,
                       rcov = ~units,
                       data = data,
                       family = "gaussian",
                       prior = priorG,
                       thin = 200,
                       burnin = 2000,
                       nitt= 251000)
summary(model1_uni)
autocorr.diag(model1_uni$Sol)
autocorr.diag(model1_uni$VCV)
plot(model1_uni$Sol)
plot(model1_uni$VCV)


###multi-response - no genetic info
model1<- MCMCglmm(cbind(rel_ARS, std_Weight)~ trait -1+ 
                      at.level(trait,1):Centered_Age + 
                      at.level(trait,1):I(Centered_Age^2)+ 
                      at.level(trait,1):Year2,
                    random = ~ us(trait):StagCode+ us(trait):Year,
                    rcov = ~ us(trait):units,
                    prior = priorGG2,
                    family = c("gaussian", "gaussian"), 
                    pedigree = pedigree1,
                    data = data,
                    thin = 200,
                    burnin = 2000,
                    nitt= 251000)



summary(model1)
autocorr.diag(model1$Sol)
autocorr.diag(model1$VCV)
plot(model1$Sol)
plot(model1$VCV)

###multi-response - genetic info included###

model2<- MCMCglmm(cbind(rel_ARS, std_Weight)~ trait -1+ 
                    at.level(trait,1):Centered_Age + 
                    at.level(trait,1):I(Centered_Age^2)+ 
                    at.level(trait,1):Year2,
                  random = ~ us(trait):StagCode + us(trait):animal+ us(trait):Year,
                  rcov = ~ us(trait):units,
                  prior = priorGG,
                  family = c("gaussian", "gaussian"), 
                  pedigree = pedigree1,
                  data = data,
                  thin = 200,
                  burnin = 2000,
                  nitt= 251000)



summary(model2)
autocorr.diag(model2$Sol)
autocorr.diag(model2$VCV)
plot(model2$Sol)
plot(model2$VCV)






saveRDS(model1_uni, "../results/Weight_gaussian_uni2")
saveRDS(model1, "../results/Weight_gaussian_multi_nogen2")
saveRDS(model2, "../results/Weight_gaussian_multi2")





model1_uni<-readRDS("../results/Weight_gaussian_uni2")
model1<- readRDS("../results/Weight_gaussian_multi_nogen2")
model2<-readRDS( "../results/Weight_gaussian_multi2")



########calculating total and decomposed selection gradients #####


#####functions####


model_fixed <- function(model) {###pulls selection gradient form univariate model
  sols <- summary(model)$solutions  ## pull out relevant info from model summary
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  ## convert to dataframes with the row.names as the first col
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  ## change the columns names to variable, so they all match
  name <- data.frame(bind_rows(fixed))
} 


antler_weight_no_gen<- function(model1){####produces posterior mean estimates and 95% CIs for the total decomposed selection gradients in the multi-response model with no animal term
model1_selection_coef_all<-{
  mean((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
          model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]+
          model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
         (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
            model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]+
            model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))
}



model1_selection_lwr_all<-{
  HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                 model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]+
                 model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]+
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[1]
}

model1_selection_upr_all<-{
  HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                 model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]+
                 model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]+
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[2]
}





model1_selection_coef_Ind<-{
  mean(model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])
}

model1_selection_lwr_Ind<-{
  HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])[1]
}


model1_selection_upr_Ind<-{
  HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])[2]
}













model1_selection_coef_year<-{
  mean(model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"])
}

model1_selection_lwr_year<-{
  HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"])[1]
}


model1_selection_upr_year<-{
  HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"])[2]
}



model1_selection_coef_units<-{
  mean(model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.units"])
}

model1_selection_lwr_units<-{
  HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.units"])[1]
}


model1_selection_upr_units<-{
  HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.units"])[2]
}






dataframe_no_gen<- data.frame(trait= "antler size",
                                post.mean = c(model1_selection_coef_all,model1_selection_coef_Ind, model1_selection_coef_year, model1_selection_coef_units),
                                l.95..CI = c(model1_selection_lwr_all,model1_selection_lwr_Ind, model1_selection_lwr_year, model1_selection_lwr_units),
                                u.95..CI = c(model1_selection_upr_all, model1_selection_upr_Ind, model1_selection_upr_year, model1_selection_upr_units),
                                Distribution = "Gaussian",
                                variance_component = c("All", "Individual", "Year", "Residual")) %>% mutate(model = "Multi+ no genetic term")



return(dataframe_no_gen)
}


antler_weight_gen<- function(model1){ ####produces posterior mean esitmates and 95% CIs for toal and decomposed selection gradients in the animal multi-response model
  model1_selection_coef_all<-{
    mean((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
            model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]+
            model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]+
            model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
           (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
              model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"]+
              model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]+
              model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))
  }
  
  
  
  model1_selection_lwr_all<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                   model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]+
                   model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]+
                   model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                  (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"]+
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]+
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[1]
  }
  
  model1_selection_upr_all<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                   model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]+
                   model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]+
                   model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                  (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"]+
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]+
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[2]
  }
  
  
  
  
  
  model1_selection_coef_StagCode<-{
    mean(model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])
  }
  
  model1_selection_lwr_StagCode<-{
    HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])[1]
  }
  
  
  model1_selection_upr_StagCode<-{
    HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])[2]
  }
  
  
  model1_selection_coef_animal<-{
    mean(model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"])
  }
  
  model1_selection_lwr_animal<-{
    HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"])[1]
  }
  
  
  model1_selection_upr_animal<-{
    HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"])[2]
  }
  
  
  

  
  
  
  model1_selection_coef_year<-{
    mean(model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"])
  }
  
  model1_selection_lwr_year<-{
    HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"])[1]
  }
  
  
  model1_selection_upr_year<-{
    HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"])[2]
  }
  
  
  
  model1_selection_coef_units<-{
    mean(model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.units"])
  }
  
  model1_selection_lwr_units<-{
    HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.units"])[1]
  }
  
  
  model1_selection_upr_units<-{
    HPDinterval(model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/model1$VCV[,"traitstd_Weight:traitstd_Weight.units"])[2]
  }
  
  
  
  
  
  
  dataframe_genetics<- data.frame(trait= "antler size",
                         post.mean = c(model1_selection_coef_all, model1_selection_coef_StagCode, model1_selection_coef_animal,  model1_selection_coef_year, model1_selection_coef_units),
                         l.95..CI = c(model1_selection_lwr_all, model1_selection_lwr_StagCode, model1_selection_lwr_animal,  model1_selection_lwr_year, model1_selection_lwr_units),
                         u.95..CI = c(model1_selection_upr_all, model1_selection_upr_StagCode, model1_selection_upr_animal, model1_selection_upr_year, model1_selection_upr_units),
                         Distribution = "Gaussian",
                         variance_component = c("All", "PE", "Va","Year", "Residual"))%>% mutate(model = "Multi+ genetic term")
  
  
  
  return(dataframe_genetics)
}






#####making dfs####



antler_uni <- model_fixed(model1_uni)%>% 
  filter(variable == "std_Weight") %>% 
  dplyr::select(post.mean, l.95..CI, u.95..CI )%>% 
  mutate(trait = "antler size",  
         variance_component = "All",
         model = "Lande and Arnold",
         Distribution = "Gaussian")


antler_weight_nogen<- antler_weight_no_gen(model1)


antler_weight_gen<-antler_weight_gen(model2)

antler_size <- rbind(antler_uni, antler_weight_nogen, antler_weight_gen)



write.csv(antler_size, "../results/antler_size_sgrads.csv")






#####Figure 5####


antler_weight<-read.csv("../results/antler_size_sgrads.csv")%>% filter(variance_component != "All" | model != "Multi+ no genetic term") %>% filter(variance_component != "All" |model != "Multi+ genetic term")




antler_weight2<- antler_weight %>% mutate(variance_component=case_when(variance_component == "All" ~ "Bz",
                                                                      variance_component == "Individual" ~ "Bu",
                                                                      variance_component == "PE"~ "Bpe",
                                                                      variance_component == "Va"~ "Ba",
                                                                      variance_component == "Year"~ "Bc",
                                                                      variance_component == "Residual"~ "Be",
                                                                      TRUE ~ variance_component))%>%
  mutate(model=case_when(model == "Lande and Arnold" ~ "Lande-Arnold",
                         model == "Multi+ no genetic term"~ "Multi-Response",
                         model == "Multi+ genetic term" ~ "Animal Multi-Response",
                         TRUE ~ model)) %>%
  mutate(variance_component = factor(variance_component, levels = c("Bz", "Bu", "Ba", "Bpe","Bc", "Be")),
         model = factor(model,levels = c("Lande-Arnold", "Multi-Response", "Animal Multi-Response")))



groupname2 <-c(expression(beta['z']),
               expression(beta['u']),
               expression(beta['a']),
               expression(beta['pe']),
               expression(beta['c']),
               expression(beta['e']))



plot1<- antler_weight2 %>% ggplot(aes(y = post.mean, x = model, colour = variance_component))+ 
  ylim(-0.5, 1.5)+geom_hline(aes(yintercept= 0), linetype = "dashed")+ geom_point(position=position_dodge(width=.5))+
  geom_pointrange(aes(ymin= l.95..CI, ymax= u.95..CI, colour = variance_component, shape = variance_component), position=position_dodge(width=.5))+
  theme_classic()+labs(x = NULL, y = "Estimate", colour = "Selection Gradient", shape = "Selection Gradient")+ 
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12))+
  scale_color_viridis(discrete = TRUE, labels = groupname2)+theme(legend.text.align = 0.5)+theme(legend.position = "none")


plot1

ggsave(filename= "Figure 5.tiff", dpi = 400, path = "../results/", plot = plot1, width= 6, height= 5, limitsize = F)

#####differences between decomposed selection gradients within each multi-response model #####

model1<- readRDS("../results/Weight_gaussian_multi_nogen2")
model2<-readRDS( "../results/Weight_gaussian_multi2")

diff_no_gen<- function(model1){ ###differences between decomposed selection gradients in the multi-response model with no genetic informaiton included
  
  
 
  
  
  diff_Bu_Bc<-{
    mean((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
            model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])-
           (model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
              model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]))
  }
  
  
  diff_Bu_Bc_upr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]))[2]
  }
  
  
  
  diff_Bu_Bc_lwr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]))[1]
  }
  
  
  diff_Bu_Be<-{
    mean((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
            model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])-
           (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
              model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))
  }
  
  
  diff_Bu_Be_upr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[2]
  }
  
  
  
  diff_Bu_Be_lwr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[1]
  }
  
  
  
  
  
  
  diff_Bc_Be<-{
    mean((model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
            model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"])-
           (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
              model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))
  }
  
  
  diff_Bc_Be_upr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[2]
  }
  
  
  
  diff_Bc_Be_lwr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[1]
  }
  
  
  
  dataframe_diff<- data.frame( post.mean = c(diff_Bu_Bc, diff_Bu_Be, diff_Bc_Be),
                               l.95..CI = c(diff_Bu_Bc_lwr, diff_Bu_Be_lwr,  diff_Bc_Be_lwr),
                               u.95..CI = c(diff_Bu_Bc_upr, diff_Bu_Be_upr, diff_Bc_Be_upr),
                               variance_component = c("Bu_Bc", "Bu_Be", "Bc_Be"))
  
  
}

diff_gen<- function(model1){###differences between decomposed selection gradients in the animal multi-response model
  
 
  diff_Ba_Bpe<-{
    mean((model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]/
            model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"])-
           (model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
              model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]))
  }
  
  
  diff_Ba_Bpe_upr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]))[2]
  }
  
  
  
  diff_Ba_Bpe_lwr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]))[1]
  }
  
  
 
  
 
  
  
  diff_Ba_Bc<-{
    mean((model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]/
            model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"])-
           (model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
              model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]))
  }
  
  
  diff_Ba_Bc_upr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]))[2]
  }
  
  
  
  diff_Ba_Bc_lwr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]))[1]
  }
  
  
  diff_Ba_Be<-{
    mean((model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]/
            model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"])-
           (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
              model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))
  }
  
  
  diff_Ba_Be_upr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[2]
  }
  
  
  
  diff_Ba_Be_lwr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[1]
  }
  
 
  
  
  
  
  
  
  diff_Bpe_Bc<-{
    mean((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
            model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])-
           (model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
              model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]))
  }
  
  
  diff_Bpe_Bc_upr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]))[2]
  }
  
  
  
  diff_Bpe_Bc_lwr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]))[1]
  }
  
  
  diff_Bpe_Be<-{
    mean((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
            model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])-
           (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
              model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))
  }
  
  
  diff_Bpe_Be_upr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[2]
  }
  
  
  
  diff_Bpe_Be_lwr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[1]
  }
  
   
  
  
  
  
  diff_Bc_Be<-{
    mean((model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
            model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"])-
           (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
              model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))
  }
  
  
  diff_Bc_Be_upr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[2]
  }
  
  
  
  diff_Bc_Be_lwr<-{
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]/
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"])-
                  (model1$VCV[,"traitstd_Weight:traitrel_ARS.units"]/
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[1]
  }
  
  
  
  dataframe_diff<- data.frame( post.mean = c(diff_Ba_Bpe, diff_Ba_Bc, diff_Ba_Be, diff_Bpe_Bc, diff_Bpe_Be, diff_Bc_Be),
                            l.95..CI = c(diff_Ba_Bpe_lwr, diff_Ba_Bc_lwr, diff_Ba_Be_lwr, diff_Bpe_Bc_lwr, diff_Bpe_Be_lwr,  diff_Bc_Be_lwr),
                            u.95..CI = c(diff_Ba_Bpe_upr, diff_Ba_Bc_upr, diff_Ba_Be_upr, diff_Bpe_Bc_upr, diff_Bpe_Be_upr, diff_Bc_Be_upr),
                            variance_component = c("Ba_Bpe","Ba_Bc", "Ba_Be", "Bpe_Bc", "Bpe_Be", "Bc_Be"))
  
  
  }


diff_nogen<-diff_no_gen(model1)


model1<- model2
diff_gen<-diff_gen(model1)

write.csv(diff_nogen, "../results/difference_bwt_sgs_nogen.csv")

write.csv(diff_gen, "../results/difference_bwt_sgs_gen.csv")
#######correcting selection gradient estimates when selection is soft######


antler_uni <- model_fixed(model1_uni)%>% 
  filter(variable == "std_Weight") %>% 
  dplyr::select(post.mean, l.95..CI, u.95..CI )%>% 
  mutate(trait = "antler size",  
         type = "Uncorrected Bz",
         Distribution = "Gaussian")
model1<-readRDS( "../results/Weight_gaussian_multi")


  
  
  corr_mm<- mean((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                            model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]+
                            model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                           (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                              model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"]+
                              model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))
  
  corr_mm_upr<-
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                   model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]+
                   model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                  (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"]+
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[2]
  

  
  
  corr_mm_lwr<-
    HPDinterval((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                   model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]+
                   model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                  (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"]+
                     model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))[1]
  
  

  
dataframe<- data.frame(trait= "antler size",
  post.mean = c(corr_mm),
  l.95..CI = c(corr_mm_lwr),
  u.95..CI = c(corr_mm_upr),
  Distribution = "Gaussian",
  type = c( "Corrected Bz"))



dataframe2<-rbind(antler_uni, dataframe)


dataframe2<- dataframe2 %>% 
  mutate(type = factor(type, levels = c( "Uncorrected Bz", "Corrected Bz")))

groupname2 <-c(expression(beta["z"]),
               expression(beta*{"`"}["z"]))


plot2<- dataframe2%>% ggplot(aes(y = post.mean, x = type, shape = type))+ 
  ylim(-0.5, 1.5)+geom_hline(aes(yintercept= 0), linetype = "dashed")+ geom_point(position=position_dodge(width=.5))+
  geom_pointrange(aes(ymin= l.95..CI, ymax= u.95..CI), position=position_dodge(width=.5))+
  theme_classic()+labs(x = NULL, y = "Estimate", shape = "Selection Gradient")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=10))+
  scale_shape_discrete(labels = groupname2)+theme(legend.text.align = 0.5)+theme(legend.position = "none")+ scale_x_discrete(labels = c('Uncorrected Bz' = expression(beta["z"]),
                                                                    'Corrected Bz' = expression(""*beta*{"`"}["z"]*"")))

plot2

ggsave(filename= "Figure6.tiff", dpi = 400, path = "../results/", plot = plot2, width= 4, height= 4, limitsize = F, bg="white")



#####difference between uncorrected and corrected selection gradients from multi-response model#s###


Lande_Arnold_diff_mm<- mean(((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                          model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]+
                          model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]+
                          model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                         (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                            model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"]+
                            model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]+
                            model1$VCV[,"traitstd_Weight:traitstd_Weight.units"])) -
                           ((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                              model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]+
                              model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                           (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                              model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"]+
                              model1$VCV[,"traitstd_Weight:traitstd_Weight.units"])))




Lande_Arnold_diff_mm_lwr<-
  HPDinterval(((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                 model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]+
                 model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]+
                 model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"]+
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]+
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))-
                ((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                    model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]+
                    model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                   (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                      model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"]+
                      model1$VCV[,"traitstd_Weight:traitstd_Weight.units"])))[1]

Lande_Arnold_diff_mm_upr<-
  HPDinterval(((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                 model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]+
                 model1$VCV[,"traitstd_Weight:traitrel_ARS.Year"]+
                 model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"]+
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.Year"]+
                   model1$VCV[,"traitstd_Weight:traitstd_Weight.units"]))- 
                ((model1$VCV[,"traitstd_Weight:traitrel_ARS.StagCode"]+
                    model1$VCV[,"traitstd_Weight:traitrel_ARS.animal"]+
                    model1$VCV[,"traitstd_Weight:traitrel_ARS.units"])/
                   (model1$VCV[,"traitstd_Weight:traitstd_Weight.StagCode"]+
                      model1$VCV[,"traitstd_Weight:traitstd_Weight.animal"]+
                      model1$VCV[,"traitstd_Weight:traitstd_Weight.units"])))[2]



dataframe3<- data.frame(trait= "antler size",
                       post.mean = c(Lande_Arnold_diff_mm),
                       l.95..CI = c(Lande_Arnold_diff_mm_lwr),
                       u.95..CI = c(Lande_Arnold_diff_mm_upr),
                       type = c( "Difference"))




#####model outputs  to go into the supplementary material####


model1_uni<-readRDS("../results/Weight_gaussian_uni")
model1<- readRDS("../results/Weight_gaussian_multi_nogen")
model2<-readRDS( "../results/Weight_gaussian_multi")



antler_uni <- model_fixed(model1_uni)

antler_nogen <- model_fixed(model1)

antler_gen <- model_fixed(model2)

summary(model1_uni)

summary(model1)

summary(model2)
