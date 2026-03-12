###Title: Theory plots.R

###description: code for plotting figures 1, 2, 3 and 4. 



rm(list=ls())


###packages required####
require(squidSim)
require(MCMCglmm)
require(tidyverse)
require(ggpubr)
require(lme4)
require(ggtext)
library(broom)
require(viridis)


####deocmposed gradient plots (Figures 1,3 and 4)#######


###functions to calculate total selection gradients from decomposed graidents (using equation 7)

calc_sg<-function(Bu, Bc, Be, Vu, Vc, Ve){ ###when there is no genetic information
  
  Bz <-(Bu*(Vu/ (Vu + Ve+ Vc))+Bc*(Vc/ (Vu + Ve+ Vc))+ Be*(Ve/ (Vu + Ve+ Vc))  )
  
  a<- data.frame(`selection gradient` = c("Bz", "Bu","Bc", "Be"),
                 estimate = c(Bz, Bu, Bc, Be))
}


calc_sg_gen<-function(Ba, Bpe, Bc, Be, Va, Vpe,Vc, Ve){ ###when there is genetic information
  
  Bz <-(Ba*(Va/ (Va +Vpe + Ve+ Vc))+Bpe*(Vpe/ (Va+Vpe + Ve+ Vc))+Bc*(Vc/ (Va+Vpe + Ve+ Vc))+ Be*(Ve/ (Va+Vpe + Ve+ Vc)) )
  
  
  
  a<- data.frame(`selection gradient` = c("Bz","Ba", "Bpe","Bc", "Be"),
                 estimate = c(Bz, Ba, Bpe, Bc, Be))
}



groupname1 <-c(expression(beta['z']),
               expression(beta['u']),
               expression(beta['c']),
               expression(beta['e']))





groupname2 <-c(expression(beta['z']),
               expression(beta['a']),
               expression(beta['pe']),
               expression(beta['c']),
               expression(beta['e']))


causal_hard <-calc_sg(0.5, 0.5, 0.5, 1, 1, 1) 

confnd <- calc_sg(0, 0.5, 0, 1, 1, 1)


causal_soft <- calc_sg(0.5, 0, 0.5, 1, 1, 1)


###simple causal vs confnd###
plot1a<-causal_hard %>%
  mutate(selection.gradient = factor(selection.gradient, levels = c("Bz", "Bu", "Bc", "Be")))%>%
  ggplot(aes(y = estimate, x = selection.gradient, colour = selection.gradient, size = selection.gradient))+ 
  geom_point()+ylim(-0.1,1)+theme_classic()+ 
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 1.5, size =0.2)+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12))+labs(y= "Estimate",
                                              subtitle = "A",
                                              x = "Selection Gradient",
                                              colour = NULL,
                                              size = NULL)+
  scale_color_viridis(discrete = TRUE, labels = groupname1)+theme(legend.text.align = 0.5)+
  scale_size_manual(values=c(10,5,5,5), labels = groupname1)+ scale_x_discrete(labels = c('Bz' = expression(beta['z']),
                                                                                          'Bu'   = expression(beta['u']),
                                                                                          'Bc'   = expression(beta['c']),
                                                                                          'Be'   = expression(beta['e'])))+theme(legend.position = "none")

plot1a

plot1b<-confnd %>%
  mutate(selection.gradient = factor(selection.gradient, levels = c("Bz", "Bu", "Bc", "Be")))%>%
  ggplot(aes(y = estimate, x = selection.gradient, colour = selection.gradient, size = selection.gradient))+  
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point()+ylim(-0.1,1)+theme_classic()+
  geom_vline(xintercept = 1.5,size =0.2)+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12))+labs(y= "Estimate",
                                              subtitle = "B",
                                              x = "Selection Gradient",
                                              colour = NULL,
                                              size = NULL)+
  scale_color_viridis(discrete = TRUE, labels = groupname1)+theme(legend.text.align = 0.5)+
  scale_size_manual(values=c(10,5,5,5), labels = groupname1)+ scale_x_discrete(labels = c('Bz' = expression(beta['z']),
                                                                                          'Bu'   = expression(beta['u']),
                                                                                          'Bc'   = expression(beta['c']),
                                                                                          'Be'   = expression(beta['e'])))+theme(legend.position = "none")
plot1b




plot1<-ggarrange(plot1a,plot1b, ncol = 2, common.legend = T, align = "hv", legend = "none" )

plot1

ggsave(filename= "Figure 1A-B.tiff", dpi = 600, path = "../results/", plot = plot1, width=8, height= 4, limitsize = F, bg="white")



#####soft vs hard selection



plot3a<-causal_hard %>%
  mutate(selection.gradient = factor(selection.gradient, levels = c("Bz", "Bu", "Bc", "Be")))%>%
  ggplot(aes(y = estimate, x = selection.gradient, colour = selection.gradient, size = selection.gradient))+ 
  geom_point()+ylim(-0.1,1)+theme_classic()+ 
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 1.5, size =0.2)+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12))+labs(y= "Estimate",
                                              subtitle = "A",
                                              x = "Selection Gradient",
                                              colour = NULL,
                                              size = NULL)+
  scale_color_viridis(discrete = TRUE, labels = groupname1)+theme(legend.text.align = 0.5)+
  scale_size_manual(values=c(10,5,5,5), labels = groupname1)+ scale_x_discrete(labels = c('Bz' = expression(beta['z']),
                                                                                          'Bu'   = expression(beta['u']),
                                                                                          'Bc'   = expression(beta['c']),
                                                                                          'Be'   = expression(beta['e'])))

plot3a




plot3b<-causal_soft %>%
  mutate(selection.gradient = factor(selection.gradient, levels = c("Bz", "Bu", "Bc", "Be")))%>%
  ggplot(aes(y = estimate, x = selection.gradient, colour = selection.gradient, size = selection.gradient))+ 
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_point()+ylim(-0.1,1)+theme_classic()+
  geom_vline(xintercept = 1.5,size = 0.2)+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12))+labs(subtitle = "B",
                                              y= "Estimate",
                                              x = "Selection Gradient",
                                              colour = NULL,
                                              size = NULL)+
  scale_color_viridis(discrete = TRUE, labels = groupname1)+theme(legend.text.align = 0.5)+
  scale_size_manual(values=c(10,5,5,5), labels = groupname1)+ scale_x_discrete(labels = c('Bz' = expression(beta['z']),
                                                                                          'Bu'   = expression(beta['u']),
                                                                                          'Bc'   = expression(beta['c']),
                                                                                          'Be'   = expression(beta['e'])))
plot3b






plot3<-ggarrange(plot3a,plot3b, ncol = 2, common.legend = T, align = "hv", legend = "none" )

plot3

ggsave(filename= "Figure 3A-B.tiff", dpi = 600, path = "../results/", plot = plot3, width=8, height= 4, limitsize = F, bg="white")






#####causality with genetic information available###

casual_gen<-calc_sg_gen(0.5, 0.5, 0.5,0.5, 1 ,1,1,1)

confnd_gen<-calc_sg_gen(0.5, 0, 0,0, 1 ,1,1,1)

confnd_perm<-calc_sg_gen(0, 0.5, 0, 0, 1 ,1,1,1)

causal_nogen<-calc_sg_gen(0, 0.5, 0.5, 0.5, 1 ,1,1,1)



plot4a<-casual_gen %>%
  mutate(selection.gradient = factor(selection.gradient, levels = c("Bz","Ba","Bpe", "Bc", "Be")))%>%
  ggplot(aes(y = estimate, x = selection.gradient, colour = selection.gradient, size = selection.gradient))+ 
  geom_hline(yintercept = 0, linetype = "dashed")+ 
  geom_point()+ylim(-0.1,1)+theme_classic()+ 
  geom_vline(xintercept = 1.5,size = 0.2)+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12))+labs(subtitle = "A",
                                              y = "Estimate",
                                              x = "Selection Gradient",
                                              colour = NULL,
                                              size = NULL)+
  scale_color_viridis(discrete = TRUE, labels = groupname2)+theme(legend.text.align = 0.5)+
  scale_size_manual(values=c(10,5,5,5,5), labels = groupname2)+ scale_x_discrete(labels = c('Bz' = expression(beta['z']),
                                                                                            'Ba'   = expression(beta['a']),
                                                                                            'Bpe'   = expression(beta['pe']),
                                                                                            'Bc'   = expression(beta['c']),
                                                                                            'Be'   = expression(beta['e'])))



plot4a


plot4b<-confnd_perm %>%
  mutate(selection.gradient = factor(selection.gradient, levels = c("Bz","Ba","Bpe", "Bc", "Be")))%>%
  ggplot(aes(y = estimate, x = selection.gradient, colour = selection.gradient, size = selection.gradient))+ 
  geom_hline(yintercept = 0, linetype = "dashed")+ 
  geom_point()+ylim(-0.1,1)+theme_classic()+ 
  geom_vline(xintercept = 1.5, size = 0.2)+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12))+labs(y= "Estimate",
                                              subtitle = "B",
                                              x = "Selection Gradient",
                                              colour = NULL,
                                              size = NULL)+
  scale_color_viridis(discrete = TRUE, labels = groupname2)+theme(legend.text.align = 0.5)+
  scale_size_manual(values=c(10,5,5,5,5), labels = groupname2)+ scale_x_discrete(labels = c('Bz' = expression(beta['z']),
                                                                                            'Ba'   = expression(beta['a']),
                                                                                            'Bpe'   = expression(beta['pe']),
                                                                                            'Bc'   = expression(beta['c']),
                                                                                            'Be'   = expression(beta['e'])))

plot4b


plot4c<-confnd_gen %>%
  mutate(selection.gradient = factor(selection.gradient, levels = c("Bz","Ba","Bpe", "Bc", "Be")))%>%
  ggplot(aes(y = estimate, x = selection.gradient, colour = selection.gradient, size = selection.gradient))+ 
  geom_hline(yintercept = 0, linetype = "dashed")+ 
  geom_point()+ylim(-0.1,1)+theme_classic()+ 
  geom_vline(xintercept = 1.5, size = 0.2)+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12))+labs(subtitle = "C",
                                              y = "Estimate",
                                              x = "Selection Gradient",
                                              colour = NULL,
                                              size = NULL)+
  scale_color_viridis(discrete = TRUE, labels = groupname2)+theme(legend.text.align = 0.5)+
  scale_size_manual(values=c(10,5,5,5,5), labels = groupname2)+ scale_x_discrete(labels = c('Bz' = expression(beta['z']),
                                                                                            'Ba'   = expression(beta['a']),
                                                                                            'Bpe'   = expression(beta['pe']),
                                                                                            'Bc'   = expression(beta['c']),
                                                                                            'Be'   = expression(beta['e'])))

plot4c


plot4d<-causal_nogen %>%
  mutate(selection.gradient = factor(selection.gradient, levels = c("Bz","Ba","Bpe", "Bc", "Be")))%>%
  ggplot(aes(y = estimate, x = selection.gradient, colour = selection.gradient, size = selection.gradient))+ 
  geom_hline(yintercept = 0, linetype = "dashed")+ 
  geom_point()+ylim(-0.1,1)+theme_classic()+ 
  geom_vline(xintercept = 1.5, size = 0.2)+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12))+labs(y= "Estimate",
                                              subtitle = "D",
                                              x = "Selection gradient",
                                              colour = NULL,
                                              size = NULL)+
  scale_color_viridis(discrete = TRUE, labels = groupname2)+theme(legend.text.align = 0.5)+
  scale_size_manual(values=c(10,5,5,5,5), labels = groupname2) + scale_x_discrete(labels = c('Bz' = expression(beta['z']),
                                                                                             'Ba'   = expression(beta['a']),
                                                                                             'Bpe'   = expression(beta['pe']),
                                                                                             'Bc'   = expression(beta['c']),
                                                                                             'Be'   = expression(beta['e'])))

plot4d



plot4<-ggarrange(plot4a,plot4b, plot4c, plot4d, ncol = 2, nrow = 2, common.legend = T, align = "hv", legend = "none" )

ggsave(filename= "Figure 4A-D.tiff", dpi = 600, path = "../results/", plot = plot4, width=10, height= 8, limitsize = F, bg = "white")




###### selection gradients Figure 2#####


soft_selection_data<- simulate_population(n = 500,
                                          response_names = "W",
                                          parameter = list(
                                            observation = list(
                                              names = c("year1", "year2"),
                                              mean = c(2.5, 7.5),
                                              beta = c(0.7,0.7)
                                            ),
                                            residual = list (
                                              vcov = 0.25)
                                          ))


###soft selection##
df<- get_population_data(soft_selection_data)%>% pivot_longer(cols= year1:year2, names_to = "year", values_to = "z")%>% mutate(w = W - mean(W))



year1<- df %>% filter(year== "year1")
year2<- df %>% filter(year == "year2")



model1<- lm(w ~ z, data = year1)
summary(model1)

model2<- lm(w ~ z, data = year2)
summary(model2)


model3<- lm(w ~ z, data = df)
summary(model3)


year1_1<- predict(model1, interval = "confidence")%>% cbind(., year1)

year2_1<- predict(model2, interval = "confidence")%>% cbind(., year2)

df_1<- predict(model3, interval = "confidence")%>% cbind(., df)




###hard selection###


hard_selection_data1<- simulate_population(n = 500,
                                          response_names = "W",
                                          parameter = list(
                                            observation = list(
                                              names = c("z"),
                                              mean = c(2.5),
                                              beta = c(0.7)
                                            ),
                                            residual = list (
                                              vcov = 0.25)
                                          ))



hard_selection_data2<- simulate_population(n = 500,
                                           response_names = "W",
                                           parameter = list(
                                             observation = list(
                                               names = c("z"),
                                               mean = c(7.5),
                                               beta = c(0.7)
                                             ),
                                             residual = list (
                                               vcov = 0.25)
                                           ))



df1<- get_population_data(hard_selection_data1)%>% mutate(year = "year1")
df2<- get_population_data(hard_selection_data2)%>% mutate(year = "year2")


df3<- full_join(df1, df2)%>% mutate(w = W - mean(W))

df1<-df3 %>% filter(year == "year1")
df2<-df3 %>% filter(year == "year2")

model1<- lm(w ~ z, data = df1)
summary(model1)

model2<- lm(w ~ z, data = df2)
summary(model2)

df3<- full_join(df1,df2)%>% mutate(w = W - mean(W))


model3<- lm(w ~ z, data = df3)
summary(model3)

year1_2<- predict(model1, interval = "confidence")%>% cbind(., df1)

year2_2<- predict(model2, interval = "confidence")%>% cbind(., df2)

df_2<- predict(model3, interval = "confidence")%>% cbind(., df3)


####plots###


group.colours <- c("blue", "orange")


plot2c<-df_2 %>% mutate(year = case_when(year== "year1" ~ "Year 1",
                                         year == "year2" ~ "Year 2"))%>% ggplot(aes(x= z, y = w))+ geom_point(aes(colour = year, shape = year), alpha = 0.5)+
  scale_colour_manual(values = group.colours) + 
  scale_y_continuous(limits = c(-5, 5))+
  geom_line(aes(y = fit, x = z), stat = "smooth", data = year1_2, size = 1)+ 
  geom_line(aes(y = fit, x = z), stat = "smooth", data = year2_2, size = 1)+ 
  geom_line(aes(y = fit), stat = "smooth", size = 1, linetype = "dashed")+
  theme_classic()+labs(y=expression('Relative Fitness '~italic(w)~''),
                       x=expression('Trait '~italic(z)~''))+ 
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12))+ labs(subtitle = "A")+ scale_x_continuous(breaks=seq(0,14,by=2))+theme(legend.title=element_blank())
plot2c



plot1c<-df_1 %>% mutate(year = case_when(year== "year1" ~ "Year 1",
                                         year == "year2" ~ "Year 2"))%>% ggplot(aes(x= z, y = w))+ geom_point(aes(colour = year, shape = year), alpha = 0.5)+
  scale_y_continuous(limits = c(-5, 5))+
  scale_colour_manual(values = group.colours) + 
   geom_line(aes(y = fit), stat = "smooth", data = year1_1, size = 1)+ 
  geom_line(aes(y = fit), stat = "smooth", data = year2_1, size = 1)+ 
  geom_line(aes(y = fit), stat = "smooth", size = 1, linetype = "dashed")+
  theme_classic()+labs(y=expression('Relative Fitness '~italic(w)~''),
                       x=expression('Trait '~italic(z)~''))+ 
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12))+ labs(subtitle = "B")+ scale_x_continuous(breaks=seq(0,14,by=2))+theme(legend.title=element_blank())

plot1c

plot1<-ggarrange(plot2c,plot1c,ncol = 2, common.legend = T, align = "hv", legend = "right" )

ggsave(filename= "Figure 2A+B.tiff", dpi = 600, path = "../results/", plot = plot1, width=10, height= 5, limitsize = F, bg = "white")














