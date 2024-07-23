# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SupFigures.R         
# - DATE:        15/09/2020
# - DESCRIPTION: Code to produce the figures.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(ggplot2)
library(broom.mixed)
library(tidybayes)
library(tidyverse)
library(dplyr)
library(rstan)
library(rstanarm)
library(brms)
library(bayesplot)
library(ggdist)
library(cowplot)
library(tidyr)
library(patchwork)
library(MetBrewer)
library(modelr)

# Set default ggplot theme

theme_set(theme_minimal()+
            theme(axis.title.x = element_text(size=12,
                                              margin = margin(t = 10, r = 0, b = 0, l = 0)), 
                  axis.title.y = element_text(size=12,
                                              margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.line.x = element_line(color="black", linewidth = 0.5),
                  axis.line.y = element_line(color="black", linewidth = 0.5),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(color="black", size = 12),
                  axis.text.y = element_text(color="black", size = 12),
                  strip.text.x = element_text(size = 12),
                  axis.ticks = element_line(color="black"),
                  plot.title = element_text(hjust = 0.5)))

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load data 
load(paste0(ResultPath, "/models.RData"))
load(paste0(DataPath,"/FinalData.RData"))

# Table S1 ---------------------------------------------------------------------

(TableS1 <- as.data.frame(tidy(m1,
                               effects="fixed", 
                               robust = TRUE,
                               conf.level=.95)) %>% 
   mutate(Parameter = gsub("b_", "", term),
          Class = gsub("Class", "", Parameter),
          Class =gsub(":", "", Class),
          Class =gsub("System", "", Class),
          Class =gsub("Freshwater", "", Class),
          Class =gsub("Marine", "", Class),
          Class =gsub("Terrestrial", "", Class),
          Class = gsub('[[:digit:]]+', '', Class),
          Class = gsub('Bonyfishes', 'Bony fishes', Class),
          Class = gsub("D", "", Class), 
          System = gsub("System", "", Parameter), 
          System = gsub("Class", "", System), 
          System = gsub("Mammalia", "", System), 
          System = gsub("Aves", "", System), 
          System = gsub("Bonyfishes", "", System), 
          System = gsub("Reptilia", "", System), 
          System = gsub("Elasmobranchii", "", System), 
          System = gsub("Amphibia", "", System), 
          System = gsub("D", "/", System), 
          System =gsub("Class", "", System),
          System =gsub(":", "", System),
          System = gsub('[[:digit:]]+', '', System), 
          Rhat=summary(m1)$fixed$Rhat)  %>%
   rename(Median = estimate, 
          CI_low=conf.low, 
          CI_high=conf.high) %>% 
   dplyr::select(Class, System, Median, CI_low, CI_high, Rhat) %>% 
   naniar::replace_with_na_all(condition= ~.x == "") %>% 
   drop_na())

# Save it 

setwd(ResultPath)
write.csv2(TableS1, "Table S1.csv")

# Figure S5: Length ------------------------------------------------------------
# Calculate mean trend

mean <- pops_data %>%
  data_grid(length=seq_range(length, n=101),
            sigma=mean(sigma,na.rm = T)) %>%
  distinct(ID, .keep_all = T) %>% 
  add_epred_draws(mle, ndraws = 250,
                  re_formula = NA) %>% 
  group_by(length) %>% 
  summarise(mean=mean(.epred))

# Calculate individual draws

(figs5 <- pops_data %>%
    data_grid(length=seq_range(length, n=101),
              sigma=mean(sigma,na.rm = T)) %>%
    distinct(ID, .keep_all = T) %>% 
    add_epred_draws(mle, ndraws = 250,
                    re_formula = NA) %>%
    ggplot(aes(x = length, y = .epred)) +
    geom_hline(yintercept = 0, linetype="dashed", 
               colour="grey25")+
    geom_line(aes(group=.draw), colour="#A6323B", alpha = .2)+
    geom_line(data = mean, aes(y=mean), linewidth=1)+
    labs(x="Length of the time series",
         y=expression(paste("Population trend (", mu, ")"))))

# Save the figure

ggsave("Figure S5.pdf", figs5,
       width = 6, height = 4, path = ResultPath)

# Table S2 ---------------------------------------------------------------------

(TableS2 <- as.data.frame(tidy(mle,
                               effects="fixed", 
                               robust = TRUE,
                               conf.level=.95))  %>% 
   mutate(Parameter = gsub("b_", "", term),
          Parameter = gsub("length", "Length of the time series", Parameter),
          Rhat=summary(mle)$fixed$Rhat)  %>%
   rename(Median = estimate, 
          CI_low=conf.low, 
          CI_high=conf.high) %>% 
   naniar::replace_with_na_all(condition= ~.x == "") %>%
   drop_na() %>% 
   dplyr::select(Parameter, Median, CI_low, CI_high, Rhat))
   
# Save it 

setwd(ResultPath)
write.csv2(TableS2, "Table S2.csv")

# Figure S6: Decade --------------------------------------------------------

# Plot

(figs6 <- mdec %>% 
   gather_draws(`b_.*`, regex = TRUE) %>%
   median_qi(.value, .width = c(.95, .8, .5)) %>% 
   mutate(.variable = gsub("b_", "", .variable),
          .variable = gsub("decade", "", .variable))  %>%
   naniar::replace_with_na_all(condition= ~.x == "") %>% 
   drop_na() %>% 
   filter(.variable!="Intercept") %>% 
   ggplot(aes(y = .variable, x = .value)) +
   geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
   geom_pointinterval(aes(xmin = .lower, 
                          xmax = .upper),
                      interval_size_range = c(0.5, 2)) +
   scale_x_continuous(
     labels = scales::number_format(accuracy = 0.01))+
   labs(x="Posterior estimate", y = "") +
   coord_flip()+
   theme(legend.position = "none"))

# Save the figure

ggsave("Figure S6.pdf", figs6,
       width = 6, height = 4, path = ResultPath)

# Figure S7: Hist years --------------------------------------------------------

(figs7 <- pops_data %>% 
   distinct(ID, .keep_all = T) %>% 
   group_by(decade) %>% 
   summarise(n=n()) %>% 
   mutate(prop=n/sum(n)) %>% 
   ggplot(aes(x=reorder(decade, prop), y=prop))+
   geom_bar(stat="identity", fill="grey90", colour="black")+
   labs(x="Decade", y="Proportion")+
   scale_y_continuous(labels = scales::percent_format(accuracy = 1)))

ggsave("Figure S7.pdf", figs7,
       width = 6, height = 4, path = ResultPath)

# Table S3 ---------------------------------------------------------------------

(TableS3 <- as.data.frame(tidy(mdec,
                               effects="fixed", 
                               robust = TRUE,
                               conf.level=.95))  %>% 
   mutate(Parameter = gsub("b_", "", term),
          Parameter = gsub("decade", "", Parameter),
          Rhat=summary(mdec)$fixed$Rhat)  %>%
   rename(Median = estimate, 
          CI_low=conf.low, 
          CI_high=conf.high) %>% 
   naniar::replace_with_na_all(condition= ~.x == "") %>%
   drop_na() %>% 
   dplyr::select(Parameter, Median, CI_low, CI_high, Rhat) %>% 
   arrange(desc(Parameter))) 

# Save it 

setwd(ResultPath)
write.csv2(TableS3, "Table S3.csv")

# Figure S8: Protection --------------------------------------------------------
# Sample size

sample_size <- pops_data %>% 
  mutate(Protected = ifelse(is.na(Protected), "Unknown", Protected),
         Protected =gsub("No", "Unprotected", Protected),
         Protected =gsub("Yes", "Protected", Protected)) %>% 
  group_by(Protected) %>% 
  summarise(median=median(mu),
            n=n(),
            max=max(mu))

med <- mp %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Protected = gsub("b_", "", .variable),
         Protected = gsub("ProtectedProtected", "Protected", Protected),
         Protected =gsub("ProtectedUnknown", "Unknown", Protected),
         Protected =gsub("ProtectedUnprotected", "Unprotected", Protected))  %>%
  group_by(Protected) %>% 
  summarise(med=quantile(.value, probs = 0.9)) %>% 
  naniar::replace_with_na_all(condition= ~.x == "") %>% 
  drop_na() 

sample_size<- sample_size %>% 
  left_join(med)

# Plot intervals 

dat <- mp %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Protected = gsub("b_", "", .variable),
         Protected = gsub("ProtectedProtected", "Protected", Protected),
         Protected =gsub("ProtectedUnknown", "Unknown", Protected),
         Protected =gsub("ProtectedUnprotected", "Unprotected", Protected))  %>%
  naniar::replace_with_na_all(condition= ~.x == "") %>% 
  drop_na()

# Plot

(figs8 <- mp %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(Protected = gsub("b_", "", .variable),
           Protected = gsub("ProtectedProtected", "Protected", Protected),
           Protected =gsub("ProtectedUnknown", "Unknown", Protected),
           Protected =gsub("ProtectedUnprotected", "Unprotected", Protected))  %>%
    naniar::replace_with_na_all(condition= ~.x == "") %>% 
    drop_na() %>% 
    ggplot(aes(y = Protected, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    geom_pointinterval(aes(xmin = .lower, 
                           xmax = .upper,
                           colour=Protected),
                       interval_size_range = c(0.5, 2)) +
    geom_text(data = sample_size, aes(x=med, y=Protected, 
                                      label = paste0("n=", n)),
              vjust   = -1)+
    scale_color_manual("", values = c("#A1D6E2","#336B87", "#334F5D")) +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    xlim(-0.025,0.025)+
    labs(x="Posterior estimate", y = "") + 
    theme(legend.position = "none"))

# Save the figure

ggsave("Figure S8.pdf", figs8,
       width = 6, height = 4, path = ResultPath)

# Table S4 ----------------------------------------------------------------------

(TableS4 <- as.data.frame(tidy(mp,
                               effects="fixed", 
                               robust = TRUE,
                               conf.level=.95))  %>%
   mutate(Parameter = gsub("b_", "", term),
          Parameter = gsub("ProtectedProtected", "Protected", Parameter),
          Parameter =gsub("ProtectedUnknown", "Unknown", Parameter),
          Parameter=gsub("ProtectedUnprotected", "Unprotected", Parameter),
          Rhat=summary(mp)$fixed$Rhat)  %>%
   rename(Median = estimate, 
          CI_low=conf.low, 
          CI_high=conf.high) %>% 
   naniar::replace_with_na_all(condition= ~.x == "") %>%
   drop_na() %>% 
   dplyr::select(Parameter, Median, CI_low, CI_high, Rhat)) 

# Save it 

setwd(ResultPath)
write.csv2(TableS4, "Table S4.csv")

# Table S5 ---------------------------------------------------------------------

(TableS5 <- as.data.frame(tidy(md,
                               effects="fixed", 
                               robust = TRUE,
                               conf.level=.95)) %>%
   mutate(Model="Database",
          Rhat=summary(md)$fixed$Rhat) %>% 
   rbind(as.data.frame(tidy(mc2,
                            effects="fixed", 
                            robust = TRUE,
                            conf.level=.95)) %>% 
           mutate(Model="IUCN",
                  Rhat=summary(mc2)$fixed$Rhat)) %>% 
   rbind(as.data.frame(tidy(ml,
                            effects="fixed", 
                            robust = TRUE,
                            conf.level=.95)) %>% 
           mutate(Model="Data source",
                  Rhat=summary(ml)$fixed$Rhat)) %>% 
   mutate(Parameter = gsub("b_threat", "", term),
          Parameter = gsub("b_", "", Parameter),
          Parameter = gsub("database", "", Parameter),
          Parameter = gsub("Data_source_type", "", Parameter),
          Parameter = gsub("IV", "IbeV", Parameter),
          Parameter = gsub("M", "-", Parameter),
          Parameter = gsub("Grey", "Grey ", Parameter))  %>%
   rename(Median = estimate, 
          CI_low=conf.low, 
          CI_high=conf.high) %>% 
   naniar::replace_with_na_all(condition= ~.x == "") %>%
   drop_na() %>% 
   dplyr::select(Model, Parameter, Median, CI_low, CI_high, Rhat)) 

# Save it 

setwd(ResultPath)
write.csv2(TableS5, "Table S5.csv")

# Table S6 ---------------------------------------------------------------------

(TableS6 <- as.data.frame(tidy(mdi,
                               effects="fixed", 
                               robust = TRUE,
                               conf.level=.95)) %>%  
   mutate(Parameter = gsub("b_", "", term),
          Parameter = gsub("database", "", Parameter),
          Parameter = gsub("Data_source_type", "", Parameter),
          Parameter = gsub("IV", "IbeV", Parameter),
          Parameter = gsub("System", "", Parameter),
          Parameter = gsub("MarineD", "Marine/", Parameter),
          Parameter = gsub("FreshwaterD", "Freshwater/", Parameter),
          Parameter = gsub("Class", "", Parameter),
          Rhat=summary(mdi)$fixed$Rhat)  %>%
   rename(Median = estimate, 
          CI_low=conf.low, 
          CI_high=conf.high) %>%  
   dplyr::select(Parameter, Median, CI_low, CI_high, Rhat) %>% 
   naniar::replace_with_na_all(condition= ~.x == "") %>% 
   drop_na())

# Save it 

setwd(ResultPath)
write.csv2(TableS6, "Table S6.csv")
