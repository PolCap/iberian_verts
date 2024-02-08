# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Figures.R         
# - DATE:        15/09/2020
# - DESCRIPTION: Code to produce the figures.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(ggplot2)
library(tidybayes)
library(tidyverse)
library(dplyr)
library(rstan)
library(rstanarm)
library(brms)
library(bayesplot)
library(bayestestR)
library(ggdist)
library(cowplot)
library(wesanderson)
library(tidyr)
library(patchwork)
library(MetBrewer)
library(ggrepel)
library(fmsb)
library(ggradar)

# Set default ggplot theme

theme_set(theme_minimal()+
            theme(axis.title.x = element_text(size=14,
                                              margin = margin(t = 10, r = 0, b = 0, l = 0)), 
                  axis.title.y = element_text(size=14,
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

# Color palette for Classs

class_pal <- wes_palette("Cavalcanti1", n = 6, type = "continuous")

# System palette

system_pal <- c("#A1D6E2", "#6aa1b5", "#336B87","#80806b", "#CC954E",
                "#b7b698", "#8b9d92")

database_pal <- c("#9A6576","#659A89")

# Data description #############################################################

pops_data %>% 
  mutate(trend=ifelse(mu>0, "Increase",
                      ifelse(mu<0, "Decrease", "Stable"))) %>%
  group_by(trend) %>% 
  summarise(n=n()) %>% 
  mutate(freq = (n / sum(n))*100)

pops_data %>% 
  group_by(database) %>% 
  summarise(n=n()) %>% 
  mutate(freq = (n / sum(n))*100)

pops_data %>% 
  ungroup() %>% 
  group_by(Species) %>%
  summarise(num_sources = n_distinct(database)) %>%
  ungroup() %>% 
  distinct(Species, .keep_all = T) %>% 
  group_by(num_sources) %>% 
  summarise(n=n()) 

pops_data %>% 
  ungroup() %>% 
  filter(database=="LPD") %>% 
  mutate(Species=as.character(as.factor(Species))) %>% 
  distinct(Species, .keep_all = T) %>% 
  nrow()

pops_data %>% 
  ungroup() %>% 
  filter(database=="IbeV") %>% 
  mutate(Species=as.character(as.factor(Species))) %>% 
  distinct(Species, .keep_all = T) %>% 
  nrow()

pops_data %>% 
  filter(database=="IbeV") %>% 
  group_by(Data_source_type) %>% 
  summarise(n=n()) %>% 
  mutate(freq = (n / sum(n))*100)

# Figure 1: Description of the database ########################################
## Proportion of LPD vs IbeV -----------------------------------------------------

(ga <- pops_data %>% 
   distinct(ID, .keep_all=T) %>% 
   mutate(database=ifelse(database=="IbeV", "IbeV", database)) %>% 
   group_by(database) %>% 
   summarise(n=n()) %>% 
   mutate(freq=(n/sum(n)),
          ymax=cumsum(freq),
          ymin=c(0, ymax[1]),
          labelPosition= (ymax+ymin)/2,
          label= paste0(database, "\n", 
                        round(freq*100, 2), "%")) %>%
   ggplot(aes(ymax=ymax, ymin=ymin, 
              xmax=4, xmin=3, 
              fill=database)) +
   geom_rect() +
   scale_color_manual("", values = datab_pal ) +
   scale_fill_manual("", values =  datab_pal ) +
   coord_polar(theta="y", clip = "off") +
   xlim(c(1, 4)) +
   theme_void() +
   geom_text(x=5.1, 
             aes(y=labelPosition, label=label, 
                 color=database), size=6) +
   theme(legend.position = "none",
         plot.margin = unit(c(0,0,0,0), units = "cm")))

## Taxon ------------------------------------------------------------------

# Set database palette

database_pal <- c("#9A6576","#659A89")

# Reshape the data so it fits into the format for the spider plot

(gb <- pops_data %>% 
    distinct(ID, .keep_all=T) %>%  
    group_by(database, Class) %>% 
    summarise(n=n()) %>% 
    mutate(freq=(n/sum(n)),
         Class=ifelse(Class=="`Bony fishes`", "Bony fishes", Class)) %>%
    ungroup() %>% 
    dplyr::select(-n) %>% 
    pivot_wider(names_from = Class,values_from = freq) %>% 
    ggradar(group.colours =   database_pal,
            group.point.size = 2,
            fill=TRUE,
            fill.alpha = .7)+
  theme(legend.position = "none")+
    coord_cartesian(clip="off"))

## System ------------------------------------------------------------------

# Reshape the data so it fits into the format for the spider plot

(gc <- pops_data %>% 
   distinct(ID, .keep_all=T) %>% 
   mutate(System, fct_relevel(System, "Marine", "Marine/Freshwater", 
                              "Freshwater", "Freshwater/Terrestrial",
                              "Terrestrial", "Marine/Terrestrial",
                              "Marine/Freshwater/Terrestrial")) %>% 
   group_by(database, System) %>% 
   summarise(n=n()) %>% 
   mutate(freq=(n/sum(n))) %>%
   ungroup() %>% 
   dplyr::select(-n) %>% 
   pivot_wider(names_from = System,
               values_from = freq) %>% 
   ggradar(group.colours =   database_pal,
           group.point.size = 2,
           fill=TRUE,
           fill.alpha = .7)+
   theme(legend.position = "none")+
   coord_cartesian(clip="off"))

## Data source -----------------------------------------------------------------

(gd <- pops_data %>% 
   distinct(ID, .keep_all=T) %>% 
   group_by(database, Data_source_type) %>% 
   summarise(n=n()) %>% 
   mutate(freq=(n/sum(n))) %>%
   ggplot(aes(x=Data_source_type,
              y=freq, group=database, fill=database)) +
   geom_bar(stat = "identity", position = position_dodge2(width = 0.9, 
                                                          preserve = "single")) +
   scale_fill_manual("", values = database_pal) +
   labs(x="Data source", y="")+ 
   scale_y_continuous(labels = scales::percent,
                      limits = c(0,1))+
   theme(legend.position = "none",
         plot.margin = unit(c(0,0,0,0), units = "cm")))


## Threat -----------------------------------------------------------------

(ge <- pops_data %>% 
   distinct(ID, .keep_all=T) %>% 
   group_by(database, threat) %>% 
   summarise(n=n()) %>% 
   mutate(freq=(n/sum(n))) %>%
   ggplot(aes(x=threat,
              y=freq, group=database, fill=database)) +
   geom_bar(stat = "identity", position = position_dodge2(width = 0.9, 
                                                          preserve = "single")) +
   scale_fill_manual("", values = database_pal) +
   labs(x="Conservation status", y="")+ 
   scale_y_continuous(labels = scales::percent,
                      limits = c(0,1))+
   theme(legend.position = "none",
         plot.margin = unit(c(0,0,0,0), units = "cm")))


# ## Protection ------------------------------------------------------------------
# 
# # LPD
# 
# (gf <- pops_data %>% 
#    distinct(ID, .keep_all=T) %>%
#    mutate(Protected=ifelse(is.na(Protected), 
#                                  "Unknown", 
#                            ifelse(Protected=="Yes", 
#                                   "Protected", "Unprotected")))%>%
#    filter(database=="LPD") %>% 
#    group_by(Protected) %>% 
#    summarise(n=n()) %>% 
#    mutate(freq=(n/sum(n)),
#           ymax=cumsum(freq),
#           ymin=c(0, ymax[1]),
#           labelPosition= (ymax+ymin)/2,
#           label= paste0(Protected, "\n", 
#                         round(freq*100, 0), "%")) %>%
#    ggplot(aes(ymax=ymax, ymin=ymin, 
#               xmax=4, xmin=3, 
#               fill=Protected)) +
#    geom_rect() +
#    scale_fill_manual(values = c("#F5C27A","#9F5B55")) +
#    scale_colour_manual(values = c("#F5C27A","#9F5B55")) +
#    coord_polar(theta="y", clip = "off") +
#    xlim(c(1, 4)) +
#    theme_void() +
#    geom_text(x=5.1,aes(y=labelPosition, label=label, 
#                      color=Protected), size=2.5) +
#    theme(legend.position = "none",
#          plot.margin = unit(c(0,0,0,0), units = "cm")))
# 
# # IbeV
# 
# (gg <- pops_data %>% 
#     distinct(ID, .keep_all=T) %>%
#     mutate(Protected=ifelse(is.na(Protected), 
#                             "Unknown", 
#                             ifelse(Protected=="Yes", 
#                                    "Protected", "Unprotected")))%>%
#     filter(database=="IbeV") %>% 
#     group_by(Protected) %>% 
#     summarise(n=n()) %>% 
#     mutate(freq=(n/sum(n)),
#            ymax=cumsum(freq),
#            ymin=c(0, ymax[-3]),
#            labelPosition= (ymax+ymin)/2,
#            label= paste0(Protected, "\n", 
#                          round(freq*100, 0), "%")) %>%
#     ggplot(aes(ymax=ymax, ymin=ymin, 
#                xmax=4, xmin=3, 
#                fill=Protected)) +
#     geom_rect() +
#     scale_fill_manual(values=c("#F5C27A","#9F9F9F", "#9F5B55")) +
#     scale_colour_manual(values=c("#F5C27A","#9F9F9F", "#9F5B55")) +
#     coord_polar(theta="y", clip = "off") +
#     xlim(c(1, 4)) +
#     theme_void() +
#     geom_text(x=5.1,aes(y=labelPosition, label=label, 
#                       color=Protected), size=2.5) +
#     theme(legend.position = "none",
#           plot.margin = unit(c(0,0,0,0), units = "cm")))

## Combine figures -------------------------------------------------------------

(fig1 <- ga+((gb+gc)/(gd+ge)) +
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 1, widths =c(1,2) ))

# Save them

ggsave("Figure 1.pdf", fig1,
       width = 35, 
       height = 20,
       units = "cm", 
       path = ResultPath)

# Figure 2: Map  ##########################################

## Map -------------------------------------------------------------------------

# Set the world map

world <- map_data("world") #, region = c("spain", "portugal"))

# Create the map

(p1 <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id=region), 
             color = "gray60", fill = "gray80", size = 0.3) +
    geom_point(data = pops_data %>% filter(mu<0.3,mu>-0.3), # I remove extreme cases for visualisation purposes
               aes(x = as.numeric(Longitude), 
                   y = as.numeric(Latitude), 
                   fill=mu), 
               alpha= 0.8, shape=21, size=3)+
    scale_fill_gradient2(expression(paste("Population trend (", mu, ")")),
                         midpoint = 0, 
                         low = "#c1666b", 
                         mid = "white",
                         high = "#4281a4",
                         guide =guide_colourbar(nbin=100,
                                                barwidth = 9,
                                                title.position="top"),
                         breaks=c(-0.3, -0.20, -0.10, 0, 
                                  0.10, 0.2, 0.30),
                         limits=c(-0.3, 0.3))+
    scale_y_continuous(limits = c(35, 45))+
    scale_x_continuous(limits = c(-10,8)) + 
    # facet_wrap(~Class)+
    theme_map()+
    theme(legend.position = c(0.6, 0.1),
          legend.direction = "horizontal",
          legend.title = element_text(size = 12, hjust =0.5),
          legend.text = element_text(size = 10),
          plot.margin = unit(c(0,0,0,0), units = , "cm")))

ggsave("Figure2.pdf", p1,
       width = 12, height = 10,
       units = "cm", 
       path = ResultPath)


# Figure 3: Trends across taxa and systems ###################

# Sample size

sample_size <- pops_data %>% 
  group_by(System, Class) %>% 
  summarise(n=n())

med <- m1 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub("b_", "", .variable),
         Class = gsub("Class", "", .variable),
         Class =gsub(":", "", Class),
         Class =gsub("System", "", Class),
         Class =gsub("Freshwater", "", Class),
         Class =gsub("Marine", "", Class),
         Class =gsub("Terrestrial", "", Class),
         Class =gsub("D", "", Class),
         Class = gsub('[[:digit:]]+', '', Class),
         Class = gsub('Bonyfishes', 'Bony fishes', Class),
         System = gsub("System", "", .variable), 
         System = gsub("D", "/", System), 
         System = gsub("Class", "", System), 
         System = gsub("Mammalia", "", System), 
         System = gsub("Aves", "", System), 
         System = gsub("Bonyfishes", "", System), 
         System = gsub("Reptilia", "", System), 
         System = gsub("Elasmobranchii", "", System), 
         System = gsub("Amphibia", "", System), 
         System =gsub(":", "", System),
         System = gsub('[[:digit:]]+', '', System))  %>%
  group_by(System, Class) %>% 
  summarise(med=quantile(.value, probs = 0.9)) %>% 
  naniar::replace_with_na_all(condition= ~.x == "") %>% 
  drop_na() 


sample_size<- sample_size %>% 
  left_join(med) %>% 
  mutate(System= factor(System, levels=c("Freshwater", "Marine", "Terrestrial",
                                         "Marine/Freshwater", "Marine/Terrestrial", 
                                         "Freshwater/Terrestrial", 
                                         "Marine/Freshwater/Terrestrial")))  


# Plot intervals 

dat <- m1 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub("b_", "", .variable),
         Class = gsub("Class", "", .variable),
         Class =gsub(":", "", Class),
         Class =gsub("System", "", Class),
         Class =gsub("Freshwater", "", Class),
         Class =gsub("Marine", "", Class),
         Class =gsub("Terrestrial", "", Class),
         Class =gsub("D", "", Class),
         Class = gsub('[[:digit:]]+', '', Class),
         Class = gsub('Bonyfishes', 'Bony fishes', Class),
         System = gsub("System", "", .variable), 
         System = gsub("D", "/", System), 
         System = gsub("Class", "", System), 
         System = gsub("Mammalia", "", System), 
         System = gsub("Aves", "", System), 
         System = gsub("Bonyfishes", "", System), 
         System = gsub("Reptilia", "", System), 
         System = gsub("Elasmobranchii", "", System), 
         System = gsub("Amphibia", "", System), 
         System =gsub(":", "", System),
         System = gsub('[[:digit:]]+', '', System))  %>%
  naniar::replace_with_na_all(condition= ~.x == "") %>% 
  drop_na() %>% 
  filter(System!="Marine"|Class!="Reptilia",
         System!="Freshwater"|Class!="Elasmobranchii",
         System!="Marine"|Class!="Amphibia",
         System!="Terrestrial"|Class!="Elasmobranchii",
         System!="Terrestrial"|Class!="Bony fishes",
         System!="Terrestrial"|Class!="Amphibia",
         System!="Freshwater"|Class!="Reptilia",
         System!="Freshwater"|Class!="Mammalia",
         System!="Freshwater"|Class!="Elasmobranchii",
         System!="Freshwater"|Class!="Amphibia",
         System!="Freshwater/Terrestrial"|Class!="Elasmobranchii",
         System!="Freshwater/Terrestrial"|Class!="Bony fishes",
         System!="Marine/Freshwater"|Class!="Reptilia",
         System!="Marine/Freshwater"|Class!="Mammalia",
         System!="Marine/Freshwater"|Class!="Elasmobranchii",
         System!="Marine/Freshwater"|Class!="Aves",
         System!="Marine/Freshwater"|Class!="Amphibia",
         System!="Marine/Freshwater/Terrestrial"|Class!="Reptilia",
         System!="Marine/Freshwater/Terrestrial"|Class!="Bony fishes",
         System!="Marine/Freshwater/Terrestrial"|Class!="Elasmobranchii",
         System!="Marine/Freshwater/Terrestrial"|Class!="Amphibia",
         System!="Marine/Terrestrial"|Class!="Reptilia",
         System!="Marine/Terrestrial"|Class!="Bony fishes",
         System!="Marine/Terrestrial"|Class!="Elasmobranchii",
         System!="Marine/Terrestrial"|Class!="Amphibia",
         System!="Marine/Terrestrial"|Class!="Mammalia") %>% 
  mutate(System= factor(System, levels=c("Freshwater", "Marine", "Terrestrial",
                                         "Marine/Freshwater", "Marine/Terrestrial", 
                                         "Freshwater/Terrestrial", 
                                         "Marine/Freshwater/Terrestrial")))  
  

# Plot

(g1 <- m1 %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(.variable = gsub("b_", "", .variable),
           Class = gsub("Class", "", .variable),
           Class =gsub(":", "", Class),
           Class =gsub("System", "", Class),
           Class =gsub("Freshwater", "", Class),
           Class =gsub("Marine", "", Class),
           Class =gsub("Terrestrial", "", Class),
           Class =gsub("D", "", Class),
           Class = gsub('[[:digit:]]+', '', Class),
           Class = gsub('Bonyfishes', 'Bony fishes', Class),
           System = gsub("System", "", .variable), 
           System = gsub("D", "/", System), 
           System = gsub("Class", "", System), 
           System = gsub("Mammalia", "", System), 
           System = gsub("Aves", "", System), 
           System = gsub("Bonyfishes", "", System), 
           System = gsub("Reptilia", "", System), 
           System = gsub("Elasmobranchii", "", System), 
           System = gsub("Amphibia", "", System), 
           System =gsub(":", "", System),
           System = gsub('[[:digit:]]+', '', System))   %>%
    naniar::replace_with_na_all(condition= ~.x == "") %>% 
    drop_na() %>% 
    filter(System!="Marine"|Class!="Reptilia",
           System!="Freshwater"|Class!="Elasmobranchii",
           System!="Marine"|Class!="Amphibia",
           System!="Terrestrial"|Class!="Elasmobranchii",
           System!="Terrestrial"|Class!="Bony fishes",
           System!="Terrestrial"|Class!="Amphibia",
           System!="Freshwater"|Class!="Reptilia",
           System!="Freshwater"|Class!="Mammalia",
           System!="Freshwater"|Class!="Elasmobrachii",
           System!="Freshwater"|Class!="Amphibia",
           System!="Freshwater/Terrestrial"|Class!="Elasmobranchii",
           System!="Freshwater/Terrestrial"|Class!="Bony fishes",
           System!="Marine/Freshwater"|Class!="Reptilia",
           System!="Marine/Freshwater"|Class!="Mammalia",
           System!="Marine/Freshwater"|Class!="Elasmobranchii",
           System!="Marine/Freshwater"|Class!="Aves",
           System!="Marine/Freshwater"|Class!="Amphibia",
           System!="Marine/Freshwater/Terrestrial"|Class!="Reptilia",
           System!="Marine/Freshwater/Terrestrial"|Class!="Bony fishes",
           System!="Marine/Freshwater/Terrestrial"|Class!="Elasmobranchii",
           System!="Marine/Freshwater/Terrestrial"|Class!="Amphibia",
           System!="Marine/Terrestrial"|Class!="Reptilia",
           System!="Marine/Terrestrial"|Class!="Bony fishes",
           System!="Marine/Terrestrial"|Class!="Elasmobranchii",
           System!="Marine/Terrestrial"|Class!="Amphibia",
           System!="Marine/Terrestrial"|Class!="Mammalia") %>% 
    mutate(System= factor(System, levels=c("Freshwater", "Marine", "Terrestrial",
                                           "Marine/Freshwater", "Marine/Terrestrial", 
                                           "Freshwater/Terrestrial", 
                                           "Marine/Freshwater/Terrestrial"))) %>% 
    ggplot(aes(y = Class, x = .value)) +
    facet_wrap(.~System, scales = "free_x")+
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper,
                           colour=Class),
                       interval_size_range = c(0.5, 2)) +
    geom_text(data = sample_size, aes(x=med, y=Class, 
                                      label = paste0("n=", n)),
              vjust   = -1)+
    scale_colour_manual(values = class_pal)+
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    xlim(-0.12,0.12)+
    labs(x="Posterior estimate", y = "Class") + 
    theme(legend.position = "none"))

# Save the figure

ggsave("Figure 3.pdf", g1,
       width = 12, height = 10, path = ResultPath)

# Figure 4: Ancillary ####################################################

## Panel a: LPD vs Iberian Vertebrates -----------------------------------------

# Sample size

sample_size <- pops_data %>% 
  drop_na(database) %>%
  group_by(database) %>% 
  summarise(median=median(mu),
            n=n(),
            max=max(mu))

med <- md %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(database = gsub("b_database", "", .variable),
         database = gsub("IV", "IbeV", database))  %>%
  group_by(database) %>% 
  summarise(med=quantile(.value, probs = 0.9)) %>% 
  naniar::replace_with_na_all(condition= ~.x == "") %>% 
  drop_na() 

sample_size<- sample_size %>% 
  left_join(med)

# Plot intervals 

dat <- md %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(database = gsub("b_database", "", .variable),
         database = gsub("IV", "IbeV", database))  %>% 
  naniar::replace_with_na_all(condition= ~.x == "") %>% 
  drop_na()

# Plot

(fig4a <- md %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(database = gsub("b_database", "", .variable),
           database = gsub("IV", "IbeV", database))  %>%
    naniar::replace_with_na_all(condition= ~.x == "") %>% 
    drop_na() %>% 
    ggplot(aes(y = database, x = .value)) +
    geom_vline(xintercept = 0, 
               linetype = "dashed", 
               colour="grey50") +
    geom_pointinterval(aes(xmin = .lower, 
                           xmax = .upper,
                           colour=database),
                       interval_size_range = c(0.5, 2)) +
    geom_text(data = sample_size, aes(x=med, y=database, 
                                      label = paste0("n=", n)),
              vjust   = -1)+
    scale_color_manual("", values = datab_pal) +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    xlim(-0.04,0.04)+
    labs(x="Posterior estimate", y = "Database") + 
    theme(legend.position = "none"))

## Panel b: IUCN ---------------------------------------------------------------

# Create a summary of the sample size

sample_size <- pops_data %>%
  drop_na(threat) %>%
  group_by(threat) %>%
  summarise(median=median(mu),
            n=n(),
            max=max(mu))

med <- mc2 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(threat = gsub("b_threat", "", .variable),
         threat = gsub("NonMth", "Non-th", threat))  %>%
  group_by(threat) %>% 
  summarise(med=quantile(.value, probs = 0.9)) %>% 
  naniar::replace_with_na_all(condition= ~.x == "") %>% 
  drop_na() 

sample_size<- sample_size %>% 
  left_join(med)

#Palette

# iucncolors <- c("grey70", "#CF000A", "#FD6A32","#E0D210","#C5E51E","#47A83E")

# Plot

(fig4b <- mc2 %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(threat = gsub("b_threat", "", .variable),
           threat = gsub("NonMth", "Non-th", threat))  %>%
    naniar::replace_with_na_all(condition= ~.x == "") %>% 
    drop_na() %>% 
    ggplot(aes(y = threat, x = .value)) +
    geom_vline(xintercept = 0, 
               linetype = "dashed", 
               colour="grey50") +
    geom_pointinterval(aes(xmin = .lower, 
                           xmax = .upper,
                           colour=threat),
                       interval_size_range = c(0.5, 2)) +
    geom_text(data = sample_size, aes(x=med, y=threat, 
                                      label = paste0("n=", n)),
              vjust   = -1)+
    scale_color_manual("", values = c("#5CC86F", "#FF8076", "grey65")) +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    xlim(-0.1,0.1)+
    labs(x="Posterior estimate", y = "IUCN category") + 
    theme(legend.position = "none"))

## Panel c: Source -------------------------------------------------------------

# Function to keep just two decimals

scaleFUN <- function(x) sprintf("%.2f", x)

# Create a summary of the sample size

sample_size <- pops_data %>%
  drop_na(Data_source_type) %>%
  rename(source=Data_source_type) %>% 
  group_by(source) %>%
  summarise(median=median(mu),
            n=n(),
            max=max(mu))

med <- ml %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(source = gsub("b_Data_source_type", "", .variable),
         source= gsub("Greyliterature", "Grey literature", source),
         source= gsub("PeerMreview", "Peer-review", source))  %>%
  group_by(source) %>% 
  summarise(med=quantile(.value, probs = 0.9)) %>% 
  naniar::replace_with_na_all(condition= ~.x == "") %>% 
  drop_na() 

sample_size<- sample_size %>% 
  left_join(med)

# Plot

(fig4c <- ml %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(source = gsub("b_Data_source_type", "", .variable),
           source= gsub("Greyliterature", "Grey literature", source),
           source= gsub("PeerMreview", "Peer-review", source))  %>%
    naniar::replace_with_na_all(condition= ~.x == "") %>% 
    drop_na() %>% 
    ggplot(aes(y = source, x = .value)) +
    geom_vline(xintercept = 0, 
               linetype = "dashed", 
               colour="grey50") +
    geom_pointinterval(aes(xmin = .lower, 
                           xmax = .upper,
                           colour=source),
                       interval_size_range = c(0.5, 2)) +
    geom_text(data = sample_size, aes(x=med, y=source, 
                                      label = paste0("n=", n)),
              vjust   = -1)+
    scale_color_manual("", values = c("#2f2c28", "#8b4513", "#cd853f")) +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
   # xlim(-0.06, 0.06)+
    labs(x="Posterior estimate", y = "Source") + 
    theme(legend.position = "none"))

## Combine plot ----------------------------------------------------------------

(fig4 <- fig4a+fig4b+fig4c+plot_annotation(tag_levels = "a"))

# Save the plot

ggsave(fig4,filename = "Figure 4.pdf",
       path = ResultPath, height = 5, width = 12)


# Figure 5: LPD system and taxa ####################################################

# Sample size

sample_size <- pops_data %>% 
  drop_na(database) %>%
  group_by(database, System, Class) %>% 
  summarise(median=median(mu),
            n=n(),
            max=max(mu))

med <- mdi %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub("b_database", "", .variable),
         database = gsub("\\:.*", "", .variable),
         System=gsub("LPD:System", "", .variable),
         System=gsub("IbeV:System", "", System),
         System = gsub("\\:.*", "", System),
         Class = gsub("Class", "", .variable),
         Class = gsub("IbeV", "", Class),
         Class = gsub("LPD", "", Class),
         Class =gsub(":", "", Class),
         Class =gsub("System", "", Class),
         Class =gsub("Freshwater", "", Class),
         Class =gsub("Marine", "", Class),
         Class =gsub("Terrestrial", "", Class),
         Class =gsub("D", "", Class),
         Class = gsub('[[:digit:]]+', '', Class),
         Class = gsub(".*:", "", Class),
         Class = gsub('Bonyfishes', 'Bony fishes', Class),
         System = gsub("D", "/", System), 
         System = gsub("Class", "", System), 
         System = gsub("Mammalia", "", System), 
         System = gsub("Aves", "", System), 
         System = gsub("Bonyfishes", "", System), 
         System = gsub("Reptilia", "", System), 
         System = gsub("Elasmobranchii", "", System), 
         System = gsub("Amphibia", "", System), 
         System =gsub(":", "", System),
         System = gsub('[[:digit:]]+', '', System)) %>% 
  drop_na() %>% 
  filter(System!="Freshwater"|Class!="Bony fishes"|database!="LPD",
         System!="Marine"|Class!="Aves"|database!="LPD",
         System!="Freshwater/Terrestrial"|Class!="Reptilia"|database!="IbeV",
         System!="Marine/Freshwater/Terrestrial"|
           Class!="Mammalia"|database!="LPD",
         System!="Marine"|Class!="Reptilia",
         System!="Freshwater"|Class!="Elasmobranchii",
         System!="Marine"|Class!="Amphibia",
         System!="Terrestrial"|Class!="Elasmobranchii",
         System!="Terrestrial"|Class!="Bony fishes",
         System!="Terrestrial"|Class!="Amphibia",
         System!="Freshwater"|Class!="Reptilia",
         System!="Freshwater"|Class!="Mammalia",
         System!="Freshwater"|Class!="Elasmobranchii",
         System!="Freshwater"|Class!="Amphibia",
         System!="Freshwater/Terrestrial"|Class!="Elasmobranchii",
         System!="Freshwater/Terrestrial"|Class!="Bony fishes",
         System!="Marine/Freshwater"|Class!="Reptilia",
         System!="Marine/Freshwater"|Class!="Mammalia",
         System!="Marine/Freshwater"|Class!="Elasmobranchii",
         System!="Marine/Freshwater"|Class!="Aves",
         System!="Marine/Freshwater"|Class!="Amphibia",
         System!="Marine/Freshwater/Terrestrial"|Class!="Reptilia",
         System!="Marine/Freshwater/Terrestrial"|Class!="Bony fishes",
         System!="Marine/Freshwater/Terrestrial"|Class!="Elasmobranchii",
         System!="Marine/Freshwater/Terrestrial"|Class!="Amphibia",
         System!="Marine/Terrestrial"|Class!="Reptilia",
         System!="Marine/Terrestrial"|Class!="Bony fishes",
         System!="Marine/Terrestrial"|Class!="Elasmobranchii",
         System!="Marine/Terrestrial"|Class!="Amphibia",
         System!="Marine/Terrestrial"|Class!="Mammalia")  %>%
  group_by(database, System, Class) %>% 
  summarise(med=quantile(.value, probs = 0.9)) %>% 
  drop_na() 

sample_size<- sample_size %>% 
  left_join(med)%>% 
  mutate(System= factor(System, levels=c("Freshwater", "Marine", "Terrestrial",
                                         "Marine/Freshwater", "Marine/Terrestrial", 
                                         "Freshwater/Terrestrial", 
                                         "Marine/Freshwater/Terrestrial")))

# Plot intervals 

dat <- mdi %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub("b_database", "", .variable),
         database = gsub("\\:.*", "", .variable),
         System=gsub("LPD:System", "", .variable),
         System=gsub("IbeV:System", "", System),
         System = gsub("\\:.*", "", System),
         Class = gsub("Class", "", .variable),
         Class = gsub("IbeV", "", Class),
         Class = gsub("LPD", "", Class),
         Class =gsub(":", "", Class),
         Class =gsub("System", "", Class),
         Class =gsub("Freshwater", "", Class),
         Class =gsub("Marine", "", Class),
         Class =gsub("Terrestrial", "", Class),
         Class =gsub("D", "", Class),
         Class = gsub('[[:digit:]]+', '', Class),
         Class = gsub(".*:", "", Class),
         Class = gsub('Bonyfishes', 'Bony fishes', Class),
         System = gsub("D", "/", System), 
         System = gsub("Class", "", System), 
         System = gsub("Mammalia", "", System), 
         System = gsub("Aves", "", System), 
         System = gsub("Bonyfishes", "", System), 
         System = gsub("Reptilia", "", System), 
         System = gsub("Elasmobranchii", "", System), 
         System = gsub("Amphibia", "", System), 
         System =gsub(":", "", System),
         System = gsub('[[:digit:]]+', '', System)) %>% 
  drop_na() %>% 
  filter(System!="Freshwater"|Class!="Bony fishes"|database!="LPD",
         System!="Marine"|Class!="Aves"|database!="LPD",
         System!="Freshwater/Terrestrial"|Class!="Reptilia"|database!="IbeV",
         System!="Marine/Freshwater/Terrestrial"|
           Class!="Mammalia"|database!="LPD",
         System!="Marine"|Class!="Reptilia",
         System!="Freshwater"|Class!="Elasmobranchii",
         System!="Marine"|Class!="Amphibia",
         System!="Terrestrial"|Class!="Elasmobranchii",
         System!="Terrestrial"|Class!="Bony fishes",
         System!="Terrestrial"|Class!="Amphibia",
         System!="Freshwater"|Class!="Reptilia",
         System!="Freshwater"|Class!="Mammalia",
         System!="Freshwater"|Class!="Elasmobranchii",
         System!="Freshwater"|Class!="Amphibia",
         System!="Freshwater/Terrestrial"|Class!="Elasmobranchii",
         System!="Freshwater/Terrestrial"|Class!="Bony fishes",
         System!="Marine/Freshwater"|Class!="Reptilia",
         System!="Marine/Freshwater"|Class!="Mammalia",
         System!="Marine/Freshwater"|Class!="Elasmobranchii",
         System!="Marine/Freshwater"|Class!="Aves",
         System!="Marine/Freshwater"|Class!="Amphibia",
         System!="Marine/Freshwater/Terrestrial"|Class!="Reptilia",
         System!="Marine/Freshwater/Terrestrial"|Class!="Bony fishes",
         System!="Marine/Freshwater/Terrestrial"|Class!="Elasmobranchii",
         System!="Marine/Freshwater/Terrestrial"|Class!="Amphibia",
         System!="Marine/Terrestrial"|Class!="Reptilia",
         System!="Marine/Terrestrial"|Class!="Bony fishes",
         System!="Marine/Terrestrial"|Class!="Elasmobranchii",
         System!="Marine/Terrestrial"|Class!="Amphibia",
         System!="Marine/Terrestrial"|Class!="Mammalia") %>% 
  mutate(System= factor(System, levels=c("Freshwater", "Marine", "Terrestrial",
                                         "Marine/Freshwater", "Marine/Terrestrial", 
                                         "Freshwater/Terrestrial", 
                                         "Marine/Freshwater/Terrestrial")))

# Plot

(fig5 <- mdi %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(.variable = gsub("b_database", "", .variable),
           database = gsub("\\:.*", "", .variable),
           System=gsub("LPD:System", "", .variable),
           System=gsub("IbeV:System", "", System),
           System = gsub("\\:.*", "", System),
           Class = gsub("Class", "", .variable),
           Class = gsub("IbeV", "", Class),
           Class = gsub("LPD", "", Class),
           Class =gsub(":", "", Class),
           Class =gsub("System", "", Class),
           Class =gsub("Freshwater", "", Class),
           Class =gsub("Marine", "", Class),
           Class =gsub("Terrestrial", "", Class),
           Class =gsub("D", "", Class),
           Class = gsub('[[:digit:]]+', '', Class),
           Class = gsub(".*:", "", Class),
           Class = gsub('Bonyfishes', 'Bony fishes', Class),
           System = gsub("D", "/", System), 
           System = gsub("Class", "", System), 
           System = gsub("Mammalia", "", System), 
           System = gsub("Aves", "", System), 
           System = gsub("Bonyfishes", "", System), 
           System = gsub("Reptilia", "", System), 
           System = gsub("Elasmobranchii", "", System), 
           System = gsub("Amphibia", "", System), 
           System =gsub(":", "", System),
           System = gsub('[[:digit:]]+', '', System)) %>% 
    drop_na() %>% 
    filter(System!="Freshwater"|Class!="Bony fishes"|database!="LPD",
           System!="Marine"|Class!="Aves"|database!="LPD",
           System!="Freshwater/Terrestrial"|Class!="Reptilia"|database!="IbeV",
           System!="Marine/Freshwater/Terrestrial"|
             Class!="Mammalia"|database!="LPD",
           System!="Marine"|Class!="Reptilia",
           System!="Freshwater"|Class!="Elasmobranchii",
           System!="Marine"|Class!="Amphibia",
           System!="Terrestrial"|Class!="Elasmobranchii",
           System!="Terrestrial"|Class!="Bony fishes",
           System!="Terrestrial"|Class!="Amphibia",
           System!="Freshwater"|Class!="Reptilia",
           System!="Freshwater"|Class!="Mammalia",
           System!="Freshwater"|Class!="Elasmobranchii",
           System!="Freshwater"|Class!="Amphibia",
           System!="Freshwater/Terrestrial"|Class!="Elasmobranchii",
           System!="Freshwater/Terrestrial"|Class!="Bony fishes",
           System!="Marine/Freshwater"|Class!="Reptilia",
           System!="Marine/Freshwater"|Class!="Mammalia",
           System!="Marine/Freshwater"|Class!="Elasmobranchii",
           System!="Marine/Freshwater"|Class!="Aves",
           System!="Marine/Freshwater"|Class!="Amphibia",
           System!="Marine/Freshwater/Terrestrial"|Class!="Reptilia",
           System!="Marine/Freshwater/Terrestrial"|Class!="Bony fishes",
           System!="Marine/Freshwater/Terrestrial"|Class!="Elasmobranchii",
           System!="Marine/Freshwater/Terrestrial"|Class!="Amphibia",
           System!="Marine/Terrestrial"|Class!="Reptilia",
           System!="Marine/Terrestrial"|Class!="Bony fishes",
           System!="Marine/Terrestrial"|Class!="Elasmobranchii",
           System!="Marine/Terrestrial"|Class!="Amphibia",
           System!="Marine/Terrestrial"|Class!="Mammalia") %>% 
    mutate(System= factor(System, levels=c("Freshwater", "Marine", "Terrestrial",
                                           "Marine/Freshwater", "Marine/Terrestrial", 
                                           "Freshwater/Terrestrial", 
                                           "Marine/Freshwater/Terrestrial"))) %>%     
    ggplot(aes(y = .value, x = Class, group=database)) +
    facet_wrap(.~System,
               nrow = 3,
               scales = "free", drop = T)+
    geom_hline(yintercept = 0, 
               linetype = "dashed", 
               colour="grey50") +
    geom_pointinterval(aes(ymin = .lower, 
                           ymax = .upper,
                           colour=database),
                       interval_size_range = c(0.5, 2),
                       position = position_dodge(width = 0.5)) +
    geom_text_repel(data = sample_size, 
              aes(x=Class, y=med, 
                  label = paste0("n=", n)),
              position = position_dodge(width = 0.5))+
    scale_color_manual("", values = database_pal) +
    scale_y_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    labs(x="Class", y = "Posterior estimate") + 
    theme(legend.position = "none",
          panel.spacing.x = unit(8, "mm"),
          axis.text.x = element_text(angle = 25, hjust = 1),
          plot.margin = unit(c(1,1,1,1), "cm")))

# Save the plot

ggsave(fig5,
       filename = "Figure 5.pdf",
       path = ResultPath, height = 10, width = 14)
