# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Descriptive analyses.R
# - DATE:        25/02/2022
# - DESCRIPTION: Code to describe the patterns observed in the Iberian populations
# - AUTHORS:
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggdist)
library(ggsci)
library(RCurl)
library(png)
library(cowplot)
library(patchwork)

# Set default ggplot theme

theme_set(theme_minimal()+
            theme(axis.title.x = element_text(size=15, margin = margin(t = 10, r = 0, b = 0, l = 0)),
                  axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.line.x = element_line(color="black", size = 0.5),
                  axis.line.y = element_line(color="black", size = 0.5),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(color="black", size = 12),
                  axis.text.y = element_text(color="black", size = 12),
                  strip.text.x = element_text(size = 12),
                  axis.ticks = element_line(color="black")))

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results")

# Read in data

load(paste0(DataPath,"/FinalData.RData"))

# Explore the data

length(unique(pops_data$ID))
length(unique(pops_data$Species))

# Map --------------------------------------------------------------------------

# Set the world map

world <- map_data("world") #, region = c("spain", "portugal"))

# Create plot1

(map <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id=region), 
             color = "gray60", fill = "gray80", size = 0.3) +
    geom_point(data = pops_data, 
               aes(x = Longitude, y = Latitude, 
                   fill= mu), 
               alpha= 0.6,shape=21,  size=5) +
    scale_fill_gradient2("Population trend",
                         midpoint = 0, 
                         low = "#c1666b", 
                         mid = "#e4dfda",
                         high = "#4281a4")+
    scale_y_continuous(limits = c(35, 45))+
    scale_x_continuous(limits = c(-10,8)) + 
    theme_map())

# Save it

ggsave("Map.pdf", map,
       width = 8, height = 6,
       path = ResultPath)

# Overall population trend -----------------------------------------------------

(p1 <- pops_data %>% 
  ggplot(aes(x=mu)) +
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_density(fill="grey70", colour=NA, alpha=.5) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x=expression(paste("Population change(", mu, ")"))))

ggsave("Dist1.pdf",p1,
       width = 8, height = 6,
       path = ResultPath)

# Population trend by taxonomic group ------------------------------------------

# Create a summary data frame 

data_summary <- pops_data %>% 
  group_by(Class) %>% 
  summarise(median=median(mu),
            n=n(),
            max=max(mu))

# Import phylopics

mammals <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/05f87521-20d4-4a05-8ac6-aa0bab7f1394.512.png"))
birds <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/a6ef5684-5683-4a46-aa2d-d39c0d134ba7.512.png"))
amphibians <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/35237d88-e323-44dc-8349-b70fe078c5e7.512.png"))
reptiles <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/a359e147-c088-4c18-a2f1-49abfb2b9325.512.png"))
fish <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/8a9bc7b5-360c-4f2a-9487-723dcd3deb8b.256.png"))
shark <-  readPNG(getURLContent("http://phylopic.org/assets/images/submissions/36d54b94-35a6-4a66-a41e-b8f28534e70c.512.png"))

# Plot 

(p2 <- pops_data %>% 
    ggplot(aes(x=mu, y=Class,
               colour=Class, group=Class))+
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    geom_point(size = 2, alpha = 0.1,
               position = position_nudge(y = -0.15) )+
    stat_halfeye(aes(color = Class,
                     fill=after_scale(colorspace::lighten(color, .3))),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8) +
    scale_color_uchicago() +
    geom_text(data=data_summary,
              aes(x = max+0.02, 
                  label = paste("n=", n)),
              stat = "unique",
              color = "black",
              fontface = "bold",
              size = 3.4,
              nudge_y = .15)+
    labs(x=expression(paste("Population change(", mu, ")")),
         y="")+
    xlim(-0.2, 0.3)+
    theme(legend.position = "none"))

# # Add the silhouettes 
# 
# (ga <- ggdraw(p2) + 
#     draw_image(reptiles, 
#                x = 0.715, 
#                y = 0.83, 
#                width = 0.1, height = 0.07)+
#     draw_image(mammals, x = 0.75, 
#                y = 0.7, 
#                width = 0.1, height = 0.07)+
#     draw_image(shark, x = 0.62, 
#                y = 0.6, 
#                width = 0.1, height = 0.07)+
#     draw_image(birds, x = 0.89, 
#                y = 0.46, 
#                width = 0.1, height = 0.07)+
#     draw_image(amphibians, x = 0.61, 
#              y =0.32, 
#              width = 0.1, height = 0.05)+
#     draw_image(fish, x = 0.71, 
#                y = 0.2, 
#                width = 0.1, height = 0.04))

  
ggsave("PopChangeByTaxa.pdf",p2,
       width = 8, height = 6,
       path = ResultPath)

# By system --------------------------------------------------------------------

# Create a summary data frame 

data_summary <- pops_data %>% 
  group_by(System) %>% 
  summarise(median=median(mu),
            n=n(),
            max=max(mu))

# Plot 

(p3 <- pops_data %>% 
    ggplot(aes(x=mu, y=System,
               colour=System, group=System))+
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    geom_point(size = 2, alpha = 0.1,
               position = position_nudge(y = -0.15) )+
    stat_halfeye(aes(color = System,
                     fill=after_scale(colorspace::lighten(color, .3))),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8) +
    scale_color_manual("", values = c("#A1D6E2","#336B87", "#CC954E")) +
    geom_text(data=data_summary,
              aes(x = max+0.02, 
                  label = paste("n=", n)),
              stat = "unique",
              color = "black",
              fontface = "bold",
              size = 3.4,
              nudge_y = .15)+
    labs(x=expression(paste("Population change(", mu, ")")),
         y="")+
    xlim(-0.2, 0.3)+
    theme(legend.position = "none"))

ggsave(p3,filename = "System.pdf", height = 6, width = 8,
       path = ResultPath)

# Protection level -------------------------------------------------------------

# Create a summary of the sample size 

data_summary <- pops_data %>% 
  drop_na(Protected) %>% 
  group_by(Protected) %>% 
  summarise(median=median(mu),
            n=n(),
            max=max(mu))

# Plot

(prot <- pops_data %>% drop_na(Protected) %>% 
   ggplot(aes(x=mu, y=Protected,
              colour=Protected, 
              group=Protected))+
   geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
   geom_point(size = 2, alpha = 0.1,
              position = position_nudge(y = -0.15) )+
   stat_halfeye(aes(color = Protected,
                    fill=after_scale(colorspace::lighten(color, .3))),
                shape = 18,
                point_size = 3,
                interval_size = 1.8) +
   scale_color_uchicago() +
   geom_text(data=data_summary,
             aes(x = max+0.02, 
                 label = paste("n=", n)),
             stat = "unique",
             color = "black",
             fontface = "bold",
             size = 3.4,
             nudge_y = .15)+
   labs(x=expression(paste("Population change(", mu, ")")),
        y="Protected area")+
   xlim(-0.2, 0.3)+
   theme(legend.position = "none")) 

# Save

ggsave(prot,filename = "Protection.pdf",
       path = ResultPath, height = 6, width = 8)

# Combined plot 

final <- p2+p3+
  prot+plot_annotation(tag_levels = "a")

# Save

ggsave(final, filename = "Combined.pdf",
       path = ResultPath, 
       height = 6, width = 12)

# Taxa vs protected ------------------------------------------------------------

# Create a summary of the sample size 

data_summary <- pops_data %>% 
  drop_na(Protected) %>% 
  group_by(Protected, Class) %>% 
  summarise(median=median(mu),
            n=n(),
            max=max(mu))

# Plot

(prot2 <- pops_data %>% drop_na(Protected) %>% 
    ggplot(aes(x=mu, y=Protected,
               colour=Protected, 
               group=Protected))+
    facet_wrap(Class~.)+
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    geom_point(size = 2, alpha = 0.1,
               position = position_nudge(y = -0.15) )+
    stat_halfeye(aes(color = Protected,
                     fill=after_scale(colorspace::lighten(color, .3))),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8) +
    scale_color_uchicago() +
    geom_text(data=data_summary,
              aes(x = max+0.02, 
                  label = paste("n=", n)),
              stat = "unique",
              color = "black",
              fontface = "bold",
              size = 3.4,
              nudge_y = .15)+
    labs(x=expression(paste("Population change(", mu, ")")),
         y="")+
    xlim(-0.2, 0.3)+
    theme(legend.position = "none")) 

# System vs protected ------------------------------------------------------------

# Create a summary of the sample size

data_summary <- pops_data %>%
  drop_na(Protected) %>%
  group_by(Protected, System) %>%
  summarise(median=median(mu),
            n=n(),
            max=max(mu))

# Plot

(prot3 <- pops_data %>% drop_na(Protected) %>%
    ggplot(aes(x=mu, y=Protected,
               colour=Protected,
               group=Protected))+
    facet_wrap(System~.)+
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    geom_point(size = 2, alpha = 0.1,
               position = position_nudge(y = -0.15) )+
    stat_halfeye(aes(color = Protected,
                     fill=after_scale(colorspace::lighten(color, .3))),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8) +
    scale_color_uchicago() +
    geom_text(data=data_summary,
              aes(x = max+0.02,
                  label = paste("n=", n)),
              stat = "unique",
              color = "black",
              fontface = "bold",
              size = 3.4,
              nudge_y = .15)+
    labs(x=expression(paste("Population change(", mu, ")")),
         y="")+
    xlim(-0.2, 0.3)+
    theme(legend.position = "none"))

# Source -----------------------------------------------------------------------

# Create a summary of the sample size

data_summary <- pops_data %>%
  drop_na(Protected) %>%
  group_by(Data_source_type) %>%
  summarise(median=median(mu),
            n=n(),
            max=max(mu))

# Plot

(source <- pops_data %>%
    ggplot(aes(x=mu, y=Data_source_type,
               colour=Data_source_type,
               group=Data_source_type))+
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    geom_point(size = 2, alpha = 0.1,
               position = position_nudge(y = -0.15) )+
    stat_halfeye(aes(color = Data_source_type,
                     fill=after_scale(colorspace::lighten(color, .3))),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8) +
    scale_color_uchicago() +
    geom_text(data=data_summary,
              aes(x = max+0.02,
                  label = paste("n=", n)),
              stat = "unique",
              color = "black",
              fontface = "bold",
              size = 3.4,
              nudge_y = .15)+
    labs(x=expression(paste("Population change(", mu, ")")),
         y="")+
    xlim(-0.2, 0.3)+
    theme(legend.position = "none"))

# Save the plot 

ggsave(source,filename = "Source.pdf",
       path = ResultPath, height = 6, width = 8)

# IUCN criteria ----------------------------------------------------------------

# Create a summary of the sample size

data_summary <- pops_data %>%
  drop_na(category) %>%
  group_by(category) %>%
  summarise(median=median(mu),
            n=n(),
            max=max(mu))

#Palette

iucncolors <- c("#47A83E","#C5E51E","#E0D210","#FD6A32","#CF000A")

# Plot

(iucn <- pops_data %>%
    drop_na(category) %>%
    mutate(category=factor(category,
                           c("CR", "EN", "VU", "NT", "LC", "DD"))) %>%
    ggplot(aes(x=mu, y=category,
               colour=category,
               group=category))+
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    geom_point(size = 2, alpha = 0.1,
               position = position_nudge(y = -0.15) )+
    stat_halfeye(aes(color = category,
                     fill=after_scale(colorspace::lighten(color, .5))),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8) +
    scale_color_manual(values = c(rev(iucncolors), "grey20")) +
    geom_text(data=data_summary,
              aes(x = max+0.02,
                  label = paste("n=", n)),
              stat = "unique",
              color = "black",
              fontface = "bold",
              size = 3.4,
              nudge_y = .15)+
    labs(x=expression(paste("Population change(", mu, ")")),
         y="")+
    xlim(-0.2, 0.3)+
    theme(legend.position = "none"))

# Save the plot

ggsave(iucn,filename = "IUCN.pdf",
       path = ResultPath, height = 6, width = 8)

# IUCN trend -------------------------------------------------------------------

# Create a summary of the sample size

data_summary <- pops_data %>%
  drop_na(population_trend) %>%
  group_by(population_trend) %>%
  summarise(median=median(mu),
            n=n(),
            max=max(mu))


# Plot

(iucn2 <- pops_data %>%
    drop_na(population_trend) %>%
    ggplot(aes(x=mu, y=population_trend,
               colour=population_trend,
               group=population_trend))+
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    geom_point(size = 2, alpha = 0.1,
               position = position_nudge(y = -0.15) )+
    stat_halfeye(aes(color = population_trend,
                     fill=after_scale(colorspace::lighten(color, .5))),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8) +
    scale_color_jco() +
    geom_text(data=data_summary,
              aes(x = max+0.02,
                  label = paste("n=", n)),
              stat = "unique",
              color = "black",
              fontface = "bold",
              size = 3.4,
              nudge_y = .15)+
    labs(x=expression(paste("Population change(", mu, ")")),
         y="")+
    xlim(-0.2, 0.3)+
    theme(legend.position = "none"))

# Save the plot

ggsave(iucn2,filename = "IUCN2.pdf",
       path = ResultPath, height = 6, width = 8)
