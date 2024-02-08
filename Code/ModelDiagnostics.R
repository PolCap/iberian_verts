# --------------------------------------------------------------------------------------- #
# - FILE NAME:   ModelDiagnostics.R         
# - DATE:        07/01/2021
# - DESCRIPTION: Basic diagnostics for the models. 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

library(tidyverse)
library(brms)
library(rstan)
library(data.table)
library(dplyr)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(cowplot)

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load the models 

load(paste0(ResultPath, "/models.RData"))

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
                  plot.title = element_text(hjust = 0.5),
                  legend.text = element_text(size = 12)))


# Model convergence ------------------------------------------------------------
# Neff for system and taxa

ratios <- neff_ratio(m1)
mcmc_neff(ratios)

# Neff for data source

ratios <- neff_ratio(md)
mcmc_neff(ratios)

# Neff for IUCN

ratios <- neff_ratio(mc2)
mcmc_neff(ratios)

# Neff for data source

ratios <- neff_ratio(ml)
mcmc_neff(ratios)

# Normality of residuals -------------------------------------------------------
# Global 

(p1a <- m1$data %>% 
   mutate(std_resid = residuals(m1)[ , "Estimate"]) %>% 
   ggplot(aes(std_resid)) + 
   geom_histogram(aes(y=..density..),
                  colour="#9B9B9B", fill="#BFBFBF")+ 
   scale_y_continuous(expand = c(0,0))+
   scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
   labs(x="Standradised residuals", y="Density"))

# Database 

(p2a <- md$data %>% 
    mutate(std_resid = residuals(md)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# IUCN category

(p3a <- mc2$data %>% 
    mutate(std_resid = residuals(mc2)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# Data source

(p4a <- ml$data %>% 
    mutate(std_resid = residuals(ml)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# Interactions data source

(p5a <- mdi$data %>% 
    mutate(std_resid = residuals(mdi)[ , "Estimate"]) %>% 
    ggplot(aes(std_resid)) + 
    geom_histogram(aes(y=..density..),
                   colour="#9B9B9B", fill="#BFBFBF")+ 
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(x="Standradised residuals", y="Density"))

# Combine the individual plots 

combined_plot <- (p1a + p2a + p3a + 
                    p4a + p5a & 
                    labs(x = NULL, y = NULL)) + 
  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 3) &
  theme(plot.tag = element_text(face = 'bold'))

# Create a separate plot for the y-axis label 

ylabel <- ggplot(data.frame(l = p1a$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size=5, angle = 90) + 
  theme_void() +
  coord_cartesian(clip = "off")

# Create a separate plot for the x-axis label 

xlabel <- ggplot(data.frame(l = p1a$labels$x, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size=5) + 
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the original figure with the y-axis label 

(top <- cowplot::plot_grid(ylabel, combined_plot,
                           rel_widths = c(1, 25)))

# Combine with the x-axis label 

(figureS2 <- cowplot::plot_grid(top, xlabel, nrow = 2,
                                rel_heights = c(25, 1)))

# Save the final plot 

ggsave("Figure S2.pdf",path = ResultPath, figureS2, 
       width = 6, height = 8)

# Homocedasticity --------------------------------------------------------------

# System

(p2a <- m1$data %>% 
   mutate(predict_y = predict(m1)[ , "Estimate"], 
          std_resid = residuals(m1)[ , "Estimate"]) %>% 
   ggplot(aes(predict_y, std_resid)) + 
   geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
   stat_smooth(se = FALSE,colour="#9B9B9B") +
   scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
   scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
   labs(y="Standardised residuals", 
        x="Predicted values"))

# Database

(p2b <- md$data %>% 
    mutate(predict_y = predict(md)[ , "Estimate"], 
           std_resid = residuals(md)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# IUCN

(p2c <- mc2$data %>% 
    mutate(predict_y = predict(mc2)[ , "Estimate"], 
           std_resid = residuals(mc2)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# Data source 

(p2d <- ml$data %>% 
    mutate(predict_y = predict(ml)[ , "Estimate"], 
           std_resid = residuals(ml)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# Interaction data source, system taxa

(p2e <- mdi$data %>% 
    mutate(predict_y = predict(mdi)[ , "Estimate"], 
           std_resid = residuals(mdi)[ , "Estimate"]) %>% 
    ggplot(aes(predict_y, std_resid)) + 
    geom_point(size = 1.5, shape=21, fill="#BFBFBF") + 
    stat_smooth(se = FALSE,colour="#9B9B9B") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
    labs(y="Standardised residuals", 
         x="Predicted values"))

# Combine the plots into a single figure

combined_plot <-  (p2a+p2b+p2c+
                     p2d+p2e & labs(x = NULL, y = NULL)) +
  plot_annotation(tag_levels = "a")+
  plot_layout(nrow = 3) &
  theme(plot.tag = element_text(face = 'bold'))

# Create a separate plot for the y-axis label 

ylabel <- ggplot(data.frame(l = p2a$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size=5, angle = 90) + 
  theme_void() +
  coord_cartesian(clip = "off")

# Create a separate plot for the x-axis label 

xlabel <- ggplot(data.frame(l = p2a$labels$x, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size=5) + 
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the original figure with the y-axis label 

(top <- cowplot::plot_grid(ylabel, combined_plot,
                           rel_widths = c(1, 25)))

# Combine with the x-axis label 

(figureS3 <- cowplot::plot_grid(top, xlabel, nrow = 2,
                                rel_heights = c(25, 1)))

# Save the combined figure 

ggsave("Figure S3.pdf", figureS3, 
       path = ResultPath,
       height = 8, width = 6)

# Posterior predictive check ---------------------------------------------------

color_scheme_set("darkgray")

# General 

(p3a <- pp_check(m1, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

pp_check(m1, type = "stat_grouped", stat = "mean", group = "System")
pp_check(m1, type = "stat_grouped", stat = "mean", group = "Class")


# Database

(p3b <- pp_check(md, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))


# IUCN

(p3c <- pp_check(mc2, ndraws = 100)+ 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

# Data source

(p3d <- pp_check(ml, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))

# Interactive data source, system and taxa

(p3e <- pp_check(mdi, ndraws =  100) + 
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits =c(-1.1, 1.1)))


# Combine  plots into a single figure

(combined_figure <- (p3a+p3b+
                       p3c+p3c+p3d & 
                        labs(x = NULL, y = NULL)) +
    plot_annotation(tag_levels = "a")+
    plot_layout(nrow = 3, guides = "collect") &
    theme(plot.tag = element_text(face = 'bold')))

# Create a separate plot for the y-axis label 

ylabel <- ggplot(data.frame(x = 1, y = 1)) +
  geom_text(aes(x, y, label = "Density"), size=5, angle = 90) + 
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the original figure with the y-axis label and x-axis label

(figureS3 <- cowplot::plot_grid(ylabel, combined_figure,
                                rel_widths = c(1, 25)))

# Save the combined figure 

ggsave("Figure S4.pdf", figureS3, 
       path = ResultPath,
       height = 8, width = 8)

