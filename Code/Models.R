# ------------------------------------------------------------------------------#
# - FILE NAME:   Models.R
# - DATE:        25/02/2022
# - DESCRIPTION: Code to perform hierarchical bayesian analyses 
# - AUTHORS: Pol Capdevila (pcapdevila.pc@gmail.com)
# ------------------------------------------------------------------------------#

rm(list=ls(all=TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

# Libraries

library(tidyverse)
library(brms)
library(data.table)
library(zoo)
library(dplyr)
library(tidybayes)

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results")

# Read in data

load(paste0(DataPath,"/FinalData.RData"))

# MODELLING #########################################################

# Set modelling parameters

iter = 8000
thin = 0.0005*iter
warmup = 0.1*iter

# Set priors

priors <- c(prior(normal(0, 1), class = b),
            prior(exponential(1), class = sigma))

# Modify the dataset to include NA in the models of protection 

pops_data <- pops_data %>%
  mutate(Protected = ifelse(is.na(Protected), "Unknown", Protected),
         Protected =gsub("No", "Unprotected", Protected),
         Protected =gsub("Yes", "Protected", Protected)) 

# Effects on mean population trend ---------------------------------------------

# System & Class

m1 <- brm(mu | se(sigma, sigma = TRUE)  ~ System:Class-1 + (1|Species), 
          iter = iter, thin = thin, warmup = warmup,
          control = list(adapt_delta = .975, max_treedepth = 12),
          data = pops_data, 
          family = gaussian, prior = priors,
          cores = 10)

# Protected

mp <- brm(mu | se(sigma, sigma = TRUE)  ~ Protected-1 + (1|Species), 
          iter = iter, thin = thin, warmup = warmup,
          control = list(adapt_delta = .975, max_treedepth = 12),
          data = pops_data, 
          family = gaussian, prior = priors,
          cores = 4)

mipd <- brm(mu | se(sigma, sigma = TRUE)  ~ Protected*Data_source_type-1 + (1|Species), 
          iter = iter, thin = thin, warmup = warmup,
          control = list(adapt_delta = .975, max_treedepth = 12),
          data = pops_data, 
          family = gaussian, prior = priors,
          cores = 10)


# Data source

ml <- brm(mu | se(sigma, sigma = TRUE)  ~ Data_source_type-1 + (1|Species), 
          iter = iter, thin = thin, warmup = warmup,
          control = list(adapt_delta = .975, max_treedepth = 12),
          data = pops_data, 
          family = gaussian, prior = priors,
          cores = 10)

mil <- brm(mu | se(sigma, sigma = TRUE)  ~ Data_source_type*System*Class-1 + (1|Species), 
          iter = iter, thin = thin, warmup = warmup,
          control = list(adapt_delta = .975, max_treedepth = 12),
          data = pops_data, 
          family = gaussian, prior = priors,
          cores = 10)

# IUCN threat

mc2 <- brm(mu | se(sigma, sigma = TRUE)  ~ threat-1 + (1|Species), 
          iter = iter, thin = thin, warmup = warmup,
          control = list(adapt_delta = .975, max_treedepth = 12),
          data = pops_data, 
          family = gaussian, prior = priors,
          cores = 10)

mci2 <- brm(mu | se(sigma, sigma = TRUE)  ~ threat*Data_source_type-1 + (1|Species), 
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 12),
           data = pops_data, 
           family = gaussian, prior = priors,
           cores = 10)

mcii <- brm(mu | se(sigma, sigma = TRUE)  ~ threat*Protected-1 + (1|Species), 
            iter = iter, thin = thin, warmup = warmup,
            control = list(adapt_delta = .975, max_treedepth = 12),
            data = pops_data, 
            family = gaussian, prior = priors,
            cores = 10)

mciii <- brm(mu | se(sigma, sigma = TRUE)  ~ threat*database-1 + (1|Species), 
            iter = iter, thin = thin, warmup = warmup,
            control = list(adapt_delta = .975, max_treedepth = 12),
            data = pops_data, 
            family = gaussian, prior = priors,
            cores = 10)


# Database

md <- brm(mu | se(sigma, sigma = TRUE)  ~ database-1 + (1|Species), 
          iter = iter, thin = thin, warmup = warmup,
          control = list(adapt_delta = .975, max_treedepth = 12),
          data = pops_data, 
          family = gaussian, prior = priors,
          cores = 10)

mdi <- brm(mu | se(sigma, sigma = TRUE)  ~ database:System:Class-1 + (1|Species), 
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 12),
           data = pops_data, 
           family = gaussian, prior = priors,
           cores = 10)

# Length of the database

mle <- brm(mu | se(sigma, sigma = TRUE)  ~ length + (1|Species),
          iter = iter, thin = thin, warmup = warmup,
          control = list(adapt_delta = .975, max_treedepth = 12),
          data = pops_data,
          family = gaussian, prior = priors,
          cores = 10)

# Decade

mdec <- brm(mu | se(sigma, sigma = TRUE)  ~ -1 + decade + (1|Species), 
           iter = iter, thin = thin, warmup = warmup,
           control = list(adapt_delta = .975, max_treedepth = 12),
           data = pops_data, 
           family = gaussian, prior = priors,
           cores = 10)

# Interactions

mip <- brm(mu | se(sigma, sigma = TRUE)  ~ Protected:System:Class-1 + (1|Species),
          iter = iter, thin = thin, warmup = warmup,
          control = list(adapt_delta = .975, max_treedepth = 15),
          data = pops_data,
          family = gaussian, prior = priors,
          cores = 10)


# Save them 

setwd(ResultPath)
save(list = ls(pattern = 'm'), 
     file = "models.RData")
