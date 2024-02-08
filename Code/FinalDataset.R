# --------------------------------------------------------------------------------------- #
# - FILE NAME:   RedList.R
# - DATE:        25/02/2022
# - DESCRIPTION: Code to match our data with the redlist categories
# - AUTHORS: Pol Capdevila (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(tidyverse)
library(dplyr)
library(cowplot)
library(rredlist)
library(foreach)
library(parallel)
library(doParallel)
library(taxize)

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

load(paste0(DataPath,"/PopChangeData.RData"))

# Resolve the taxonomy --------------------------------------------------------

# Remove some mistakes in the supplied names

pops_data <- pops_data %>% 
  mutate(Species= sub("\\([^)]+\\)", "", Species))

# Using data source from gbif

sources <- gnr_datasources()
taxonomy <- gnr_resolve(pops_data$Species,
                        data_source_ids = sources$id[sources$title=="GBIF Backbone Taxonomy"])

# Remove double spaces and correct those with only genus

taxonomy <- taxonomy %>% 
  mutate(Species=sub("^([[:alnum:]]+ [[:alnum:]]+) .*", "\\1", matched_name),
         Species=ifelse(grepl("sp\\.|spp\\.", submitted_name), sub(" .*", "", Species),
                        Species))


#Match the IUCN criteria ------------------------------------------------------

IUCN_REDLIST_KEY = "8d9c556bb4aa2efeacdc49775c0d0514b42139c96d5ed0f6b9ca5a495ffa49a4"

# Get Red List version

rl_version(key=IUCN_REDLIST_KEY)

# # Correct names with a space
# 
# pops_data <- 
#   pops_data %>% 
#   mutate(Species=sub("^(\\S+\\s+\\S+)\\s+", "\\1", Species))

# Register the clusters

cl <- makeCluster(detectCores())
registerDoParallel(cl)

# Run the for each

red_list <- foreach(i = 1:length(taxonomy$Species), 
                    .combine="rbind", .packages="rredlist") %dopar%
  {tryCatch(rl_search(taxonomy$Species[i],
                          key = IUCN_REDLIST_KEY,
                          parse = T)$result,
                error=function(e) NULL)}

# Stop the cluster

stopCluster(cl)

# Distinct redlist

red_list <- red_list %>% 
  distinct(scientific_name, .keep_all = T)  

# We match it with the original data set

pops_data <- pops_data %>% 
  left_join(red_list[,c("scientific_name", "category", "population_trend", 
  "marine_system", "freshwater_system","terrestrial_system")],
            by=c("Species"="scientific_name")) %>%
  mutate(category=ifelse(category=="LR/nt", "NT",
                         category),
         System2=ifelse(marine_system==TRUE & 
                         freshwater_system==TRUE &
                         terrestrial_system==TRUE, 
                       "Marine/Freshwater/Terrestrial", 
                       ifelse(marine_system==TRUE & 
                                freshwater_system==FALSE &
                                terrestrial_system==TRUE, 
                              "Marine/Terrestrial",
                              ifelse(marine_system==TRUE & 
                                       freshwater_system==TRUE &
                                       terrestrial_system==FALSE, 
                                     "Marine/Freshwater",
                                     ifelse(marine_system==FALSE & 
                                              freshwater_system==TRUE &
                                              terrestrial_system==TRUE, 
                                            "Freshwater/Terrestrial",
                                            ifelse(marine_system==FALSE & 
                                                     freshwater_system==FALSE &
                                                     terrestrial_system==TRUE, 
                                                   "Terrestrial",
                                                   ifelse(marine_system==FALSE & 
                                                            freshwater_system==TRUE &
                                                            terrestrial_system==FALSE, 
                                                          "Freshwater",
                                                          ifelse(marine_system==TRUE & 
                                                                   freshwater_system==FALSE &
                                                                   terrestrial_system==FALSE,
                                                                 "Marine", 
                                                                 System))))))),
         System=ifelse(is.na(System2), System, System2)) %>% 
  dplyr::select(-System2)


# Factor of data base 

pops_data <- pops_data %>% 
  mutate(database=gsub(" .*", "", ID), 
         database=gsub("Vertebichos", "IbeV", database))

# Length of the dataset 

pops_data <- pops_data %>% 
  group_by(ID) %>% 
  mutate(End = as.numeric(End),
         Start = as.numeric(Start),
         length=End-Start)

# Create a decade 

floor_decade <- function(value){ return(value - value %% 10) }#Function to round 
# to decades

pops_data <- pops_data %>% 
  mutate(decade=as.factor(floor_decade(Start)))

# Conservation category

pops_data <- pops_data %>% 
  mutate(threat= ifelse(category=="LC"|category=="NT", "Non-threatened",
                        ifelse(category=="DD", 
                               "Unknown", 
                               "Threatened")),
         threat=ifelse(is.na(category), "Unknown", threat))

# Save it 

setwd(DataPath)
save(pops_data, file = "FinalData.RData")
