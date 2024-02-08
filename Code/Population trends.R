# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Population trends.R         
# - DATE:        25/02/2022
# - DESCRIPTION: Code to estimate the population trends 
# - AUTHORS:     
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(piecewiseSEM)
library(tidyverse)
library(data.table)
library(zoo)
library(dplyr)
library(expss)
library(MASS)
library(VGAM)
library(goeveg)
library(RcppRoll)
library(tidyr)
library(foreach)
library(doParallel)

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Read in data

lpi_data <- read.csv(paste0(DataPath,"/LPDIberia.csv"))
iberian_data <- read.csv2(paste0(DataPath,"/data_all.csv"),
                         dec = ",", encoding = "UTF-8")
iberian_metadata <- read.csv2(paste0(DataPath,"/metadata_new.csv"))
data_clean <- read.csv(paste0(DataPath,"/data_clean.csv"))

# Join the iberian data with the metadata

iberian_metadata<- iberian_metadata %>%
  filter(Series_ID!="gris3.7",
         Series_ID!="") 

iberian_pops <- iberian_data %>%
  filter(Series_ID!="gris13.6",
         Series_ID!="") %>% 
  left_join(iberian_metadata) %>%
  mutate(Year=as.numeric(Year),
         Count=as.numeric(Abundance),
         ID=paste(StudyID, Series_ID)) 

# Data cleaning

lpi_data <- lpi_data %>% 
  mutate(Year=as.numeric(Year), # Make count and year numeric
         Count=as.numeric(Count),
         # Change species names
         Species= gsub("_", " ", lpi_data$Binomial)) %>% 
  filter(!ID%in%c(9135, 9126, 23363, 23364, 
                  23365, 23366, 23367, 1505, 
                  1701, 4636)) # Filter out the IDs which were not inside the Iberian Peninsula 

# Estimate population change

lpi_data <- lpi_data %>%
  group_by(ID) %>%  
  # Calculate population change
  mutate(popchange= log(Count+(max(Count, na.rm=T) /100))) %>% 
  # Remove any groupings we have created in the pipe
  ungroup() %>% 
  # Create a factor for year
  mutate(year=Year,
         Dataset="LPD", 
         ID=paste(Dataset, ID)) %>% 
  rename(Protected=Protected_status) %>% 
  dplyr::select(ID, Species, Class, Order,
                System, Latitude, Longitude, 
                Data_source_type, Protected,
                Year, year, popchange, Dataset) 

# Now for our data

iberian_pops <- iberian_pops %>%
  group_by(ID) %>% 
  filter(Count>0) %>% 
  # Calculate population change
  mutate(popchange= (Count)/(max(Count,na.rm = T))) %>% 
  # Remove any groupings we have created in the pipe
  ungroup() %>% 
  # Create a factor for year
  mutate(year=Year,
         Dataset="Vertebichos",
         ID=as.factor(paste(Dataset, ID))) %>% 
  dplyr::select(ID, Species, Class, Order,
                System, Latitude, Longitude,
                Data_source_type, Protected,
                Year, year, popchange, Dataset) 
  
# Join the two datasets 

total_data <- rbind(lpi_data, iberian_pops) %>% 
  filter(ID%in%data_clean$ID)

# Spread the data 

pops <- total_data %>%
  group_by(ID) %>%
  drop_na(popchange) %>% 
  filter(length(unique(Year)) > 4) %>%
  dplyr::select(-Year) %>%
  arrange(year) %>%
  pivot_wider(id_cols = c(ID, Species, Class, Order,
                          System, Latitude, Longitude,
                          Data_source_type, Protected,
                          Dataset), 
              names_from = year, values_from=popchange)

# Create a separate dataset with metadata

pops_data <- pops %>%
  dplyr::select(ID, Species, Class, Order,
                System, Latitude, Longitude,
                Data_source_type, Protected) 

# Load the functions to estimate state-space trends from Humbert et al. 2009

source(paste0(CodePath,"/Humbert_function.R"))

# Compile all time-series into a list to estimate state-space population trends##

pop <- list()
for(i in 1:nrow(pops)){
  data <- pops[i, 11:106]
  Year <- colnames(data)[is.na(data) == FALSE]
  data <- data[Year]
  N <- as.vector(t(data))
  pop[[i]] <- data.frame(pops$ID[i], Year, N)
}

 
# Loop through every time-series, estimating mu from state-space (REML) models ----
 
for(i in 1:nrow(pops)) {
    Observed.t <- as.numeric(pop[[i]]$N)
    Time.t <- as.numeric(pop[[i]]$Year)
    print(i)

    T.t = Time.t - Time.t[1] # Time starts at zero

    for(z in 1:length(Observed.t)){
      if(Observed.t[z] == 0){
        if(max(Observed.t %% 1) == 0){
          Observed.t <- Observed.t + 1}
        if(max(Observed.t %% 1) > 0){
          Observed.t <- Observed.t + max(Observed.t %% 1)
        }
      }
    }
    Y.t <-  Observed.t # The observations
    q <-  length(Y.t) - 1 # Number of time series transitions
    qp1 <-  q + 1 # q+1
    S.t <-  T.t[2:qp1] - T.t[1:q] # Time intervals
    m <-  rep(1, qp1) # Room for Kalman means
    v <-  rep(1, qp1) # Room form variances for Kalman calculations

    # Calculating EGOE AND EGPN estimates for use as initial values ----

    # The EGOE estimates
    Ybar <-  mean(Y.t)
    Tbar <-  mean(T.t)
    mu.egoe <-  sum((T.t - Tbar)*(Y.t - Ybar))/sum((T.t - Tbar)*(T.t - Tbar))
    x0.egoe <-  Ybar - mu.egoe*Tbar
    ssq.egoe <-  0
    Yhat.egoe <-  x0.egoe + mu.egoe*T.t
    tsq.egoe <-  sum((Y.t - Yhat.egoe)*(Y.t - Yhat.egoe))/(q - 1)

    # The EGPN estimates
    Ttr <-  sqrt(S.t)
    Ytr <-  (Y.t[2:qp1] - Y.t[1:q])/Ttr
    mu.egpn <-  sum(Ttr*Ytr)/sum(Ttr*Ttr)
    Ytrhat <-  mu.egpn*Ttr
    ssq.egpn <-  sum((Ytr - Ytrhat)*(Ytr - Ytrhat))/(q - 1)
    tsq.egpn <-  0
    x0.egpn <-  Y.t[1]

    # Initial values for EGSS are averages of EGOE and EGPN values
    ssq0 <-  ssq.egpn/2 # For ML and REML
    tsq0 <-  tsq.egoe/2 # For ML and REML

    # Calculate ML & REML parameter estimates

    # The REML estimates.
    if(ssq0==0 & tsq0==0){
      EGSSreml <-  optim(par =  c(ssq0, tsq0),
                         negloglike.reml, NULL, method = "Nelder-Mead", yt = Y.t, tt = T.t)
    }else{
      EGSSreml <-  optim(par =  c(log(ssq0), log(tsq0)),
                         negloglike.reml, NULL, method = "Nelder-Mead", yt = Y.t, tt = T.t)
    }
      
    params.reml <-  c(exp(EGSSreml$par[1]), exp(EGSSreml$par[2]))
    ssq.reml <-  params.reml[1]
    #  These are the REML estimates
    tsq.reml <-  params.reml[2]
    vx <-  matrix(0, qp1, qp1)
    for (ti in 1:q){
      vx[(ti + 1):qp1, (ti + 1):qp1] <-  matrix(1, 1, (qp1 - ti))*T.t[ti + 1]
    }
    Sigma.mat <-  ssq.reml*vx
    Itausq <-  matrix(rep(0,(qp1*qp1)), nrow = q + 1, ncol = q + 1)
    diag(Itausq) <-  rep(tsq.reml, q + 1)
    V <-  Sigma.mat + Itausq
    D1mat <-  cbind(-diag(1/S.t), matrix(0, q, 1)) + cbind(matrix(0, q, 1), diag(1/S.t))
    V1mat <-  D1mat%*%V%*%t(D1mat)
    W.t <-  (Y.t[2:qp1] - Y.t[1:q])/S.t
    j1 <-  matrix(1, q, 1)
    V1inv <-  ginv(V1mat)
    mu.reml <-  (t(j1)%*%V1inv%*%W.t)/(t(j1)%*%V1inv%*%j1)
    j <-  matrix(1, qp1, 1)
    Vinv <-  ginv(V)
    x0.reml <-  (t(j)%*%Vinv%*%(Y.t-as.vector(mu.reml)*T.t))/(t(j)%*%Vinv%*%j)
    Var_mu.reml <-  1/(t(j1)%*%V1inv%*%j1) # Variance of mu
    mu_hi.reml <-  mu.reml + 1.96*sqrt(Var_mu.reml) # 95% CI for mu
    mu_lo.reml <-  mu.reml - 1.96*sqrt(Var_mu.reml)

    #  Calculate estimated population sizes for EGSS model with Kalman filter, for plotting.
    #  Choose REML estimates here for calculating model values
    mu <-  mu.reml
    ssq <-  ssq.reml
    tsq <-  tsq.reml
    x0 <-  x0.reml
    m[1] <-  x0
    #  Initial mean of Y(t)
    v[1] <-  tsq
    #  Initial variance of Y(t)
    for (ti in 1:q){
      # Loop to generate estimated population abundances
      # using Kalman filter (see equations 6 & 7, Dennis et al. (2006)).
      m[ti + 1] <-  mu + (m[ti] + ((v[ti] - tsq)/v[ti])*(Y.t[ti] - m[ti]))
      v[ti + 1] <-  tsq*((v[ti] - tsq)/v[ti]) + ssq + tsq
    }
    Predict.t <-  exp(m + ((v - tsq)/v)*(Y.t - m))
    pop[[i]]$Pred.N <- Predict.t

    #  The following statement calculates exp{E[X(t) | Y(t), Y(t-1),...,Y(0)]};
    #  Print the parameter estimates
    parms.egoe <-  c(mu.egoe, ssq.egoe, tsq.egoe, x0.egoe) #  Collect for printing
    parms.egpn <-  c(mu.egpn, ssq.egpn, tsq.egpn, x0.egpn)
    parms.reml <-  c(mu.reml, ssq.reml, tsq.reml, x0.reml)
    names <-  c("mu", "ssq", "tsq", "x0")
    types <-  c("EGOE","EGPN","EGSS-ML","EGSS-REML")

  # Add to dataframe
  pops_data[i,10] <- parms.reml[1]
  pops_data[i,11] <- mu_lo.reml[1,1]
  pops_data[i,12] <- mu_hi.reml[1,1]
  pops_data[i,13] <- Var_mu.reml[1,1]
  pops_data[i,14] <- parms.reml[2]
  pops_data[i,15] <- parms.reml[3]
  pops_data[i,16] <- parms.reml[4]
  pops_data[i,17] <- min(Time.t)
  pops_data[i,18] <- max(Time.t)

}

colnames(pops_data)[10:18]<- c("mu", "lCI", "uCI", "var", "sigma",
                               "tau", "x0", "Start", "End")

# Correct the dataset

pops_data <- pops_data %>%
  drop_na(Species) %>%
  mutate(System=ifelse(System=="terrestrial",
                       "Terrestrial",
                       System),
         Protected=gsub("Both ", "Yes", Protected),
         Protected=gsub("Both", "Yes", Protected),
         Protected=gsub("Yes ", "Yes", Protected),
         Protected=gsub(" (large survey area)", "", Protected, fixed = T),
         Protected=gsub("-|\\s+|^$", NA, Protected),
         Data_source_type=ifelse(Data_source_type=="Unpublished report"|
                                   Data_source_type=="Government report"|
                                   Data_source_type=="Secondary source"|
                                   Data_source_type=="grey literature",
                                 "Grey literature", Data_source_type),
         Data_source_type=ifelse(Data_source_type=="Book"|
                                   Data_source_type=="Symposium"|
                                   Data_source_type=="Journal",
                                 "Peer-review", Data_source_type),
         Class=ifelse(Class=="Osteichthyes"|Class=="Actinopterygii",
                      "Bony fishes", Class),
         Class=ifelse(Class=="Chondrichthyes",
                      "Elasmobranchii", Class)) %>%
  filter(Class!="Malacostraca",
         Class!="Cephalopoda",
         Class!="Bivalvia")

# Save the data ----

setwd(DataPath)
save(pops_data, file = "PopChangeData.RData")

# # Load the Dennies and Ponciano 2014 functions 
# setwd(CodePath)
# source("State-space modelling.R")
# 
# # Run a parallel loop to estimate the EGSS and the OU models and their fit
# # Detect the number of cores
# 
# n.cores <- parallel::detectCores() - 2
# 
# # Create the cluster
# 
# my.cluster <- parallel::makeCluster(n.cores, 
#                                     type = "PSOCK")
# 
# #register it to be used by %dopar%
# 
# doParallel::registerDoParallel(cl = my.cluster)
# 
# # Run the foreach loop
# 
# ss_data <- foreach(i = 1:length(pop), 
#                    .combine='rbind', 
#                    .multicombine=TRUE,
#                    .packages = "MASS") %dopar% {
#                      PBLRT.ddpcont(B=10,
#                                    Yobs=pop[[i]][[3]],
#                                    Tvec=as.numeric(pop[[i]][[2]]),
#                                    alpha=0.05,
#                                    plot.PBLR=FALSE)
#                    }
# 
# # Stop the cluster 
# 
# stopCluster(my.cluster)
# 
# # Keep the egss models
# 
# trends <- ss_data %>% 
#   mutate(mu=egss.mu,
#          tau=egss.tau,
#          sigma= egss.sigma) %>% 
#   dplyr::select(mu, tau, sigma) 
# 
# # Join it
# 
# pops_data <- pops_data %>% 
#   cbind(trends) %>% 
#   drop_na(mu)
# 
# # Estimate the start and end of the time series
# 
# dates <- pops %>% 
#   pivot_longer(cols = 11:106,names_to = "Year", values_to = "Count") %>% 
#   drop_na(Count) %>% 
#   group_by(ID) %>% 
#   mutate(Start=min(Year), 
#          End=max(Year)) %>%
#   ungroup() %>% 
#   distinct(ID,.keep_all = T) %>% 
#   dplyr::select(ID, Start, End)
# 
# # Join with the previous dataset 
# 
# pops_data <- pops_data %>% 
#   left_join(dates)
