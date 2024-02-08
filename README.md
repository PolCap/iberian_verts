# The importance of locally sourced data in identifying population trends: insights from Iberian vertebrates 

Roberto C. Rodríguez-Caro<sup>1</sup>,<sup>2</sup>, Zebensui Morales-Reyes<sup>3</sup>, Alba Aguión<sup>1</sup>, Rebeca Arias-Real, Eneko Arrondo, Eneko Aspillaga, Jordi Boada, Andrea Campos-Candela, Mónica Expósito-Granados, Aitor Forcada, Robin Freeman, Miguel Ángel Gómez-Serrano, Cayetano Gutiérrez-Cánovas, Roberto Pascual-Rico, Louise McRae, Maria Montseny, Andreu Rotger, Graciel·la Rovira, Amalia Segura, Iván Sola, Carlos Valle, Pol Capdevila

<sup>1</sup>School of Biological Sciences, University of Bristol, 24 Tyndall Ave, BS8 1TQ, Bristol, UK. 

<sup>2</sup>Institute of Zoology, Zoological Society of London, Regent’s Park, London NW1 4RY, UK.

<sup>3</sup>Ecology and Evolutionary Biology, School of Biosciences, University of Sheffield, Sheffield, UK.

#### Contact: pcapdevila.pc[at]gmail.com

---

## Abstract

_Anthropogenic threats are reshaping Earth’s biodiversity at an unprecedented pace and scale. Conservation policies aiming to halt the loss of biodiversity often rely on global rankings of threats based on their prevalence, where habitat loss and exploitation are often considered the most important. However, these global assessments rarely quantify the impacts of single or multiple threats, which could mask the true effects of the Anthropocene. Here, we quantify the impacts of threats by analysing 3,415 vertebrate populations worldwide, where information about the threats affecting these populations is known. We show that, despite being the most prevalent threats, habitat loss and exploitation are not causing the most rapid population declines. Instead, less prevalent threats, such as invasive species, diseases, and climate change have the strongest negative impacts. However, habitat loss and exploitation have the strongest impacts when acting in conjunction with other threats, even if their cumulative effects are mostly additive (equal to the sum of their individual effects). Finally, we show that mitigating the effects of climate change is as important as mitigating exploitation and/or habitat loss in slowing population declines. These results imply that the most prevalent threats from previous reports are not necessarily the most impactful at local scales and that consequences of climate change on vertebrate populations may have been underestimated._

---

## Data

- __`FiveData.RData`__: contains the population trends of 2,948 vertebrate time-series from the state-space models. 

---

# Code

To run the statistical analyses we used different R scripts: 

- __`AnalysisSSData.R`__: code to clean and prepare the state-space data.
- __`CompareYears.R`__: code to compare the results from the models using 5, 10 and 20 data points.
- __`Counterfactual tests.R`__: code to simulate the counterfactual scenarios, as well as to generate figures 4 and S12.
- __`Figures.R`__: code to create the figures 1-3 and tables S1-S6 of the study. 
- __`ModelDiagnostics.R`__: code to perform the diagnostics of the multilevel Bayesian models. 
- __`MultiEffects.R`__: code to calculate the different proportions of the cummulative effects. 
- __`SensitivityInteractions.R`__: code to calculate the sensitivity of the cummulative effects analyses to the null model choice. 
- __`State-SpaceModels.R`__: code to estimate the state-space models and their diagnostics.
- __`StepTwoModels.R`__: code to apply multilevel Bayesian models to estimate the influence of threats, System and Taxa on population trends. These analyses require some computing power.
- __`SupFigures.R`__: code to create the supplementary figures and tables.
- __`TwentyAnalyses.R`__: code to replicate the multilevel analyses but with 10 and 20 data points. 

---

# Software

_R version 4.0.2 or greater_

To download `R`, go to https://www.r-project.org and for `RStudio`, go to https://www.rstudio.com/products/rstudio/download/ .
 
