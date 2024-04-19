# The importance of locally sourced data in identifying population trends: insights from Iberian vertebrates 

Roberto C. Rodríguez-Caro<sup>1</sup>,<sup>2</sup>, Zebensui Morales-Reyes<sup>3</sup>, Alba Aguión<sup>4</sup>, Rebeca Arias-Real<sup>5</sup>, Eneko Arrondo<sup>6</sup>, Eneko Aspillaga<sup>7</sup>, Jordi Boada<sup>8</sup>, Andrea Campos-Candela<sup>9</sup>, Mónica Expósito-Granados<sup>10</sup>, Aitor Forcada<sup>11</sup>, Robin Freeman<sup>12</sup>, Miguel Ángel Gómez-Serrano<sup>13</sup><sup>14</sup>, Cayetano Gutiérrez-Cánovas<sup>15</sup>, Roberto Pascual-Rico<sup>16</sup>, Louise McRae<sup>12</sup>, Maria Montseny<sup>17</sup><sup>18</sup>, Andreu Rotger<sup>19</sup>, Graciel·la Rovira<sup>17</sup><sup>18</sup>, Amalia Segura<sup>20</sup>, Iván Sola<sup>21</sup>, Carlos Valle<sup>22</sup>, Pol Capdevila<sup>17</sup><sup>18</sup>

<sup>1</sup> Departamento de Ecología, Universidad de Alicante, San Vicent del Raspeig, 03690, Alicante, Spain
<sup>2</sup> Department of Biology, University of Oxford, 11a Mansfield Road, OX1 3SZ, Oxford, UK., 
<sup>3</sup> Instituto de Estudios Sociales Avanzados (IESA), CSIC, Plaza Campo Santo de los Mártires, 7, 14004 Córdoba, Spain.
<sup>4</sup> Coasts and Commons Co-Lab, Nicholas School of the Environment, Duke University, USA.
<sup>5</sup> Museo Nacional de Ciencias Naturales, Consejo Superior de Investigaciones Científicas, Serrano 115 bis, 28006, Madrid, Spain.
<sup>6</sup> Department of Zoology, University of Granada, Spain.
<sup>7</sup> Instituto Mediterráneo de Estudios Avanzados (IMEDEA, CSIC-UIB), 07190 Esporles (Illes Balears), Spain. 
<sup>8</sup> Centre d’Estudis Avançats de Blanes (CEAB-CSIC), 17300 Blanes, Girona.
<sup>9</sup> Department of Biological Sciences, University of Bergen, Bergen 5020, Norway.
<sup>10</sup> Department of Biology, University of Cádiz, Cádiz, Spain. 
<sup>11</sup> Department of Marine Sciences and Applied Biology, University of Alicante, Alicante, Spain.
<sup>12</sup> Institute of Zoology, Zoological Society of London, Regent's Park, London, NW1 4RY UK
<sup>13</sup> Department of Microbiology and Ecology. Faculty of Biological Sciences. University of Valencia. E-46100 Burjassot, Valencia, Spain.
<sup>14</sup> Servicio de Vida Silvestre y Red Natura 2000-Generalitat Valenciana, Spain.
<sup>15</sup> Biodiversity and Conservation Area, Rey Juan Carlos University, Tulipán s/n, 28933 Madrid, Spain.
<sup>16</sup> Grupo SaBio, Instituto de Investigación en Recursos Cinegéticos (IREC), UCLM-CSIC-JCCM, Ciudad Real, Spain.
<sup>17</sup> Departament de Biologia Evolutiva, Ecologia i Ciències Ambientals, Facultat de Biologia, Universitat de Barcelona (UB), Barcelona 08028, Spain.
<sup>18</sup> Institut de Recerca de la Biodiversitat (IRBio), Universitat de Barcelona (UB), Barcelona 08028, Spain.
<sup>19</sup> Animal Demography and Ecology Unit, GEDA – IMEDEA (CSIC/UIB), c. M Marques, 21, 01790, Esporles, Spain.

---

## Abstract

_The Anthropocene is driving profound changes in global biodiversity. However, studies are reporting contrasting biodiversity trends over time, including increases, declines and static temporal trends. These discrepancies can be partly attributed to biases in global datasets, which might not capture the representativeness of regional and local processes, preventing a real identification of population trends. In this study, we aimed to address this gap of knowledge by complementing data included in the Living Planet Database (LPD), one of the largest repositories of vertebrate population time-series, with locally sourced data from the Iberian Peninsula. The study had two main objectives: (i) to assess the state of wildlife populations in the Iberian Peninsula and (ii) to determine if population trends in locally sourced data could differ from those found in LPD data. To supplement the population time-series from the LPD, we conducted an extensive review, analysing over 5,000 peer-reviewed manuscripts and grey literature documents, including population databases from vertebrate monitoring programs. After applying quality checks, we obtained 999 population time-series for 294 vertebrate species compiled in an Iberian Vertebrate (IbeV) database, twice compared to the previous information in the LPD in the Iberian peninsula (~500). Our results indicate contrasting population trends across taxonomic groups, with freshwater amphibians and bony fishes showing the steepest declines, and birds and mammals exhibiting positive trends. Both datasets show contrasting overall trends of Iberian vertebrate populations, where the LPD shows a positive trend and IbeV indicates no net change over time. In addition, the population trends from IbeV differ from the LPD for some taxonomic groups and systems. For both datasets, peer-reviewed manuscripts and grey literature yielded similar trends. Threatened species did not exhibit net changes in population trends, while non-threatened species showed positive population trends. We showed that local databases can provide distinct population trends compared to global databases, making evident the necessity of incorporating local data. This approach highlights the need to bridge the gap between global and local datasets, to support context-specific management and conservation programmes._ 
_

---

## Data

- __`PopChange.RData`__: contains the population trends of 1,543 time series for 430 vertebrate species calculated with the state-space models. 

---

# Code

To run the statistical analyses we used different R scripts: 

- __`FinalDataset.R`__: code to clean and prepare the state-space data.
- __`Figures.R`__: code to create the figures 1-3 and tables S1-S6 of the study. 
- __`ModelDiagnostics.R`__: code to perform the diagnostics of the multilevel Bayesian models. 
- __`Population trends.R`__: code to estimate the state-space models and their diagnostics.
- __`Models.R`__: code to apply multilevel Bayesian models. These analyses require some computing power.
- __`SupFigus.R`__: code to create the supplementary figures and tables.

---

# Software

_R version 4.2.3 _

To download `R`, go to https://www.r-project.org and for `RStudio`, go to https://www.rstudio.com/products/rstudio/download/.
 
