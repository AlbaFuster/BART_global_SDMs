# Machine learning applied to global scale species distribution models (SDMs)

<p align="justify">
In this repository we present the code to reproduce our manuscript "Machine learning applied to global scale species distribution models (SDMs)". Our goal here is to apply BART on a global scale for the estimation and prediction of spatial-temporal distributions of different species and their relationship with environmental variables as a first step to eventually driving a global MEM with species distributions predicted by BART. Also, with these models users are able to obtain results for their target species in a simple and up-to-date way.
<p align="justify">
The available code allows to 1) download presence data from GBIF and environmental variables data from ISIMIP: 2) apply BART models for historical scenarios and future climate change scenarios; and 3) simulate space and time scenarios of presence/absence data to test the properties of the model. 
  
## Repository's structure

<p align="justify">

1. **Case Study:**
This folder with all the code to reproduce the case study for marine turtles, follow these steps:

  1.1 **functions.R:** This script contains all the functions necessary for the study. **Warning: Load the 'functions.R' script before any other script to ensure proper functionality.**

  1.2. **GBIF_data:** Utilize the script `gbif_data.R` to download and clean the presence data. An accompanying HTML file illustrates the functionality of this script.

  1.3. **ISIMIP_data:** Employ the script `isimip_data.R` to download environmental data. An associated HTML file provides an illustration of how this script works.

  1.4. **models_bart.R:** This script includes all the code for applying BART models on a global scale, covering model fitting, prediction, and validation. 


  2. **Simulation Study:**
This folder contains all the code necessary to reproduce the simulation study for spatio-temporal scenarios:

  2.1. **cosmopolitan:** This folder includes scripts required to replicate the simulation of a cosmopolitan species' behavior. The process involves:
       
       1) Running the script script_simulation.R.
       
       2) Executing the models using the script model_bart.R.
       
       3) Visualizing the results by running simulation.rmd.

  2.2. **persistent:** This folder contains scripts essential for reproducing the simulation of a persistent species' behavior. The procedure is identical to the previous one explained.


3. Session information
   
 ```
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)
 ```


