This repository contains the code and data used in the manuscript: 

_"Artificial light at night and invasive signal crayfish alter aquatic-terrestrial food webs"_

The repository include R function and visualization scripts adapted from `MixSIAR` and `SIBER` packages for
Bayesian mixing models and isotopic niche width estimation using carbon and nitrogen stable isotope data.

**Code description**

To ensure code display and runs as originally configured, download and extract the `.zip` file, then open the project by double-clicking the `ALAN_crayfish_SIA.Rproj` file in `RStudio`. Opening the `.Rproj` file will automatically set the correct working directory and project environment.

The `code` folder include the main scripts and should be run in the following order:

1. `00_Septup.R` - Initialize the environment, loads required packages, and run the function
2. `01_SIA_model_visualization.R` – Prepare the data, runs the Bayesian mixing models and produces visualization 

Additional script are included for niche width estimation, plotting and summarizing results
including other code and output for emergence models and visualization, macroinvertebrate, environmental parameter, and light intensity visualization.


**Folder Structure **

`Code/`

## Contains:
- Setup and model scripts
- Niche widht estimation and visualization code
- Other analysis script and output

`Data/`
Contains input and supporting datasets:
- Mixture data for spider and crayfish
- Source data for spider and crayfish
- Trophic enrichment factor data for spider and crayfish
- Other data_output subfolder: Aquatic insect emergence, macroinvertebrate community, environmental parameter, and light intensity data (.rds) and output with generated plots 

`Output/`
- Spider polygon and diet proportion plots
- Crayfish polygon and diet proportion plots
- Spider trace plots for model diagnostics
- Crayfish trace plots for model diagnostics

`rds data/`
Contains the pre-run model file (.rds) for spider and crayfish that can be loaded directly without 
rerunning the full model

**Running the Code**
1. Download the repository and open it using the .Rproj file.
2. Restore the required R environment by running:
   ```
   renv::restore()
   
   ```
   




