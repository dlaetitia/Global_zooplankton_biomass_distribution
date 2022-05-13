#
# Master script
#
# (c) 2022 L Drago, T Panaïotis, J-O Irisson, R Kiko, GNU General Public License v3

#!/usr/bin/Rscript

# This is your master script. 
# It will allow you to run the scripts according to the arguments you provide
# In the script, "outside_40S_40N" and "inside_40S_40N" models correspond to 
# the models refered as "high latitude" and "low latitude" in the article

# Config files
config_ready <- readline(prompt="Are the config files ready ? (yes or no) : ")
paste0("Option chosen : ", config_ready)
if (config_ready == "no") {
  paste0("Please create the config file : more information in config_parameters_instructions.txt")}

if (config_ready == "yes") {
  paste0("Let's start !")
  
  # Run the model or download and prepare the data ? 
  start_or_prepare <- readline(prompt="Are you ready to run the model or do you need to download and/or prepare the data ? (run or prepare) : ")
  paste0("Option chosen : ", start_or_prepare)
  
  if (start_or_prepare == "prepare") {
    
    ### First, choose some main parameters in a .txt file
    # open config_parameters_instructions.txt for more information
    library("readr")
    group_study <- readline(prompt="What group do you want to study ? (All_C was used in Drago et al., 2022) : ")
    paste0("Option chosen : ", group_study)
    config <- read_delim(paste0("Scripts/config_parameters_",group_study,".txt"), 
                         ";", escape_double = FALSE, trim_ws = TRUE)
    # Extract data
    region <- config$region
    n_cores <- unique(config$n_cores)
    beginning_WOA <- unique(config$beginning_year_WOA)
    end_WOA <- unique(config$end_year_WOA)
    max_depth_WOA <- config$max_depth

    # we want to be able to do the model for one layer and one group 
    # and make the script loop on the depths
    # let's try first for only 0-100m
    layers <- read_delim(paste0("Scripts/config_parameters_",group_study,".txt", ";", 
                                escape_double = FALSE, trim_ws = TRUE))
    
    # Create directory for figures
    dir.create(paste0("Figures/", config$element, "/"), showWarnings=F, recursive=T)
    figure_dir <- paste0("Figures/", config$element, "/")
    
    
    # Step 1: fit gradient boosted trees for all taxa in epi and mesopelagic layers
    config_region <- read_delim(paste0("Scripts/config_parameters_regions_", group_study,".txt"), ";", escape_double = FALSE, trim_ws = TRUE)
    for (r in 1:nrow(config_region)) {
      source("Scripts/1.model-fit_habitat_models.R")
    }
    
    ## Step 2: Explore results and compute biomass at global scale
    
    # 1: Extract data from the model (R2 and total biomass)
    source("Scripts/2.global-1.biomass_table_and_maps.R")
    
    # 2: Extract the predicted biomass 
    source("Scripts/2.global-2.biomass_seen_vs_predicted.R")
    
    # 3: Correlation
    source("Scripts/2.global-3.correlation.R")
    
    # 4: Look at the participation in biomass of each taxa
    source("Scripts/2.global-4.barplots.R")
    
    # 5: Compute the global biomass
    source("Scripts/2.global-5.global_biomass.R")
    
    # 6: Plot environmental data maps
    source("Scripts/2.global-6.env_data.R")
    
    # 6: Compare the latitudinal distribution of our Copepoda 
    # and Rhizaria models with the data from Tara Oceans
    source("Scripts/2.global-7.Tara_data.R")
    
  }
}