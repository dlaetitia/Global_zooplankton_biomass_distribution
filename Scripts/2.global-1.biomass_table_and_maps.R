#
# Explore the result of models
#
# Creation of biomass_table containing 
# - the R2 from the different types of models for each zone, layer and taxa
# - the total biomass for each zone, layer and taxa
#
# Takes around 30min  
#
# (c) 2022 L Drago, T Panaïotis, J-O Irisson, R Kiko, GNU General Public License v3

source("Scripts/0.2.config.R")
library("chroma")
library("grid")
library("strex")
library("sf")
library("RColorBrewer")
library("gridExtra")
library("patchwork")
library("ggpubr")

# Create an empty table with all the layers and taxa to be filled later
list_taxa <- read_csv("Data/",config$element,"/taxon_list_models.csv")
l_taxa <- list_taxa %>%  rename("taxa" = "taxon")
l_taxa <- l_taxa[rep(seq_len(nrow(l_taxa)), times = nrow(layers)), ]
l_layers <- layers 
l_layers <- l_layers[rep(seq_len(nrow(l_layers)), each = nrow(list_taxa)),]
taxa_table <- cbind(l_taxa, l_layers)


# Create an empty table to add the metrics and total biomass values for each situation
biomass_table <- matrix(ncol = 26, 
                        dimnames=list(c(), c(
                          #zones, layer, taxa
                          "region","start_layer","end_layer", "taxa",
                          #R2 for the world model
                          "R2_world_mean","R2_world_sd","R2_world_se","R2_world_cv",
                          # R2 for the regional model
                          "R2_region_mean","R2_region_sd","R2_region_se","R2_region_cv",
                          # R2 for outside 40-40 from the world model
                          "R2_world_outside_mean","R2_world_outside_sd","R2_world_outside_se","R2_world_outside_cv",
                          # R2 for inside 40-40 from the world model
                          "R2_world_inside_mean","R2_world_inside_sd","R2_world_inside_se","R2_world_inside_cv",
                          # Biomass tot for each zone, layer, taxa
                          "biomass_tot_tonC_mean","biomass_tot_tonC_sd","biomass_tot_tonC_se","biomass_tot_tonC_cv", "pred_mean", "pred_sd")))
biomass_table=data.frame(biomass_table)
counter = 0

# preparation for computation of total biomass for each group
# compute the area of a 1º pixel at each latitude
lats <- -89.5:89.5
# Earth radius
R = 6378.137 # km # https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
# Compute area in m2, Cf.
# https://fr.mathworks.com/help/map/ref/areaquad.html
# https://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2004/msg00023.html
lat1 <- (lats - 0.5) * pi / 180
lat2 <- (lats + 0.5) * pi / 180
area <- ((pi/180)*R^2 * abs(sin(lat1)-sin(lat2)))*10^6
lat_area <- tibble(lat=lats, area)

dir.create(paste0("results/",config$region,"/biomass_tot_1x1"), showWarnings=F, recursive=T)
dir.create(paste0("results/",config$region,"/predicted_vs_true"), showWarnings=F, recursive=T)
dir.create(paste0("results/",config$region,"/importance"), showWarnings=F, recursive=T)

# get the variables names from the models
load(paste0("Data/",config$region,"/3.biomass_env_world_0_200m_",config$layers,".Rdata"))
data_matrix <- d %>% select(-lon, -lat, -taxon, -c(contains("sd")), -c(contains("biomass")), -c(contains("AOU")))
vars <- colnames(data_matrix)


# read all model results
for (r in 1:nrow(config_region)) {
  
  # create a folder to store the global biomass on 1°x1° grid
  dir.create(paste0("results/",config$region,"/biomass_tot_1x1/",config_region$zone[r]), showWarnings=F, recursive=T)
  dir.create(paste0("results/",config$region,"/predicted_vs_true/",config_region$zone[r]), showWarnings=F, recursive=T)
  
  for (lay in 1:nrow(layers)) {
    start_layer = layers$start_layer[lay]
    end_layer = layers$end_layer[lay]
    message(paste0("Global maps and partial dependance plots for: ", config_region$zone[r], " ", start_layer, "-", end_layer, "m"))
    
    dir.create(paste0("results/",config$region,"/biomass_tot_1x1/",config_region$zone[r], "/", start_layer,"_",end_layer, "_", config$sd_or_no_sd), showWarnings=F, recursive=T)
    dir.create(paste0("results/",config$region,"/predicted_vs_true/",config_region$zone[r], "/", start_layer,"_",end_layer,  "_", config$sd_or_no_sd), showWarnings=F, recursive=T)
    results_dir <- paste0("results/",config$element,"/", config$element, "_", config_region$zone[r], "/", start_layer,"_",end_layer,  "_", config$sd_or_no_sd)
    
    f <- list.files(results_dir, full.names=T) 
    models <- map(f, readRDS)
    names(models) <- basename(f)
    names_groups <- str_before_nth(basename(f), pattern = paste0("_",config_region$zone[r]), n = 1)
    
    ## Extract prediction statistics on test set ----
    message("Extract prediction statistics")
    # fits of the models
    fits <- map_dfr(models, function(l) {data.frame(taxa = l$t, l$cv_metrics)})
    # their predictions
    prediction_all <- map_dfr(models, function(l) {data.frame(taxa = l$t, l$prediction)})
    
    # For the world we need
    if (config_region$zone[r] == "world") {
      # the R2 for the whole world
      fits_summary <- fits %>%  group_by(taxa) %>%  
        summarise(R2_world_mean = mean(R2_correl), R2_world_sd = sd(R2_correl), 
                  R2_world_se = joml::se(R2_correl), R2_world_cv = R2_world_sd/R2_world_mean)
      # the R2 computed only inside or outside 40°N-40°S band
      fits_summary_inside <- map_dfr(models, function(l) {data.frame(taxa = l$t, l$cv_metrics_inside)}) %>% 
        group_by(taxa) %>%  
        summarise(R2_world_inside_mean = mean(R2_correl), R2_world_inside_sd = sd(R2_correl), 
                  R2_world_inside_se = joml::se(R2_correl), R2_world_inside_cv = R2_world_inside_sd/R2_world_inside_mean)
      fits_summary_outside <- map_dfr(models, function(l) {data.frame(taxa = l$t, l$cv_metrics_outside)}) %>% 
        group_by(taxa)%>%  
        summarise(R2_world_outside_mean = mean(R2_correl), R2_world_outside_sd = sd(R2_correl), 
                  R2_world_outside_se = joml::se(R2_correl), R2_world_outside_cv = R2_world_outside_sd/R2_world_outside_mean)
    }
    
    # For the regional models, we only need the R2 values from the regional models
    if (config_region$zone[r] != "world") {
      fits_summary <- fits %>%  group_by(taxa) %>%  
        summarise(R2_region_mean = mean(R2_correl), R2_region_sd = sd(R2_correl), 
                  R2_region_se = joml::se(R2_correl), R2_region_cv = R2_region_sd/R2_region_mean)
    }
    
    ## Plot report sheet for all models ----
    message(paste0("Plot report sheets for ",config_region$zone[r], ": ", start_layer, "-",end_layer,"m"))
    # Save in a PDF file
    dir.create(paste0("Figures/",config$element,"/", config$element, "_", config_region$zone[r], "/"), showWarnings=FALSE)
    figure_dir <- paste0("Figures/",config$element,"/", config$element, "_", config_region$zone[r], "/")
    openfig(paste0("models_plots_", config_region$zone[r], "_", start_layer,"_",end_layer,"m", "_", config$sd_or_no_sd,".pdf"), width=18, height=20)
    
    for (t in 1:length(names_groups)) {
      counter = counter + 1
      prediction_t <- prediction_all %>%  filter(taxa == names_groups[t])
      models_t <- fits %>%  filter(taxa == names_groups[t])
      
      # get R2 values
      value_R2 = fits_summary %>%  filter(taxa == names_groups[t]) 
      value_R2 <- round(value_R2[2], digits = 1)
      
      map_mean <- ggplot(prediction_t) + 
        geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
        geom_raster(aes(lon, lat, fill=pred_mean)) +
        scale_fill_viridis_c(trans="sqrt") +
        coord_quickmap() + scale_xy_map() +
        theme_dark() +
        labs(fill="Predicted\nbiomass\n(mgC/m3)\n", title=paste0(names_groups[t], " ", start_layer, "-", end_layer, "m (R2=",round(value_R2/100, digits = 2),")"))
      
      map_cv <- ggplot(prediction_t) + 
        geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
        geom_raster(aes(lon, lat, fill=pred_cv)) +
        scale_fill_distiller(palette="RdYlGn") +
        coord_quickmap() + scale_xy_map() +
        theme_dark() +
        labs(fill="Coefficient\nof variation")
      
      # Get the models for the group t
      models_t <- fits %>%  filter (taxa == names_groups[t])
      f_t <- list.files(results_dir, pattern= paste0(names_groups[t]),full.names=T)
      models_t <- map(f_t, readRDS)
      m_best <- map_dfr(models_t, function(l) {l$model_full})
      
      # compute the importance of variables for each of the 100 models
      imp <- xgbr.importance(m_best)
      imp_data <- plot.xgbr.importance(imp)
      # save this information
      imp_info <- data.frame(imp_data[[1]])
      write_excel_csv(imp_info, file = paste0("results/",config$region,"/importance/Importance_variables_", names_groups[t],"_",config_region$zone[r],"_",start_layer,"_",end_layer,".csv"))
      # get the order of variables
      imp_order <- as.character(data.frame(imp_data[[1]])$Feature)
      # plot the according to the order of variables' importance
      imp_plot <- imp_data[[2]]

      plot(map_mean/map_cv/imp_plot)
      
      # choose the name you want
      openfig(paste0("One_group_",names_groups[t], "_plots_", config_region$zone[r], "_", start_layer,"_",end_layer,"m", "_", config$sd_or_no_sd,".pdf"), width=18, height=20)
      # pdp
      system.time(
        m_best <- partials(m_best, expl = c(imp_order[1:3]), cores = 50, quantiles=TRUE, probs=0:10/10)
      )
      partial_plot <- plot_pdp(m_best)
      grid.arrange(map_mean, map_cv,arrangeGrob(imp_plot,partial_plot,widths=c(3/10,7/10),ncol=2),ncol=1)
      dev.off()
      
      # fill the table of information
      ### context
      biomass_table[counter,"region"] <- config_region$zone[r]
      biomass_table[counter,"start_layer"] <- start_layer
      biomass_table[counter,"end_layer"] <- end_layer
      biomass_table[counter,"taxa"] <- names_groups[t]
      
      ### R2 values
      value_R2 = fits_summary %>%  filter(taxa == names_groups[t]) 
      value_R2 <- round(value_R2[2], digits = 1)
      num <- which(fits_summary$taxa == names_groups[t])
      
      if (config_region$zone[r] != "world") {
        biomass_table[counter,"R2_region_mean"] <- fits_summary[num,"R2_region_mean"] 
        biomass_table[counter,"R2_region_sd"] <- fits_summary[num,"R2_region_sd"] 
        biomass_table[counter,"R2_region_se"] <- fits_summary[num,"R2_region_se"] 
        biomass_table[counter,"R2_region_cv"] <- fits_summary[num,"R2_region_cv"] 
      }
      
      if (config_region$zone[r] == "world") {
        biomass_table[counter,"R2_world_mean"] <- fits_summary[num,"R2_world_mean"] 
        biomass_table[counter,"R2_world_sd"] <- fits_summary[num,"R2_world_sd"] 
        biomass_table[counter,"R2_world_se"] <- fits_summary[num,"R2_world_se"] 
        biomass_table[counter,"R2_world_cv"] <- fits_summary[num,"R2_world_cv"] 
        
        # R2 from cookie cutting inside 40°N-40°S
        biomass_table[counter,"R2_world_inside_mean"] <- fits_summary_inside[num,"R2_world_inside_mean"] 
        biomass_table[counter,"R2_world_inside_sd"] <- fits_summary_inside[num,"R2_world_inside_sd"] 
        biomass_table[counter,"R2_world_inside_se"] <- fits_summary_inside[num,"R2_world_inside_se"] 
        biomass_table[counter,"R2_world_inside_cv"] <- fits_summary_inside[num,"R2_world_inside_cv"] 
        
        # R2 from cookie cutting outside 40°N-40°S
        biomass_table[counter,"R2_world_outside_mean"] <- fits_summary_outside[num,"R2_world_outside_mean"] 
        biomass_table[counter,"R2_world_outside_sd"] <- fits_summary_outside[num,"R2_world_outside_sd"] 
        biomass_table[counter,"R2_world_outside_se"] <- fits_summary_outside[num,"R2_world_outside_se"] 
        biomass_table[counter,"R2_world_outside_cv"] <- fits_summary_outside[num,"R2_world_outside_cv"] 
      }
      
      ### total biomass in tonC and biomass in mgC/m3
      # add area to the data
      prediction_t <- left_join(prediction_t, lat_area,by="lat")
      biomass_tot <- prediction_t %>% 
        group_by(lon, lat) %>% 
        summarise(         # biomass    # depth.                  # scale by area # in tons
          biomass_tonC_mean=(pred_mean * (end_layer - start_layer)) * area          / 10^9,
          biomass_tonC_sd=(pred_sd * (end_layer - start_layer)) * area          / 10^9,
          biomass_tonC_se=(pred_se * (end_layer - start_layer)) * area          / 10^9,
          biomass_PgC_mean = biomass_tonC_mean / 10^9,
          biomass_PgC_sd = biomass_tonC_sd / 10^9,
          biomass_PgC_se = biomass_tonC_se / 10^9,
          pred_mean = pred_mean,
          pred_sd = pred_sd) %>% 
        ungroup()
      
      write_excel_csv(biomass_tot, paste0("results/", config$element,"/biomass_tot_1x1/",config_region$zone[r], "/",start_layer,"_",end_layer , "_", config$sd_or_no_sd,"/", names_groups[t],"_",config_region$zone[r],"_",start_layer,"_",end_layer,"m_biomass_tot_1x1.csv"))

      biomass_tot_summary <- biomass_tot %>% 
        summarise(biomass_tonC_mean = sum(biomass_tonC_mean),
                  biomass_tonC_sd = sum(biomass_tonC_sd),
                  biomass_tonC_se = sum(biomass_tonC_se),
                  biomass_tonC_cv = biomass_tonC_sd/biomass_tonC_mean)
      biomass_table[counter,"biomass_tot_tonC_mean"] <- biomass_tot_summary[1,"biomass_tonC_mean"]
      biomass_table[counter,"biomass_tot_tonC_sd"] <- biomass_tot_summary[1,"biomass_tonC_sd"]
      biomass_table[counter,"biomass_tot_tonC_se"] <- biomass_tot_summary[1,"biomass_tonC_se"]
      biomass_table[counter,"biomass_tot_tonC_cv"] <- biomass_tot_summary[1,"biomass_tonC_cv"]
      
      # put everything in one line for one taxa, one layer and all regions
      biomass_world <- biomass_table %>%  filter(region == "world") %>% select(-contains("region")) 
      biomass_outside <- biomass_table %>% filter(region == "outside_40S_40N") %>% select(-contains("world")) %>% select(-region)
      biomass_inside <- biomass_table %>% filter(region == "inside_40S_40N") %>% select(-contains("world")) %>% select(-region)
      
      colnames(biomass_outside) <- gsub("region", "region_outside", colnames(biomass_outside))
      colnames(biomass_outside) <- gsub("tot", "outside", colnames(biomass_outside))
      
      colnames(biomass_inside) <- gsub("region", "region_inside", colnames(biomass_inside))
      colnames(biomass_inside) <- gsub("tot", "inside", colnames(biomass_inside))
      
      colnames(biomass_world) <- gsub("tot", "world", colnames(biomass_world))
      
      R2_biomass_table <- taxa_table %>% 
        left_join(biomass_world) %>% 
        left_join(biomass_outside) %>% 
        left_join(biomass_inside)
    }
    dev.off()
  }
}

# save it !
write_excel_csv(R2_biomass_table, paste0("results/",config$region,"/R2_biomass_table.csv"))
