#
# Explore the result of models
#
# Creation of cor_data table containing 
# - the correlation results from the different types of models for each zone, layer and taxa
#
# (c) 2022 L Drago, T Panaïotis, J-O Irisson, R Kiko, GNU General Public License v3

source("Scripts/0.2.config.R")

# Charge the results from the script before containing the 
d_test_all <- read_csv(paste0("results/",config$region,"/d_test_all.csv"))

## Correlation table
cor_data <- matrix(ncol = 7, 
                   dimnames=list(c(), c(
                     #taxa, zones, layer
                     "taxa","start_layer","end_layer", "Region",
                     #correlation info for the world model
                     "p_value","cor","stars")))
cor_data = data.frame(cor_data)

for (r in 1:nrow(config_region)) {
  for (lay in 1:nrow(layers)) {
    start_layer = layers$start_layer[lay]
    end_layer = layers$end_layer[lay] 
    message(paste0(config_region$zone[r], ": ", start_layer, "-", end_layer, "m"))
    
    # get the correlation results for 
    # the models: world, high latitude, low latitude
    cor_data_temp <- read_csv(paste0("results/",config$region,"/cor_table/",config_region$zone[r],"_",start_layer,"_",end_layer,"m.csv")) %>% 
      select_if(!names(.) %in% c("X1", "p_value_spearman", "cor_spearman")) %>%  # remove the column X1 if it exists
      rename(p_value = p_value_pearson, cor = cor_pearson) %>% 
      filter(!is.na(taxa)) 
    cor_data_temp$stars <- cut(cor_data_temp$p_value, breaks=c(-Inf, 0.05, Inf), label=c("*", ""))
    if (config_region$zone[r] == "outside_40S_40N") {cor_data_temp$Region = "region_outside"}
    if (config_region$zone[r] == "inside_40S_40N") {cor_data_temp$Region = "region_inside"}
    
    
    # the part of the model world called "world_outside" and "world_inside"
    if (config_region$zone[r] == "world") {
      # Get the d_test data
      d_test_lay <- d_test_all %>%  filter(start_layer == layers$start_layer[lay] & end_layer == layers$end_layer[lay] )
      
      # get the list of taxa on which to compute the correlations
      results_dir <- paste0("results/",config$element,"/", config$element, "_world/", start_layer,"_",end_layer,  "_", config$sd_or_no_sd)
      f <- list.files(results_dir, full.names=T) 
      models <- map(f, readRDS)
      names(models) <- basename(f)
      names_groups <- str_before_nth(basename(f), pattern = paste0("_",config_region$zone[r]), n = 1)
      
      # create a temporary table for this zone and depth layer
      cor_temp <- matrix(ncol = 7, 
                         dimnames=list(c(), c("taxa","start_layer","end_layer", "Region","p_value","cor","stars")))
      cor_temp = data.frame(cor_temp)
      counter_temp = 0
      
      # now fill the cor_temp with the test for all groups 
      # first for inside 40N-40S
      for (t in 1:length(names_groups)) {
        counter_temp = counter_temp + 1
        predict_test_set_inside <- d_test_lay %>%  filter(taxa== names_groups[t] & Region == "world_inside")
        test_pearson_inside <- cor.test(predict_test_set_inside$pred_mean, predict_test_set_inside$biomass_count, method = "pearson")
        cor_temp[counter_temp,"taxa"] = names_groups[t]
        cor_temp[counter_temp,"p_value"] = test_pearson_inside$p.value
        cor_temp[counter_temp,"cor"] = test_pearson_inside$estimate
        cor_temp[counter_temp,"Region"] = "world_inside"
        
      }
      # then for outside 40N-40S
      for (t in 1:length(names_groups)) {
        counter_temp = counter_temp + 1
        predict_test_set_outside <- d_test_lay %>%  filter(taxa== names_groups[t] & Region == "world_outside")
        test_pearson_outside <- cor.test(predict_test_set_outside$pred_mean, predict_test_set_outside$biomass_count, method = "pearson")
        cor_temp[counter_temp,"taxa"] = names_groups[t]
        cor_temp[counter_temp,"p_value"] = test_pearson_outside$p.value
        cor_temp[counter_temp,"cor"] = test_pearson_outside$estimate
        cor_temp[counter_temp,"Region"] = "world_outside"
      }
      # If the result is significant, add in the column star, a "*"
      cor_temp$stars <- cut(cor_temp$p_value, breaks=c(-Inf, 0.05, Inf), label=c("*", ""))
      # fill out the rest of the information
      cor_temp$start_layer = start_layer
      cor_temp$end_layer = end_layer
    }
    
    # Add these temporary tables to cor_data
    if (config_region$zone[r] == "world") {cor_data <- rbind(cor_data, cor_data_temp, cor_temp)}
    if (config_region$zone[r] != "world") {cor_data <- rbind(cor_data, cor_data_temp)}
  }
}

# Remove the first line containing NA values
cor_data <- cor_data %>% filter(!is.na(start_layer))

# save this table
write_excel_csv(cor_data, paste0("results/",config$region,"/cor_data.csv"))
