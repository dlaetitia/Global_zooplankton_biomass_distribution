#
# Create barplots figures of the groups' biomass
#
# (c) 2022 L Drago, T Panaïotis, J-O Irisson, R Kiko, GNU General Public License v3

source("Scripts/0.2.config.R")
library("chroma")
library("grid")
library("strex")
library("sf")
library("joml")
library("RColorBrewer")
library("patchwork")
library("ggpubr")

R2_biomass_table <- read_csv("results/",config$region,"/R2_biomass_table.csv") 
cor_data <- read_csv("results/",config$region,"/cor_data.csv")

# Global biomass 
biomass_table <- R2_biomass_table %>% select(-c(contains("R2"))) 

for (lay in 1:nrow(layers)) {
  start_layer_lay = layers$start_layer[lay]
  end_layer_lay = layers$end_layer[lay] 
  
  # get the correlation results 
  cor_data_world <- cor_data %>%  filter(Region == "world")
  cor_data_outside <- cor_data %>%  filter(Region == "region_outside")
  cor_data_inside <- cor_data %>%  filter(Region == "region_inside")
  
  # world
  biom <- biomass_table %>% filter(start_layer == start_layer_lay & end_layer == end_layer_lay) 
  biom <- biom %>% group_by(taxa) %>% 
    mutate(biomass_world_tonC_mean = round(biomass_world_tonC_mean, digits = 2),
           biomass_world_tonC_sd = round(biomass_world_tonC_sd, digits = 2),
           
           biomass_world_PgC_mean = biomass_world_tonC_mean/(10^9),
           biomass_world_PgC_sd = biomass_world_tonC_sd/(10^9),
           biomass_outside_PgC_mean = round(biomass_outside_tonC_mean, digits = 2)/(10^9),
           biomass_outside_PgC_sd = round(biomass_outside_tonC_sd, digits = 2)/(10^9),
           biomass_inside_PgC_mean = round(biomass_inside_tonC_mean, digits = 2)/(10^9),
           biomass_inside_PgC_sd = round(biomass_inside_tonC_sd, digits = 2)/(10^9))
  # Change names of columns
  colnames(biom) <- gsub("biomass_", "", colnames(biom))

  R2_world <- R2_biomass_table %>%  filter(start_layer == start_layer_lay & end_layer == end_layer_lay) %>% 
    select(taxa, start_layer, end_layer, R2_world_mean) %>% 
    rename(R2 = R2_world_mean) %>% 
    mutate(R2=R2/100) %>% 
    left_join(cor_data_world, by = c("taxa", "start_layer","end_layer"))
  
  # get the data for the world models outputs
  biom_world <- biom %>% select(taxa, start_layer, end_layer, world_PgC_mean, world_PgC_sd) %>% 
    filter(!is.na(world_PgC_mean)) %>% 
    arrange(desc(world_PgC_mean)) %>% 
    mutate(world_PgC_mean_cum = cumsum(world_PgC_mean)) %>% 
    arrange(world_PgC_mean_cum) %>% 
    mutate(rn = row_number())  %>% 
    left_join(R2_world)
  write_excel_csv(biom_world, paste0("results/",config$region,"/biomass_tot_1x1/biom_world_",start_layer_lay,"_",end_layer_lay,"m.csv"))
  
  biom_layer <- biom_world %>%  filter(stars == "*")
  biom_layer = sum(biom_layer$world_PgC_mean)
  biom_layer_boundaries <- biom_world  %>%  filter(stars == "*") %>%
    mutate(min_world_PgC = world_PgC_mean - world_PgC_sd,
           max_world_PgC = world_PgC_mean + world_PgC_sd)
  biom_layer_mean = sum(biom_layer_boundaries$world_PgC_mean)
  biom_layer_mean
  biom_layer_bound_low = sum(biom_layer_boundaries$min_world_PgC)
  biom_layer_bound_low
  biom_layer_bound_high = sum(biom_layer_boundaries$max_world_PgC)
  biom_layer_bound_high
  
  if (start_layer_lay == 0 & end_layer_lay == 200) {
    biomass_order_world_0_200 <- biom_world
    write_excel_csv(biomass_order_world_0_200, paste0("results/",config$region,"/biomass_tot_1x1/biomass_order_world_0_200.csv"))
  }
  
  # Save the names of the taxa to plot the good order on x axis
  taxa_order <- as.character(biom_world$taxa)
  
  # choose the sizes for the figures
  s_t = 8 # Size text title
  s_a = 5 # Size text axis
  s_lt = 8 # Size text title legend
  s_la = 6 # Size text axis legend
  s_l = 0.4 # Size legend in cm
  
  #rename some groups for the figure
  biom_world['taxa'][biom_world['taxa'] == 'Collodaria_others'] <- 'Solitary Collodaria'
  biom_world['taxa'][biom_world['taxa'] == 'Collodaria_colonial'] <- 'Colonial Collodaria'
  biom_world['taxa'][biom_world['taxa'] == 'Hydrozoa_others'] <- 'other Hydrozoa'
  biom_world['taxa'][biom_world['taxa'] == 'Mollusca_others'] <- 'other Mollusca'
  biom_world['taxa'][biom_world['taxa'] == 'Rhizaria_others'] <- 'other Rhizaria'
  biom_world['taxa'][biom_world['taxa'] == 'Cnidaria_others'] <- 'other Cnidaria'
  biom_world['taxa'][biom_world['taxa'] == 'Crustacea_others'] <- 'other Crustacea'
  biom_world['taxa'][biom_world['taxa'] == 'Thecosomata_cavo_or_creseis'] <- 'Thecosomata'
  
  barplot_world <- ggplot(biom_world, aes(x=fct_reorder(taxa, world_PgC_mean), y=world_PgC_mean, fill = R2)) +  
    scale_y_continuous(limits = c(0, 0.12)) +
    geom_bar(stat="identity", alpha=1, width=.8) + coord_flip()  + 
    scale_fill_gradientn(colours = brewer.pal(9,"GnBu"), limits = c(0,1)) +
    geom_errorbar(aes(x=taxa, ymin=world_PgC_mean, ymax=world_PgC_mean+world_PgC_sd), width=0.2, size = 0.2) +
    labs(x= "", y = "Mean biomass PgC") + 
    theme(axis.text = element_text(size=s_a),
          axis.title = element_text(size=s_a),
          plot.title = element_text(size = s_t),
          legend.title = element_text(size=s_lt),
          legend.text = element_text(size=s_la),
          legend.key.size = unit(s_l, 'cm')) +
    geom_text(aes(label=stars), vjust = 0.4, hjust = -0.2, size=2) +
    ggtitle("World")

  # high latitude
  R2_outside <- R2_biomass_table %>%  filter(start_layer == start_layer_lay & end_layer == end_layer_lay) %>% 
    select(taxa, start_layer, end_layer, R2_region_outside_mean) %>% 
    rename(R2 = R2_region_outside_mean) %>% 
    mutate(R2=R2/100) %>% 
    left_join(cor_data_outside, by = c("taxa", "start_layer","end_layer"))
  
  biom_outside <- biom %>% select(taxa, start_layer, end_layer, outside_PgC_mean, outside_PgC_sd) %>% 
    filter(!is.na(outside_PgC_mean)) %>% 
    arrange(desc(outside_PgC_mean)) %>% 
    mutate(outside_PgC_mean_cum = cumsum(outside_PgC_mean)) %>% 
    arrange(outside_PgC_mean_cum) %>% 
    mutate(rn = row_number()) %>% 
    left_join(R2_outside)
  
  # Save the names of the taxa to plot the good order on x axis
  taxa_order <- as.character(biom_outside$taxa)
  
  #rename some groups for the figure
  biom_outside['taxa'][biom_outside['taxa'] == 'Collodaria_others'] <- 'Solitary Collodaria'
  biom_outside['taxa'][biom_outside['taxa'] == 'Collodaria_colonial'] <- 'Colonial Collodaria'
  biom_outside['taxa'][biom_outside['taxa'] == 'Hydrozoa_others'] <- 'other Hydrozoa'
  biom_outside['taxa'][biom_outside['taxa'] == 'Mollusca_others'] <- 'other Mollusca'
  biom_outside['taxa'][biom_outside['taxa'] == 'Rhizaria_others'] <- 'other Rhizaria'
  biom_outside['taxa'][biom_outside['taxa'] == 'Cnidaria_others'] <- 'other Cnidaria'
  biom_outside['taxa'][biom_outside['taxa'] == 'Crustacea_others'] <- 'other Crustacea'
  biom_outside['taxa'][biom_outside['taxa'] == 'Thecosomata_cavo_or_creseis'] <- 'Thecosomata'
  
  barplot_outside <- ggplot(biom_outside, aes(x=fct_reorder(taxa, outside_PgC_mean), y=outside_PgC_mean,fill = R2)) + 
    scale_y_continuous(limits = c(0, 0.12)) +
    geom_bar(stat="identity", alpha=1, width=.8) + 
    scale_fill_gradientn(colours = brewer.pal(9,"GnBu"), limits = c(0,1)) +
    geom_errorbar(aes(x=taxa, ymin=outside_PgC_mean, ymax=outside_PgC_mean+outside_PgC_sd), width=0.2, size = 0.2) +
    labs(x= "", y = "Mean biomass PgC") + coord_flip()  +
    theme(axis.text = element_text(size=s_a),
          axis.title = element_text(size=s_a),
          plot.title = element_text(size = s_t),
          legend.title = element_text(size=s_lt),
          legend.text = element_text(size=s_la),
          legend.key.size = unit(s_l, 'cm')) +
    geom_text(aes(label=stars), vjust = 0.4, hjust = -0.2, size=2) +
    ggtitle("High latitudes") 
  
  
  # low latitude
  R2_inside <- R2_biomass_table %>%  filter(start_layer == start_layer_lay & end_layer == end_layer_lay) %>% 
    select(taxa, start_layer, end_layer, R2_region_inside_mean) %>% 
    rename(R2 = R2_region_inside_mean) %>% 
    mutate(R2=R2/100) %>% 
    left_join(cor_data_inside, by = c("taxa", "start_layer","end_layer"))
  
  biom_inside <- biom %>% select(taxa, start_layer, end_layer, inside_PgC_mean, inside_PgC_sd) %>% 
    # Remove NA data
    filter(!is.na(inside_PgC_mean)) %>% 
    arrange(desc(inside_PgC_mean)) %>% 
    mutate(inside_PgC_mean_cum = cumsum(inside_PgC_mean)) %>% 
    arrange(inside_PgC_mean_cum) %>% 
    mutate(rn = row_number()) %>% 
    left_join(R2_inside)
  
  # Save the names of the taxa to plot the good order on x axis
  taxa_order <- as.character(biom_inside$taxa)
  
  #rename some groups for the figure
  biom_inside['taxa'][biom_inside['taxa'] == 'Collodaria_others'] <- 'Solitary Collodaria'
  biom_inside['taxa'][biom_inside['taxa'] == 'Collodaria_colonial'] <- 'Colonial Collodaria'
  biom_inside['taxa'][biom_inside['taxa'] == 'Hydrozoa_others'] <- 'other Hydrozoa'
  biom_inside['taxa'][biom_inside['taxa'] == 'Mollusca_others'] <- 'other Mollusca'
  biom_inside['taxa'][biom_inside['taxa'] == 'Rhizaria_others'] <- 'other Rhizaria'
  biom_inside['taxa'][biom_inside['taxa'] == 'Cnidaria_others'] <- 'other Cnidaria'
  biom_inside['taxa'][biom_inside['taxa'] == 'Crustacea_others'] <- 'other Crustacea'
  biom_inside['taxa'][biom_inside['taxa'] == 'Thecosomata_cavo_or_creseis'] <- 'Thecosomata'
  
  barplot_inside <- ggplot(biom_inside, aes(x=fct_reorder(taxa, inside_PgC_mean), y=inside_PgC_mean, fill = R2)) +     
    scale_y_continuous(limits = c(0, 0.12)) +
    geom_bar(stat="identity", alpha=1, width=.8) +
    geom_errorbar(aes(x=taxa, ymin=inside_PgC_mean, ymax=inside_PgC_mean+inside_PgC_sd), width=0.2, size = 0.2) +
    scale_fill_gradientn(colours = brewer.pal(9,"GnBu"), limits = c(0,1)) +
    labs(x= " ", y = "Mean biomass PgC") + coord_flip()  +
    theme(axis.text = element_text(size=s_a),
          axis.title = element_text(size=s_a), 
          plot.title = element_text(size = s_t),
          legend.title = element_text(size=s_lt),
          legend.text = element_text(size=s_la),
          legend.key.size = unit(s_l, 'cm')) +
    geom_text(aes(label=stars), vjust = 0.4, hjust = -0.2, size=2) +
    ggtitle("Low latitudes") 
  
  # Barplots
  print(barplot_world) / (barplot_outside) / (barplot_inside)
  ggsave(file=paste0("Figures/",config$region,"/barplot_",start_layer_lay,"_",end_layer_lay,".pdf"), 
         plot = last_plot(),
         width=9.5, height = 16, units = "cm")
}