#
# Compare the models outputs to the ones of TARA Ocean 
#
# (c) 2022 L Drago, T Panaïotis, J-O Irisson, R Kiko, GNU General Public License v3

source("Scripts/0.2.config.R")
library("patchwork")
library("readxl")

# 1. Get the data
# 1.1 Tara data : 
table_Tara <- read_excel("Data/",config$region,"/table_Tara.xls")
table_Tara <- table_Tara %>% rename(Copepoda_biom = 'Copepod Biomass') %>% 
  rename(Protist_biom = 'Protist Biomass')
table_Tara$Copepoda_biom <- as.numeric(table_Tara$Copepoda_biom)
table_Tara$Protist_biom <- as.numeric(table_Tara$Protist_biom)
table_Tara$Latitude <- as.numeric(table_Tara$Latitude)

# 1.2. Charge the predictions and correlation data from our models
load("results/",config$region,"/every_pred.Rdata")
cor_data <- read_csv("results/",config$region,"/cor_data.csv")

# Now add the info about the correlation for our 3 regions 
cor_d <- cor_data %>%  filter(Region == "world" | Region == "region_outside" | Region == "region_inside") 
# low latitude models results
cor_d['Region'][cor_d['Region'] == 'region_inside'] <- 'inside_40S_40N'
# high latitude models results
cor_d['Region'][cor_d['Region'] == 'region_outside'] <- 'outside_40S_40N'

every_pred <- every_pred %>%  left_join(cor_d) %>% 
  # keep only world data with significant results
  filter(Region == "world" & stars == "*") 

biomass_mgC_m3 <- every_pred %>%  left_join(cor_d)

# 1.3. Raw data as seen by the UVP5
load(paste0("Data/",config$region,"/3.biomass_env_0_200m.Rdata"))
d_0_200 <- d 
remove(d)
load(paste0("Data/",config$region,"/3.biomass_env_200_500m.Rdata"))
d_200_500 <- d
remove(d)

# 2. Plot the data for different groups
## 2.1. COPEPODA
### 2.1.1. Data from UVP5 raw data
d_0_200_cop <- d_0_200 %>%  filter(taxon == "Copepoda") %>%  
  select(lon, lat, biomass_sph) %>% 
  # convert biomass from mgC/m3 -> mgC/m2
  mutate(biomass_200 = biomass_sph*200) %>%  
  select(-biomass_sph)
d_200_500_cop <- d_200_500 %>%  filter(taxon == "Copepoda") %>%  
  select(lon, lat, biomass_sph) %>% 
  mutate(biomass_500 = biomass_sph *300) %>%  
  select(-biomass_sph)
d_0_500_cop <- d_0_200_cop %>%  left_join(d_200_500_cop, by = c("lon", "lat")) %>% 
  group_by(lon, lat) %>% 
  # Add both layers 
  mutate(biomass = biomass_200 + biomass_500) %>% 
  ungroup()

### 2.1.2. Data from models
# Compute the biomass on 200m depth 
biomass_mgC_m3_cop_200 <- biomass_mgC_m3 %>%  
  filter(taxa == "Copepoda", Region == "world") %>% 
  filter(start_layer == 0 & end_layer == 200) %>% 
  mutate(pred_mean_200 = 200*pred_mean) %>% select(-pred_mean)

# Between 200-500
biomass_mgC_m3_cop_500 <- biomass_mgC_m3 %>%  
  filter(taxa == "Copepoda", Region == "world") %>% 
  filter(start_layer ==200 & end_layer == 500) %>% 
  mutate(pred_mean_500 = 300*pred_mean) %>% 
  select(-pred_mean)

# Add them up
biomass_mgC_m3_cop <- biomass_mgC_m3_cop_200 %>%  
  left_join(biomass_mgC_m3_cop_500, by = c("lat", "lon")) %>% 
  group_by(lon, lat) %>% 
  mutate(pred_mean = pred_mean_200 + pred_mean_500) %>% 
  ungroup()

# Plot the comparison
s = 8# size of text in figure
t = 10 # size of title

plot_cop <- ggplot()+
  # Models
  geom_smooth(data = biomass_mgC_m3_cop, aes(y=log10(pred_mean+1), x=lat, col = "#332288"), se=T, alpha=1, size=0.6) +
  # UVP5
  geom_smooth(data = d_0_500_cop, aes(y=log10(biomass+1), x=lat, col = "#44A998"), se=T, alpha=0.5, size=0.6) +
  # Tara
  geom_smooth(data = table_Tara, aes(y=log10(Copepoda_biom+1), x=Latitude, col = "#CC6677"), se=T, alpha=0.5, size=0.6) +
  coord_flip() +
  scale_y_continuous("Biomass 0-500m (mgC.m-2)", breaks=c(0,1,2,3,4,5,6,7), limits = c(0,4.5),
                     labels=c("0","10","100","1 000","10 000","100 000","1 000 000", "10 000 000"))+
  scale_x_continuous("Latitude", limits=c(-80, 90), breaks=seq(-80, 90, 10),
                     labels=c("","","60°S","","","30°S","","","0°","","","30°N","","","60°N","","","90°N"))+
  scale_color_identity(guide = "legend",
                       name = "Trends", 
                       breaks = c("#332288", "#44A998","#CC6677"),
                       labels = c("Our models", "UVP5", "TARA Ocean")) +
  ggtitle("Copepoda") +
  theme(legend.position="none", axis.text = element_text(size=s), 
        axis.title = element_text(size=s), plot.title = element_text(size = t))


## 2.2. PROTISTS : Acantharea, Phaeodaria, Foraminifera
### 2.1.1. Data from UVP5 raw data
d_0_200_rhiz <- d_0_200 %>% 
  filter(taxon == "Acantharea" | taxon == "Phaeodaria" | taxon == "Foraminifera"|
           taxon == "Collodaria_others" | taxon == "Rhizaria_others" | taxon == "Collodaria_colonial") %>% 
  # convert biomass from mgC/m3 -> mgC/m2
  group_by(lon, lat) %>% 
  summarise(biomass_200 = 200*sum(biomass_sph)) %>% 
  ungroup() 

d_200_500_rhiz <- d_200_500 %>%  
  filter(taxon == "Acantharea" | taxon == "Phaeodaria" | taxon == "Foraminifera"|
           taxon == "Collodaria_others" | taxon == "Rhizaria_others" | taxon == "Collodaria_colonial") %>% 
  group_by(lon, lat) %>% 
  summarise(biomass_500 = 300*sum(biomass_sph)) %>%  
  ungroup()

d_0_500_rhiz <- d_0_200_rhiz %>%  left_join(d_200_500_rhiz, by = c("lon", "lat")) %>% 
  group_by(lon, lat) %>% 
  # Add both layers 
  mutate(biomass = biomass_200 + biomass_500) %>% 
  ungroup()

### 2.1.2. Data from models
# Compute the biomass on 200m depth
biomass_mgC_m3_rhiz_200 <- biomass_mgC_m3 %>%  
  filter(Region == "world") %>% 
  filter(taxa == "Acantharea" | taxa == "Phaeodaria" | taxa == "Foraminifera"|
           taxa == "Collodaria_others" | taxa == "Rhizaria_others" | taxa == "Collodaria_colonial") %>% 
  filter(start_layer == 0 & end_layer == 200) %>% 
  # compute the biomass on 200m depth for all 3 groups for each 1°x1° cell
  group_by(lon, lat) %>% 
  summarise(pred_mean_200 = 200*sum(pred_mean)) %>% 
  ungroup() 

# Between 200-500
biomass_mgC_m3_rhiz_500 <- biomass_mgC_m3 %>%  
  filter(Region == "world") %>% 
  filter(taxa == "Acantharea" | taxa == "Phaeodaria" | taxa == "Foraminifera" |
           taxa == "Collodaria_others" | taxa == "Rhizaria_others" | taxa == "Collodaria_colonial") %>% 
  filter(start_layer == 200 & end_layer == 500) %>% 
  # compute the biomass on 200m depth for all 3 groups for each 1°x1° cell
  group_by(lon, lat) %>% 
  summarise(pred_mean_500 = 300*sum(pred_mean)) %>% 
  ungroup() 

# Add them up
biomass_mgC_m3_rhiz <- biomass_mgC_m3_rhiz_200 %>%  
  left_join(biomass_mgC_m3_rhiz_500, by = c("lat", "lon")) %>% 
  group_by(lon, lat) %>% 
  summarise(pred_mean = pred_mean_200 + pred_mean_500) %>% 
  ungroup()

plot_rhiz <- ggplot()+
  # Model
  geom_smooth(data = biomass_mgC_m3_rhiz, aes(y=log10(pred_mean+1), x=lat, col = "#332288"), se=T, alpha=1, size=0.6) +
  # UVP5
  geom_smooth(data = d_0_500_rhiz, aes(y=log10(biomass+1), x=lat, col = "#44A998"), se=T, alpha=0.5, size=0.6) +
  # Tara
  geom_smooth(data = table_Tara, aes(y=log10(Protist_biom+1), x=Latitude, col = "#CC6677"), se=T, alpha=0.5, size=0.6) +
  coord_flip() + 
  scale_y_continuous("Biomass 0-500m (mgC.m-2)", breaks=c(0,1,2,3,4,5,6,7), limits = c(-0.3,4.5),
                     labels=c("0","10","100","1 000","10 000","100 000","1 000 000", "10 000 000"))+
  scale_x_continuous("Latitude", limits=c(-80, 90), breaks=seq(-80, 90, 10),
                     labels=c("","","60°S","","","30°S","","","0°","","","30°N","","","60°N","","","90°N"))+
  scale_color_identity(guide = "legend",
                       name = "Trends", 
                       breaks = c("#332288", "44A998","#CC6677"),
                       labels = c("Our models", "UVP5", "TARA Ocean")) +
  ggtitle("Rhizaria") +
  theme(legend.position="none", axis.text = element_text(size=s), 
        axis.title = element_text(size=s), plot.title = element_text(size = t))

#Get legend for the geom_smooth
plot_smooth <- ggplot()+
  # models
  geom_smooth(data = biomass_mgC_m3_rhiz, aes(y=log10(pred_mean+1), x=lat, col = "#332288"), se=T, alpha=1, size=0.6) +
  # UVP5
  geom_smooth(data = d_0_500_rhiz, aes(y=log10(biomass+1), x=lat, col = "#44A998"), se=T, alpha=0.5, size=0.6) +
  # Tara
  geom_smooth(data = table_Tara, aes(y=log10(Protist_biom+1), x=Latitude, col = "#CC6677"), se=T, alpha=0.5, size=0.6) +
  coord_flip() +
  scale_y_continuous("Biomass 0-500m (mgC.m-2)", breaks=c(0,1,2,3,4,5,6,7),
                     labels=c("0","10","100","1 000","10 000","100 000","1 000 000", "10 000 000"))+
  scale_x_continuous("Latitude", limits=c(-80, 90), breaks=seq(-80, 90, 10),
                     labels=c("","","60°S","","","30°S","","","0°","","","30°N","","","60°N","","","90°N"))+
  scale_color_identity(guide = "legend",
                       name = "Trends", 
                       breaks = c("#332288", "#44A998","#CC6677"),
                       labels = c("BRT models", "UVP5", "TARA Ocean nets")) 
leg_smooth <- cowplot::get_legend(plot_smooth + theme(legend.position = "bottom"))

## Plot only copepoda and rhizaria
## Plot all together
layout_TARA <- "
AAABBB
AAABBB
AAABBB
AAABBB
##CC##
"
plot_cop + plot_rhiz + cowplot::plot_grid(leg_smooth) +
  plot_layout(design = layout_TARA) +
  plot_annotation(title = paste0(),
                  theme = theme(plot.title = element_text(size = 10)))

ggsave(file=paste0("Figures/",config$region,"/TARA_comparison.pdf"), 
       plot = last_plot(),
       width=18, height = 13, units = "cm")
