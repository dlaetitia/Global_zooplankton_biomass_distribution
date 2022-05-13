#
# Explore the result of models
#
# Look at the predicted VS seen biomass on the vertical layers
# Compute the global biomass for each group and for all groups in a layer
#
# (c) 2022 L Drago, T Panaïotis, J-O Irisson, R Kiko, GNU General Public License v3

source("Scripts/0.2.config.R")

## Get data
# table for correlation pearson
cor_data <- read_csv("results/",config$region,"/cor_data.csv") 

# table of prediction data (pred_mean) and data as seen by UVP (biomass_count)
d_test_all <- read_csv("results/",config$region,"/d_test_all.csv") 

all <- left_join(d_test_all, cor_data)
cor_test_all <- matrix(ncol = 6, dimnames=list(c(), c("start_layer","end_layer","p_value","cor", "R2","star")))
cor_test_all = data.frame(cor_test_all)


### 1. Regression between seen biomass and predicted biomass for only groups with pvalue < 0.05 -----
all_p <- all %>%  filter(stars == "*")
cor_test_all[4:6,"star"] = "*"

## 1.A. for 0-200 ------
all_0_200_star  <- all_p %>%  filter(start_layer == 0 & end_layer == 200)

test_0_200_star  <- cor.test(all_0_200_star$pred_mean, all_0_200_star$biomass_count, method = "pearson", exact = FALSE)
cor_test_all[4,"p_value"] = test_0_200_star$p.value
cor_test_all[4,"cor"] = test_0_200_star$estimate
cor_test_all[4,"R2"] = stats::cor(all_0_200_star$pred_mean, all_0_200_star$biomass_count)^2 * 100
cor_test_all[4,"start_layer"] = 0
cor_test_all[4,"end_layer"] = 200

plot_0_200_star <- ggplot(all_0_200_star) + coord_fixed() +  geom_abline(intercept=0, slope=1, colour="#3B87B3") +
  geom_point(aes(biomass_count, pred_mean), size=0.5, alpha=0.5) +
  scale_x_continuous(trans="log1p") + scale_y_continuous(trans="log1p") + # nicer aspect
  labs(x="True biomass", y="Predicted biomass") +
  theme(plot.title = element_text(size = 5), axis.text=element_text(size=4), axis.title=element_text(size=5)) +
  ggtitle(paste0("*0-200m: cor=",round(test_0_200_star$estimate, digits = 2)," pvalue=", test_0_200_star$p.value))


## 1.B. for 200-500 ------
all_200_500_star  <- all_p %>%  filter(start_layer == 200 & end_layer == 500)

test_200_500_star <- cor.test(all_200_500_star$pred_mean, all_200_500_star$biomass_count, method = "pearson", exact = FALSE)
cor_test_all[5,"p_value"] = test_200_500_star$p.value
cor_test_all[5,"cor"] = test_200_500_star$estimate
cor_test_all[5,"R2"] = stats::cor(all_200_500_star$pred_mean, all_200_500_star$biomass_count)^2 * 100
cor_test_all[5,"start_layer"] = 200
cor_test_all[5,"end_layer"] = 500

plot_200_500_star <- ggplot(all_200_500_star) + coord_fixed() + geom_abline(intercept=0, slope=1, colour="#3B87B3") +
  geom_point(aes(biomass_count, pred_mean), size=0.5, alpha=0.5) +
  scale_x_continuous(trans="log1p") + scale_y_continuous(trans="log1p") +   # nicer aspect
  labs(x="True biomass", y="Predicted biomass") +
  theme(plot.title = element_text(size = 5), axis.text=element_text(size=4), axis.title=element_text(size=5)) +
  ggtitle(paste0("*200-500m: cor=",round(test_200_500_star$estimate, digits = 2)," pvalue=", test_200_500_star$p.value))


## 1.C. for 0-500 ------
test_0_500_star <- cor.test(all_p$pred_mean, all_p$biomass_count, method = "pearson", exact = FALSE)
cor_test_all[6,"p_value"] = test_0_500_star$p.value
cor_test_all[6,"cor"] = test_0_500_star$estimate
cor_test_all[6,"R2"] = stats::cor(all_p$pred_mean, all_p$biomass_count)^2 * 100
cor_test_all[6,"start_layer"] = 0
cor_test_all[6,"end_layer"] = 500

plot_0_500_star  <- ggplot(all_p) + coord_fixed() + geom_abline(intercept=0, slope=1, colour="#3B87B3") +
  geom_point(aes(biomass_count, pred_mean), size=0.5, alpha=0.5) +
  scale_x_continuous(trans="log1p") + scale_y_continuous(trans="log1p") +
  labs(x="True biomass", y="Predicted biomass") +
  theme(plot.title = element_text(size = 5), axis.text=element_text(size=4), axis.title=element_text(size=5)) +
  ggtitle(paste0("*0-500m: cor=",round(test_0_500_star$estimate, digits = 2)," pvalue=", test_0_500_star$p.value))

cor_test_all
print(plot_0_200_star + plot_200_500_star + plot_0_500_star)


# print all
plot_all_0_500 <- print((plot_0_200 + plot_200_500 + plot_0_500)/ (plot_0_200_star + plot_200_500_star + plot_0_500_star))
ggsave(file=paste0("results/",config$region,"/predicted_vs_true/comparison_p_value.pdf"),
               plot = plot_all_0_500,
              width=18, height=10, units = "cm")


### 2. Compute global biomass numbers for each 1°x1° cell in PgC-----

# Get values for world model for the correlation stars from the barplot script
biom_world_0_200m <- read_csv("results/",config$region,"/biomass_tot_1x1/biom_world_0_200m.csv")
biom_world_200_500m <- read_csv("results/",config$region,"/biomass_tot_1x1/biom_world_200_500m.csv")

biom_world_all <- rbind(biom_world_0_200m,biom_world_200_500m) %>% 
  group_by(start_layer, end_layer) %>% 
  summarise(world_PgC_mean = sum(world_PgC_mean))
biom_world_stars <- rbind(biom_world_0_200m,biom_world_200_500m)

# Get the predicted biomass for all groups, regions, layers and store in a giant table
every_pred <- matrix(ncol =11, 
                     dimnames=list(c(), c(
                       #taxa, zones, layer
                       "taxa","start_layer","end_layer", "Region","lon","lat",
                       #biomass count (seen by UVP5) and predicted biomass
                       "biomass_tonC_mean","biomass_tonC_sd", "biomass_tonC_se", "pred_mean", "pred_sd")))
every_pred = data.frame(every_pred)

for (r in 1:nrow(config_region)) {
  for (lay in 1:nrow(layers)) {
    start_layer = layers$start_layer[lay]
    end_layer = layers$end_layer[lay] 
    # Get the groups names
    results_dir <- paste0("results/",config$element,"/", config$element, "_", config_region$zone[r], "/", start_layer,"_",end_layer,  "_", config$sd_or_no_sd)
    f <- list.files(results_dir, full.names=T) 
    models <- map(f, readRDS)
    names(models) <- basename(f)
    names_groups <- str_before_nth(basename(f), pattern = paste0("_",config_region$zone[r]), n = 1)
    
    for (t in 1:length(names_groups)) {
      message(paste0(config_region$zone[r], ": ", names_groups[t], " in ", start_layer, "-", end_layer, "m"))
      # get the predictions for each group
      biom_t <- read_csv(paste0("results/",config$region,"/biomass_tot_1x1/",config_region$zone[r], "/",start_layer,"_",end_layer , "_", config$sd_or_no_sd,"/", names_groups[t],"_",config_region$zone[r],"_",start_layer,"_",end_layer,"m_biomass_tot_1x1.csv")) %>% 
        mutate(taxa = names_groups[t], start_layer = start_layer, end_layer = end_layer, Region= config_region$zone[r]) %>% 
        select(taxa, start_layer, end_layer, Region, lon, lat, biomass_tonC_mean,biomass_tonC_sd, biomass_tonC_se, pred_mean, pred_sd)
      every_pred <- rbind(every_pred,biom_t)
    }
  }
}
every_pred <- every_pred %>% filter(!is.na(start_layer))
save(every_pred, file = "results/",config$region,"/every_pred.Rdata")
load("results/",config$region,"/every_pred.Rdata")

### 2.1. Numbers -------
### 2.1.1. All ------

# Now add the info about the correlation for our 3 regions 
cor_d <- cor_data %>%  
  filter(Region == "world" | Region == "region_outside" | Region == "region_inside") 
cor_d['Region'][cor_d['Region'] == 'region_inside'] <- 'inside_40S_40N'
cor_d['Region'][cor_d['Region'] == 'region_outside'] <- 'outside_40S_40N'

every_pred <- every_pred %>%  left_join(cor_d) %>%
  # keep only world data 
  filter(Region == "world") 

## 0-200m
pred_0_200_all_star <- every_pred %>%  filter(start_layer == 0 & end_layer == 200) %>% filter(stars =="*") %>% 
  group_by(lon, lat) %>%  
  summarise(biomass_tonC_m = sum(biomass_tonC_mean), biomass_tonC_sd = sd(biomass_tonC_mean)) %>% 
  ungroup()
biomass_tot_0_200_star <- pred_0_200_all_star %>% 
  summarise(biomass_tonC_mean = sum(biomass_tonC_m),
            biomass_tonC_sd = sd(biomass_tonC_m),
            biomass_PgC_mean = biomass_tonC_mean/10^9,
            biomass_PgC_sd = biomass_tonC_sd/10^9) %>% 
  mutate(start_layer = 0, end_layer = 200, star = "*")


## 200-500m
pred_200_500_all_star <- every_pred %>%  filter(start_layer == 200 & end_layer == 500) %>% filter(stars =="*") %>% 
  group_by(lon, lat) %>%  
  summarise(biomass_tonC_m = sum(biomass_tonC_mean), biomass_tonC_sd = sd(biomass_tonC_mean)) %>% 
  ungroup()
biomass_tot_200_500_star <- pred_200_500_all_star %>% 
  summarise(biomass_tonC_mean = sum(biomass_tonC_m),
            biomass_tonC_sd = sd(biomass_tonC_m),
            biomass_PgC_mean = biomass_tonC_mean/10^9,
            biomass_PgC_sd = biomass_tonC_sd/10^9) %>% 
  mutate(start_layer = 200, end_layer = 500, star = "*")


## 0-500m
pred_0_500_all_star <- every_pred %>% filter(stars =="*") %>%  
  filter(start_layer == 0 & end_layer == 200 | start_layer == 200 & end_layer == 500) %>% 
  group_by(lon, lat) %>%  
  summarise(biomass_tonC_m = sum(biomass_tonC_mean), biomass_tonC_sd = sd(biomass_tonC_mean)) %>% 
  ungroup()
biomass_tot_0_500_star <- pred_0_500_all_star %>% 
  summarise(biomass_tonC_mean = sum(biomass_tonC_m),
            biomass_tonC_sd = sd(biomass_tonC_m),
            biomass_PgC_mean = biomass_tonC_mean/10^9,
            biomass_PgC_sd = biomass_tonC_sd/10^9) %>% 
  mutate(start_layer = 0, end_layer = 500, star = "*")

biomass_tot <- rbind(biomass_tot_0_200_star,biomass_tot_200_500_star,biomass_tot_0_500_star) %>% 
  select(start_layer, end_layer, star, biomass_PgC_mean, biomass_PgC_sd, biomass_tonC_mean, biomass_tonC_sd)

write_excel_csv(biomass_tot, paste0("results/",config$region,"/biomass_tot_1x1/biomass_tot.csv"))


### 2.1.2 Per group -------
# use what was done early on in 4.global-3.barplots.R and regrouped earlier in this script
biomass_world_0_500 <- rbind(biom_world_0_200m, biom_world_200_500m)

# Compute % of biomass 
biomass_tot_PgC_0_200 <- biomass_tot %>%  filter(start_layer == 0 & end_layer == 200 & star == "*")
biom_world_0_200m_final <- biom_world_0_200m %>% filter(stars == "*") %>% 
  mutate(percent = world_PgC_mean / biomass_tot_PgC_0_200$biomass_PgC_mean * 100) %>% 
  select(-R2, -stars, - world_PgC_mean_cum, -rn, -Region, -p_value, -cor) %>% 
  mutate()
biom_world_0_200m_final <- biom_world_0_200m_final[order(biom_world_0_200m_final$percent, decreasing = TRUE), ]


biomass_tot_PgC_200_500 <- biomass_tot %>%  filter(start_layer == 200 & end_layer == 500 & star == "*")
biom_world_200_500m_final <- biom_world_200_500m %>% filter(stars == "*") %>% 
  mutate(percent = world_PgC_mean / biomass_tot_PgC_200_500$biomass_PgC_mean * 100) %>% 
  select(-R2, -stars, - world_PgC_mean_cum, -rn, -Region, -p_value, -cor)
biom_world_200_500m_final <- biom_world_200_500m_final[order(biom_world_200_500m_final$percent, decreasing = TRUE), ]

biom_world_0_200_500_final <- rbind(biom_world_0_200m_final, biom_world_200_500m_final)
biom_world_0_200_500_final <- biom_world_0_200_500_final[order(biom_world_0_200_500_final$taxa), ]

biomass_tot_PgC_0_500 <- biomass_tot %>%  filter(start_layer == 0 & end_layer == 500 & star == "*")
biomass_world_0_500_final <- biomass_world_0_500 %>% filter(stars == "*")  %>% 
  # compute biomass per group 
  group_by(taxa) %>% 
  summarise(world_PgC_m = sum(world_PgC_mean)) %>% 
  mutate(percent = (world_PgC_m / biomass_tot_PgC_0_500$biomass_PgC_mean * 100))

biomass_world_0_500_final <- biomass_world_0_500_final[order(biomass_world_0_500_final$percent, decreasing = TRUE), ]


## Final table 
biom_world_0_200_500_final
biomass_world_0_500_final
write_excel_csv(biom_world_0_200_500_final, paste0("results/",config$region,"/biom_world_0__200_500_final.csv"))
write_excel_csv(biomass_world_0_500_final, paste0("results/",config$region,"/biomass_world_0_500_final.csv"))


## 3. Plots of global biomass in mgC/m3----
biomass_mgC_m3 <- every_pred %>%  left_join(cor_d) 

biomass_mgC_m3_0_200 <- biomass_mgC_m3 %>% filter(start_layer == 0 & end_layer == 200) %>% filter(stars =="*") %>% 
  group_by(lon, lat) %>%  
  summarise(biomass_mgC_m3_mean = sum(pred_mean), biomass_mgC_m3_sd = sd(pred_mean)) %>% 
  ungroup()
biomass_mgC_m3_200_500 <- biomass_mgC_m3 %>% filter(start_layer ==200 & end_layer == 500) %>% filter(stars =="*") %>% 
  group_by(lon, lat) %>%  
  summarise(biomass_mgC_m3_mean = sum(pred_mean), biomass_mgC_m3_sd = sd(pred_mean)) %>% 
  ungroup()
biomass_mgC_m3_0_500 <- biomass_mgC_m3 %>% 
  filter(start_layer == 0 & end_layer == 200 & stars =="*" |
         start_layer == 200 & end_layer == 500 & stars =="*") %>%  
  group_by(lon, lat) %>%  
  summarise(biomass_mgC_m3_mean = sum(pred_mean), biomass_mgC_m3_sd = sd(pred_mean)) %>% 
  ungroup()

cor_test_all # look at the very top for the correlation tests

text_size = 6
global_biomass_0_200 <- ggplot(biomass_mgC_m3_0_200) + 
  geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
  geom_raster(aes(lon, lat, fill=biomass_mgC_m3_mean)) +
  scale_fill_viridis_c(trans="sqrt", lim = c(0,200)) +
  coord_quickmap() + scale_xy_map() +
  theme_dark() +
  theme(axis.text = element_text(size=text_size), axis.title = element_text(size=text_size), plot.title = element_text(size = 8), 
        legend.text = element_text(size = text_size),legend.title = element_text(size = text_size)) +
  labs(fill="Predicted\nbiomass\n(mgC/m3)\n", 
       title=paste0("0-200m: R2 = ", round(cor_test_all$R2[4], digits = 2), ", cor = ", round(cor_test_all$cor[4], digits = 2), ", p_value = ", cor_test_all$p_value[4]))

global_biomass_200_500 <- ggplot(biomass_mgC_m3_200_500) + 
  geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
  geom_raster(aes(lon, lat, fill=biomass_mgC_m3_mean)) +
  scale_fill_viridis_c(trans="sqrt") +
  coord_quickmap() + scale_xy_map() +
  theme_dark() +
  theme(axis.text = element_text(size=text_size), axis.title = element_text(size=text_size), plot.title = element_text(size = 8), 
        legend.text = element_text(size = text_size),legend.title = element_text(size = text_size)) +
  labs(fill="Predicted\nbiomass\n(mgC/m3)\n", 
       title=paste0("200-500: R2 = ", round(cor_test_all$R2[5], digits = 2), ", cor = ", round(cor_test_all$cor[5], digits = 2), ", p_value = ", cor_test_all$p_value[5]))

global_biomass_0_500 <- 
  ggplot(biomass_mgC_m3_0_500) + 
  geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
  geom_raster(aes(lon, lat, fill=biomass_mgC_m3_mean)) +
  scale_fill_viridis_c(trans="sqrt", lim = c(0,100)) +
  coord_quickmap() + scale_xy_map() +
  theme_dark() +
  labs(fill="Predicted\nbiomass\n(mgC/m3)\n", 
       title=paste0("0-500m: R2 = ", round(cor_test_all$R2[6], digits = 2), ", cor = ", round(cor_test_all$cor[6], digits = 2), ", p_value = ", cor_test_all$p_value[6]))

print((global_biomass_0_200) / (global_biomass_200_500) / (global_biomass_0_500))
ggsave(file=paste0("results/",config$region,"/figures_dvpmt/global_maps.pdf"),
       plot = last_plot(),
       width=18, height=27, units = "cm")


## 4. Latitudinal patterns of global biomass 
## 4.1. in PgC ------

# Use these 
pred_0_200_all_star
pred_200_500_all_star
pred_0_500_all_star

# 0-500m
ggplot()+
  #0-500m
  geom_point(data = pred_0_500_all_star, aes(y=log10(biomass_tonC_m+1), x=lat, col = "#9ecae1"), alpha=1, size=0.6) +   
  geom_smooth(data = pred_0_500_all_star, aes(y=log10(biomass_tonC_m+1), x=lat, col = "#08519c"), se=T, alpha=1, size=0.6) + 
  coord_flip() +
  scale_y_continuous("Predicted integrated biomass 0-500m (PgC)", breaks=c(0,1,2,3,4,5,6,7),
                     labels=c("0","10","100","1 000","10 000","100 000","1 000 000", "10 000 000"))+
  scale_x_continuous("Latitude", limits=c(-80, 90), breaks=seq(-80, 90, 10),
                     labels=c("","","60°S","","","30°S","","","0°","","","30°N","","","60°N","","","90°N"))+
  theme(legend.position="bottom")  +
  scale_color_identity(guide = "legend",
                       name = "Data",
                       breaks = c("#08519c"),
                       labels = c("0-500m")) 
ggsave(file=paste0("results/",config$region,"/figures_dvpmt/distribution_lat_biomass.pdf"), 
       plot = last_plot(),
       width=18, height=18, units = "cm") 

ggplot()+
  #0-500m
  geom_point(data = pred_200_500_all_star, aes(y=log10(biomass_tonC_m+1), x=lat, col = "#fdd0a2"), alpha=0.3, size=0.6) +   
  geom_point(data = pred_0_500_all_star, aes(y=log10(biomass_tonC_m+1), x=lat, col = "#9ecae1"), alpha=0.3, size=0.6) +   
  
  geom_smooth(data = pred_200_500_all_star, aes(y=log10(biomass_tonC_m+1), x=lat, col = "#d94801"), se=T, alpha=1, size=0.6) + 
  geom_smooth(data = pred_0_500_all_star, aes(y=log10(biomass_tonC_m+1), x=lat, col = "#4292c6"), se=T, alpha=1, size=0.6) + 
  coord_flip() +
  scale_y_continuous("Integrated biomass 0-500m (PgC)", breaks=c(0,1,2,3,4,5,6,7),
                     labels=c("0","10","100","1 000","10 000","100 000","1 000 000", "10 000 000"))+
  scale_x_continuous("Latitude", limits=c(-80, 90), breaks=seq(-80, 90, 10),
                     labels=c("","","60°S","","","30°S","","","0°","","","30°N","","","60°N","","","90°N"))+
  theme(legend.position="bottom")  +
   scale_color_identity(guide = "legend",
                        name = "Data",
                        breaks = c("#4292c6", "#f16913"),
                        labels = c("0-500m", "200m-500m")) +
  ggtitle("Predicted biomass")


## 4.2. in mgC/m3 ------
# Use these 
biomass_mgC_m3_0_200 <- biomass_mgC_m3 %>% filter(start_layer == 0 & end_layer == 200) %>% filter(stars =="*") %>% 
  group_by(lon, lat) %>%  
  summarise(pred_mean_200 = sum(200 * pred_mean)) %>% 
  ungroup()
biomass_mgC_m3_200_500 <- biomass_mgC_m3 %>% filter(start_layer ==200 & end_layer == 500) %>% filter(stars =="*") %>% 
  group_by(lon, lat) %>%  
  summarise(pred_mean_500 = sum(300 * pred_mean)) %>% 
  ungroup()

biomass_mgC_m3_0_500 <- biomass_mgC_m3_0_200 %>%  
  left_join(biomass_mgC_m3_200_500, by = c("lat", "lon")) %>% 
  group_by(lon, lat) %>% 
  summarise(pred_mean = pred_mean_200 + pred_mean_500) %>% 
  ungroup()

ggplot()+
  # geom_point(data = biomass_mgC_m3_0_500, aes(y=log10(pred_mean+1), x=lat, col = "#9ecae1"), alpha=0.8, size=0.6)+ 
  geom_smooth(data = biomass_mgC_m3_0_500, aes(y=log10(pred_mean+1), x=lat, col = "#332288"), se=T, alpha=1, size=0.6) +
  coord_flip() +
  scale_y_continuous("Biomass 0-500m (mgC.m-2)", breaks=c(0,1,2,3,4,5,6,7),
                     labels=c("0","10","100","1 000","10 000","100 000","1 000 000", "10 000 000"))+
  scale_x_continuous("Latitude", limits=c(-80, 90), breaks=seq(-80, 90, 10),
                     labels=c("","","60°S","","","30°S","","","0°","","","30°N","","","60°N","","","90°N"))+
  scale_color_identity(guide = "legend",
                       name = "Data", 
                       breaks = c("#332288"),
                       labels = c("Our models")) +
  ggtitle("Biomass 0-500m") +
  theme(legend.position="none") 

ggsave(file=paste0("results/",config$region,"/figures_dvpmt/distribution_lat_biomass.pdf"), 
       plot = last_plot(),
       width=18, height=18, units = "cm") 
