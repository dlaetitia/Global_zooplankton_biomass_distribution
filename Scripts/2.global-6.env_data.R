#
# Explore the environmental data
#
# (c) 2022 L Drago, T Panaïotis, J-O Irisson, R Kiko, GNU General Public License v3

source("Scripts/0.2.config.R")
library("chroma")
library("RColorBrewer")
library("patchwork")
library("ggpubr")
library("cmocean")


# ### 1. Datapoints -----
load(paste0("Data/", config$element, "/1.uvp_samples.Rdata"))

ggplot() +
  geom_polygon(aes(lon, lat), data=coast, fill="grey75")+
  coord_quickmap() + scale_xy_map() +
  geom_point(aes(lon, lat), size = 0.75, data=samples, alpha=0.5) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave(file=paste0("Figures/",config$region,"/campaigns_metadata_black.pdf"),
       plot = last_plot(), width=18, units = "cm")

ggplot() +
  geom_polygon(aes(lon, lat), data=coast, fill="grey75")+
  coord_quickmap() +scale_xy_map() +
  geom_point(aes(lon, lat, color = project), size = 0.75, data=samples, alpha=0.5) +
  theme(legend.position="bottom", legend.text = element_text(size=4), legend.title = element_text(size=8)) +
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave(file=paste0("Figures/",config$region,"/campaigns_metadata_colour.pdf"),
       plot = last_plot(), width=18, units = "cm")


### 2. Environmental data -----
start_layer = 0
end_layer = 200

load(paste0("Data/",config$region,"/2.env_world_",start_layer,"_",end_layer,"m.Rdata"))

m = 0 # size of marge of figure
s = 6 # size of text figure
t = 8 # size of title figure
size_leg = 0.5 # size of legend

temp <- ggplot(env_y) +
  geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
  geom_raster(aes(lon, lat, fill=mean.temperature)) +
  scale_fill_gradientn(colours = cmocean('thermal')(256), name = "°C") +
  coord_quickmap() + scale_xy_map() + ggtitle("Temperature")  +
  theme(plot.margin = margin(t = m, r = m, b = m, l =m, unit = "pt"),
        axis.text = element_text(size=s), axis.title = element_text(size=s),
        plot.title = element_text(size = t),  
        legend.text = element_text(size = s),legend.title = element_text(size = s), legend.key.size = unit(size_leg, 'cm'))

sal <- ggplot(env_y) +
  geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
  geom_raster(aes(lon, lat, fill=mean.salinity)) +
  scale_fill_gradientn(colours = cmocean('haline')(256), name = " ") +
  coord_quickmap() + scale_xy_map() + ggtitle("Salinity") +
  theme(plot.margin = margin(t = m, r = m, b = m, l =m, unit = "pt"),
        axis.text = element_text(size=s), axis.title = element_text(size=s),
        plot.title = element_text(size = t),  
        legend.text = element_text(size = s),legend.title = element_text(size = s), legend.key.size = unit(size_leg, 'cm'))

oxy <- ggplot(env_y) +
  geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
  geom_raster(aes(lon, lat, fill=mean.oxygen)) +
  scale_fill_gradientn(colours = cmocean('oxy')(256), name = "kPa") +
  coord_quickmap() + scale_xy_map() + ggtitle("Oxygen") +
  theme(plot.margin = margin(t = m, r = m, b = m, l =m, unit = "pt"),
        axis.text = element_text(size=s), axis.title = element_text(size=s),
        plot.title = element_text(size = t),  
        legend.text = element_text(size = s),legend.title = element_text(size = s), legend.key.size = unit(size_leg, 'cm'))

silicate <- ggplot(env_y) +
  geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
  geom_raster(aes(lon, lat, fill=mean.silicate)) +
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9,'YlOrBr'))(64), name = "µmol.kg-1") +
  coord_quickmap() + scale_xy_map() + ggtitle("Silicate") +
  theme(plot.margin = margin(t = m, r = m, b = m, l =m, unit = "pt"),
        axis.text = element_text(size=s), axis.title = element_text(size=s),
        plot.title = element_text(size = t),  
        legend.text = element_text(size = s),legend.title = element_text(size = s), legend.key.size = unit(size_leg, 'cm'))

nitrate <- ggplot(env_y) +
  geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
  geom_raster(aes(lon, lat, fill=mean.nitrate)) +
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9,'PuBu'))(64), name = "µmol.kg-1") +
  coord_quickmap() + scale_xy_map() + ggtitle("Nitrate") +
  theme(plot.margin = margin(t = m, r = m, b = m, l =m, unit = "pt"),
        axis.text = element_text(size=s), axis.title = element_text(size=s),
        plot.title = element_text(size = t),  
        legend.text = element_text(size = s),legend.title = element_text(size = s), legend.key.size = unit(size_leg, 'cm'))

phosphate <- ggplot(env_y) +
  geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
  geom_raster(aes(lon, lat, fill=mean.phosphate)) +
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9,'BuPu'))(64), name = "µmol.kg-1") +
  coord_quickmap() + scale_xy_map() + ggtitle("Phosphate") +
  theme(plot.margin = margin(t = m, r = m, b = m, l =m, unit = "pt"),
        axis.text = element_text(size=s), axis.title = element_text(size=s),
        plot.title = element_text(size = t),  
        legend.text = element_text(size = s),legend.title = element_text(size = s), legend.key.size = unit(size_leg, 'cm'))

# Arrrange everything
layout_env <- "
AABB
AABB
CCDD
CCDD
EEFF
EEFF
"
temp + sal + oxy +
  silicate + nitrate + phosphate +
  plot_layout(design = layout_env) +
  plot_annotation(title = paste0(start_layer, "-", end_layer, "m"),
                  theme = theme(plot.title = element_text(size = 10)))

ggsave(file=paste0("Figures/",config$region,"/env_maps_", start_layer, "_", end_layer, ".pdf"),
       plot = last_plot(),
       width=18, height=14.5, units = "cm", dpi = 300)


## Other variables in 0-200m
load(paste0("Data/",config$region,"/2.env_world_0_200m.Rdata"))
chla <- ggplot(env_y) +
  geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
  geom_raster(aes(lon, lat, fill=mean.chla)) +
  scale_fill_gradientn(colours = cmocean('algae')(256), name = "mg.m-3", trans = 'log10') +
  coord_quickmap() + scale_xy_map() + ggtitle("Chlorophyll a") +
  theme(plot.margin = margin(t = m, r = m, b = m, l =m, unit = "pt"),
        axis.text = element_text(size=s), axis.title = element_text(size=s),
        plot.title = element_text(size = t),  
        legend.text = element_text(size = s),legend.title = element_text(size = s), legend.key.size = unit(size_leg, 'cm'))

bathymetry <- ggplot(env_y) +
  geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
  geom_raster(aes(lon, lat, fill=bathymetry)) +
  scale_fill_gradientn(colours = cmocean('deep')(256), name = "m") +
  coord_quickmap() + scale_xy_map() + ggtitle("Bathymetry") +
  theme(plot.margin = margin(t = m, r = m, b = m, l =m, unit = "pt"),
        axis.text = element_text(size=s), axis.title = element_text(size=s),
        plot.title = element_text(size = t),  
        legend.text = element_text(size = s),legend.title = element_text(size = s), legend.key.size = unit(size_leg, 'cm'))

dist <- ggplot(env_y) +
  geom_raster(aes(lon, lat, fill=dist)) +
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9,'Blues'))(64), name = "km") +
  geom_polygon(aes(lon, lat), data=coast, fill="grey40") +
  coord_quickmap() + scale_xy_map() + ggtitle("Distance to coast") +
  theme(plot.margin = margin(t = m, r = m, b = m, l =m, unit = "pt"),
        axis.text = element_text(size=s), axis.title = element_text(size=s),
        plot.title = element_text(size = t),  
        legend.text = element_text(size = s),legend.title = element_text(size = s), legend.key.size = unit(size_leg, 'cm'))

# Arrrange everything
layout_env2 <- "
AABB
AABB
CC##
CC##
"
dist + chla + bathymetry +
  plot_layout(design = layout_env2) +
  plot_annotation(title = paste0(" "),
                  theme = theme(plot.title = element_text(size = 10)))

ggsave(file=paste0("Figures/",config$region,"/env_maps_others.pdf"),
       plot = last_plot(),
       width=18, height=10, units = "cm", dpi = 300)


