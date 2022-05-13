#
# Explore the result of models
#
# Creation of table d_test_all containing all the predicted data on test set "d_test"
#
# Plots of biomass as seen by the UVP5 VS predicted by our models ~ 45min for all layers, regions, groups
#
# (c) 2022 L Drago, T Panaïotis, J-O Irisson, R Kiko, GNU General Public License v3

source("Scripts/0.2.config.R")
library("chroma")
library("grid")
library("strex")
library("sf")
library("joml")
library("RColorBrewer")
library("gridExtra")
library("patchwork")
library("ggpubr")

# create a table containing all the predicted data on test set so that we can do it for all groups later
d_test_all <- matrix(ncol = 6, 
                     dimnames=list(c(), c(
                       #taxa, zones, layer
                       "taxa","start_layer","end_layer", "Region",
                       #biomass count (seen by UVP5) and predicted biomass
                       "biomass_count","pred_mean")))
d_test_all = data.frame(d_test_all)

for (r in 1:nrow(config_region)) {
  for (lay in 1:nrow(layers)) {
    
    start_layer = layers$start_layer[lay]
    end_layer = layers$end_layer[lay]
    
    # get the list of files and taxa for this region and layer 
    results_dir <- paste0("results/",config$element,"/", config$element, "_", config_region$zone[r], "/", start_layer,"_",end_layer,  "_", config$sd_or_no_sd)
    f <- list.files(results_dir, full.names=T) 
    models <- map(f, readRDS)
    names(models) <- basename(f)
    names_groups <- str_before_nth(basename(f), pattern = paste0("_",config_region$zone[r]), n = 1)
    
    # extract from the models the fits and the predictions for all groups
    fits <- map_dfr(models, function(l) {data.frame(taxa = l$t, l$cv_metrics)})
    prediction_all <- map_dfr(models, function(l) {data.frame(taxa = l$t, l$prediction)})
    
    # Now group by group
    for (t in 1:length(names_groups)) {
      message(paste0(config_region$zone[r], ": ", names_groups[t], " in ", start_layer, "-", end_layer, "m"))
      # select the predictions and models associated
      prediction_t <- prediction_all %>%  filter (taxa == names_groups[t])
      models_t <- fits %>%  filter (taxa == names_groups[t])
      f_t <- list.files(results_dir, pattern= paste0(names_groups[t]),full.names=T)
      if (names_groups[t] == "Copepoda") {f_t <- f_t[2]}
      if (names_groups[t] == "Copepoda_refit") {f_t <- f_t[1]}
      
      models_t <- map(f_t, readRDS)
      
      # extract the test set, training set, best 100 models and lowest error combination
      d_test <- map_dfr(models_t, function(l) {data.frame(taxa = l$t, l$d_test)})
      d_train <- map_dfr(models_t, function(l) {data.frame(taxa = l$t, l$d_train)})
      m_best <- map_dfr(models_t, function(l) {l$model_full})
      lowest_error <- map_dfr(models_t, function(l) {l$lowest_error})
      # for the world model, get the part of the training set from inside 40S_40N and outside
      if (config_region$zone[r] == "world") {
        d_test_inside <- map_dfr(models_t, function(l) {data.frame(taxa = l$t, l$d_test_inside)})
        d_test_outside <- map_dfr(models_t, function(l) {data.frame(taxa = l$t, l$d_test_outside)})
      }
      
      # predict on training and testing data
      predict_test_set <- m_best %>%
        xgb_predict(ntrees=lowest_error$iter, newdata=d_test, fns=list(mean=mean, sd=sd, se=se))
      d_test_set <- d_test %>%  left_join(predict_test_set)
      # save it in d_test_all
      d_test_set_save <- d_test_set %>%  select(taxa, biomass_count, pred_mean) %>%  
        mutate(start_layer = start_layer, end_layer= end_layer, Region = config_region$zone[r])
      d_test_all <- rbind(d_test_all, d_test_set_save)
      
      if (config_region$zone[r] == "world") {
        # for the world model, do the same but for the d_test_inside and d_test_outside 
        predict_test_set <- m_best %>%
          xgb_predict(ntrees=lowest_error$iter, newdata=d_test_inside, fns=list(mean=mean, sd=sd, se=se))
        d_test_set <- d_test_inside %>%  left_join(predict_test_set)
        # save it in d_test_all
        d_test_set_save <- d_test_set %>%  select(taxa, biomass_count, pred_mean) %>%  
          mutate(start_layer = start_layer, end_layer= end_layer, Region = "world_inside")
        d_test_all <- rbind(d_test_all, d_test_set_save)
        
        predict_test_set <- m_best %>%
          xgb_predict(ntrees=lowest_error$iter, newdata=d_test_outside, fns=list(mean=mean, sd=sd, se=se))
        d_test_set <- d_test_outside %>%  left_join(predict_test_set)
        # save it in d_test_all
        d_test_set_save <- d_test_set %>%  select(taxa, biomass_count, pred_mean) %>%  
          mutate(start_layer = start_layer, end_layer= end_layer, Region = "world_outside")
        d_test_all <- rbind(d_test_all, d_test_set_save)
      }
      
      # plot it
      fig_test <- ggplot(d_test_set) + coord_fixed() +
        geom_abline(intercept=0, slope=1, colour="#3B87B3") +  # 1:1 line
        geom_point(aes(biomass_count, pred_mean), size=0.5, alpha=0.5) + # data points
        scale_x_continuous(trans="log1p") + scale_y_continuous(trans="log1p") + # nicer aspect
        labs(x="True biomass", y="Predicted biomass") +
        ggtitle("test") +
        theme(plot.title = element_text(size = 5),
              axis.text=element_text(size=4),
              axis.title=element_text(size=5))
      fig_train <- ggplot(d_train_set) + coord_fixed() + # 1:1 line
        geom_abline(intercept=0, slope=1, colour="#3B87B3") + # data points
        geom_point(aes(biomass_count, pred_mean), size=0.5, alpha=0.5) + # nicer aspect
        scale_x_continuous(trans="log1p") + scale_y_continuous(trans="log1p") +
        labs(x="True biomass", y="Predicted biomass") +
        ggtitle("train") +
        theme(plot.title = element_text(size = 5),
              axis.text=element_text(size=4),
              axis.title=element_text(size=5))
      
      plot_train_test <- (fig_train + fig_test) + plot_annotation(
        title = paste0(names_groups[t], ": ",config_region$zone[r]," ", start_layer, "-", end_layer, "m"),
        tag_levels = c('A')) & theme(plot.tag = element_text(size = 5))
    }
  }
}

d_test_all <- d_test_all %>% filter(!is.na(start_layer))

# save it
write_excel_csv(d_test_all, paste0("results/",config$region,"/d_test_all.csv"))