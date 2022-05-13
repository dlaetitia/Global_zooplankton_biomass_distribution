#
# Fit model for all taxa x layers x biovolume measurement
#
# (c) 2022 L Drago, T Panaïotis, J-O Irisson, R Kiko, GNU General Public License v3

# NB: this can be run as an RStudio "job"

source("Scripts/0.2.config.R")
library("devtools")
library("multiscales")
# devtools::install_github("clauswilke/multiscales")
library("geometry")
library("gridExtra")
library("jihoml")
# remotes::install_github("jiho/jihoml")


# Loop through the layers
for (layer_num in 1:nrow(layers)) { 
  start_layer = layers$start_layer[layer_num]
  end_layer = layers$end_layer[layer_num]
  
  ## Fit model for all taxa,layers,biomass measurements ----
  
  # get biomass+matched env
  load(paste0("Data/", config$element, "/2.biomass_env_",start_layer,"_",end_layer,"m.Rdata"))
  
  # Now, let's only select the elements which are in the zone we want   
  if (config_region$zone[r] == "outside_40S_40N") {
    d <- d %>% 
      # Keep only the data above 40°N and below 40°S
      filter(lat > 40 | lat < -40) 
  }
  if (config_region$zone[r] == "inside_40S_40N") {
    d <- d %>% 
      # Keep only the data above 40°N and below 40°S
      filter(between(lat,-40,40))
  }
  
  ### 3. Handle the environmental data  -----
  ### 3.1. Remove AOU -----
  data_d <- d %>% 
    # remove AOU
    select(-contains("AOU")) %>% 
    rename(biomass = biomass_sph) %>% select(-biomass_ell) 
  if (config$sd_or_no_sd == "no_sd") {data_d <- data_d %>% select(-contains("sd"))}
  
  load(paste0("Data/",config$element,"/2.env_world_",start_layer,"_",end_layer,"m.Rdata"))
  env <- env_y %>%  
    # remove sd values in env dataset
    select(-c(contains("sd"))) %>% 
    # remove AOU 
    select(-c(contains("AOU")))
  remove(env_m, env_y)
  
  # select the data corresponding to the region of interest
  if (config_region$zone[r] != "world") {
    if (config_region$zone[r] == "outside_40S_40N") {
      env <- env %>% 
        # Keep only the data above 40°N and below 40°S
        filter(lat > 40 | lat < -40) 
    }
    if (config_region$zone[r] == "inside_40S_40N") {
      env <- env %>% 
        # Keep only the data between 40°N and 40°S
        filter(between(lat,-40,40))
    }
  }
  # reduce to actually defined pixels
  env <- drop_na(env, mean.temperature)
  
  # Prepare storage
  dir.create(paste0("results/",config$element,"/", config$element, "_", config_region$zone[r]), showWarnings=FALSE)
  dir.create(paste0("results/",config$element,"/", config$element, "_", config_region$zone[r], "/", start_layer,"_",end_layer, "_", config$sd_or_no_sd), showWarnings=FALSE)
  dir.create(paste0("results/",config$element,"/predicted_vs_true"), showWarnings=FALSE)
  dir.create(paste0("results/",config$element,"/predicted_vs_true/", config_region$zone[r]), showWarnings=FALSE)
  dir.create(paste0("results/",config$element,"/predicted_vs_true/", config_region$zone[r],"/",start_layer,"_",end_layer, "_", config$sd_or_no_sd), showWarnings=FALSE)
  
  results_dir <- paste0("results/",config$element,"/", config$element, "_", config_region$zone[r], "/", start_layer,"_",end_layer, "/")
  
  # Load stats of groups to know the names of the groups kept in this layer
  stats <- read.csv(paste0("Data/",config$element,"/stats_biom_percent_0_in_",start_layer,"_",end_layer,"m_",config$layers,".csv")) %>% 
    spread(key = "taxon", value = "mean_biom") %>% 
    select(-X, -percent_0)
  taxa <- names(stats)
  
  # prepare a table for the stats
  dir.create(paste0("results/",config$element,"/cor_table"), showWarnings=FALSE)
  cor_table <- matrix(ncol = 8, 
                      nrow = length(taxa),
                      dimnames=list(c(), c("taxa","start_layer","end_layer","Region", 
                                           "p_value_pearson","cor_pearson", 
                                           "p_value_spearman", "cor_spearman")))
  cor_table=data.frame(cor_table)
  counter_cor = 0
  cor_table[,"start_layer"] = start_layer
  cor_table[,"end_layer"] = end_layer
  cor_table[,"Region"] = config_region$zone[r]
  
  start_layer 
  end_layer

  for (t in taxa) {
    message(t)
    counter_cor = counter_cor + 1
    # Clean data
    d <- data_d %>%
      # pick a taxon
      filter(taxon==t) 
    
    # determine the scaling level
    small_quant <- quantile(d$biomass[d$biomass>0], 0.01)
    # make sure this becomes 1 after scaling (and round to the nearest 100, for cleanliness)
    scale_f <- ceiling(1/small_quant/100)*1000
    # Transform biomass 
    d <- d %>% 
      mutate(
        # for Poisson fitting
        biomass_count = round(biomass*scale_f),
        # for stratification by biomass
        biomass_quant=cut(biomass,
                          breaks=c(0, quantile(d$biomass[d$biomass > 0], seq(0, 1, length=9))),
                          include.lowest=TRUE, labels=FALSE)
      )
    
    set.seed(1)
    # split the data according to the distribution of biomass with 80% in train set and 20% in test set
    train_test <- resample_split(d, p=0.8, biomass_quant)
    # identify the train and test sets
    d_train <- as.data.frame(train_test$train) %>% select(-(lon:taxon))
    d_test  <- as.data.frame(train_test$val)  %>% select(-(lon:taxon))
    # for the world model, extract regional data to be able to compute correlations on that further down
    if (config_region$zone[r] == "world") {
      d_test_inside  <- as.data.frame(train_test$val) %>%  filter(between(lat,-40,40))
      d_test_outside  <- as.data.frame(train_test$val) %>%  filter(lat > 40 | lat < -40) 
    }
    
    # 2. Split the training set for stratified cross-val
    # in addition, to be more robust, perform few repetitions of the cross-val splits
    rs <- resample_cv(d_train, k=5, n=20, biomass_quant)
    rs
    # = now we have 20 times 5 folds = 100 resamples
    
    # prepare a grid of parameters for tuning
    eta_values=c(0.05, 0.075, 0.1)
    max_depth_values=c(2, 4, 6)
    min_child_weight_values=c(1,3,5)
    rs_grid <- param_grid(
      rs,                # the CV resamples we will run for each parameter combination
      eta=eta_values,   # learning rate
      max_depth=max_depth_values,  # maximum depth of trees
      min_child_weight=min_child_weight_values
    )
    # we now have the parameter values in our grid
    head(rs_grid)
    
    # we have 100 resamples per combination and 27 combinations => 2700 resamples to fit
    nrow(rs_grid)
    
    # 3. Fit the model on the data
    # for each of the resamples, fit on each train part, validate on the val part
    var <- names(select(d_train, -contains(c("biomass"))))
    
    system.time(
      m <- xgb_fit(rs_grid, resp="biomass_count", expl=var,
                   nrounds=600,
                   subsample=0.5, objective="count:poisson", cores=20)
    )
    
    # 4. Explore the model for different hyperparameter combination
    # now the statistics will be computed for each combination of parameters separately
    # (because m, above, is still grouped by values of eta and max_depth)
    head(m)
    # so we will plot those combinations separately
    fit_sum <- xgb_summarise_fit(m, fns=c(mean=mean, se=se))
    head(fit_sum)
    ggplot(fit_sum) + facet_wrap(~max_depth+eta, labeller=label_context) +
      # NB: facets are used here but colour could be too, just turn values into factors first
      geom_ribbon(aes(x=iter,
                      ymin=val_poisson_nloglik_mean-val_poisson_nloglik_se,
                      ymax=val_poisson_nloglik_mean+val_poisson_nloglik_se), alpha=0.3) +
      geom_path(aes(x=iter, y=val_poisson_nloglik_mean)) + 
      geom_hline(aes(yintercept=min(val_poisson_nloglik_mean)), colour="red")
    
    # plot this with colours
    xgb_summarise_fit(m) %>% 
      mutate(
        eta=factor(eta),
        max_depth = factor(max_depth),
        min_child_weight=factor(min_child_weight)
      ) %>% 
      ggplot() + facet_wrap(~max_depth+min_child_weight, labeller=label_context) +
      # NB: facets are used here but colour could be too, just turn values into factors first
      geom_ribbon(aes(x=iter,
                      ymin=val_poisson_nloglik_mean-val_poisson_nloglik_se,
                      ymax=val_poisson_nloglik_mean+val_poisson_nloglik_se,
                      fill=eta
      ), alpha=0.3) +
      geom_path(aes(x=iter, y=val_poisson_nloglik_mean, colour=eta)) + 
      geom_hline(aes(yintercept=min(val_poisson_nloglik_mean)), colour="red", linetype = "dotted")
    
    
    # 5. extract the best hyperparameter combination
    eval_table <- matrix(ncol = length(fit_sum), 
                         nrow = length(eta_values) + length(max_depth_values) + length (min_child_weight_values),
                         dimnames=list(c(), c(names(eval))))
    eval_table=data.frame(eval_table)
    counter = 0
    for (i in 1:length(eta_values)) {
      for (j in 1:length(max_depth_values)) {
        for (k in 1:length(min_child_weight_values)) {
          sub_table <- fit_sum %>%  filter(eta == eta_values[i] 
                                           & max_depth == max_depth_values[j] 
                                           & min_child_weight==min_child_weight_values[k])
          counter = counter + 1
          # create an empty vector of the length sub_table
          eval_error <- rep(NA, nrow(sub_table)) 
          # find the lowest error and stop if the error hasn't lowered by 1% for every 10 trees step
          for (l in 11:nrow(sub_table)) {
            # if the value of error is between 
            if (between(sub_table$val_poisson_nloglik_mean[l], 
                        # the error 10 steps before - 1% of this error
                        sub_table$val_poisson_nloglik_mean[l-10] - 0.01*sub_table$val_poisson_nloglik_mean[l-10],
                        # and the error 10 steps before
                        sub_table$val_poisson_nloglik_mean[l-10])) {
              # then give the value 1 to eval_error vector at position l
              eval_error[l] = 1
            } else {
              eval_error[l] = NA
            }
          }
          # find the first value equal to 1 if we have one
          if ( min(which(eval_error == 1)) != Inf) {
            min_val <- min(which(eval_error == 1))
            # now that we have the minimum value of error, find the value Z of 95% upper confidence boundary
            # reached for the same number of trees
            sub_table[min_val,] 
            Z95_value = sub_table[min_val,"val_poisson_nloglik_mean"] + sub_table[min_val,"val_poisson_nloglik_se"]
            
            # find the closest value of the mean error, for the same loss as Z for values before
            sub_table <- sub_table %>%
              mutate(Z95 = Z95_value,
                     Diff = val_poisson_nloglik_mean - Z95) 
            
            # Function to find the lowest positive value 
            minpositive = function(x) min(x[x > 0])
            
            if ( minpositive(sub_table$Diff[0:min_val,]) != Inf) {
              # If the model is working for this combination, 
              # get the row corresponding to the lowest error
              min_diff <- which(sub_table$Diff == minpositive(sub_table$Diff[0:min_val,])) 
              min_diff
              sub_table[min_diff,]
              # and get the hyperparameter combination
              sub_table <- sub_table %>% select(-c(Z95, Diff))
              colnames(eval_table) =  colnames(sub_table)
              eval_table[counter,] <- sub_table[min_diff,]
            } else {
              # If the model is not working for this combination:
              min_diff <- 99999
              # and get the hyperparameter combination
              sub_table <- sub_table %>% select(-c(Z95, Diff))
              colnames(eval_table) =  colnames(sub_table)
              eval_table[counter,"eta"] <- eta_values[i]
              eval_table[counter,"max_depth"] <- max_depth_values[i]
              eval_table[counter,"min_child_weight"] <- min_child_weight_values[i]
              eval_table[counter,"iter"] <- 99999
              eval_table[counter,"val_poisson_nloglik_mean"] <- 99999
              eval_table[counter,"val_poisson_nloglik_se"] <- 99999
            }
          } else {
            # if we don't have a value equal to 1, then write a document 
            # containing the name of the group for which it didn't work
            min_val = NA
          }
        }
      }
    }
    # Look at the evaluation table containing the lowest error for each combination of hyperparameters
    is_it_working <- eval_table %>% 
      drop_na(val_poisson_nloglik_mean)
    # If we have only NA or 99999 values for val_poisson_nloglik_mean
    if (any(is_it_working$val_poisson_nloglik_mean != 99999) == FALSE) {
      # tell that it's not working
      indicator = "not_working"
      # create a document to have the info somewhere
      didnt_work <- data.frame(
        "start_layer" = start_layer,
        "end_layer" = end_layer,
        "group" = t,
        "region" = config_region$zone[r])
      write.csv(didnt_work, file=paste0("results/",config$element,"/",config$element,"_",config_region$zone[r],"/", t, "_didnt_work_", config_region$zone[r], "_",start_layer,"_",end_layer,"m_", config$layers,".csv"))
    } else {indicator = "working"}
    
    if (indicator == "working") {
      # What is the lowest error value in all this table ?
      eval_table <- eval_table %>% filter(!is.na(eta))
      lowest_error <- eval_table[eval_table$val_poisson_nloglik_mean == min(eval_table$val_poisson_nloglik_mean, na.rm = TRUE),]
      lowest_error
      
      if (any(lowest_error$val_poisson_nloglik_mean != 99999) == FALSE) {
        # tell that it's not working
        indicator = "not_working"
        # create a document to have the info somewhere
        didnt_work <- data.frame(
          "start_layer" = start_layer,
          "end_layer" = end_layer,
          "group" = t,
          "region" = config_region$zone[r])
        write.csv(didnt_work, file=paste0("results/",config$element,"/",config$element,"_",config_region$zone[r],"/", t, "_didnt_work_", config_region$zone[r], "_",start_layer,"_",end_layer,"m_", config$layers,".csv"))
      } else {
        
        # Look at the distribution of R2 in all the bootstraps
        m_best <- filter(m, eta==lowest_error$eta, max_depth==lowest_error$max_depth,
                         min_child_weight==lowest_error$min_child_weight) 
        m_best$nb_cv <- 1:nrow(m_best) 
        cv_metrics <- matrix(ncol = 8, 
                             nrow = max(m_best$nb_cv),
                             dimnames=list(c(), c("eta", "max_depth","min_child_weight","RMSE", "MAE","R2","R2_correl","R2_correl_log")))
        cv_metrics=data.frame(cv_metrics)
        if (config_region$zone[r] == "world") {
          cv_metrics_inside <- cv_metrics
          cv_metrics_outside <- cv_metrics
        }
        counter = 0
        for (i in 1:max(m_best$nb_cv)) {
          counter = counter + 1
          sub_table <- m_best[i,] %>% 
            xgb_predict(niter=lowest_error$iter, newdata=d_test, fns=list(mean=mean, sd=sd, se=se)) %>%
            summarise(pred_metrics(pred_mean, biomass_count))
          test <-m_best[i,] %>% 
            xgb_predict(niter=lowest_error$iter, newdata=d_test, fns=list(mean=mean, sd=sd, se=se)) 
          test_results <- cor.test(test$pred_mean, test$biomass_count, method = "pearson", exact = FALSE)
          cv_metrics[i,] = sub_table
          
          # Predict on the regions 
          if (config_region$zone[r] == "world") {
            sub_table_inside <- m_best[i,] %>% 
              xgb_predict(niter=lowest_error$iter, newdata=d_test_inside, fns=list(mean=mean, sd=sd, se=se)) %>%
              summarise(pred_metrics(pred_mean, biomass_count))
            cv_metrics_inside[i,] = sub_table_inside
            
            sub_table_outside <- m_best[i,] %>% 
              xgb_predict(niter=lowest_error$iter, newdata=d_test_outside, fns=list(mean=mean, sd=sd, se=se)) %>%
              summarise(pred_metrics(pred_mean, biomass_count))
            cv_metrics_outside[i,] = sub_table_outside
          }
        }
        summary(cv_metrics)
        
        cv_metrics_mean <- cv_metrics %>% summarise(RMSE = mean(RMSE), R2_correl=mean(R2_correl))
        if (config_region$zone[r] == "world") {
          cv_metrics_inside_mean <- cv_metrics_inside %>% summarise(RMSE = mean(RMSE), R2_correl=mean(R2_correl))
          cv_metrics_outside_mean <- cv_metrics_outside %>% summarise(RMSE = mean(RMSE), R2_correl=mean(R2_correl))
        }  
        
        # do and save some stats on the mean model from the 100 CV models
        predict_test_set <- m_best %>% 
          xgb_predict(niter=lowest_error$iter, newdata=d_test, fns=list(mean=mean, sd=sd, se=se)) 
        results_test_pearson <- cor.test(predict_test_set$pred_mean, predict_test_set$biomass_count, method = "pearson")
        results_test_spearman <- cor.test(predict_test_set$pred_mean, predict_test_set$biomass_count, method = "spearman", exact=FALSE)
        cor_table[counter_cor,"taxa"] = t
        cor_table[counter_cor,"p_value_pearson"] = results_test_pearson$p.value
        cor_table[counter_cor,"cor_pearson"] = results_test_pearson$estimate
        cor_table[counter_cor,"p_value_spearman"] = results_test_spearman$p.value
        cor_table[counter_cor,"cor_spearman"] = results_test_spearman$estimate
        
        predict_train_set <- m_best %>% 
          xgb_predict(niter=lowest_error$iter, newdata=d_train, fns=list(mean=mean, sd=sd, se=se)) 
        
        ## Figure of biomass seen vs predicted
        d_test_set <- d_test %>%  left_join(predict_test_set)
        d_train_set <- d_train %>%  left_join(predict_train_set)
        fig_test <- ggplot(d_test_set) + coord_fixed() +
          geom_abline(intercept=0, slope=1, colour="#3B87B3") +  # 1:1 line
          geom_point(aes(biomass_count, pred_mean), size=0.5, alpha=0.5) +  # data points
          scale_x_continuous(trans="log1p") + scale_y_continuous(trans="log1p") +  # nicer aspect
          labs(x="True biomass", y="Predicted biomass") +
          ggtitle("test") +
          theme(plot.title = element_text(size = 5),
                axis.text=element_text(size=4),
                axis.title=element_text(size=5))
        fig_train <- ggplot(d_train_set) + coord_fixed() + 
          geom_abline(intercept=0, slope=1, colour="#3B87B3") +  # 1:1 line
          geom_point(aes(biomass_count, pred_mean), size=0.5, alpha=0.5) + # data points
          scale_x_continuous(trans="log1p") + scale_y_continuous(trans="log1p") +  # nicer aspect
          labs(x="True biomass", y="Predicted biomass") +
          ggtitle("train") +
          theme(plot.title = element_text(size = 5),
                axis.text=element_text(size=4),
                axis.title=element_text(size=5))
        plot_train_test <- grid.arrange(fig_test, fig_train, ncol=2)
        ggsave(file=paste0("results/",config$element,"/predicted_vs_true/",config_region$zone[r], "/",start_layer,"_",end_layer, "_", config$sd_or_no_sd,"/",t,"_",start_layer,"_",end_layer, "m_predicted_vs_true.pdf"), 
               plot = plot_train_test, 
               width=18, height=10, units = "cm")
        
        # predict over the whole world
        prediction <- filter(m, eta==lowest_error$eta, max_depth==lowest_error$max_depth,
                             min_child_weight==lowest_error$min_child_weight) %>% 
          ungroup() %>% 
          xgb_predict(niter=lowest_error$iter, newdata=env, fns=list(mean=mean, sd=sd, se=se)) %>% 
          # NB: remember that the biomass was multiplied by scale_f to be compatible with
          #     Poisson counts, at the very beginning
          mutate(pred_mean = pred_mean/scale_f,
                 pred_sd = pred_sd/scale_f,
                 pred_se = pred_se/scale_f,
                 pred_cv = pred_sd/pred_mean)
        # map it
        ggplot() + coord_quickmap() +
          geom_raster(aes(lon, lat, fill=pred_mean), data=prediction) +
          scale_fill_viridis_c(trans="sqrt") +
          geom_polygon(aes(lon, lat), data=coast)
        
        # store results in a list
        if (config_region$zone[r] != "world") {
          res <- list(
            # identifier
            taxon = t,
            # data
            scale_factor=scale_f,
            d_train = d_train,
            d_test = d_test,
            # models
            model_full=m_best,
            lowest_error=lowest_error,
            # predictions
            cv_metrics=cv_metrics,
            prediction=prediction
          )
        }
        
        if (config_region$zone[r] == "world") {
          res <- list(
            # identifier
            taxon = t,
            # data
            scale_factor=scale_f,
            d_train = d_train,
            d_test = d_test,
            d_test_inside = d_test_inside,
            d_test_outside = d_test_outside,
            # models
            model_full=m_best,
            lowest_error=lowest_error,
            # predictions
            cv_metrics=cv_metrics,
            cv_metrics_inside = cv_metrics_inside,
            cv_metrics_outside = cv_metrics_outside,
            prediction=prediction
          )
        }
        
        # and save it
        if (config$region == "region") {
          saveRDS(res, file=paste0(results_dir, t, "_",config_region$zone[r], "_",start_layer,"_",end_layer,"m.rds"))
        }
      }
    }
  }
  write.csv(cor_table, file=paste0("results/", config$element,"/cor_table/", config_region$zone[r], "_",start_layer,"_",end_layer,"m.csv"))
}