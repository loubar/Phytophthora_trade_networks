
rm(list=ls())
library(brms)
#require(dplyr)
#require(reshape2)
library(tidybayes)
library(ggplot2)
library(sf)
library(ggpubr)
library(viridis)

minRE_spp <- "arrivals_minREspp0"
com_type <- "trade_hort"

rank_models_all <- readRDS(paste0(minRE_spp,"/rank_models_all.rds"))
rank_models_all$model <- gsub("_rslurm_trade_hort_2023/", "", rank_models_all$model)
model_info <- readRDS(paste0(minRE_spp,"/model_info.rds"))
rownames(model_info) <- rank_models_all$model[rank_models_all$diffic >=2]
#model_summary_df <- readRDS(paste0(minRE_spp,"/model_summary_df.rds"))
arrivals_unscaled <- readRDS(paste0(minRE_spp, "/", "arrivals_unscaled.rds"))
arrivals <- readRDS(paste0(minRE_spp, "/", "arrivals.rds"))
included_traits <- readRDS("included_traits.rds")

best_models <- rank_models_all[rank_models_all$diffic >=-2,]


top_subset <- paste0(minRE_spp, "/", "best_models", "/", best_models$model)
top_models <- lapply(top_subset, readRDS)

####################### map predictions of emergence risk and export risk ######################

# using mnodel averaging across the top 13 best performing models 
# to produce the predictions

hort_tradeMatrix <- readRDS("hort_tradeMatrix.rds")
sourceMatrix <- readRDS("phytophthora_source_matrix.rds")
arrivalMatrix <- readRDS("phytophthora_arrival_matrix.rds")
# Phytophthora distributions
arrivalMatrix[is.na(arrivalMatrix)] <- 0
presenceMatrix <- sourceMatrix + arrivalMatrix
table(presenceMatrix, useNA = "always") # should be present or absent or NA (e.g. no countries should be both source and recipient == 2)
climate_dist <- readRDS("MD_mean.rds")
EPPO_reporting_service <- read.csv("biosecurity_EPPORS.csv", stringsAsFactors = FALSE)


# countries135 <- readRDS("countries135.rds")
# sourceMatrix <- sourceMatrix[countries135,]
# climate_dist <- climate_dist[countries135,countries135]
# hort_tradeMatrix <- hort_tradeMatrix[countries135, countries135]

# scaling
trade_scale <- 2 * sd(log(1 + arrivals_unscaled$trade_hort), na.rm = TRUE)
trade_center <- mean(log(1 + arrivals_unscaled$trade_hort), na.rm = TRUE)

climate_scale <- 2 * sd(-sqrt(arrivals_unscaled$climate_distance), na.rm = TRUE)
climate_center <- mean(-sqrt(arrivals_unscaled$climate_distance), na.rm = TRUE)

surv_scale <- 2 * sd(log(arrivals_unscaled$EPPO_reporting_service), na.rm = TRUE)
surv_center <- mean(log(arrivals_unscaled$EPPO_reporting_service), na.rm = TRUE)

CH_OO_scale <- 2*sd(arrivals_unscaled$CH_OO)
CH_OO_center <- mean(arrivals_unscaled$CH_OO)

temp_range_scale <- 2* sd(arrivals_unscaled$temp_range)
temp_range_center <- mean(arrivals_unscaled$temp_range)

date_scale <- 2*sd(log(1+(2022 - arrivals_unscaled$date)))
date_center <- mean(log(1 + (2022 - arrivals_unscaled$date)))

### map emergence risk by country: risk of new plant pathogen arrivals/existing non-detections ###
# sum of imports from all regions

#create newdata to predict composite risk for all 135 countries with trade and climate-matching data
# we will adjust for biosecurity effort (fixing this value at the maximum 
# observed to try and estimate relative risk of new Phytophthora detections in the 
# absence of surveillance bias



## marginal risk from trade connectivity (fixing all other risk factors at mean)
## species are unknown (e.g. drawn from the estimated species-level random effect distributions)

newdata_trade_risk <- data.frame(importer_iso3 = names(colSums(hort_tradeMatrix,
                                                               na.rm = TRUE)),
                                 trade_hort = (log(1+colSums(hort_tradeMatrix,
                                                             na.rm = TRUE))-trade_center)/trade_scale, 
                                 phyto = "unknown_species",
                                 phylo = "unknown_species",
                                 CH_OO = median(arrivals$CH_OO),
                                 temp_range = median(arrivals$temp_range),
                                 climate_match = median(arrivals$climate_match),
                                 EPPO_reporting_service = max(arrivals$EPPO_reporting_service),
                                 date = median(arrivals$date),
                                 stringsAsFactors = FALSE)

newdata_climate_risk <- data.frame(importer_iso3 = colnames(climate_dist),
                                   trade_hort = median(arrivals$trade_hort), 
                                   phyto = "unknown_species",
                                   phylo = "unknown_species",
                                   CH_OO = median(arrivals$CH_OO),
                                   temp_range = median(arrivals$temp_range),
                                   climate_match = ((-sqrt(colMeans(climate_dist, na.rm = TRUE)))-climate_center)/climate_scale,
                                   EPPO_reporting_service = max(arrivals$EPPO_reporting_service),
                                   date = max(arrivals$date),
                                   stringsAsFactors = FALSE)

newdata_surveillance <- data.frame(importer_iso3 = EPPO_reporting_service$importer_iso3,
                                   trade_hort = median(arrivals$trade_hort), 
                                   phyto = "unknown_species",#paste0("unknown_species", 1:length(EPPO_reporting_service$importer_iso3)),
                                   phylo = "unknown_species",
                                   CH_OO = median(arrivals$CH_OO),
                                   temp_range = median(arrivals$temp_range),
                                   climate_match = median(arrivals$climate_match),
                                   EPPO_reporting_service = (log(EPPO_reporting_service$EPPO_reporting_service)-surv_center)/surv_scale,
                                   date = max(arrivals$date),
                                   stringsAsFactors = FALSE)


# here we are fixing recording effort a high level to attempt to remove this effect when predicting risk
newdata_overall_risk <- left_join(newdata_trade_risk[,grep("climate_match", names(newdata_trade_risk), invert = TRUE)], 
                                  newdata_climate_risk[,c("importer_iso3", "climate_match")])

newdata_overall_risk <- inner_join(newdata_overall_risk[,grep("EPPO_reporting_service", 
                                                              names(newdata_overall_risk), 
                                                              invert = TRUE)], 
                                   newdata_surveillance[,c("importer_iso3","EPPO_reporting_service")])


newdata_emergence_risk <- newdata_overall_risk
newdata_emergence_risk$EPPO_reporting_service <- median(arrivals$EPPO_reporting_service)
newdata_emergence_risk$date <- median(arrivals$date)
# risk map of non-dectection
# predict zero-inflation (non-detection)

newdata_zi <- data.frame(importer_iso3 = EPPO_reporting_service$importer_iso3,
                         trade_hort = median(arrivals$trade_hort), 
                         phyto = "unknown_species",#paste0("unknown_species", 1:length(EPPO_reporting_service$importer_iso3)),
                         phylo = "unknown_species",
                         CH_OO = median(arrivals$CH_OO),
                         temp_range = median(arrivals$temp_range),
                         climate_match = median(arrivals$climate_match),
                         EPPO_reporting_service = (log(EPPO_reporting_service$EPPO_reporting_service)-surv_center)/surv_scale,
                         date = max(arrivals$date),
                         stringsAsFactors = FALSE)

#newdata_zi <- newdata_zi[newdata_zi$EPPO_reporting_service > quantile(newdata_zi$EPPO_reporting_service, c(0.30)),]

# based on advice in
# https://www.andrewheiss.com/blog/2021/11/08/beta-regression-guide/#zero-inflated-beta-regression-bayesian-style

# model_av includes only the predictors that are present in all 13 best-perfroming models
# these are
colnames(posterior_average(top_models[[1]],
                           top_models[[2]], 
                           top_models[[3]],
                           top_models[[4]],
                           top_models[[5]], 
                           top_models[[6]], 
                           top_models[[7]],
                           top_models[[8]],
                           top_models[[9]], 
                           top_models[[10]], 
                           top_models[[11]],
                           top_models[[12]],
                           top_models[[13]],
                           weights = "loo",
                           missing = 0))

# not adjusting for species- and country-level variance/uncertainty
# just overall risk based on deterministic drivers
pred_draws <- pp_average(top_models[[1]],
                         top_models[[2]], 
                         top_models[[3]],
                         top_models[[4]],
                         top_models[[5]], 
                         top_models[[6]], 
                         top_models[[7]],
                         top_models[[8]],
                         top_models[[9]], 
                         top_models[[10]], 
                         top_models[[11]],
                         top_models[[12]],
                         top_models[[13]],
                         weights = "loo", method = "posterior_epred",
                         newdata = newdata_emergence_risk[complete.cases(newdata_emergence_risk),],
                         dpar = "mu", # excluding zero-inflation/non-detection 
                         ndraws = 1000,
                         re_formula = NA,
                         #allow_new_levels = TRUE,
                         #sample_new_levels = "gaussian",
                         summary = FALSE,
                         missing = 0) %>% as_draws_df() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  # group_by(importer_iso3) %>% 
  # summarize(`50%` = median(zi),
  #           `2.5%` = quantile(zi, c(0.025)),
  #           `97.5` = quantile(zi, c(0.975))) %>%
  mutate(GU_A3 = newdata_emergence_risk[complete.cases(newdata_emergence_risk),]$importer_iso3)

# pred_draws <- best_model %>% 
#   epred_draws(newdata = newdata_zi,
#             dpar = "zi",
#             allow_new_levels = TRUE,
#             sample_new_levels = "gaussian"#,
#             #re_formula = NA
#             ) %>%
#   group_by(importer_iso3) %>% 
#   summarize(`50%` = median(zi),
#             `2.5%` = quantile(zi, c(0.025)),
#             `97.5` = quantile(zi, c(0.975))) %>%
#   rename(GU_A3 = importer_iso3)


world_sf <- st_read(dsn = "ne_50m_admin_0_countries", layer = "ne_50m_admin_0_countries",
                    stringsAsFactors = FALSE)
#setdiff(pred_trade_risk$GU_A3, world_sf$GU_A3) # check all 135 countries present on map
risk_map_overall <- left_join(world_sf, pred_draws)
risk_map_overall$risk_factor <- "Composite risk"
risk_map_overall$uncertainty <- risk_map_overall$`97.5%` - risk_map_overall$`2.5%`

# predict zero-inflation/non-detection
# fixing all variables except recording effort at 0
pred_draws_zi <- pp_average(top_models[[1]],
                         top_models[[2]], 
                         top_models[[3]],
                         top_models[[4]],
                         top_models[[5]], 
                         top_models[[6]], 
                         top_models[[7]],
                         top_models[[8]],
                         top_models[[9]], 
                         top_models[[10]], 
                         top_models[[11]],
                         top_models[[12]],
                         top_models[[13]],
                         weights = "loo", method = "posterior_epred",
                         newdata = newdata_zi,
                         dpar = "zi", # excluding zero-inflation/non-detection 
                         ndraws = 1000,
                         allow_new_levels = TRUE,
                         sample_new_levels = "gaussian",
                         summary = FALSE,
                         missing = 0) %>% as_draws_df() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  # group_by(importer_iso3) %>% 
  # summarize(`50%` = median(zi),
  #           `2.5%` = quantile(zi, c(0.025)),
  #           `97.5` = quantile(zi, c(0.975))) %>%
  mutate(GU_A3 = newdata_zi$importer_iso3)

world_sf <- st_read(dsn = "ne_50m_admin_0_countries", layer = "ne_50m_admin_0_countries",
                    stringsAsFactors = FALSE)
#setdiff(pred_trade_risk$GU_A3, world_sf$GU_A3) # check all 135 countries present on map
risk_map_zi <- left_join(world_sf, pred_draws_zi)
risk_map_zi$risk_factor <- "Non-detection"

key_size <- 0.25
position <- c(0.06,0.3)
legend_text <- 5
# plot the risk map of overall emergence risk 
overall_risk_map <- ggplot(data = risk_map_overall, aes(fill = `50%`)) + # create a ggplot object and 
  # change its fill colour according to median_age
  geom_sf(size = 0.5) + scale_fill_viridis_c("Probability", option = "magma") +
  theme(legend.text = element_text(size = legend_text),
        legend.title = (element_text(size = legend_text)),
    legend.position.inside  = position,
    #legend.background = element_rect(fill ="white", colour = "black"),
    legend.key.size = unit(key_size, units = "cm") ,
    #legend.position = "bottom",    
    axis.text = element_text(size = 7)) + coord_sf(expand = FALSE) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = c( -120, -60, 0, 60, 120)) +
  scale_y_continuous(limits = c(-60, 80),
                     breaks = c(-60, -40, -20, 0, 20, 40, 60, 80)) +
  facet_wrap(~risk_factor, nrow = 1)



# risk map of known Phytophthora diversity
sf_use_s2(FALSE)
known_diversity <- readRDS("phyto_diversity.rds")
#setdiff()
#known_diversity <- left_join(newdata_zi[,"importer_iso3", drop = FALSE], known_diversity)
#known_diversity$known_diversity[is.na(known_diversity$known_diversity)] <- 0
#known_diversity$`50%` <- known_diversity$known_diversity / ncol(presenceMatrix)

# risk_map_diversity <- epred_draws(best_model,
#                                   newdata = newdata_zi,
#                                   dpar = "zi",
#                                   #allow_new_levels = TRUE,
#                                   #sample_new_levels = "gaussian",
#                                   re_formula = NA) %>%
#   mutate(predicted_diversity = .epred / )
world_sf <- st_read(dsn = "ne_50m_admin_0_countries", 
                    layer = "ne_50m_admin_0_countries", stringsAsFactors = FALSE)
risk_map_diversity <- left_join(world_sf, known_diversity %>% rename(GU_A3 = importer_iso3))
risk_map_diversity$risk_factor <- "'Reported'~italic('Phytophthora')~'diversity'"





lenUnique <- function(x){length(unique(x))}
sf_use_s2(FALSE)
library(countrycode)
map_risk_factors <- function(newdata_marginal, type = "map", pred = "marginal"){
  # check which variable is being conditioned on (defined in newdata_)
  marg_eff <- apply(newdata_marginal[, c("trade_hort", 
                                         "climate_match",
                                         "EPPO_reporting_service")], 2, lenUnique)
  marg_eff <- names(which(marg_eff >1))
  
  if(length(marg_eff) == 3){
    risk_factor <- "Composite risk of detection"
  }
  if(length(marg_eff) == 2){
    risk_factor <- "Composite risk of arrival"
  }
  if(length(marg_eff) == 1 & marg_eff[1] == "climate_match"){
    risk_factor <- "Climate matching"
  }
  if(length(marg_eff) == 1 & marg_eff[1] == "EPPO_reporting_service"){
    risk_factor <- "Plant health surveillance"
  }
  if(length(marg_eff) == 1 & marg_eff[1] == "trade_hort"){
    risk_factor <- "Horticultural imports"
  }
  
  if(pred == "marginal"){
    pred_draws <- pp_average(top_models[[1]],
                             top_models[[2]], 
                             top_models[[3]],
                             top_models[[4]],
                             top_models[[5]], 
                             top_models[[6]], 
                             top_models[[7]],
                             top_models[[8]],
                             top_models[[9]], 
                             top_models[[10]], 
                             top_models[[11]],
                             top_models[[12]],
                             top_models[[13]],
                             weights = "loo", method = "posterior_epred",
                             newdata = newdata_marginal,
                             dpar = "mu",
                             missing = 0,
                             ndraws = 1000,
                             re_formula = NA,
                             summary = FALSE) %>% as_draws_df() %>% 
      summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
      # group_by(importer_iso3) %>% 
      # summarize(`50%` = median(zi),
      #           `2.5%` = quantile(zi, c(0.025)),
      #           `97.5` = quantile(zi, c(0.975))) %>%
      mutate(GU_A3 = newdata_marginal$importer_iso3)
  }
  if(pred == "conditional"){
    pred_draws <- pp_average(top_models[[1]],
                             top_models[[2]], 
                             top_models[[3]],
                             top_models[[4]],
                             top_models[[5]], 
                             top_models[[6]], 
                             top_models[[7]],
                             top_models[[8]],
                             top_models[[9]], 
                             top_models[[10]], 
                             top_models[[11]],
                             top_models[[12]],
                             top_models[[13]],
                             weights = "loo", method = "posterior_epred",
                             newdata = newdata_marginal,
                             ndraws = 1000,
                             missing = 0,
                             dpar = "mu", # predict based on the process model only (factoring out non-detections)
                             allow_new_levels = TRUE,
                             sample_new_levels = "gaussian",
                             summary = FALSE) %>% as_draws_df() %>% 
      summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
      # group_by(importer_iso3) %>% 
      # summarize(`50%` = median(zi),
      #           `2.5%` = quantile(zi, c(0.025)),
      #           `97.5` = quantile(zi, c(0.975))) %>%
      mutate(GU_A3 = newdata_marginal$importer_iso3)
  }
  
  if(type == "map"){
    world_sf <- st_read(dsn = "ne_50m_admin_0_countries",
                        layer ="ne_50m_admin_0_countries", stringsAsFactors = FALSE)
    #world_sf <- ne_countries(returnclass = "sf", type = "map_units")
    #setdiff(pred_trade_risk$GU_A3, world_sf$GU_A3) # all 135 countries present on map
    risk_map <- left_join(world_sf, pred_draws)
    risk_map$risk_factor <- risk_factor
    return(risk_map)
  } 
  if(type == "table"){
    pred_draws$risk_factor <- risk_factor
    return(pred_draws)
  }
}

newdata_list <- list(newdata_trade_risk,
                     newdata_climate_risk,
                     newdata_surveillance)

# risk_table <- lapply(1:length(newdata_list), function(i)
#   map_risk_factors(newdata_marginal = newdata_list[[i]],
#                    type = "table"))

marginal_risk_maps <- lapply(newdata_list,
                             map_risk_factors, type = "map", 
                             pred = "marginal")



marginal_risk_map_panels <- bind_rows(marginal_risk_maps)
#marginal_risk_map_panels <- bind_rows(marginal_risk_map_panels, risk_map_zi)
#marginal_risk_map_panels <- bind_rows(marginal_risk_map_panels, risk_map_diversity)
#marginal_risk_map_panels <- bind_rows(risk_map_zi, marginal_risk_map_panels)
marginal_risk_map_panels$risk_factor <- gsub(" ", "~", marginal_risk_map_panels$risk_factor)
marginal_risk_map_panels$risk_factor <- factor(marginal_risk_map_panels$risk_factor, 
                                               levels = c(#"'Reported'~italic('Phytophthora')~'diversity'",
                                                          "Horticultural~imports",
                                                          "Climate~matching",
                                                          "Plant~health~surveillance"#,
                                                          #"Non-detection",
                                                          #"Composite~risk~of~detection"
                                               ))



# sf_use_s2(FALSE) # 
# 
# marginal_risk_map_plot <- ggplot(data = marginal_risk_map_panels, aes(fill = `50%`) ) + # create a ggplot object and 
#   # change its fill colour according to median_age
#   geom_sf(size = 0.1) +scale_fill_viridis("Probability", option = "magma") +
#   theme(#legend.title = element_text(size = 5),
#     #legend.text = element_text(size = 4),
#     #legend.position = c(0.1,0.4),
#     #legend.background = element_rect(fill ="white", colour = "black"),
#     #legend.key.size = unit(0.2, units = "cm") ,
#     legend.position = "right",    
#     axis.text = element_text(size = 7)) + coord_sf(expand = FALSE) +
#   scale_x_continuous(limits = c(-180, 180),
#                      breaks = c( -120, -60, 0, 60, 120)) +
#   scale_y_continuous(limits = c(-60, 80),
#                      breaks = c(-60, -40, -20, 0, 20, 40, 60, 80)) +
#   facet_wrap(~risk_factor, labeller = label_parsed, ncol = 2)


# using unique scales for each plot (to highlight differences in climate matching, despite weaker effect)

climate_risk_map_plot <- ggplot(data = marginal_risk_maps[[2]], aes(fill = `50%`) ) + # create a ggplot object and 
  # change its fill colour according to median_age
  geom_sf(size = 0.1) +scale_fill_viridis("Probability\n", direction = -1) +
  theme(legend.text = element_text(size = legend_text),
        legend.title = (element_text(size = legend_text)),
        legend.position = position,
        #legend.background = element_rect(fill ="white", colour = "black"),
        legend.key.size = unit(key_size, units = "cm") ,
        #legend.position = "bottom",    
        axis.text = element_text(size = 7)) + coord_sf(expand = FALSE) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = c( -120, -60, 0, 60, 120)) +
  scale_y_continuous(limits = c(-60, 80),
                     breaks = c(-60, -40, -20, 0, 20, 40, 60, 80)) +
  facet_wrap(~risk_factor, nrow = 1)

# risk_factor_climate <- left_join(world_sf, newdata_climate_risk %>% rename(GU_A3 = importer_iso3) )
# risk_factor_climate_map <- ggplot(data = risk_factor_climate, aes(fill = climate_match * climate_scale + climate_center) ) + # create a ggplot object and 
#   # change its fill colour according to median_age
#   geom_sf(size = 0.1) +scale_fill_viridis(expression(Climate~matching~(-sqrt(Mahalanobis~distance))), 
#                                           option = "magma",
#                                           breaks = c(-5:-1),
#                                           labels = c(-5:-1)) +
#   theme(#legend.title = element_text(size = 5),
#     #legend.text = element_text(size = 4),
#     #legend.position = c(0.1,0.4),
#     #legend.background = element_rect(fill ="white", colour = "black"),
#     #legend.key.size = unit(0.2, units = "cm") ,
#     legend.position = "bottom",    
#     axis.text = element_text(size = 7)) + coord_sf(expand = FALSE) +
#   scale_x_continuous(limits = c(-180, 180),
#                      breaks = c( -120, -60, 0, 60, 120)) +
#   scale_y_continuous(limits = c(-60, 80),
#                      breaks = c(-60, -40, -20, 0, 20, 40, 60, 80))+
#   facet_wrap(~risk_factor, nrow = 1)

trade_risk_map_plot <- ggplot(data = marginal_risk_maps[[1]], aes(fill = `50%`) ) + # create a ggplot object and 
  # change its fill colour according to median_age
  geom_sf(size = 0.1) +scale_fill_viridis("Probability\n ", direction = -1) +
  theme(legend.text = element_text(size = legend_text),
        legend.title = (element_text(size = legend_text)),
        legend.position = position,
        #legend.background = element_rect(fill ="white", colour = "black"),
        legend.key.size = unit(key_size, units = "cm") ,
        #legend.position = "bottom",    
        axis.text = element_text(size = 7)) + coord_sf(expand = FALSE) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = c( -120, -60, 0, 60, 120)) +
  scale_y_continuous(limits = c(-60, 80),
                     breaks = c(-60, -40, -20, 0, 20, 40, 60, 80)) +
  facet_wrap(~risk_factor, nrow = 1)


surv_risk_map_plot <- ggplot(data = marginal_risk_maps[[3]], aes(fill = `50%`) ) + # create a ggplot object and 
  # change its fill colour according to median_age
  geom_sf(size = 0.1) +scale_fill_viridis("Probability\n", direction = -1) +
  theme(legend.text = element_text(size = legend_text),
        legend.title = (element_text(size = legend_text)),
        legend.position = position,
        #legend.background = element_rect(fill ="white", colour = "black"),
        legend.key.size = unit(key_size, units = "cm") ,
        #legend.position = "bottom",    
        axis.text = element_text(size = 7)) + coord_sf(expand = FALSE) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = c( -120, -60, 0, 60, 120)) +
  scale_y_continuous(limits = c(-60, 80),
                     breaks = c(-60, -40, -20, 0, 20, 40, 60, 80)) +
  facet_wrap(~risk_factor, nrow = 1)


nondetection_risk_map <- ggplot(data = risk_map_zi, aes(fill = `50%`) ) + # create a ggplot object and 
  # change its fill colour according to median_age
  geom_sf(size = 0.1) + scale_fill_viridis("Probability\n ", direction = -1) +
  theme(legend.text = element_text(size = legend_text),
        legend.title = (element_text(size = legend_text)),
        legend.position = position,
        #legend.background = element_rect(fill ="white", colour = "black"),
        legend.key.size = unit(key_size, units = "cm") ,
        #legend.position = "bottom",    
        axis.text = element_text(size = 7)) + coord_sf(expand = FALSE) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = c( -120, -60, 0, 60, 120)) +
  scale_y_continuous(limits = c(-60, 80),
                     breaks = c(-60, -40, -20, 0, 20, 40, 60, 80)) +
  facet_wrap(~risk_factor, nrow = 1)

diversity_risk_map <- ggplot(data = risk_map_diversity, aes(fill = phyto_diversity) ) + # create a ggplot object and 
  # change its fill colour according to median_age
  geom_sf(size = 0.1) + scale_fill_distiller("Detected\nspecies", palette = "Reds" ,
                                             direction = -1) +
  theme(legend.text = element_text(size = legend_text),
        legend.title = (element_text(size = legend_text)),
        legend.position = position,
        #legend.background = element_rect(fill ="white", colour = "black"),
        legend.key.size = unit(key_size, units = "cm") ,
        #legend.position = "bottom",    
        axis.text = element_text(size = 7)) + coord_sf(expand = FALSE) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = c( -120, -60, 0, 60, 120)) +
  scale_y_continuous(limits = c(-60, 80),
                     breaks = c(-60, -40, -20, 0, 20, 40, 60, 80)) +
  facet_wrap(~risk_factor, nrow = 1, labeller = label_parsed)

ggarrange(overall_risk_map,
  diversity_risk_map,
          trade_risk_map_plot,
          climate_risk_map_plot,
          nondetection_risk_map,
          surv_risk_map_plot,
          nrow = 3, ncol = 2)
ggsave(paste0(minRE_spp, "/", "risk_factor_maps_unique_prob_scales.png"),
       height = 6,
       width = 8,
       units = "in",
       dpi = 500)

# use the models to produce a ranked list of Phytophthora's likely to emerge in the UK

newdata_rank <- top_models[[1]]$data[top_models[[1]]$data$importer_iso3 == "GBR",]
# 76 species
# for those without a source region prior to 2005, 
# can we use the average of trade connectivity and climate matching with
# all source regions?


pred_threats_uk <- pp_average(top_models[[1]],
                         top_models[[2]], 
                         top_models[[3]],
                         top_models[[4]],
                         top_models[[5]], 
                         top_models[[6]], 
                         top_models[[7]],
                         top_models[[8]],
                         top_models[[9]], 
                         top_models[[10]], 
                         top_models[[11]],
                         top_models[[12]],
                         top_models[[13]],
                         weights = "loo", method = "posterior_epred",
                         newdata = newdata_rank,
                         ndraws = 1000,
                         missing = 0,
                         #re_formula = NA, # remove species- and country-level variation?
                         dpar = "mu", # predict based on the process model only (factoring out non-detections)
                         #allow_new_levels = TRUE,
                         #sample_new_levels = "gaussian",
                         summary = FALSE) %>% as_draws_df() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  bind_cols(newdata_rank[,c("phyto", "importer_iso3", "arrival")])

rank_threats_uk <- pred_threats_uk[order(pred_threats_uk$`50%`, decreasing = TRUE),] 

# this gives predictions for 76 Phytophthora species with known source region(s)
# prior to 2005
all_threats <-
  rank_threats_uk %>% 
  mutate(risk = paste0(sprintf(`50%`, fmt = '%#.2f'),
                       " (",
                       sprintf(`2.5%`, fmt = '%#.2f'),
                       ", ",
                       sprintf(`97.5%`, fmt = '%#.2f'),
                       ")"),
         threat = ifelse(arrival == 1, yes = "Current (present)",
                         no = "Anticipated (not yet reported)")) %>%
  select(phyto, threat, risk) 
write.csv(all_threats, file = paste0(minRE_spp, "/", "threats_uk_76spp.csv"),
          row.names = FALSE)


## current threats (already present)
current_threats_uk <- rank_threats_uk[rank_threats_uk$arrival == 1,][1:20,] %>% 
  mutate(threat = "Current")

#### future threats (e.g. ranking predictions for species that have yet to arrive)
future_threats_uk <- rank_threats_uk[rank_threats_uk$arrival == 0,][1:20,] %>% 
  mutate(threat = "Anticipated")

threats_uk <- bind_rows(current_threats_uk, future_threats_uk) %>%
  mutate(risk = paste0(sprintf(`50%`, fmt = '%#.2f'),
                       " (",
                       sprintf(`2.5%`, fmt = '%#.2f'),
                       ", ",
                       sprintf(`97.5%`, fmt = '%#.2f'),
                       ")")) %>%
  select(phyto, threat, risk) 

# join with impact metrics
impact_metrics <- readRDS("impact_metrics.rds") %>%
  mutate(phyto = gsub("×", "x", phyto)) %>%
  mutate(phyto = gsub("phytophthora", "Phytophthora", phyto))
# check all species are represented in the impacts data
setdiff(threats_uk$phyto, impact_metrics$phyto)

threats_uk_impact <- left_join(threats_uk, impact_metrics)
threats_uk_impact$phyto <- gsub("Phytophthora",
                                "P.",
                                threats_uk_impact$phyto)
threats_uk_impact$lat_limit <- sprintf(threats_uk_impact$lat_limit, fmt = '%#.2f')

# export the table for the manuscript
kbl(threats_uk_impact[,-2]) %>%
  kable_classic(full_width = T) %>%
  column_spec(1, italic = TRUE) %>%
  add_header_above(c(" " = 5, "Ecosystems impacted" = 6)) %>%
  #add_header_above(c(" " = 2, "Impact metrics" = 9)) %>%
  pack_rows(index = c("Current threats (present)" = 20,
                      "Anticipated threats (absent)" = 20)) %>%
  save_kable(file = paste0(minRE_spp, "/risky_species.html"),
             self_contained = TRUE)



# and other countries, too?


# map export risk to focal importers (GBR and China, Brazil) low, medium, high risk, respectively
importer_choice <- "GBR"





library(tibble)
map_export_risk <- function(importer){
  export_newdata <- data.frame(trade_hort = (log(1 + hort_tradeMatrix[,importer]) - trade_center) / trade_scale,
                               exporter_iso3 = names(hort_tradeMatrix[,importer]),
                               row.names = NULL,
                               stringsAsFactors = FALSE) 
  export_newdata <- left_join(export_newdata,
                              data.frame(exporter_iso3 = names(climate_dist[,importer]),
                                         climate_match = (-sqrt(climate_dist[,importer])-climate_center)/climate_scale,
                                         row.names = NULL,
                                         stringsAsFactors = FALSE)) %>%
    #left_join(EPPO_reporting_service %>% rename(exporter_iso3 = importer_iso3)) %>%
    mutate(EPPO_reporting_service = median(arrivals$EPPO_reporting_service),#(log(1+EPPO_reporting_service)-surv_center)/surv_scale,
           importer_iso3 = importer,
           date = median(arrivals$date),
           temp_range = median(arrivals$temp_range),
           CH_OO = median(arrivals$CH_OO),
           phylo = "unknown_species",
           phyto = "unknown_species")
    
  export_newdata <- export_newdata[complete.cases(export_newdata),]
  
  
  export_risk <- pp_average(top_models[[1]],
                                top_models[[2]], 
                                top_models[[3]],
                                top_models[[4]],
                                top_models[[5]], 
                                top_models[[6]], 
                                top_models[[7]],
                                top_models[[8]],
                                top_models[[9]], 
                                top_models[[10]], 
                                top_models[[11]],
                                top_models[[12]],
                                top_models[[13]],
                                weights = "loo", method = "posterior_epred",
                                newdata = export_newdata,
                                ndraws = 1000,
                                missing = 0,
                                re_formula = NA, # remove species- and country-level variation?
                                dpar = "mu", # predict based on the process model only (factoring out non-detections)
                                # allow_new_levels = TRUE,
                                # sample_new_levels = "gaussian",
                                summary = FALSE) %>% as_draws_df() %>% 
    summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
    bind_cols(export_newdata[,c("exporter_iso3", "importer_iso3")]) %>%
    rename(GU_A3 = exporter_iso3)
  # for species and countries that are not in the original dataset, posterior samples can be drawn from
  # the multivariate normal distribution defined by the group-level standard deviations and correlations
  # 
  
  
  high_risk_exporters <- export_risk[order(export_risk$`50%`, decreasing = TRUE),]
  high_risk_exporters$rank <- 1:nrow(high_risk_exporters)
  high_risk_exporters$country <- countrycode(sourcevar = high_risk_exporters$GU_A3,
                                             origin = "iso3c",
                                             destination = "country.name")
  
  
  
  
  world_sf <- st_read(dsn = "ne_50m_admin_0_countries", 
                      layer = "ne_50m_admin_0_countries", 
                      stringsAsFactors = FALSE)
  #world_sf <- ne_countries(returnclass = "sf", type = "map_units")
  #setdiff(high_risk_exporters$GU_A3, world_sf$GU_A3) # all 135 countries present on map
  exporter_risk_map <- left_join(world_sf, high_risk_exporters)
  #exporter_risk_map$importer_iso3 <- importer
  return(exporter_risk_map)
  
  
}
GB_source_risk <- map_export_risk(importer = "GBR") 
GB_source_risk[which.max(GB_source_risk$`50%`),"country"]
# by continent
GB_source_risk[,c("country", "CONTINENT", "50%")] %>% 
  st_drop_geometry() %>% group_by(CONTINENT) %>% summarise(greatest = max(`50%`, na.rm = TRUE))  

GB_source_risk_summary <- GB_source_risk_summary[!is.infinite(GB_source_risk_summary$greatest),]
GB_source_risk_summary[order(GB_source_risk_summary$greatest, decreasing = TRUE), ]

library(viridis)
### map the high risk source regions for the UK
ggplot(data = GB_source_risk, aes(fill = `50%`) ) + # create a ggplot object and 
  # change its fill colour according to median
  geom_sf(size = 0.1) + scale_fill_viridis("Relative risk of arrival by source region") +
  # geom_sf(data = GB_source_risk[GB_source_risk$GU_A3 %in% "GBR"],
  #         size = 0.1, fill = "black", col = "black") +
  theme(#legend.title = element_text(size = 5),
    #legend.text = element_text(size = 4),
    #legend.position = c(0.1,0.4),
    #legend.background = element_rect(fill ="white", colour = "black"),
    #legend.key.size = unit(0.2, units = "cm") ,
    legend.position = "top",    
    axis.text = element_text(size = 7)) + coord_sf(expand = FALSE) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = c( -120, -60, 0, 60, 120)) +
  scale_y_continuous(limits = c(-60, 80),
                     breaks = c(-60, -40, -20, 0, 20, 40, 60, 80))
ggsave("UK_importation_risk.png",
       dpi = 500, units = "in")

## ranked list of Phytophthora threats to the UK
#newdata_UK_Phytophthora <- data.frame(phyto = unique(arrivals$))


# exp_draws <- epred_draws(best_model, newdata = export_newdata[complete.cases(export_newdata),], 
#                          allow_new_levels = TRUE,
#                          sample_new_levels = "gaussian") %>%
#   group_by(exporter_iso3) %>% 
#   summarize(`50%` = median(.epred),
#             `2.5%` = quantile(.epred, c(0.25)),
#             `97.5%` = quantile(.epred, c(0.975))) %>%
#   rename(GU_A3 = exporter_iso3)

# choose exmaple countries based on the predictions in Fig. 4c and 4d

# rank_importers <- conditional_risk_maps[order(conditional_risk_maps$`50%`, decreasing = TRUE),c("GU_A3", "50%", "2.5%", "97.5")]
# rank_importers$rank <- 1:nrow(rank_importers)

export_risk_maps <- lapply(c("ESP", "BRA", "USA", "CHN", "ZAF", "GBR"),  
                           map_export_risk)

export_risk_maps_panels <- bind_rows(export_risk_maps)
export_risk_maps_panels$country <- countrycode(sourcevar = export_risk_maps_panels$importer_iso3,
                                               origin = "iso3c",
                                               destination = "country.name")

export_risk_maps_panels$country <- factor(export_risk_maps_panels$country, 
       levels = countrycode(sourcevar = c("ESP", "BRA", "USA", "CHN", "ZAF", "GBR"),
                            origin = "iso3c",
                            destination = "country.name"))

# high_risk_countries$text <- paste0(high_risk_countries$country,
#                                    " ",
#                                    round(high_risk_countries$`50%`, 1),
#                                    " (",
#                                    round(high_risk_countries$`2.5%`, 1),
#                                    ", ",
#                                    round(high_risk_countries$`97.5%`, 1),
#                                    ")")





plot_importer <- world_sf[world_sf$GU_A3 %in% c("ESP", "BRA", "USA", "CHN", "ZAF", "GBR"),]

#plot_importer <- bind_rows(plot_importer)
plot_importer$country <- countrycode(sourcevar = plot_importer$GU_A3,
                                     origin = "iso3c",
                                     destination = "country.name")
plot_importer$country <- factor(plot_importer$country, 
       levels = countrycode(sourcevar = c("ESP", "BRA", "USA", "CHN", "ZAF", "GBR"),
                            origin = "iso3c",
                            destination = "country.name"))
export_risk_maps_panels2 <- export_risk_maps_panels[!is.na(export_risk_maps_panels$country),]

sf_use_s2(FALSE) # 
#export_risk <- st_crop(export_risk, xmin = -180, xmax = 180, ymin = -60, ymax = 80)
export_risk_maps <- ggplot(data = export_risk_maps_panels2, aes(fill = `50%`) ) + # create a ggplot object and 
  # change its fill colour according to median
  geom_sf(size = 0.1) + scale_fill_viridis("Exporter risk") +
  geom_sf(data = plot_importer,
          size = 0.1, fill = "black", col = "black") +
  theme(#legend.title = element_text(size = 5),
    #legend.text = element_text(size = 4),
    #legend.position = c(0.1,0.4),
    #legend.background = element_rect(fill ="white", colour = "black"),
    #legend.key.size = unit(0.2, units = "cm") ,
    legend.position = "bottom",    
    axis.text = element_text(size = 7)) + coord_sf(expand = FALSE) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = c( -120, -60, 0, 60, 120)) +
  scale_y_continuous(limits = c(-60, 80),
                     breaks = c(-60, -40, -20, 0, 20, 40, 60, 80)) +
  facet_wrap(~country, ncol = 2) 
#export_risk_maps
ggsave(export_risk_maps, 
       filename = paste0(minRE_spp, "/export_risk_maps.png"),
       height = 6, width = 6, units= "in",
       dpi = 500)  

export_risk_maps_panels2 %>% split(.$country) %>%
  map(~ ggplot(data = export_risk_maps_panels2, aes(fill = `50%`) ) + # create a ggplot object and 
  # change its fill colour according to median
  geom_sf(size = 0.1) + scale_fill_viridis("Exporter risk") +
  geom_sf(data = plot_importer,
          size = 0.1, fill = "black", col = "black") +
  theme(#legend.title = element_text(size = 5),
    #legend.text = element_text(size = 4),
    #legend.position = c(0.1,0.4),
    #legend.background = element_rect(fill ="white", colour = "black"),
    #legend.key.size = unit(0.2, units = "cm") ,
    legend.position = "bottom",    
    axis.text = element_text(size = 7)) + coord_sf(expand = FALSE) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = c( -120, -60, 0, 60, 120)) +
  scale_y_continuous(limits = c(-60, 80),
                     breaks = c(-60, -40, -20, 0, 20, 40, 60, 80)) +
  facet_wrap(~country, ncol = 1)) %>%
  cowplot::plot_grid(plotlist = .)

export_panels3 <- export_risk_maps_panels2 %>% split(.$country)
map_panels <- function(x){
  map_export_risk <- ggplot(data = x, aes(fill = `50%`) ) + # create a ggplot object and 
    # change its fill colour according to median
    geom_sf(size = 0.1) + scale_fill_viridis("Exporter risk") +
    geom_sf(data = plot_importer[plot_importer$country %in% 
                                   intersect(plot_importer$country, x$country),],
            size = 0.1, fill = "black", col = "black") +
    theme(legend.text = element_text(size = legend_text),
          legend.title = (element_text(size = legend_text)),
          legend.position = position,
          #legend.background = element_rect(fill ="white", colour = "black"),
          legend.key.size = unit(key_size, units = "cm") ,
          #legend.position = "bottom",    
          axis.text = element_text(size = 7)) + coord_sf(expand = FALSE) +
    scale_x_continuous(limits = c(-180, 180),
                       breaks = c( -120, -60, 0, 60, 120)) +
    scale_y_continuous(limits = c(-60, 80),
                       breaks = c(-60, -40, -20, 0, 20, 40, 60, 80)) +
    facet_wrap(~country)
  return(map_export_risk)
}

exporter_maps <- lapply(export_panels3, map_panels)

ggarrange(exporter_maps[[1]],
          exporter_maps[[2]],
          exporter_maps[[3]],
          exporter_maps[[4]],
          exporter_maps[[5]],
          exporter_maps[[6]],
          nrow = 3, ncol = 2)
ggsave(paste0(minRE_spp, "/", "export_risk_maps.png"),
       height = 6,
       width = 8,
       units = "in",
       dpi = 500)


# map the country-level uncertainty in predictions
pred_draws <- pp_average(top_models[[1]],
                         top_models[[2]], 
                         top_models[[3]],
                         top_models[[4]],
                         top_models[[5]], 
                         top_models[[6]], 
                         top_models[[7]],
                         top_models[[8]],
                         top_models[[9]], 
                         top_models[[10]], 
                         top_models[[11]],
                         top_models[[12]],
                         top_models[[13]],
                         weights = "loo", method = "posterior_epred",
                         newdata = newdata_emergence_risk[complete.cases(newdata_emergence_risk),],
                         #dpar = "mu", # excluding zero-inflation/non-detection 
                         ndraws = 1000,
                         re_formula = ~(1|importer_iso3),
                         allow_new_levels = TRUE,
                         sample_new_levels = "gaussian",
                         summary = FALSE,
                         missing = 0) %>% as_draws_df() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  # group_by(importer_iso3) %>% 
  # summarize(`50%` = median(zi),
  #           `2.5%` = quantile(zi, c(0.025)),
  #           `97.5` = quantile(zi, c(0.975))) %>%
  mutate(GU_A3 = newdata_emergence_risk[complete.cases(newdata_emergence_risk),]$importer_iso3)

world_sf <- st_read(dsn = "ne_50m_admin_0_countries", layer = "ne_50m_admin_0_countries",
                    stringsAsFactors = FALSE)
#setdiff(pred_trade_risk$GU_A3, world_sf$GU_A3) # check all 135 countries present on map
risk_map_overall <- left_join(world_sf, pred_draws)
risk_map_overall$risk_factor <- "Composite risk"
risk_map_overall$uncertainty <- risk_map_overall$`97.5%` - risk_map_overall$`2.5%`
#plot uncertainty in predictions highlights regions with relatively well-documented 
# Phytophthora diversity: this is slightly counter-intuitive: reflects substantial species-level variation?
# it is reflecting species-level variation.
# we probably need to plot country-level uncertainty only here?
ggplot(data = risk_map_overall, aes(fill = uncertainty)) + # create a ggplot object and 
  # change its fill colour according to median_age
  geom_sf(size = 0.5) + scale_fill_viridis_c("Uncertainty (95% credible interval)", option = "magma") +
  theme(#legend.title = element_text(size = 5),
    #legend.text = element_text(size = 4),
    #legend.position = c(0.1,0.4),
    #legend.background = element_rect(fill ="white", colour = "black"),
    #legend.key.size = unit(0.2, units = "cm") ,
    legend.position = "bottom",    
    axis.text = element_text(size = 7)) + coord_sf(expand = FALSE)

#plot observed against predicted for UK
pred_draws <- pp_average(top_models[[1]],
                         top_models[[2]], 
                         top_models[[3]],
                         top_models[[4]],
                         top_models[[5]], 
                         top_models[[6]], 
                         top_models[[7]],
                         top_models[[8]],
                         top_models[[9]], 
                         top_models[[10]], 
                         top_models[[11]],
                         top_models[[12]],
                         top_models[[13]],
                         weights = "loo", method = "posterior_epred",
                         newdata = top_models[[1]]$data,
                         dpar = "mu", # excluding zero-inflation/non-detection 
                         ndraws = 1000,
                         #re_formula = NA,
                         #allow_new_levels = TRUE,
                         #sample_new_levels = "gaussian",
                         summary = FALSE,
                         missing = 0) %>% as_draws_df() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  # group_by(importer_iso3) %>% 
  # summarize(`50%` = median(zi),
  #           `2.5%` = quantile(zi, c(0.025)),
  #           `97.5` = quantile(zi, c(0.975))) %>%
  mutate(importer_iso3 = top_models[[1]]$data$importer_iso3,
         phyto = top_models[[1]]$data$phyto,
         arrival = top_models[[1]]$data$arrival)

# look at the non-arrivals
# which of these have the greatest predicted probability of arrival?
#pred_draws0 <- pred_draws[pred_draws$arrival == 0,]
pred_draws$phyto <- factor(pred_draws$phyto, 
                           levels = names(sort(tapply(pred_draws$`50%`, 
                                                      pred_draws$phyto, FUN = median), 
                                               decreasing = TRUE)))
num_detections <- data.frame(phyto = names(tapply(pred_draws$arrival, 
                                                  pred_draws$phyto, FUN = sum)),
                             num_detections = tapply(pred_draws$arrival, 
                                                     pred_draws$phyto, FUN = sum),
                             stringsAsFactors = FALSE)

pred_detections <- data.frame(phyto = names(tapply(pred_draws$`50%`, 
                                                  pred_draws$phyto, FUN = sum)),
                             pred_detections = tapply(pred_draws$`50%`, 
                                                     pred_draws$phyto, FUN = sum),
                             stringsAsFactors = FALSE)

ob_pred_det <- left_join(num_detections,
                         pred_detections)
ob_pred_det$diff <- ob_pred_det$pred_detections - ob_pred_det$num_detections
top10underpredicted <- ob_pred_det[order(ob_pred_det$diff, decreasing = TRUE),][1:15,]$phyto
ob_pred_det$underpred <- ob_pred_det$phyto %in% top10underpredicted

ggplot() + geom_boxplot(data = pred_draws, aes(x = phyto, y = `50%`))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_point(data = num_detections, aes(x = phyto, y = num_detections/20), col = "red")



ggplot(data = ob_pred_det, aes(y = pred_detections, 
                               x = num_detections, 
                               fill = diff,
                               col = diff))  +
  scale_fill_gradient2("Predicted - Observed", mid = "white" , low = "steelblue", high = "darkred" ) +
  scale_colour_gradient2("Predicted - Observed", mid = "white" , low = "steelblue", high = "darkred" ) +
  geom_point(shape = 21, size = 2) +
  geom_point(shape = 1, aes(y = pred_detections, 
                             x = num_detections), col = "grey") +
  geom_text(data = ob_pred_det[ob_pred_det$underpred,], aes(y = pred_detections, 
                 x = num_detections, 
                 label = gsub("Phytophthora", "P.", phyto)),
             hjust = 1, vjust = 0, check_overlap = TRUE, fontface = "italic" , angle = 0) +
  geom_abline(slope = 1, linetype = 3) + 
  theme(panel.grid = element_blank()) + scale_x_continuous(limits = c(-6, 20)) +
  ylab("Predicted detections") + xlab("Observed detections")
ggsave("Predicted_observed.png",
       height = 6, width = 6,
       units = "in", dpi = 500)
