rm(list=ls())
library(dplyr)
library(ggplot2)
library(brmstools)
library(brms)
library(ggpubr)
library(tidybayes)
library(boot)
library(tidyverse)

minRE_spp <- "arrivals_minREspp0"
com_type <- "trade_hort"

dredge_arrival <- readRDS(paste0(minRE_spp,"/dredge_arrival.rds"))
arrivals_unscaled <- readRDS(paste0(minRE_spp, "/arrivals_unscaled.rds"))

rank_models <- list()
models <- list.files(paste0(minRE_spp,"/", com_type), 
                     pattern = "fit", full = TRUE)
# check for missing models
# models <- setdiff(list.files(paste0("arrivals_minREspp", minRE_spp, "/", com_type), 
#                              pattern = "fit", full = TRUE),
#                   rank_models_all$model)


for(i in 1:length(models)){
  print(models[i])
  fit <- readRDS(models[i])
  rank_models[[i]] <- data.frame(model = models[i], 
                                 looic = fit$criteria$loo$estimates["looic",
                                                                    "Estimate"],
                                 #rhat_gt1.05 = length(which(rhat(fit)[grep("b_", names(rhat(fit)))] > 1.05)),
                                 stringsAsFactors = FALSE)
}
rank_models_all <- bind_rows(rank_models)
rank_models_all$diffic <- min(rank_models_all$looic) - rank_models_all$looic
# missing <- setdiff(1:nrow(dredge_arrival), 
#         as.integer(gsub("P:/NEC05617_PHYTOTHREATS/trade_hort_minREspp2/|\\.rds|fit", 
#                         "", rank_models_all$model)))
# saveRDS(missing, file = "missing.rds")
rank_models_all<- rank_models_all[order(rank_models_all$diffic, decreasing = TRUE),]
rownames(rank_models_all) <- 1:nrow(rank_models_all)
saveRDS(rank_models_all, file = paste0(minRE_spp, "/", "rank_models_all.rds"))
rank_models_all <- readRDS(paste0(minRE_spp, "/", "rank_models_all.rds"))





#rank_models[order(rank_models$diffic, decreasing = TRUE),]



best_models <- rank_models_all[rank_models_all$diffic >=-2,]
best_models$model<- gsub("_rslurm_trade_hort_2023\\/", paste0(minRE_spp, "/", com_type, "/"), best_models$model)

#best_models <- rank_models[1:5,]
# plot IC against rank - is there a drop off somewhere?
plot(rownames(rank_models_all), rank_models_all$diffic, pch = 16,
     xlab = "model rank", ylab = "difference LOO-IC", las = 1, cex = 0.6)
abline(h = -5, lty = 4)
# 16 best models (diff LOO-IC >= -5, within 5IC untis of best model)

png(paste0(minRE_spp, "/model_rank_selection.png"),
    heigh = 6, width = 6,
    units = "in", res = 500)
plot(rownames(rank_models_all), rank_models_all$diffic, pch = 16,
     xlab = "model rank", ylab = "difference LOO-IC", las = 1, cex = 0.6)
abline(h = -2, lty = 4)
dev.off()
# 16 best models (diff LOO-IC >= -2, within 2IC untis of best model)

for(i in 1:nrow(best_models)){
  print(best_models$model[i])
  assign(gsub("arrivals_minREspp0\\/trade_hort\\/|\\.rds", "", best_models$model[i]),
         readRDS(file = best_models$model[i]))
}


########### get variance components manually ###############
# performance package not working with zero inflated models
# null model with no zero inflation predictors
# prior <-  prior(cauchy(0, 1), "sd")
# Phytophthora_covariance <- readRDS("PhytophthoraA.rds")
# fit_null <-
#  brm(arrival | trials(1) ~ (1 | importer_iso3) + (1 | gr(phylo, cov = Phytophthora_covariance)) + (1 | phyto),
#      data = arrivals, family = zero_inflated_binomial() ,
#      data2 = list(Phytophthora_covariance = Phytophthora_covariance),
#      prior = prior,
#      control = list(adapt_delta = 0.999,
#                     max_treedepth = 15),
#     iter = 6000,
#      chains = 3
#  )
#fit_null <- add_criterion(fit_null, criterion = c("loo", "loo_R2", "bayes_R2"))
#saveRDS(fit_null, file = "fit_null.rds")

#saveRDS(fit, file = paste0(index, ".rds"))

#model_null <- readRDS(paste0(minRE_spp, "/fit_null.rds"))
model_null <- readRDS(paste0(minRE_spp,"/", com_type, "/fit0001.rds")) # including zero-inflation component
full_model <- readRDS(paste0(minRE_spp,"/", com_type, "/fit0178.rds"))

# check for differences in conditional and marginal R squared using a null model
# where
# a) covariates are used to predict zero-inflation
# b) no covariates are used to predict zero-inflation

partition_variance <- function(model_f, model_null){
  model_f <- readRDS(model_f)
  get_draws_null <- as_draws_df(model_null)
  get_draws_fixed <- as_draws_df(model_f)

  # fixed effects
  # sigma2_fixed <- apply(posterior_linpred(model_f,
  #                                       #scale = "linear", # linear predictor scale
  #                                       robust = TRUE,
  #                                       re_formula = NA), 1, var)
  sigma2_fixed <- linpred_draws(model_f,
                           newdata = model_f$data) %>% group_by(.draw) %>% 
    summarise(var = var(.linpred)) %>% pull(var) 
  
  # total structured variance of the null model (random effects only)
  Vt <- as.numeric(rowSums(get_draws_fixed[,c(grep("sd", names(get_draws_fixed)))]^2))
  
  
  
  # Intercept of the null model: expected grand mean on the linear predictor scale
  # accounting for phylogenetically structured random intercept
  b0 <- get_draws_null$b_Intercept
  n = 1 # trials = 1: detection(1), non-detection (0)
  pmean<-plogis(as.numeric(get_draws_null$b_Intercept)-0.5*Vt*tanh(as.numeric(get_draws_null$b_Intercept)* (1+2*exp(-0.5*Vt))/6))
  # convert to the response scale using Eqn. 6.8 (Nakagawa et al. 2017)
  VarOL<-1/(pmean*(1-pmean)) # residual variance (Nakagawa et al. 2017)
  #VarDS <- pi^2/3

  # the variances of the random effects (including the random slopes)
  sigma2_RE <- get_draws_fixed[,c(grep("sd", names(get_draws_fixed)))]^2
  
  # random intercept variances
  sigma2_phylo <- as.numeric(sigma2_RE[,c(grep("sd_phylo__Intercept", names(sigma2_RE)))])# phylogenetically structured variance (of candidate model - model_f)
  sigma2_phyto <- as.numeric(sigma2_RE[,c(grep("sd_phyto__Intercept", names(sigma2_RE)))]) # species-level variance
  sigma2_importer_iso3 <- as.numeric(sigma2_RE[,c(grep("sd_importer_iso3__Intercept", names(sigma2_RE)))]) # country-level variance
  # random slope variances
  sigma2_sd_phyto__climate_match <- as.numeric(sigma2_RE[,c(grep("sd_phyto__climate_match", names(sigma2_RE)))]) #
  sigma2_sd_phyto__EPPO_reporting_service <- as.numeric(sigma2_RE[,c(grep("sd_phyto__EPPO_reporting_service", names(sigma2_RE)))]) #
  sigma2_sd_phyto__trade_hort <- as.numeric(sigma2_RE[,c(grep("sd_phyto__trade_hort", names(sigma2_RE)))]) #
  
  parests <- posterior_summary(model_f, variable = grep("b_", variables(model_f), value = TRUE))
  
  results <-
    array(dim = c(1,24,4),
          dimnames = list(NULL, 
                          c(grep("b_", variables(full_model), value = TRUE),
                            "R2_marginal",
                            "R2_conditional",
                            "ICC_phylo",
                            "ICC_phyto",
                            "ICC_importer_iso3",
                            "ICC_phyto__climate_match",
                            "ICC_phyto__EPPO_reporting_service",
                            "ICC_phyto__trade_hort",
                            "ICC_intercept", # all intercept variance
                            "ICC_slope"), # all slope variance),
                          c("est_same_draw", "lower", "est", "upper" )))
  
  R2glmmM<-sigma2_fixed/(sigma2_fixed + rowSums(sigma2_RE) + VarOL)
  R2glmmC<-(sigma2_fixed + rowSums(sigma2_RE)) / 
  (sigma2_fixed + rowSums(sigma2_RE) + VarOL)
  
  # get the draw representative of the median: the varinace components will only add up to the conditional R2 within the same sample/draw
  draw_aprox_median <- which.min(abs(R2glmmM-quantile(R2glmmM, probs = 0.5)))
  # using the same draw to extract the median estimates of the marginal, 
  # conditional and ICCs to ensure they exactly add up to the total variance (there will be some variation among draws otherwise)
  results[,"R2_marginal",]       <-  c(((sigma2_fixed)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL))[draw_aprox_median],
                                       ((sigma2_fixed)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL)) %>% quantile(probs = c(0.025, 0.5, 0.975)))
  results[,"R2_conditional",]       <-  c(((sigma2_fixed + rowSums(sigma2_RE))/(sigma2_fixed + rowSums(sigma2_RE) + VarOL))[draw_aprox_median],
                                       ((sigma2_fixed + rowSums(sigma2_RE))/(sigma2_fixed + rowSums(sigma2_RE) + VarOL)) %>% quantile(probs = c(0.025, 0.5, 0.975)))
  results[,"ICC_phylo",]         <-  c(((sigma2_phylo)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL))[draw_aprox_median],
                                       ((sigma2_phylo)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL)) %>% quantile(probs = c(0.025, 0.5, 0.975)))
  results[,"ICC_phyto",]         <-  c(((sigma2_phyto)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL))[draw_aprox_median],
                                       ((sigma2_phyto)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL)) %>% quantile(probs = c(0.025, 0.5, 0.975)))
  results[,"ICC_importer_iso3",]         <-  c(((sigma2_importer_iso3)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL))[draw_aprox_median],
                                       ((sigma2_importer_iso3)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL)) %>% quantile(probs = c(0.025, 0.5, 0.975)))
  results[,"ICC_phyto__climate_match",]         <-  c(((sigma2_sd_phyto__climate_match)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL))[draw_aprox_median],
                                       ((sigma2_sd_phyto__climate_match)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL)) %>% quantile(probs = c(0.025, 0.5, 0.975)))
  results[,"ICC_phyto__EPPO_reporting_service",]         <-  c(((sigma2_sd_phyto__EPPO_reporting_service)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL))[draw_aprox_median],
                                       ((sigma2_sd_phyto__EPPO_reporting_service)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL)) %>% quantile(probs = c(0.025, 0.5, 0.975))) 
  results[,"ICC_phyto__trade_hort",]         <-  c(((sigma2_sd_phyto__trade_hort)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL))[draw_aprox_median],
                                       ((sigma2_sd_phyto__trade_hort)/(sigma2_fixed + rowSums(sigma2_RE) + VarOL)) %>% quantile(probs = c(0.025, 0.5, 0.975)))                                      
  results[,"ICC_intercept",]           <-  c(((rowSums(sigma2_RE[,grep("Intercept", names(sigma2_RE))]))/(sigma2_fixed + rowSums(sigma2_RE) + VarOL))[draw_aprox_median],
                                       ((rowSums(sigma2_RE[,grep("Intercept", names(sigma2_RE))]))/(sigma2_fixed + rowSums(sigma2_RE) + VarOL)) %>% quantile(probs = c(0.025, 0.5, 0.975)))
  results[,"ICC_slope",]           <-  c(((rowSums(sigma2_RE[,grep("Intercept", names(sigma2_RE), invert = TRUE)]))/(sigma2_fixed + rowSums(sigma2_RE) + VarOL))[draw_aprox_median],
                                       ((rowSums(sigma2_RE[,grep("Intercept", names(sigma2_RE), invert = TRUE)]))/(sigma2_fixed + rowSums(sigma2_RE) + VarOL)) %>% quantile(probs = c(0.025, 0.5, 0.975)))
  
  results[,grep("b_", variables(model_f), value = TRUE),c("est", "lower", "upper")] <- parests[,c(1,3,4)]
  
  
  #print(model_f)
  rm(list = c("model_f"))
  
  return(results)
}


# extract variance components from the top subset of models (13 models) 
model_summary <- lapply(paste0(minRE_spp, "/", com_type, "/", best_models$model), 
                        partition_variance, model_null)
# names(model_summary) <- gsub("_rslurm_trade_hort_2023/", "", best_models$model)
# saveRDS(model_summary, file = paste0(minRE_spp, "/model_summary.rds"))

model_info <- abind(model_summary, along = 1)
#model_info <- bind_rows(model_summary, .id = "")
saveRDS(model_info, file = paste0(minRE_spp, "/model_info.rds"))
model_info <-readRDS(paste0(minRE_spp, "/model_info.rds"))



#### report top-model subset variance components ####
round(model_info[order(model_info[,"R2_conditional","est"]),"R2_conditional",
           c("lower", "est", "upper")]*100, 1)

round(model_info[order(model_info[,"R2_marginal","est"]),"R2_marginal",
                 c("lower", "est", "upper")]*100, 1)

# random intercept variance
round(model_info[order(model_info[,"ICC_phyto","est"]),"ICC_phyto",
                 c("lower", "est", "upper")]*100, 1)

round(model_info[order(model_info[,"ICC_importer_iso3","est"]),"ICC_importer_iso3",
                 c("lower", "est", "upper")]*100, 1)

round(model_info[order(model_info[,"ICC_phylo","est"]),"ICC_phylo",
                 c("lower", "est", "upper")]*100, 1)

# random slope variance
round(model_info[order(model_info[,"ICC_phyto__trade_hort","est"]),"ICC_phyto__trade_hort",
                 c("lower", "est", "upper")]*100, 1)

round(model_info[order(model_info[,"ICC_phyto__climate_match","est"]),"ICC_phyto__climate_match",
                 c("lower", "est", "upper")]*100, 1)

round(model_info[order(model_info[,"ICC_phyto__EPPO_reporting_service","est"]),"ICC_phyto__EPPO_reporting_service",
                 c("lower", "est", "upper")]*100, 1)

# which of the best-performing models contain interaction terms
# get the rank of the best-performing model with one or more interaction terms
int_terms <- model_info[,grep(":", colnames(model_info)),"est"]
anynotNA <- function(x){
  any(!is.na(x))
}
allNA <- function(x){
  all(is.na(x))
}

sum(apply(int_terms, 1, anynotNA)) # 9/13 top-performing models cont

################### report effect sizes using model averaging #####################
# the modelled data (identical for all models)
arrivals <- readRDS(paste0(minRE_spp, "/",
                             com_type, "/",
                             gsub("_rslurm_trade_hort_2023\\/", "",
                                  rank_models_all$model[1])))$data


top_subset <- best_models$model
top_models <- lapply(top_subset, readRDS)

# to reverse the scaling, get the mean(centre) and scaling factor (2 SDs)  of the original data
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

date_scale <- 2*sd(log(1+(arrivals_unscaled$date)))
date_center <- mean(log(1 + (arrivals_unscaled$date)))

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#rank_models_all <- readRDS(paste0(minRE_spp, "/rank_models_all.rds"))

# overall grand mean: global probability fo a new detection
# spread_draws(best_model, b_Intercept) %>% 
#   mutate(b_Intercept = inv.logit(b_Intercept)) %>%
#   summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) 

averaged_model_all <- posterior_average(top_models[[1]],
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
                                    weights = "loo") %>% as_draws()

#global detection probability
averaged_model_all %>% select(b_Intercept) %>%
  mutate(logit2prob(b_Intercept)) %>% # convert to probability
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))



# overall zero-inflation: probability of a non-detection 
# spread_draws(best_model, b_zi_Intercept) %>% 
#   mutate(b_zi_Intercept = inv.logit(b_zi_Intercept)) %>%
#   summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))



averaged_model_all %>% select(b_zi_Intercept) %>%
  mutate(logit2prob(b_zi_Intercept)) %>% # convert to probability
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))


# a 10% increase in EPPO reporting reduces the odds of a non-detection to
averaged_model_all %>% select(b_zi_EPPO_reporting_service) %>%
  mutate(b_zi_EPPO_reporting_service = 1.1^b_zi_EPPO_reporting_service) %>%
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))  

quantile(arrivals_unscaled$EPPO_reporting_service, c(0.25,0.5, 0.75))
# in absolute terms an increase in historical EPPO Reporting Service activity 
# from 26 (25th percentile) to 161 (75th percentile) reports reduced the 
# predicted probability of a non-detection from... to...


newdata_zi <- data.frame(importer_iso3 = "any",
                         trade_hort = mean(arrivals$trade_hort), 
                         phyto = "unknown_species",
                         phylo = "unknown_species",
                         CH_OO = mean(arrivals$CH_OO),
                         temp_range = mean(arrivals$temp_range),
                         climate_match = mean(arrivals$climate_match),
                         EPPO_reporting_service = (log(quantile(arrivals_unscaled$EPPO_reporting_service, c(0.25,0.5, 0.75)))-surv_center)/surv_scale,
                         date = mean(arrivals$date),
                         stringsAsFactors = FALSE)

pp_average(top_models[[1]],
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
           dpar = "zi",
           allow_new_levels = TRUE,
           sample_new_levels = "gaussian",
           summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))

# epred_draws(best_model,
#             newdata = newdata_zi,
#             dpar = "zi",
#             allow_new_levels = TRUE,
#             sample_new_levels = "gaussian") %>% 
#   group_by(EPPO_reporting_service) %>%
#   summarise(zi = quantile(zi, c(0.025, 0.5, 0.975))) %>% ungroup()%>%
#   mutate(quantile = rep(c("2.5%", "50%", "97.5%"),3)) %>%
#   pivot_wider(names_from = quantile, values_from = zi ) %>%
#   mutate(EPPO_reporting_service_raw = quantile(arrivals_unscaled$EPPO_reporting_service, c(0.25,0.5, 0.75)))
# 
# 

# a 10% increase in years since pathogen was described increased the odds of a 
# non-detection to...
averaged_model_all %>% select(b_zi_date) %>%
  mutate(b_zi_date = 1.1^b_zi_date) %>%
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))  

quantile(arrivals_unscaled$date, c(0.25,0.5, 0.75))
# in absolute terms an increase in years since pathogen was described 
# from 12 (25th percentile) to 46 (75th percentile) reports reduced the 
# predicted probability of a non-detection from... to...
newdata_zi <- data.frame(importer_iso3 = "any",
                         trade_hort = mean(arrivals$trade_hort), 
                         phyto = "unknown_species",
                         phylo = "unknown_species",
                         CH_OO = mean(arrivals$CH_OO),
                         temp_range = mean(arrivals$temp_range),
                         climate_match = mean(arrivals$climate_match),
                         EPPO_reporting_service = mean(arrivals$EPPO_reporting_service),
                         date = (log(1+quantile(unique(arrivals_unscaled$date), c(0.25,0.5, 0.75)))-date_center)/date_scale,
                         stringsAsFactors = FALSE)

pp_average(top_models[[1]],
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
           dpar = "zi",
           allow_new_levels = TRUE,
           sample_new_levels = "gaussian",
           summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))
# epred_draws(best_model,
#             newdata = newdata_zi,
#             dpar = "zi",
#             allow_new_levels = TRUE,
#             sample_new_levels = "gaussian") %>% 
#   group_by(date) %>%
#   summarise(zi = quantile(zi, c(0.025, 0.5, 0.975))) %>% ungroup()%>%
#   mutate(quantile = rep(c("2.5%", "50%", "97.5%"),3)) %>%
#   pivot_wider(names_from = quantile, values_from = zi ) %>%
#   mutate(date_raw = quantile(arrivals_unscaled$date, c(0.25,0.5, 0.75)))


options(pillar.sigfig = 4)
# a 10% increase in imports increases the odds of a new detection by 
#round(1.1^posterior_summary(best_model)["b_trade_hort",]*100, 1)
# 121.1%
averaged_model_all %>% select(b_trade_hort) %>%
  mutate(b_trade_hort = (1.1^b_trade_hort)*100) %>%
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) 
# spread_draws(best_model, b_trade_hort) %>% 
#   mutate(b_trade_hort = (1.1^b_trade_hort)*100) %>%
#   summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))  

quantile(arrivals_unscaled$trade_hort, c(0.25,0.5, 0.75))
#An increase in horticultural trade flows from 446 (25th percentile) to 114,397
# (75th percentile) metric tons increased the probability of a new detection 
# from... to...
newdata <- data.frame(importer_iso3 = "any",
                         trade_hort = (log(1+quantile(arrivals_unscaled$trade_hort, c(0.25,0.5, 0.75)))-trade_center)/trade_scale, 
                         phyto = "unknown_species",
                         phylo = "unknown_species",
                         CH_OO = mean(arrivals$CH_OO),
                         temp_range = mean(arrivals$temp_range),
                         climate_match = mean(arrivals$climate_match),
                         EPPO_reporting_service = mean(arrivals$EPPO_reporting_service),
                         date = mean(arrivals$date),
                         stringsAsFactors = FALSE)

pp_average(top_models[[1]],
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
           newdata = newdata,
           allow_new_levels = TRUE,
           sample_new_levels = "gaussian",
           summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))
# epred_draws(best_model,
#             newdata = newdata,
#             allow_new_levels = TRUE,
#             sample_new_levels = "gaussian") %>% 
#   group_by(trade_hort) %>%
#   summarise(.epred = quantile(.epred, c(0.025, 0.5, 0.975))) %>% ungroup()%>%
#   mutate(quantile = rep(c("2.5%", "50%", "97.5%"),4)) %>%
#   pivot_wider(names_from = quantile, values_from = .epred ) %>%
#   mutate(trade_raw = quantile(arrivals_unscaled$trade_hort, c(0.25,0.5, 0.75, 0.90)))


# climate matching effects

#round(1.1^posterior_summary(best_model)["b_climate_match",]*100, 1)
# 109.8%
# spread_draws(best_model, b_climate_match) %>% 
#   mutate(b_climate_match = (1.1^b_climate_match)*100) %>%
#   summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))  
averaged_model_all %>% select(b_climate_match) %>%
  mutate(b_climate_match = (1.1^b_climate_match)*100) %>%
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) 

-sqrt(quantile(arrivals_unscaled$climate_distance, c(0.25,0.5, 0.75)))
newdata <- data.frame(importer_iso3 = "any",
                      trade_hort = mean(arrivals$trade_hort), 
                      phyto = "unknown_species",
                      phylo = "unknown_species",
                      CH_OO = mean(arrivals$CH_OO),
                      temp_range = mean(arrivals$temp_range),
                      climate_match = ((-sqrt(quantile(arrivals_unscaled$climate_distance, c(0.25,0.5, 0.75))))-climate_center)/climate_scale,
                      EPPO_reporting_service = mean(arrivals$EPPO_reporting_service),
                      date = mean(arrivals$date),
                      stringsAsFactors = FALSE)
# epred_draws(best_model,
#             newdata = newdata,
#             allow_new_levels = TRUE,
#             sample_new_levels = "gaussian") %>% 
#   group_by(climate_match) %>%
#   summarise(.epred = quantile(.epred, c(0.025, 0.5, 0.975))) %>% ungroup()%>%
#   mutate(quantile = rep(c("2.5%", "50%", "97.5%"),4)) %>%
#   pivot_wider(names_from = quantile, values_from = .epred ) %>%
#   mutate(climate_raw = quantile(arrivals_unscaled$climate_distance, c(0.25,0.5, 0.75, 0.90)))
pp_average(top_models[[1]],
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
           newdata = newdata,
           allow_new_levels = TRUE,
           sample_new_levels = "gaussian",
           summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))




# surveillance effects
#round(1.1^posterior_summary(best_model)["b_EPPO_reporting_service",]*100, 1)

# spread_draws(best_model, b_EPPO_reporting_service) %>% 
#   mutate(b_EPPO_reporting_service = (1.1^b_EPPO_reporting_service)*100) %>%
#   summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))  
averaged_model_all %>% select(b_EPPO_reporting_service) %>%
  mutate(b_EPPO_reporting_service = (1.1^b_EPPO_reporting_service)*100) %>%
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) 


newdata <- data.frame(importer_iso3 = "any",
                      trade_hort = mean(arrivals$trade_hort), 
                      phyto = "unknown_species",
                      phylo = "unknown_species",
                      CH_OO = mean(arrivals$CH_OO),
                      temp_range = mean(arrivals$temp_range),
                      climate_match = mean(arrivals$climate_match),
                      EPPO_reporting_service = (log(quantile(arrivals_unscaled$EPPO_reporting_service, c(0.25,0.5, 0.75)))-surv_center)/surv_scale,
                      date = mean(arrivals$date),
                      stringsAsFactors = FALSE)
pp_average(top_models[[1]],
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
           newdata = newdata,
           allow_new_levels = TRUE,
           sample_new_levels = "gaussian",
           summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))
# epred_draws(best_model,
#             newdata = newdata,
#             allow_new_levels = TRUE,
#             sample_new_levels = "gaussian") %>% 
#   group_by(EPPO_reporting_service) %>%
#   summarise(.epred = quantile(.epred, c(0.025, 0.5, 0.975))) %>% ungroup()%>%
#   mutate(quantile = rep(c("2.5%", "50%", "97.5%"),4)) %>%
#   pivot_wider(names_from = quantile, values_from = .epred ) %>%
#   mutate(EPPO_reporting_service_raw = quantile(arrivals_unscaled$EPPO_reporting_service, c(0.25,0.5, 0.75, 0.90)))

# thermal tolerance range effects
#round(1.1^posterior_summary(best_model)["b_EPPO_reporting_service",]*100, 1)

# for consistency witht he other predictors, also report effect of incfrease in effect size of 10%
# this predictor is not log transformed so we can't use the *1.1 for odds ratio
# median is 24.5, median plus 10% is 27: extract predictions using these figures
# spread_draws(best_model, b_temp_range) %>% 
#   mutate(b_temp_range = (1.1^b_temp_range)*100) %>%
#   summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))  
# temp range is not present in all of the top models
which(!is.na(model_info[,"b_temp_range","est"]))
averaged_models_incl_temp_range <- posterior_average(top_models[[1]],
                                                     #top_models[[2]], 
                                                     #top_models[[3]],
                                                     top_models[[4]],
                                                     top_models[[5]], 
                                                     top_models[[6]], 
                                                     top_models[[7]],
                                                     #top_models[[8]],
                                                     top_models[[9]], 
                                                     top_models[[10]], 
                                                     #top_models[[11]],
                                                     #top_models[[12]],
                                                     #top_models[[13]],
                                                     weights = "loo") %>% as_draws()


averaged_models_incl_temp_range %>% select(b_temp_range ) %>%
  mutate(b_temp_range = (1.1^b_temp_range)*100) %>%
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) 

quantile(arrivals_unscaled$temp_range, c(0.25, 0.5, 0.75))
newdata <- data.frame(importer_iso3 = "any",
                      trade_hort = mean(arrivals$trade_hort), 
                      phyto = "unknown_species",
                      phylo = "unknown_species",
                      CH_OO = mean(arrivals$CH_OO),
                      temp_range = (quantile(arrivals_unscaled$temp_range, c(0.25, 0.5, 0.75))-temp_range_center)/temp_range_scale,
                      climate_match = mean(arrivals$climate_match),
                      EPPO_reporting_service = mean(arrivals$EPPO_reporting_service),
                      date = mean(arrivals$date),
                      stringsAsFactors = FALSE)
pp_average(top_models[[1]],
           top_models[[2]], 
           #top_models[[3]],
           top_models[[4]],
           top_models[[5]], 
           top_models[[6]], 
           top_models[[7]],
           #top_models[[8]],
           top_models[[9]], 
           top_models[[10]], 
           #top_models[[11]],
           #top_models[[12]],
           #top_models[[13]],
           weights = "loo", method = "posterior_epred",
           newdata = newdata,
           allow_new_levels = TRUE,
           sample_new_levels = "gaussian",
           summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))
# epred_draws(best_model,
#             newdata = newdata,
#             allow_new_levels = TRUE,
#             sample_new_levels = "gaussian") %>% 
#   group_by(temp_range) %>%
#   summarise(.epred = quantile(.epred, c(0.025, 0.5, 0.975))) %>% ungroup()%>%
#   mutate(quantile = rep(c("2.5%", "50%", "97.5%"),4)) %>%
#   pivot_wider(names_from = quantile, values_from = .epred ) %>%
#   mutate(temp_range_raw = c(quantile(arrivals_unscaled$temp_range, c(0.25, 0.5, 0.75)), 27))
# # predict odds ratio

              

  
# highest and lowest sensitivity to trade
r_draws <- averaged_model_all %>% 
  spread_draws(b_trade_hort, r_phyto[phyto,risk_factor]) %>% 
  filter(risk_factor == "trade_hort" | is.na(risk_factor)) %>%
  mutate(.value = b_trade_hort + r_phyto) %>%
  select(risk_factor, .chain, .iteration, .draw, phyto, .value) 

mean_responses <- r_draws %>% group_by(phyto) %>% summarise(median = quantile(.value, 0.5)) 
low_sens <- as.character(mean_responses[which.min(mean_responses$median),"phyto"])
high_sens <- as.character(mean_responses[which.max(mean_responses$median),"phyto"])


# weakest/most negative
low_sens
quantile(1.1^r_draws[r_draws$phyto %in% low_sens,]$.value, c(0.025, 0.5, 0.975))*100
#strongest/most positive
high_sens
quantile(1.1^r_draws[r_draws$phyto %in% high_sens,]$.value, c(0.025, 0.5, 0.975))*100

newdata <- expand.grid((log(1+quantile(arrivals_unscaled$trade_hort, c(0.25,0.5, 0.75)))-trade_center)/trade_scale,
                       gsub("\\.", " ", c(high_sens, low_sens))) %>% 
  dplyr::rename(trade_hort = 1,
                phyto = 2) %>%
  mutate(importer_iso3 = "any",
         phylo = phyto,
         CH_OO = mean(arrivals$CH_OO),
         temp_range = mean(arrivals$temp_range),
         climate_match = mean(arrivals$climate_match),
         EPPO_reporting_service = mean(arrivals$EPPO_reporting_service),
         date = mean(arrivals$date))
pp_average(top_models[[1]],
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
           newdata = newdata,
           allow_new_levels = TRUE,
           sample_new_levels = "gaussian",
           summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975)))  
# epred_draws(best_model,
#             newdata = newdata,
#             allow_new_levels = TRUE,
#             sample_new_levels = "gaussian") %>% 
#   group_by(trade_hort, phyto) %>%
#   summarise(.epred = quantile(.epred, c(0.025, 0.5, 0.975))) %>% ungroup()%>%
#   mutate(quantile = rep(c("2.5%", "50%", "97.5%"),8)) %>%
#   pivot_wider(names_from = quantile, values_from = .epred ) %>%
#   mutate(trade_hort_raw = exp(trade_hort * trade_scale + trade_center) - 1)

# highest and lowest sensitivity to surveillance
r_draws <- averaged_model_all %>% 
  spread_draws(b_EPPO_reporting_service, r_phyto[phyto,risk_factor]) %>% 
  filter(risk_factor == "EPPO_reporting_service" | is.na(risk_factor)) %>%
  mutate(.value = b_EPPO_reporting_service + r_phyto) %>%
  select(risk_factor, .chain, .iteration, .draw, phyto, .value) 

mean_responses <- r_draws %>% group_by(phyto) %>% summarise(median = quantile(.value, 0.5)) 
low_sens <- as.character(mean_responses[which.min(mean_responses$median),"phyto"])
high_sens <- as.character(mean_responses[which.max(mean_responses$median),"phyto"])


low_sens
quantile(1.1^r_draws[r_draws$phyto %in% low_sens,]$.value, c(0.025, 0.5, 0.975))*100
high_sens
quantile(1.1^r_draws[r_draws$phyto %in% high_sens,]$.value, c(0.025, 0.5, 0.975))*100

newdata <- expand.grid((log(quantile(arrivals_unscaled$EPPO_reporting_service, c(0.25,0.5, 0.75)))-surv_center)/surv_scale,
                       gsub("\\.", " ", c(high_sens, low_sens))) %>% 
  dplyr::rename(EPPO_reporting_service = 1,
                phyto = 2) %>%
  mutate(importer_iso3 = "any",
         phylo = phyto,
         CH_OO = mean(arrivals$CH_OO),
         temp_range = mean(arrivals$temp_range),
         climate_match = mean(arrivals$climate_match),
         trade_hort = mean(arrivals$trade_hort),
         date = mean(arrivals$date))
pp_average(top_models[[1]],
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
           newdata = newdata,
           allow_new_levels = TRUE,
           sample_new_levels = "gaussian",
           summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  mutate(EPPO_reporting_service_raw = exp(newdata$EPPO_reporting_service * surv_scale + surv_center) - 1)
# epred_draws(best_model,
#             newdata = newdata,
#             allow_new_levels = TRUE,
#             sample_new_levels = "gaussian") %>% 
#   group_by(EPPO_reporting_service, phyto) %>%
#   summarise(.epred = quantile(.epred, c(0.025, 0.5, 0.975))) %>% ungroup()%>%
#   mutate(quantile = rep(c("2.5%", "50%", "97.5%"),8)) %>%
#   pivot_wider(names_from = quantile, values_from = .epred ) %>%
#   mutate(EPPO_reporting_service_raw = exp(EPPO_reporting_service * surv_scale + surv_center) - 1)

# highest and lowest sensitivity to climate matching
r_draws <- averaged_model_all %>% 
  spread_draws(b_climate_match, r_phyto[phyto,risk_factor]) %>% 
  filter(risk_factor == "climate_match" | is.na(risk_factor)) %>%
  mutate(.value = b_climate_match + r_phyto) %>%
  select(risk_factor, .chain, .iteration, .draw, phyto, .value) 

mean_responses <- r_draws %>% group_by(phyto) %>% summarise(median = quantile(.value, 0.5),
                                                            lower  = quantile(.value, 0.025),
                                                            upper = quantile (.value, 0.975)) %>%
  mutate(uncertainty = upper-lower)

mean_responses[order(mean_responses$median, decreasing = TRUE),]
low_sens <- as.character(mean_responses[which.min(mean_responses$median),"phyto"])
high_sens <- as.character(mean_responses[which.max(mean_responses$median),"phyto"])

# weakest and strongest responses to climate
low_sens
quantile(1.1^r_draws[r_draws$phyto %in% low_sens,]$.value, c(0.025, 0.5, 0.975))*100
high_sens
quantile(1.1^r_draws[r_draws$phyto %in% high_sens,]$.value, c(0.025, 0.5, 0.975))*100

newdata <- expand.grid(((-sqrt(quantile(arrivals_unscaled$climate_distance, c(0.25,0.5, 0.75))))-climate_center)/climate_scale,
                       gsub("\\.", " ", c(high_sens, low_sens))) %>% 
  dplyr::rename(climate_match = 1,
                phyto = 2) %>%
  mutate(importer_iso3 = "any",
         phylo = phyto,
         CH_OO = mean(arrivals$CH_OO),
         temp_range = mean(arrivals$temp_range),
         EPPO_reporting_service = mean(arrivals$EPPO_reporting_service),
         trade_hort = mean(arrivals$trade_hort),
         date = mean(arrivals$date))
pp_average(top_models[[1]],
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
           newdata = newdata,
           allow_new_levels = TRUE,
           sample_new_levels = "gaussian",
           summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  mutate(climate_match_raw = newdata$climate_match)
# epred_draws(best_model,
#             newdata = newdata,
#             allow_new_levels = TRUE,
#             sample_new_levels = "gaussian") %>% 
#   group_by(climate_match, phyto) %>%
#   summarise(.epred = quantile(.epred, c(0.025, 0.5, 0.975))) %>% ungroup()%>%
#   mutate(quantile = rep(c("2.5%", "50%", "97.5%"),8)) %>%
#   pivot_wider(names_from = quantile, values_from = .epred ) %>%
#   mutate(climate_match_raw = (climate_match * climate_scale + climate_center))


###### reporting interactions between risk factors and traits (trait-mediated effects) ####
# specify the absolute change in probability at specific values of trade and other risk factors
# as well as odds ratios at the mean
# first the odds ratio (this is the slope at the mean of the curve)
which(!is.na(model_info[,"b_CH_OO:EPPO_reporting_service","est"]))

averaged_models_incl_surv_CHOO <- posterior_average(top_models[[1]],
                                                     top_models[[2]], 
                                                     #top_models[[3]],
                                                     #top_models[[4]],
                                                     top_models[[5]], 
                                                     top_models[[6]], 
                                                     #top_models[[7]],
                                                     #top_models[[8]],
                                                     #top_models[[9]], 
                                                     #top_models[[10]], 
                                                     top_models[[11]],
                                                     #top_models[[12]],
                                                     top_models[[13]],
                                                     weights = "loo") %>% as_draws()
# odds ratios/slopes for each trait group

averaged_models_incl_surv_CHOO %>% select(`b_CH_OO:EPPO_reporting_service`) %>%
  mutate(`b_CH_OO:EPPO_reporting_service` = (1.1^`b_CH_OO:EPPO_reporting_service`)*100) %>%
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) 

quantile(arrivals_unscaled$EPPO_reporting_service, c(0.25,0.5, 0.75))
newdata <- expand.grid((log(quantile(arrivals_unscaled$EPPO_reporting_service, c(0.25,0.5, 0.75)))-surv_center)/surv_scale,
                       sort(unique(arrivals$CH_OO))) %>%
  rename(EPPO_reporting_service = Var1, CH_OO = Var2) %>%
  mutate(importer_iso3 = "any",
         trade_hort = mean(arrivals$trade_hort), 
         phyto = "unknown_species",
         phylo = "unknown_species",
         temp_range = mean(arrivals$temp_range),
         climate_match = mean(arrivals$climate_match),
         date = mean(arrivals$date))
pp_average(top_models[[1]],
           top_models[[2]], 
           #top_models[[3]],
           #top_models[[4]],
           top_models[[5]], 
           top_models[[6]], 
           #top_models[[7]],
           #top_models[[8]],
           #top_models[[9]], 
           #top_models[[10]], 
           top_models[[11]],
           #top_models[[12]],
           top_models[[13]],
           weights = "loo", method = "posterior_epred",
           newdata = newdata,
           allow_new_levels = TRUE,
           sample_new_levels = "gaussian",
           summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  mutate(EPPO_reporting_service = newdata$EPPO_reporting_service,
         CH_OO = newdata$CH_OO)


# epred_draws(best_model,
#             newdata = newdata,
#             allow_new_levels = TRUE,
#             sample_new_levels = "gaussian") %>%
#   group_by(CH_OO, EPPO_reporting_service) %>% 
#   summarise(lower = round(quantile(.epred, 0.025), 2),
#             est = round(quantile(.epred, 0.5),2),
#             upper = round(quantile(.epred, 0.975), 2))

#### interactions between trade and survival structures ###
which(!is.na(model_info[,"b_CH_OO:trade_hort","est"]))

averaged_models_incl_trade_CHOO <- posterior_average(#top_models[[1]],
                                                    top_models[[2]], 
                                                    top_models[[3]],
                                                    top_models[[4]],
                                                    top_models[[5]], 
                                                    top_models[[6]], 
                                                    top_models[[7]],
                                                    #top_models[[8]],
                                                    #top_models[[9]], 
                                                    #top_models[[10]], 
                                                    #top_models[[11]],
                                                    #top_models[[12]],
                                                    #top_models[[13]],
                                                    weights = "loo") %>% as_draws()
# odds ratios/slopes for each trait group

averaged_models_incl_trade_CHOO %>% select(`b_CH_OO:trade_hort`) %>%
  mutate(`b_CH_OO:trade_hort` = (1.1^`b_CH_OO:trade_hort`)*100) %>%
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) 

scale2sd(log(1+quantile(arrivals_unscaled$trade_hort, c(0.25,0.5, 0.75))))

quantile(arrivals_unscaled$trade_hort, c(0.25,0.5, 0.75))
newdata <- expand.grid((log(1+quantile(arrivals_unscaled$trade_hort, c(0.25,0.5, 0.75)))-trade_center)/trade_scale,
                       unique(sort(arrivals$CH_OO))) %>%
  rename(trade_hort = Var1, CH_OO = Var2) %>%
  mutate(importer_iso3 = "any",
         EPPO_reporting_service = mean(arrivals$EPPO_reporting_service), 
         phyto = "unknown_species",
         phylo = "unknown_species",
         temp_range = mean(arrivals$temp_range),
         climate_match = mean(arrivals$climate_match),
         date = mean(arrivals$date))
pp_average(#top_models[[1]],
  top_models[[2]], 
  top_models[[3]],
  top_models[[4]],
  top_models[[5]], 
  top_models[[6]], 
  top_models[[7]],
  #top_models[[8]],
  #top_models[[9]], 
  #top_models[[10]], 
  #top_models[[11]],
  #top_models[[12]],
  #top_models[[13]],
  weights = "loo", method = "posterior_epred",
           newdata = newdata,
           allow_new_levels = TRUE,
           sample_new_levels = "gaussian",
           summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  mutate(trade_hort = newdata$trade_hort,
         CH_OO = newdata$CH_OO)

# epred_draws(best_model,
#             newdata = newdata,
#             allow_new_levels = TRUE,
#             sample_new_levels = "gaussian") %>%
#   group_by(CH_OO, trade_hort) %>% 
#   summarise(lower = round(quantile(.epred, 0.025), 2),
#             est = round(quantile(.epred, 0.5),2),
#             upper = round(quantile(.epred, 0.975), 2))

  


