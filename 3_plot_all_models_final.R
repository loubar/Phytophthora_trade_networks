
rm(list=ls())
require(brms)
require(dplyr)
require(reshape2)
library(tidybayes)
library(ggplot2)
library(sf)
library(ggpubr)
library(abind)
library(DHARMa)
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



countries_per_species <- aggregate(arrival~phyto, data = arrivals,
                                   FUN = sum)
range(countries_per_species$arrival)
median(countries_per_species$arrival)
species_per_country <- aggregate(arrival~importer_iso3, data = arrivals,
                                   FUN = sum)
range(species_per_country$arrival)
median(species_per_country$arrival)
# number of new countries invded per per species since 2005
# countries by continent: have we captured sufficient climate, trade and surveillance across these regions?
add_continent <- left_join(arrivals_unscaled |> 
                 rename(iso3c = importer_iso3), countrycode::codelist) 


options(scipen = 10)
# plot variation in trade, climate and surveillance among countries
ggplot(data = arrivals_unscaled |> 
         pivot_longer(cols = c(trade_hort,
                               climate_distance, 
                               EPPO_reporting_service ),
                      names_to = "covariate", values_to = "value" ) |>
         mutate(covariate = factor(covariate,
                                   levels = c("trade_hort",
                                              "climate_distance", 
                                              "EPPO_reporting_service" ),
                                   labels = c("Horticultural trade connectivity to source country (1 + 1000 kgs)",
                                              "Climatic distance from source country (Mahalanobis distance)",
                                              "Surveillance effort (Number of items in EPPO Reporting Service archive)")))) + 
  geom_histogram(aes(x = 1+value)) +
  scale_x_continuous(trans = "log10") +
  facet_wrap(~covariate, scales = "free_x", labeller = label_wrap_gen(width = 25, multi_line = TRUE)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0))
ggsave(filename = "covariate_variance.png",
       height = 5, width = 6, units = "in",
       dpi = 500)

# plot correlations among trade, climate matching and surveillance
library(GGally)
ggpairs(arrivals |> select(trade_hort,
                           climate_distance,
                           EPPO_reporting_service,
                           date,
                           temp_range,
                           CH_OO) |>
        mutate(CH_OO = factor(CH_OO,
                                 labels = c(0,1,2))),
        upper = list(continuous = GGally::wrap(ggally_cor, stars = F)),
        columnLabels = c("Horticultural trade connectivity",
                         "Climate distance",
                         "Surveillance",
                         "year species described",
                         "Thermal tolerance\nrange",
                         "Number of survival structures"),
        labeller = label_wrap_gen(width = 12, 
                                  multi_line = TRUE)) +
  theme(axis.text.x = element_text(size = 8,
                                 angle = 45, 
                                 vjust = 1, 
                                 hjust = 1))
ggsave(filename = "variable_correlations.png",
       height = 6.5, width = 6.5, units= "in", dpi = 500)  

# VIF to check for collinearity
# using full model (all covariates)
fit0178 <- readRDS(paste0(minRE_spp, "/", com_type, "/", "fit0178.rds"))
check_collinearity(fit0178)

rownames(model_info)<- rank_models_all$model[rank_models_all$diffic >=-2]


#model <- rank_models_all$model[9]
get_sig <- function(model){ 
  b_pars <- c( "b_CH_OO",
               "b_climate_match" ,
               "b_date",
               "b_EPPO_reporting_service",
               "b_temp_range" ,
               "b_trade_hort",
               "b_CH_OO:EPPO_reporting_service",
               "b_CH_OO:trade_hort",
               "b_climate_match:temp_range",
               #"b_temp_range:trade_hort",
               "b_zi_EPPO_reporting_service",
               "b_zi_date")
  pos_sig <- model_info[model,b_pars,"upper"] > 0 & model_info[model,b_pars,"lower"] > 0
  neg_sig <- model_info[model,b_pars,"upper"] < 0 & model_info[model,b_pars,"lower"] < 0
  sig <- ifelse(pos_sig|neg_sig, yes = "*", no = " ")
  return(data.frame(model = model, variable = names(sig), sig = sig ,
                    stringsAsFactors = FALSE,
                    row.names = NULL))
  
}

model_summary_df <- as.data.frame(model_info[,,"est"])
model_summary_df$rank <- 1:nrow(model_info)
model_summary_df$model <- rank_models_all$model[rank_models_all$diffic >=-2]




model_summary_df <- left_join(rank_models_all[rank_models_all$diffic >=-2,], 
                              model_summary_df)
model_summary_df$rank <- 1:nrow(model_summary_df)
#saveRDS(model_summary_df, file = paste0(minRE_spp ,"/model_summary_df.rds"))

# prepare data for plotting heat map of effect sizes
model_summary_df <- model_summary_df[,c("model" ,
                                        "rank",
                    "diffic",
                    #"b_Intercept" ,
                    #"b_zi_Intercept",
                    "b_CH_OO",
                    "b_climate_match" ,
                    "b_date",
                    "b_EPPO_reporting_service",
                    "b_temp_range" ,
                    "b_trade_hort",
                    "b_CH_OO:EPPO_reporting_service",
                    "b_CH_OO:trade_hort",
                    "b_climate_match:temp_range",
                    "b_temp_range:trade_hort",
                    "b_zi_EPPO_reporting_service",
                    "b_zi_date",
                    "R2_marginal",
                    "R2_conditional",
                    "ICC_phylo",
                    "ICC_phyto",
                    "ICC_importer_iso3",
                    "ICC_phyto__trade_hort",
                    "ICC_phyto__climate_match",
                    "ICC_phyto__EPPO_reporting_service")]

# plot the top model subset

# 13 models

#check the variance component add up to the conditional R2
partitions <- c("R2_marginal",  "ICC_phylo", "ICC_phyto", "ICC_importer_iso3", "ICC_phyto__climate_match", "ICC_phyto__EPPO_reporting_service", "ICC_phyto__trade_hort")
rowSums(model_info[,partitions,"est_same_draw"])



model_summary_df$rank <- as.integer(model_summary_df$rank)
model_performance <- model_summary_df[,c("rank",
                                         "diffic",
                                         "R2_marginal",
                                         "R2_conditional",
                                         "ICC_phylo",
                                         "ICC_phyto",
                                         "ICC_importer_iso3",
                                         "ICC_phyto__trade_hort",
                                         "ICC_phyto__climate_match",
                                         "ICC_phyto__EPPO_reporting_service")]
model_performance <- melt(model_performance, id.vars = "rank")


model_summary_df[,c("diffic",
                    "R2_marginal",
                    "R2_conditional",
                    "ICC_phylo",
                    "ICC_phyto",
                    "ICC_importer_iso3",
                    "ICC_phyto__trade_hort",
                    "ICC_phyto__climate_match",
                    "ICC_phyto__EPPO_reporting_service")] <- NA

parest <- melt(model_summary_df, 
               id.vars = c("model", "rank"))
#parest$rank <- as.integer(parest$rank)

sig_df <- lapply(rank_models_all$model[rank_models_all$diffic  >=-2 ],
                 get_sig)
sig_df <- bind_rows(sig_df)

nrow(sig_df <- left_join(sig_df, unique(parest[, c("model", "rank")])))


#order variables by number of times they appear in the top subset (n=14)
#top <- max(which(rank_models$diffic>=-5))
#times_in_top_subset <- function(x){sum(!is.na(x))}
#mean_parest <- function(x){mean(x, na.rm = TRUE)}

# order_vars <- sort(abs(tapply(parest$value, 
#                               parest$variable,
#                               mean_parest)) + tapply(parest$value, 
#                                                      parest$variable,
#                                                      times_in_top_subset))

# group variable as a) zero-inflation (non-detection), b) arrival risk factors, c) pathogen traits d) trait-mediated effects

order_vars <- rev(c("b_zi_EPPO_reporting_service", "b_EPPO_reporting_service", # importer surveillance
                    "b_zi_date", "b_date", # knowledge of pathogen
                    "b_temp_range", "b_CH_OO", # pathogen traits
                    "b_trade_hort", "b_climate_match", # arrival risk factors
                    "b_CH_OO:EPPO_reporting_service", "b_CH_OO:trade_hort", "b_temp_range:trade_hort", "b_climate_match:temp_range")) # trait-mediated effects



parest$variable <- factor(parest$variable,
                          c(rev(c("diffic",
                                  "R2_marginal",
                                  "R2_conditional",
                                  "ICC_phylo",
                                  "ICC_phyto",
                                  "ICC_importer_iso3",
                                  "ICC_phyto__trade_hort",
                                  "ICC_phyto__climate_match",
                                  "ICC_phyto__EPPO_reporting_service")),
                          #names(order_vars)
                          order_vars))                            

#model_summary_df <-array(model_summary, dim = c(178,21,3))



gg <- ggplot() 
gg <- gg + geom_tile(data = parest,aes(rank, variable, fill = value)) 
gg <- gg + scale_fill_gradient2(low = "blue",
                                high = "red",
                                mid = "white",
                                name = expression(bold(atop(standardised, effect~size~(beta)))),
                                midpoint = 0)

# gg <- gg + facet_grid(~impact_metric,
#                       drop = TRUE,
#                       scales = "free",
#                       space = "free"
# )
#make.bold <- function(x) as.expression(lapply(x, function(y) bquote(bold(.(y)))))

gg <- gg + scale_y_discrete(labels=c("b_zi_EPPO_reporting_service" = expression(bold(zero-inflation%~%surveillance~effort)),
                                     "b_trade_hort" = expression(bold(horticultural~imports)),
                                     "b_EPPO_reporting_service" = expression(bold(surveillance~effort)),
                                     "b_CH_OO:EPPO_reporting_service" = expression(bold(surveillance~effort~x~survival~structures)),
                                     "b_zi_date" = expression(bold(zero-inflation%~%years~since~described)),
                                     "b_CH_OO:trade_hort" = expression(bold(horticultural~imports~x~survival~structures)),
                                     "b_climate_match" = expression(bold(climate~matching)),
                                     "b_CH_OO" = expression(bold(survival~structures)),
                                     "b_temp_range" = expression(bold(thermal~tolerance~range)),
                                     "b_date" = expression(bold(years~since~described)),
                                     "b_temp_range:trade_hort" = expression(bold(horticultural~imports~x~thermal~tolerance~range)),
                                     "b_climate_match:temp_range" = expression(bold(climate~matching~x~thermal~tolerance~range)),
                                     "R2_conditional" = expression(bold(conditional~R^2)),
                                     "R2_marginal" = expression(bold(marginal~R^2)),
                                     #"ICCadj" = expression(bold(adjusted~phylogenetic~intra-class~correlation)),
                                     "ICC_phylo" = expression(bold(phylogenetic~intra-class~correlation)),
                                     "ICC_phyto" = expression(bold(species-level~intra-class~correlation)),
                                     "ICC_importer_iso3" = expression(bold(country-level~intra-class~correlation)),
                                     "ICC_phyto__trade_hort" = expression(bold(species-specific~responses~to~trade)),
                                     "ICC_phyto__climate_match" = expression(bold(species-specific~responses~to~climate-matching)),
                                     "ICC_phyto__EPPO_reporting_service" = expression(bold(species-specific~responses~to~surveillance~effort)),
                                     "diffic" = expression(bold(Delta~IC~(LOO~cross-validation)))),
                            
                            
                            # "R2_conditional" = expression(conditional~R^2),
                            # "R2_traits" = expression(marginal~R^2),
                            # "diff_IC" = expression(Delta~IC~(10-fold~cross-validation)),
                            # "lambda" = expression(phylogenetic~signal~(lambda))),  
                            
                            expand = c(0, 0)) + ylab("") 
gg <- gg + geom_text(data = sig_df,
                               aes(x = rank,
                                   y = variable,
                                   label = sig))
gg <- gg + geom_tile(data = model_performance,aes(rank, variable), fill = "white",
                     show.legend = FALSE) 
gg <- gg + scale_x_continuous("model rank\n(approximate leave-one-out cross-validation)",
                              breaks = c(1:13),
                              expand = c(0,0))

#gg <- gg + scale_x_discrete(expand = c(0,0))
gg <- gg + theme(#axis.text.x=element_text(angle=90, vjust = 0.5, hjust = #1),
  #text = element_text(size = 5),               
  axis.text = element_text(colour = "black", size = 10, face = "bold"),
  axis.title.x = element_text(size = 16, face = "bold"),
  legend.title = element_text(size = 12),
  panel.grid = element_blank(),
  panel.background = element_rect(fill = "white"),
  panel.border = element_rect(colour = "grey", fill = NA),
  plot.margin=unit(c(0,0,0,2), "cm"))

gg <- gg + geom_text(data = model_performance,
                     aes(x = rank,
                         y = variable, 
                         label  = sprintf("%.2f", round(value, 2)),
                         angle = 0),
                     size = 3,
                     show.legend = FALSE)


gg <- gg + annotate("rect", xmin = -14, xmax = 0.5, ymin = 19.5, ymax = 21.5,
              alpha = .2) + 
  coord_cartesian(xlim = c(0.5, 13.5), ylim = c(0.5, 21.5), clip = "off") + 
  annotate("text", x = -10, y = 20.5, label = "importer\nsurveillance", angle = 90)
gg<- gg + annotate("text", x = -10, y = 18.5, label = "knowledge of\npathogen", angle = 90)
gg <- gg + annotate("rect", xmin = -12, xmax = 0.5, ymin = 15.5, ymax = 17.5,
              alpha = .2) + 
  annotate("text", x = -10, y = 16.5, label = "pathogen\ntraits", angle = 90) 
gg <- gg +  annotate("text", x = -10, y = 14.5, label = "connectivity\nrisk factors", angle = 90)  
gg <- gg + annotate("rect", xmin = -12, xmax = 0.5, ymin = 9.5, ymax = 13.5,
              alpha = .2) + 
  annotate("text", x = -10, y = 11.5, label = "trait-mediated\neffects", angle = 90) 
gg <- gg + annotate("rect", xmin = -12, xmax = 13.4, ymin = 0.59, ymax = 8.5,
              alpha = .2, fill = NA, col = "black") + 
  annotate("text", x = -10, y = 5, label = "variance components", angle = 90) 
gg
ggsave(paste0(minRE_spp, "/model_comparison_top_subset.png"),
       height = 10.4, width = 10.2,
       units = "in", dpi = 500
       )


############ report and plot the parameter estimates averaged across the top model subset #################

best_models <- rank_models_all[rank_models_all$diffic >=-2,]


top_subset <- paste0(minRE_spp, "/best_models/", best_models$model)
top_models <- lapply(top_subset, readRDS)



order_vars <- rev(c("b_zi_EPPO_reporting_service", "b_EPPO_reporting_service", # importer surveillance
                    "b_zi_date", "b_date", # knowledge of pathogen
                    "b_temp_range", "b_CH_OO", # pathogen traits
                    "b_trade_hort", "b_climate_match", # arrival risk factors
                    "b_CH_OO:EPPO_reporting_service", "b_CH_OO:trade_hort",  "b_climate_match:temp_range")) # trait-mediated effects

variable_present <- t(!apply(model_info[,grep("^b_", colnames(model_info)),"est"], 1, is.na))
# which variables are absent from at least one of the top models
missing_variable_names <- names(which(apply(variable_present, 2, any)))
# for each get the models where they are present
average_over <- apply(variable_present[,order_vars], 2, which)

present_in_all <- which(lapply(average_over, length) == 13)
present_in_any <- which(lapply(average_over, length) >0)
average_over <- average_over[present_in_any]

# for each parameter estimate, only average over the variables that are present
# which variables are present in all top models?
model_names_to_av <- lapply(1:length(average_over),
                       function(i)
                         gsub("\\.rds", "", best_models$model[average_over[[i]]]))
names(model_names_to_av) <- names(average_over)



# model_av <- posterior_average(top_models[[1]],
#                               top_models[[2]],
#                               top_models[[3]],
#                               top_models[[4]],
#                               top_models[[5]],
#                               top_models[[6]],
#                               top_models[[7]],
#                               top_models[[8]],
#                               top_models[[9]],
#                               top_models[[10]],
#                               top_models[[11]],
#                               top_models[[12]],
#                               top_models[[13]],
#                               weights = "loo",
#                               variable = order_vars,
#                               missing = 0) # assumes a parameter estimate of 0 for variables that are absent from a model

model_av <- posterior_average(top_models[[1]],
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
                              variable = names(present_in_all))


setdiff(names(average_over), names(model_av))
average_over[["b_climate_match:temp_range"]]
model_av$`b_climate_match:temp_range` <- posterior_average(top_models[[5]],
                                                           top_models[[7]],
                                                           weights = "loo",
                                                           variable = "b_climate_match:temp_range")
average_over[["b_CH_OO:trade_hort"]]
model_av$`b_CH_OO:trade_hort` <- posterior_average(top_models[[2]],
                                                           top_models[[3]],
                                                           top_models[[4]],
                                                           top_models[[5]],
                                                           top_models[[6]],
                                                           top_models[[7]],
                                                           weights = "loo",
                                                           variable = "b_CH_OO:trade_hort")
average_over[["b_CH_OO"]]
model_av$b_CH_OO <- posterior_average(top_models[[1]],
                                           top_models[[2]],
                                           top_models[[3]],
                                           top_models[[4]],
                                           top_models[[5]],
                                           top_models[[6]],
                                           top_models[[7]],
                                           top_models[[9]],
                                           top_models[[11]],
                                           top_models[[12]],
                                           top_models[[13]],
                                           weights = "loo",
                                           variable = "b_CH_OO")

average_over[["b_CH_OO:EPPO_reporting_service"]]
model_av$`b_CH_OO:EPPO_reporting_service` <- posterior_average(top_models[[1]],
                                      top_models[[2]],
                                      top_models[[5]],
                                      top_models[[6]],
                                      top_models[[11]],
                                      top_models[[13]],
                                      weights = "loo",
                                      variable = "b_CH_OO:EPPO_reporting_service")

average_over[["b_temp_range"]]
model_av$b_temp_range <- posterior_average(top_models[[1]],
                                           top_models[[4]],
                                           top_models[[5]],
                                           top_models[[6]],
                                           top_models[[7]],
                                           top_models[[9]],
                                           top_models[[10]],
                                           weights = "loo",
                                           variable = "b_temp_range"
                                                               )
average_over[["b_date"]]
model_av$b_date <- posterior_average(top_models[[1]],
                                           top_models[[2]],
                                           top_models[[3]],
                                           top_models[[4]],
                                           top_models[[5]],
                                           top_models[[6]],
                                           top_models[[7]],
                                           top_models[[8]],
                                           top_models[[9]],
                                           top_models[[10]],
                                           top_models[[12]],
                                           top_models[[13]],
                                           weights = "loo",
                                           variable = "b_date"
)



# draws for model averaging across top performing models
plot_model_av <- model_av %>% as_draws() %>% gather_draws(`b_.*`, regex = TRUE) %>%
  filter(.variable %in% intersect(order_vars, .variable))  #%>%
  #mutate(model = "best model with interactions")

plot_model_av$effect <- ifelse(plot_model_av$.value > 0,
                                      yes = "beta >= 0",
                                      no = "beta < 0")
plot_model_av$submodel <- ifelse(plot_model_av$.variable %in% c("b_zi_EPPO_reporting_service",
                                                                              "b_zi_date"),
                                        yes = "non-detection\nsub-model",
                                        no = "emergence\nmodel")
plot_model_av$submodel <- factor(plot_model_av$submodel, 
                                 levels = c("non-detection\nsub-model",
                                            "emergence\nmodel"))

plot_model_av$.variable <- factor(plot_model_av$.variable, 
                                  levels = order_vars)

#plot_model_av$model <- "Model averaging"

# plot_model_best <- best_model  %>% gather_draws(`b_.*`, regex = TRUE) %>%
#   filter(.variable %in% intersect(order_vars, .variable))  %>%
#   mutate(model = "Top-ranked model")
# plot_model_best$submodel <- ifelse(plot_model_best$.variable %in% c("b_zi_EPPO_reporting_service",
#                                                                                "b_zi_date"),
#                                               yes = "non-detection\nsub-model",
#                                               no = "emergence\nmodel")
# 
# plot_models <- bind_rows(plot_model_best, plot_model_av)
# 
# plot_models$submodel <- factor(plot_models$submodel, levels = c("non-detection\nsub-model",
#                                                                     "emergence\nmodel"))
# 
# plot_models$.variable <- factor(plot_models$.variable,
#                                   levels = order_vars[order_vars %in% plot_models$.variable])
# 
# plot_models$model <- factor(plot_models$model, 
#                             levels = c("Top-ranked model",
#                                        "Model averaging"))

gg1 <- ggplot() +
  stat_halfeye(data = plot_model_av, aes(y = .variable, x = .value, fill = stat(x > 0)), size = 0.2) + 
  scale_fill_manual("", values = c("steelblue", "salmon"), labels = c(expression(~beta < 0), expression(~beta > 0))) +
  xlab(expression(Posterior~distributions~of~standardised~effects~on~probability~of~new~italic(Phytophthora)~detection))

gg2 <- gg1 + scale_y_discrete(labels=c("b_zi_EPPO_reporting_service" = expression(bold(zero-inflation%~%surveillance~effort)),
                                       "b_trade_hort" = expression(bold(horticultural~imports)),
                                       "b_EPPO_reporting_service" = expression(bold(surveillance~effort)),
                                       "b_CH_OO:EPPO_reporting_service" = expression(bold(surveillance~effort~x~survival~structures)),
                                       "b_zi_date" = expression(bold(zero-inflation%~%years~since~described)),
                                       "b_CH_OO:trade_hort" = expression(bold(horticultural~imports~x~survival~structures)),
                                       "b_climate_match" = expression(bold(climate~matching)),
                                       "b_CH_OO" = expression(bold(survival~structures)),
                                       "b_temp_range" = expression(bold(thermal~tolerance~range)),
                                       "b_date" = expression(bold(years~since~described)),
                                       #"b_temp_range:trade_hort" = expression(bold(horticultural~imports~x~thermal~tolerance~range)),
                                       "b_climate_match:temp_range" = expression(bold(climate~matching~x~thermal~tolerance~range)))) + ylab("") 

gg3 <- gg2 + theme(#axis.text.x=element_text(angle=90, vjust = 0.5, hjust = #1),
  #text = element_text(size = 5),               
  axis.text.x = element_text(colour = "black", size = 12),
  axis.text.y = element_text(colour = "black", size = 9),
  axis.title.x = element_text(size = 11, face = "bold"),
  legend.text = element_text(size = 12),
  legend.position = "bottom",
  strip.text = element_text(size = 9) , 
  #panel.grid = element_blank(),
  panel.background = element_rect(fill = "white"),
  panel.border = element_rect(colour = "grey", fill = NA),
  plot.margin=unit(c(0,0,0,3), "cm")) + geom_vline(xintercept = 0) +
  facet_grid(rows = vars(submodel), scales = "free_y", space = "free")

ann_boxes <- data.frame(.value = c(-32, -32),
                        .variable = grep("zi", order_vars, value = TRUE, invert = TRUE)[c(1, 10)],
                        xmax = c(-11,  -11),
                        ymin = c(8.5, 0.5),
                        ymax = c(5.5, 3.5),
                        submodel = factor("emergence\nmodel", levels = c("non-detection\nsub-model",
                                                                         "emergence\nmodel")),
                        model = factor("Top-ranked model", 
                                       levels = c("Top-ranked model",
                                                  "Model averaging")))

gg4 <- gg3 + geom_rect(data = ann_boxes, aes(xmin = .value, xmax = xmax, 
                                             ymin = ymin, ymax= ymax),
                       fill = "lightblue",
                       col = NA,
                       alpha = 0.2) +
  coord_cartesian(xlim = c(-10, 8.5), clip = "off")

ann_text <- data.frame(.value = c(-23, -23, -23, -23),
                       .variable = grep("zi", order_vars, value = TRUE, invert = TRUE)[c(1,4,6, 10)],
                       y_coord = c(9, 7, 4.5, 2 ),
                       y_label = c("importer surveillance\n(country-level)",
                                   "pathogen traits\n(species-level)",
                                   "connectivity\nrisk factors\n(species x country level)",
                                   "trait-mediated effects\n(interaction terms)"),
                       submodel = factor("emergence\nmodel", levels = c("non-detection\nsub-model",
                                                                        "emergence\nmodel")), 
                       model = factor("Top-ranked model", 
                                      levels = c("Top-ranked model",
                                                 "Model averaging")))

gg5 <- gg4 + geom_text(data = ann_text, aes(x = .value, y = y_coord,
                                            label = y_label ), size = 3)

gg5

ggsave(paste0(minRE_spp, "/best_model_effects.png"),
       height = 19, width = 25,
       units = "cm", dpi = 500
)



# 
# # best model with interactions
# best_model_int_draws <- best_model_int  %>% gather_draws(`b_.*`, regex = TRUE) %>%
#   filter(.variable %in% intersect(order_vars, .variable))  %>%
#   mutate(model = "best model with interactions")
# 
# # best model overall
# best_model_draws <- best_model  %>% gather_draws(`b_.*`, regex = TRUE) %>%
#   filter(.variable %in% intersect(order_vars, .variable)) %>%
#   mutate(model = "best model overall")
# 
# 
# best_model_int_draws <- bind_rows(best_model_int_draws, best_model_draws)
# 
# best_model_int_draws$effect <- ifelse(best_model_int_draws$.value > 0,
#                                       yes = "beta >= 0",
#                                       no = "beta < 0")
# best_model_int_draws$submodel <- ifelse(best_model_int_draws$.variable %in% c("b_zi_EPPO_reporting_service",
#                                                                               "b_zi_date"),
#                                         yes = "non-detection\nsub-model",
#                                         no = "emergence\nmodel")
# best_model_int_draws$submodel <- factor(best_model_int_draws$submodel, levels = c("non-detection\nsub-model",
#                                                                                   "emergence\nmodel"))
# best_model_int_draws$.variable <- factor(best_model_int_draws$.variable,
#                                          levels = order_vars[order_vars %in% best_model_int_draws$.variable])
# 
# 
# best_model_int_draws$model <- factor(best_model_int_draws$model,
#                                      levels = c("best model overall",
#                                                 "best model with interactions"))
# 
# library(scales)

# gg1 <- ggplot() +
#   stat_halfeye(data = best_model_int_draws, aes(y = .variable, x = .value, fill = stat(x > 0)), size = 2) + 
#   scale_fill_manual("", values = c("steelblue", "salmon"), labels = c(expression(~beta < 0), expression(~beta > 0))) +
#   xlab(expression(Posterior~distributions~of~standardised~effects~on~probability~of~new~italic(Phytophthora)~detection))
# gg2 <- gg1 + scale_y_discrete(labels=c("b_zi_EPPO_reporting_service" = expression(bold(zero-inflation%~%surveillance~effort)),
#                                        "b_trade_hort" = expression(bold(horticultural~imports)),
#                                        "b_EPPO_reporting_service" = expression(bold(surveillance~effort)),
#                                        "b_CH_OO:EPPO_reporting_service" = expression(bold(surveillance~effort~x~survival~structures)),
#                                        "b_zi_date" = expression(bold(zero-inflation%~%years~since~described)),
#                                        "b_CH_OO:trade_hort" = expression(bold(horticultural~imports~x~survival~structures)),
#                                        "b_climate_match" = expression(bold(climate~matching)),
#                                        "b_CH_OO" = expression(bold(survival~structures)),
#                                        "b_temp_range" = expression(bold(thermal~tolerance~range)),
#                                        "b_date" = expression(bold(years~since~described)),
#                                        #"b_temp_range:trade_hort" = expression(bold(horticultural~imports~x~thermal~tolerance~range)),
#                                        "b_climate_match:temp_range" = expression(bold(climate~matching~x~thermal~tolerance~range)))) + ylab("")
# 
# gg3 <- gg2 + theme(#axis.text.x=element_text(angle=90, vjust = 0.5, hjust = #1),
#   #text = element_text(size = 5),               
#   axis.text.x = element_text(colour = "black", size = 12),
#   axis.text.y = element_text(colour = "black", size = 9),
#   axis.title.x = element_text(size = 11, face = "bold"),
#   legend.text = element_text(size = 12),
#   legend.position = "bottom",
#   strip.text = element_text(size = 9) , 
#   #panel.grid = element_blank(),
#   panel.background = element_rect(fill = "white"),
#   panel.border = element_rect(colour = "grey", fill = NA),
#   plot.margin=unit(c(0,0,0,3), "cm")) + geom_vline(xintercept = 0) +
#   facet_grid(submodel~model , scales = "free_y", space = "free")
# 
# ann_boxes <- data.frame(.value = c(-32, -32),
#                         .variable = grep("zi", order_vars, value = TRUE, invert = TRUE)[c(1, 10)],
#                         xmax = c(-8.25,  -8.25),
#                         ymin = c(8.5, 0.5),
#                         ymax = c(5.5, 3.5),
#                         submodel = factor("emergence\nmodel", levels = c("non-detection\nsub-model",
#                                                                         "emergence\nmodel")),
#                         model = factor("best model overall", levels = c("best model overall",
#                                                                         "best model with interactions")))
# 
# gg4 <- gg3 + geom_rect(data = ann_boxes, aes(xmin = .value, xmax = xmax, 
#                                              ymin = ymin, ymax= ymax),
#                        fill = "lightblue",
#                        col = NA,
#                        alpha = 0.2) +
#   coord_cartesian(xlim = c(-7.5, 7.5), clip = "off")
# 
# ann_text <- data.frame(.value = c(-28, -28, -28, -28),
#                        .variable = grep("zi", order_vars, value = TRUE, invert = TRUE)[c(1,4,6, 10)],
#                        y_coord = c(9, 7, 4.5, 2 ),
#                        y_label = c("importer surveillance\n(country-level)",
#                                    "pathogen traits\n(species-level)",
#                                    "connectivity\nrisk factors\n(species x country level)",
#                                    "trait-mediated effects\n(interaction terms)"),
#                        submodel = factor("emergence\nmodel", levels = c("non-detection\nsub-model",
#                                                                        "emergence\nmodel")),
#                        model = factor("best model overall", levels = c("best model overall",
#                                                                        "best model with interactions")))
# 
# gg5 <- gg4 + geom_text(data = ann_text, aes(x = .value, y = y_coord,
#                                             label = y_label ), size = 3)
# 
# gg5
# 
# ggsave(paste0(minRE_spp, "/best_model_effects.png"),
#        height = 19, width = 25,
#        units = "cm", dpi = 500
# )

# species-level variation in responses to trade, climate matching and
# 1) trade
# we present the results of the model including only the main effects of the risk factors
# to show the extent of the unexplained variability between species in the absence of the interaction terms
# find the model with main effects of risk factors only (e.g. no species-level effects)

# this time also extract the random effects from the model averaging
model_av <- posterior_average(top_models[[1]],
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
                              weights = "loo")

b_draws <- model_av %>% as_draws() %>% spread_draws(b_trade_hort) %>%
  rename(.value = b_trade_hort) %>%
  mutate(phyto = "Average effect") %>%
  select(.value, .chain, .iteration, .draw, phyto) %>%
  mutate(risk_factor = "trade_hort")


r_draws <- model_av %>% as_draws() %>%
  spread_draws(b_trade_hort, r_phyto[phyto,risk_factor]) %>% 
  filter(risk_factor == "trade_hort" | is.na(risk_factor)) %>%
  mutate(.value = b_trade_hort + r_phyto) %>%
  select(risk_factor, .chain, .iteration, .draw, phyto, .value)



species_sensitivity <- as_tibble(bind_rows(b_draws, r_draws))

species_sensitivity$phyto <- gsub("\\.", " ", species_sensitivity$phyto)

# highlight the species with both survival structures with an asterisk
# spp_surv_str <- left_join(data.frame(phyto = unique(species_sensitivity$phyto)[2:71]),
#                           included_traits[,c("phyto", "CH_OO")])

# str2 <- spp_surv_str$phyto[spp_surv_str$CH_OO == 2]
# 
# species_sensitivity$phyto <- ifelse(species_sensitivity$phyto %in% str2,
#        yes = paste0("* ", species_sensitivity$phyto),
#        no = species_sensitivity$phyto)

species_sensitivity$phyto <- gsub("Phytophthora", "P.", species_sensitivity$phyto)

species_sensitivity$phyto <- factor(species_sensitivity$phyto, 
                        levels = c("Average effect",
                                   grep("Average effect", 
                                        names(sort(tapply(species_sensitivity$.value, 
                                                          species_sensitivity$phyto, 
                                                          FUN = "median"),
                                                   decreasing = TRUE)),
                                        invert = TRUE, value = TRUE)))




g <- ggplot() + stat_halfeye(data = species_sensitivity,
                      aes(x = .value, y = phyto, fill = stat(x>0), col = stat(x>0)), 
                      #scale = scale,
                      alpha = 0.6,
                      size = 1) + 
  scale_fill_manual("", values = c("steelblue", "salmon"), labels = c(expression(~~~beta<0), expression(~~~beta>0)))+
  scale_colour_manual("", values = c("steelblue", "salmon"), labels = c(expression(~~~beta<0), expression(~~~beta>0)))
g <- g + geom_vline(xintercept = 0, linetype = 1) +
  geom_vline(xintercept = median(species_sensitivity$.value[species_sensitivity$phyto == "Average effect"]),
             linetype = 2)

# g <- g + geom_text(data = samples_sum,
#                    aes_string(label = "phyto",x = -2.5, col = "associated"),
#                    hjust = "outward", size = 4,
#                    show.legend = FALSE)
#g <- g + theme_forest()

h <- g + scale_x_continuous(breaks = -2:7, 
                            labels = -2:7) + 
  xlab("Posterior distributions for effect of\nhorticultural imports on new detections")

col_spp_by_surv_str <- left_join(data.frame(phyto = gsub("P\\.", "Phytophthora",
                                                         levels(species_sensitivity$phyto))[2:71]),
                                 included_traits[,c("phyto", "CH_OO")])
col_spp_by_surv_str$col_spp <- "grey40"
col_spp_by_surv_str$col_spp[col_spp_by_surv_str$CH_OO >1] <- "black"

h <- h + theme(axis.line.x.bottom = element_line(),
               panel.grid.major.x = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.ticks.y = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(face = c("bold", rep("italic", 70) ),
                                          colour = c("black", col_spp_by_surv_str$col_spp),
                                          size = 7),
               axis.title.x = element_text(size = 8),
               axis.text.x = element_text(size = 10),
               legend.text = element_text(size = 14),
               legend.key = element_rect(size = 6, fill = "white", colour = "white")) 
trade_sens <- h + coord_cartesian(xlim = c(-2,7))

trade_sens


# 2) sensitivity to surveillance

b_draws <- model_av %>% as_draws() %>% spread_draws(b_EPPO_reporting_service) %>%
  rename(.value = b_EPPO_reporting_service) %>%
  mutate(phyto = "Average effect") %>%
  select(.value, .chain, .iteration, .draw, phyto) %>%
  mutate(risk_factor = "EPPO_reporting_service")


r_draws <- model_av %>% as_draws() %>%
  spread_draws(b_EPPO_reporting_service, r_phyto[phyto,risk_factor]) %>% 
  filter(risk_factor == "EPPO_reporting_service" | is.na(risk_factor)) %>%
  mutate(.value = b_EPPO_reporting_service + r_phyto) %>%
  select(risk_factor, .chain, .iteration, .draw, phyto, .value)

species_sensitivity <- as_tibble(bind_rows(b_draws, r_draws))

species_sensitivity$phyto <- gsub("\\.", " ", species_sensitivity$phyto)


species_sensitivity$phyto <- gsub("Phytophthora", "P.", species_sensitivity$phyto)

species_sensitivity$phyto <- factor(species_sensitivity$phyto, 
                                    levels = c("Average effect",
                                               grep("Average effect", 
                                                    names(sort(tapply(species_sensitivity$.value, 
                                                                      species_sensitivity$phyto, 
                                                                      FUN = "median"),
                                                               decreasing = TRUE)),
                                                    invert = TRUE, value = TRUE)))

g <- ggplot() + stat_halfeye(data = species_sensitivity,
                             aes(x = .value, y = phyto, fill = stat(x>0), col = stat(x>0)), 
                             #scale = scale,
                             alpha = 0.6,
                             size = 1) + 
  scale_fill_manual("", values = c("steelblue", "salmon"), labels = c(expression(~~~beta<0), expression(~~~beta>0)))+
  scale_colour_manual("", values = c("steelblue", "salmon"), labels = c(expression(~~~beta<0), expression(~~~beta>0)))
g <- g + geom_vline(xintercept = 0, linetype = 1) +
  geom_vline(xintercept = median(species_sensitivity$.value[species_sensitivity$phyto == "Average effect"]),
             linetype = 2)

h <- g + scale_x_continuous(breaks = c(-6, -4, -2, 0, 2,4,6,8), 
                            labels = c(-6, -4, -2, 0, 2,4,6,8)) + 
  xlab("Posterior distributions for effect of\nEPPO reporting activity on new detections")

col_spp_by_surv_str <- left_join(data.frame(phyto = gsub("P\\.", "Phytophthora",
                                                         levels(species_sensitivity$phyto))[2:71]),
                                 included_traits[,c("phyto", "CH_OO")])
col_spp_by_surv_str$col_spp <- "grey40"
col_spp_by_surv_str$col_spp[col_spp_by_surv_str$CH_OO == 2] <- "black"

col_spp_by_temp <- left_join(data.frame(phyto = gsub("P\\.", "Phytophthora",
                                                         levels(species_sensitivity$phyto))[2:71]),
                                 included_traits[,c("phyto", "temp_range")])
col_spp_by_temp$col_spp <- ifelse(col_spp_by_temp$temp_range > 25, 
                                  yes = "black",
                                  no = "grey40")


h <- h + theme(axis.line.x.bottom = element_line(),
               panel.grid.major.x = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.ticks.y = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(face = c("bold", rep("italic", 70) ),
                                          vjust = 0,
                                          colour = c("black", col_spp_by_surv_str$col_spp),
                                          size = 7),
               axis.title.x = element_text(size = 8),
               axis.text.x = element_text(size = 10),
               legend.text = element_text(size = 14),
               legend.key = element_rect(size = 6, fill = "white", colour = "white")) 
surv_sens <- h +coord_cartesian(xlim = c(-6,8))
surv_sens


# 3) sensitivity to climate matching

b_draws <- model_av %>% as_draws() %>% spread_draws(b_climate_match) %>%
  rename(.value = b_climate_match) %>%
  mutate(phyto = "Average effect") %>%
  select(.value, .chain, .iteration, .draw, phyto) %>%
  mutate(risk_factor = "climate_match")


r_draws <- model_av %>% as_draws() %>%
  spread_draws(b_climate_match, r_phyto[phyto,risk_factor]) %>% 
  filter(risk_factor == "climate_match" | is.na(risk_factor)) %>%
  mutate(.value = b_climate_match + r_phyto) %>%
  select(risk_factor, .chain, .iteration, .draw, phyto, .value)



species_sensitivity <- as_tibble(bind_rows(b_draws, r_draws))

species_sensitivity$phyto <- gsub("Phytophthora", "P.", 
                                  gsub("\\.", " ", species_sensitivity$phyto))

species_sensitivity$phyto <- factor(species_sensitivity$phyto, 
                                    levels = c("Average effect",
                                               grep("Average effect", 
                                                    names(sort(tapply(species_sensitivity$.value, 
                                                                      species_sensitivity$phyto, 
                                                                      FUN = "median"),
                                                               decreasing = TRUE)),
                                                    invert = TRUE, value = TRUE)))


g <- ggplot() + stat_halfeye(data = species_sensitivity,
                             aes(x = .value, y = phyto, fill = stat(x>0), col = stat(x>0)), 
                             #scale = scale,
                             alpha = 0.6,
                             size = 1) + 
  scale_fill_manual("", values = c("steelblue", "salmon"), labels = c(expression(~~~beta<0), expression(~~~beta>0)))+
  scale_colour_manual("", values = c("steelblue", "salmon"), labels = c(expression(~~~beta<0), expression(~~~beta>0)))
g <- g + geom_vline(xintercept = 0, linetype = 1) +
  geom_vline(xintercept = median(species_sensitivity$.value[species_sensitivity$phyto == "Average effect"]),
             linetype = 2)

h <- g + scale_x_continuous(breaks = c(-6, -4, -2, 0, 2,4,6,8), 
                            labels = c(-6, -4, -2, 0, 2,4,6,8)) + 
  xlab("Posterior distributions for effect of\nclimate matching on new detections")

# col_spp_by_surv_str <- left_join(data.frame(phyto = gsub("P\\.", "Phytophthora",
#                                                          levels(species_sensitivity$phyto))[2:71]),
#                                  included_traits[,c("phyto", "CH_OO")])
# col_spp_by_surv_str$col_spp <- "grey40"
# col_spp_by_surv_str$col_spp[col_spp_by_surv_str$CH_OO == 2] <- "black"

h <- h + theme(axis.line.x.bottom = element_line(),
               panel.grid.major.x = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.ticks.y = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(face = c("bold", rep("italic", 70) ),
                                          vjust = 0#,
                                          #colour = c("black", col_spp_by_surv_str$col_spp)),
               ),
               axis.title.x = element_text(size = 8),
               axis.text.x = element_text(size = 10),
               legend.text = element_text(size = 14),
               legend.key = element_rect(size = 6)) 
climate_sens <- h +coord_cartesian(xlim = c(-6,8))
climate_sens

# # highest and lowest sensitivity to surveillance
# quantile(1.1^r_draws[r_draws$phyto == "Phytophthora.pseudocryptogea",]$.value, c(0.025, 0.5, 0.975))
# quantile(1.1^r_draws[r_draws$phyto == "Phytophthora.uniformis",]$.value, c(0.025, 0.5, 0.975))


# # plot differences in species-level intercept (overall prevalnce of new detections)


r_draws <- model_av %>% as_draws() %>%
  spread_draws(b_Intercept, r_phyto[phyto,Intercept]) %>%
  filter(Intercept == "Intercept")%>%
  mutate(.value = b_Intercept + r_phyto) %>%
  rename(.variable = Intercept) %>%
  select(!c(b_Intercept, r_phyto))
  #select(risk_factor, .chain, .iteration, .draw, phyto, .value)

b_draws <-model_av %>% as_draws() %>% spread_draws(b_Intercept) %>%
  mutate(phyto = "Average effect") 

species_prevalence <- bind_rows(b_draws, r_draws)

species_prevalence$phyto <- gsub("Phytophthora", "P.",
                                  gsub("\\.", " ", species_prevalence$phyto))

species_prevalence$phyto <- factor(species_prevalence$phyto,
                                    levels = c("Average effect",
                                               grep("Average effect",
                                                    names(sort(tapply(species_prevalence$.value,
                                                                      species_prevalence$phyto,
                                                                      FUN = "median"),
                                                               decreasing = TRUE)),
                                                    invert = TRUE, value = TRUE)))


g <- ggplot() + stat_halfeye(data = species_prevalence,
                             aes(x = .value, y = phyto, fill = stat(x>0), col = stat(x>0)),
                             #scale = scale,
                             alpha = 0.6,
                             size = 1) +
  scale_fill_manual("", values = c("steelblue", "salmon"), labels = c(expression(~~~beta<0), expression(~~~beta>0)))+
  scale_colour_manual("", values = c("steelblue", "salmon"), labels = c(expression(~~~beta<0), expression(~~~beta>0)))
g <- g + geom_vline(xintercept = 0, linetype = 1) +
  geom_vline(xintercept = median(species_prevalence$b_Intercept[species_prevalence$phyto == "Average effect"]),
             linetype = 2) +

# g <- g + geom_text(data = samples_sum,
#                    aes_string(label = "phyto",x = -2.5, col = "associated"),
#                    hjust = "outward", size = 4,
#                    show.legend = FALSE)
#g <- g + theme_forest()

# h <- g + scale_x_continuous(breaks = -2:7,
#                             limits = c(-2,7),
#                             labels = -2:7) +
  xlab("Posterior distributions for effect of\nhorticultural imports on new detections")


# col_spp_by_surv_str <- left_join(data.frame(phyto = gsub("P\\.", "Phytophthora", levels(species_prevalence$phyto))[2:71]), included_traits)
# col_spp_by_surv_str$col_spp <- ifelse(col_spp_by_surv_str$CH_OO >0, yes = "red4",
#                                       no = "black")
# col_spp_by_surv_str$col_spp[col_spp_by_surv_str$CH_OO >0 & col_spp_by_surv_str$CH_OO <1 ] <- "orange"


h <- g + theme(axis.line.x.bottom = element_line(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.ticks.y = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(face = c("bold", rep("italic", 70) )),
               axis.title.x = element_text(size = 8),
               axis.text.x = element_text(size = 10),
               legend.text = element_text(size = 14),
               legend.key = element_rect(size = 6, fill = "white"))
new_det_prevalence <- h
new_det_prevalence
ggsave(new_det_prevalence,
       filename = "species_new_detections_intercepts.png", 
       height = 22, width = 9,
       dpi = 500,
       inits = "cm")


#library(ggpubr)
ggarrange(trade_sens, 
          surv_sens, 
          ncol = 2, 
          common.legend = TRUE, legend = "top",
          labels = "auto")
ggsave(paste0(minRE_spp, "/species_sensitivity_risk_factors.png"),
       height = 24, width = 17.5, units = "cm",
       dpi = 500, bg = "white")



############## plot significant interactions between traits and risk factors ##########
# use the brms defaults to choose the predicted values
trade_int_plot <- conditional_effects(best_model, effects = "trade_hort:CH_OO",
                                      int_conditions = list(CH_OO = setNames(sort(unique(arrivals$CH_OO), decreasing = TRUE), 
                                                                             c("2", "1", "0"))),
                                      #conditions = cond,
                                      type = "fitted",
                                      re_formula = NA)


newdata <- expand.grid(trade_int_plot$`trade_hort:CH_OO`$trade_hort,
                       unique(sort(arrivals$CH_OO))) %>%
  rename(trade_hort = Var1, CH_OO = Var2) %>%
  mutate(importer_iso3 = "any",
         EPPO_reporting_service = mean(arrivals$EPPO_reporting_service), 
         phyto = "unknown_species",
         phylo = "unknown_species",
         temp_range = mean(arrivals$temp_range),
         climate_match = mean(arrivals$climate_match),
         date = mean(arrivals$date))
# which models contain the interation term?
which(!is.na(model_info[,"b_CH_OO:trade_hort","est"]))

# 
plot_trade_int <- pp_average(#top_models[[1]],
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
  re_formula  = NA,
  summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  mutate(trade_hort = newdata$trade_hort,
         CH_OO = newdata$CH_OO,
         est = `50%`,
         lower = `2.5%`,
         upper = `97.5%`) %>%
  mutate(trade_hort_raw = exp(trade_hort * trade_scale + trade_center),
         CH_OO_raw = factor(newdata$CH_OO, labels = c(0,1,2)))

plot_int_with_trade <- 
ggplot() + geom_ribbon(data = plot_trade_int,
                       aes(x = trade_hort_raw, ymin = lower, ymax = upper, 
                           fill = CH_OO_raw), alpha = 0.4) +
  geom_line(data = plot_trade_int,
                            aes(x = trade_hort_raw, y = est, col = CH_OO_raw),
            size = 1.2) +
  
  scale_fill_manual("survival structures",
                      values = alpha(c("0" = "steelblue",
                                       "1" = "orange",
                                       "2" = "red4"), 0.5)) +
  scale_color_manual("survival structures",
                     values = c("0" = "steelblue",
                                "1" = "orange",
                                "2" = "red4")) + 
  scale_x_continuous("1 + horticultural imports from source regions\n(1000 kg)", 
                                                                    trans = "log2",
                                                                    breaks = c(1, 100, 10000, 1000000),
                                                                    labels = comma) +
  ylim(c(0, 0.6)) + theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        axis.title = element_text(size = 16),
                        axis.text = element_text(size = 14),#
                        legend.position = c(0.35,0.75),
                        legend.key.size = unit(0.75, 'cm'),
                        legend.title = element_text(size = 14),
                        legend.text = element_text(size = 14)) + ylab("") +
  geom_vline(xintercept = quantile(arrivals_unscaled$trade_hort, prob = 0.25), linetype = 2) +
  geom_vline(xintercept = quantile(arrivals_unscaled$trade_hort, prob = 0.75), linetype = 2)


plot_int_with_trade


# best_model <- readRDS(paste0(minRE_spp, "/", com_type, "/", best_models$model[1]))
# 
# trade_int_plot <- conditional_effects(best_model, effects = "trade_hort:CH_OO",
#                                       int_conditions = list(CH_OO = setNames(sort(unique(arrivals$CH_OO), decreasing = TRUE), 
#                                                                              c("2", "1", "0"))),
#                                    #conditions = cond,
#                                    type = "fitted",
#                                    re_formula = NA)
# 
# 
# trade_plot <- plot(trade_int_plot, plot = FALSE)[[1]] 
# 
# 
# 
# # to reverse scaling and centering:
# trade_scale <- 2 * sd(log(1 + arrivals_unscaled$trade_hort), na.rm = TRUE)
# trade_center <- mean(log(1 + arrivals_unscaled$trade_hort), na.rm = TRUE)
# # unscale the axes and labels for plotting
# trade_plot$data$trade_hort <- exp(trade_plot$data$trade_hort * trade_scale + trade_center)
# trade_plot$data$effect1__ <- exp(trade_plot$data$effect1__ * trade_scale + trade_center)
# library(scales) # avoid scientific notation of numbers
# trade_plot <- trade_plot + scale_x_continuous("1 + horticultural imports from source regions\n(1000 kg)", 
#                                 trans = "log2",
#                                 breaks = c(1, 100, 10000, 1000000),
#                                 labels = comma) +
#   scale_y_continuous("") + ylab("")
# # trade_plot <- trade_plot + scale_fill_colorblind("survival structures") +
# #   scale_color_colorblind("survival structures")
# trade_plot <- trade_plot + scale_fill_manual("survival structures",
#                                             values = alpha(c("0" = "steelblue",
#                                               "1" = "orange",
#                                               "2" = "red4"), 0.5)) +
#  scale_color_manual("survival structures",
#                       values = c("0" = "steelblue",
#                                    "1" = "orange",
#                                    "2" = "red4"))
# trade_plot <- trade_plot + theme(panel.grid.major = element_blank(),
#                                  panel.grid.minor = element_blank(),
#                                  panel.background = element_blank(),
#                                  axis.line = element_line(colour = "black"),
#                                  axis.title = element_text(size = 16),
#                                  axis.text = element_text(size = 14),#
#                                  legend.position = c(0.35,0.75),
#                                  legend.key.size = unit(0.75, 'cm'),
#                                  legend.title = element_text(size = 14),
#                                  legend.text = element_text(size = 14))
# trade_CH_OO_interaction <- trade_plot + ylim(c(0,1)) #+ geom_point(data = arrivals_unscaled, # add observed data if needed
#                         #aes(x = 1+trade_hort, y = arrival),
#                         #inherit.aes = FALSE)
# trade_CH_OO_interaction
ggsave(filename = "trade_CHOO_interaction.png",
       height = 6, width = 5,
       dpi = 600,
       units = "in")

########### plot main effect of thermal tolerance ########
therm_tol_effect <- conditional_effects(best_model, effects = "temp_range", type = "fitted")
#therm_tol_plot <- plot(therm_tol_effect, plot = FALSE, line_args = c("colour" = "grey30")  )[[1]]

newdata <- data.frame(importer_iso3 = "any",
                      trade_hort = mean(arrivals$trade_hort), 
                      phyto = "unknown_species",
                      phylo = "unknown_species",
                      CH_OO = mean(arrivals$CH_OO),
                      temp_range = therm_tol_effect$temp_range$temp_range,
                      climate_match = mean(arrivals$climate_match),
                      EPPO_reporting_service = mean(arrivals$EPPO_reporting_service),
                      date = mean(arrivals$date),
                      stringsAsFactors = FALSE)
which(!is.na(model_info[,"b_temp_range","est"])) 
plot_thermal_tol <- pp_average(
  top_models[[1]],
  # top_models[[2]], 
  # top_models[[3]],
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
  re_formula  = NA,
  summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  mutate(temp_range = newdata$temp_range,
         est = `50%`,
         lower = `2.5%`,
         upper = `97.5%`) %>%
  mutate(temp_range_raw = temp_range * temp_range_scale + temp_range_center)

plot_temp_range_main <- 
  ggplot() + geom_ribbon(data = plot_thermal_tol,
                         aes(x = temp_range_raw, ymin = lower, ymax = upper), alpha = 0.4) +
  geom_line(data = plot_thermal_tol,
            aes(x = temp_range_raw, y = est),
            size = 1.2) +
  ylim(c(0, 0.6)) + theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        axis.title = element_text(size = 16),
                        axis.text = element_text(size = 14),#
                        legend.position = c(0.35,0.75),
                        legend.key.size = unit(0.75, 'cm'),
                        legend.title = element_text(size = 14),
                        legend.text = element_text(size = 14)) +
  xlab(expression(atop("Thermal tolerance range",
                       paste((T[max] - T[min]))))) + ylim(c(0,0.6)) +
  ylab("Probability of arrival") +
  geom_vline(xintercept = quantile(arrivals_unscaled$temp_range, prob = 0.25), linetype = 2) +
  geom_vline(xintercept = quantile(arrivals_unscaled$temp_range, prob = 0.75), linetype = 2)


plot_temp_range_main



# temp_range_scale <- 2* sd(arrivals_unscaled$temp_range)
# temp_range_center <- mean(arrivals_unscaled$temp_range)
# 
# therm_tol_plot$data$temp_range <- therm_tol_plot$data$temp_range * temp_range_scale + temp_range_center
# therm_tol_plot$data$effect1__ <- therm_tol_plot$data$effect1__ * temp_range_scale + temp_range_center
# 
# therm_tol_plot <- therm_tol_plot + theme(panel.grid.major = element_blank(),
#                          panel.grid.minor = element_blank(),
#                          panel.background = element_blank(),
#                          axis.line = element_line(colour = "black"),
#                          axis.title = element_text(size = 16),
#                          axis.text = element_text(size = 14),#
#                          legend.position = c(0.35,0.75),
#                          legend.key.size = unit(0.75, 'cm'),
#                          legend.title = element_text(size = 14),
#                          legend.text = element_text(size = 14)) +
#   xlab(expression(Thermal~tolerance~range~(T[max] - T[min]))) + ylim(c(0,1)) +
#   ylab("Probability of arrival") 

#climate and temp_range
# climate_int_plot <- conditional_effects(fit0129, effects = "climate_match:temp_range",
#                                         int_conditions = list(temp_range = setNames(quantile(fit0129$data$temp_range, 
#                                                                                              c(0.10, 0.5, 0.90), na.rm = TRUE), 
#                                                                                     c("narrow", "median", "broad"))),
#                                       #conditions = cond,
#                                       type = "fitted",
#                                       re_formula = NA)
# 
# 
# climate_plot <- plot(climate_int_plot, plot = FALSE)[[1]] 
# 
# 
# 
# # to reverse scaling and centering:
# climate_scale <- 2 * sd(-sqrt(arrivals_unscaled$climate_distance), na.rm = TRUE)
# climate_center <- mean(-sqrt(arrivals_unscaled$climate_distance), na.rm = TRUE)
# # unscale the axes and labels for plotting
# climate_plot$data$climate_match <- climate_plot$data$climate_match * climate_scale + climate_center
# climate_plot$data$effect1__ <- climate_plot$data$effect1__ * climate_scale + climate_center
# library(scales) # avoid scientific notation of numbers
# 
# climate_plot <- climate_plot + scale_x_continuous(expression(Climate~matching~(-sqrt(Mahalanobis~distance)))) #+
#   #scale_y_continuous("")
# climate_plot <- climate_plot + scale_fill_discrete("thermal tolerance range") +
#   scale_color_discrete("thermal tolerance range")
# # climate_plot <- climate_plot + scale_fill_manual("thermal tolerance range" ,
# #                                                  values = c("broad" = "salmon",
# #                                                             "median" = "darkgreen",
# #                                                             "narrow" = "blue")) +
# #   scale_color_manual( "thermal tolerance range" ,
# #                      values = c("broad" = "salmon",
# #                                 "median" = "darkgreen",
# #                                 "narrow" = "blue"))
# climate_plot <- climate_plot + theme(panel.grid.major = element_blank(),
#                                      panel.grid.minor = element_blank(),
#                                      panel.background = element_blank(),
#                                      axis.line = element_line(colour = "black"),
#                                      axis.title.y = element_blank(),
#                                      axis.title.x = element_text(size = 16),
#                                      axis.text = element_text(size = 14),#
#                                      legend.position = c(0.4,0.75),
#                                      legend.key.size = unit(0.75, 'cm'),
#                                      legend.title = element_text(size = 14),
#                                      legend.text = element_text(size = 14))
# climate_temp_interaction <- climate_plot + geom_point(data = arrivals_unscaled,
#                                              aes(x = -sqrt(climate_distance), y = arrival),
#                                              inherit.aes = FALSE)
# climate_temp_interaction
# ggsave(filename = "climate_temp_interaction.png",
#        height = 6, width = 5,
#        dpi = 600,
#        units = "in")

# EPPO reporting service CH_OO
surv_int_plot <- conditional_effects(best_model, effects = "EPPO_reporting_service:CH_OO",
                                     int_conditions = list(CH_OO = setNames(sort(unique(arrivals$CH_OO), decreasing = TRUE), 
                                                                            c("2", "1", "0"))),
                                      #conditions = cond,
                                      type = "fitted",
                                      re_formula = NA)

newdata <- expand.grid(surv_int_plot$`EPPO_reporting_service:CH_OO`$EPPO_reporting_service,
                       unique(sort(arrivals$CH_OO))) %>%
  rename(EPPO_reporting_service = Var1, CH_OO = Var2) %>%
  mutate(importer_iso3 = "any",
         trade_hort = mean(arrivals$trade_hort),
         #EPPO_reporting_service = mean(arrivals$EPPO_reporting_service), 
         phyto = "unknown_species",
         phylo = "unknown_species",
         temp_range = mean(arrivals$temp_range),
         climate_match = mean(arrivals$climate_match),
         date = mean(arrivals$date))
# which models contain the interation term?
which(!is.na(model_info[,"b_CH_OO:EPPO_reporting_service","est"]))

# 
plot_surv_int <- pp_average(
  top_models[[1]],
  top_models[[2]], 
  top_models[[3]],
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
  re_formula  = NA,
  summary = FALSE) %>% as_draws() %>% 
  summarise_draws(~quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  mutate(EPPO_reporting_service = newdata$EPPO_reporting_service,
         CH_OO = newdata$CH_OO,
         est = `50%`,
         lower = `2.5%`,
         upper = `97.5%`) %>%
  mutate(EPPO_reporting_service_raw = exp(EPPO_reporting_service * surv_scale + surv_center),
         CH_OO_raw = factor(newdata$CH_OO, labels = c(0,1,2)))

plot_int_with_surv <- 
  ggplot() + geom_ribbon(data = plot_surv_int,
                         aes(x = EPPO_reporting_service_raw, ymin = lower, ymax = upper, 
                             fill = CH_OO_raw), alpha = 0.4) +
  geom_line(data = plot_surv_int,
            aes(x = EPPO_reporting_service_raw, y = est, col = CH_OO_raw),
            size = 1.2) +
  
  scale_fill_manual("survival structures",
                    values = alpha(c("0" = "steelblue",
                                     "1" = "orange",
                                     "2" = "red4"), 0.5)) +
  scale_color_manual("survival structures",
                     values = c("0" = "steelblue",
                                "1" = "orange",
                                "2" = "red4")) + 
  scale_x_continuous("Importer plant health surveillance effort\n(EPPO Reporting Service entries 1967-2021)",
                     trans = "log2",
                     breaks = c(1, 5, 10, 50,  100, 500, 1000),
                     labels = comma) +
  ylim(c(0, 0.6)) + theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black"),
                        axis.title = element_text(size = 16),
                        axis.text = element_text(size = 14),#
                        legend.position = c(0.35,0.75),
                        legend.key.size = unit(0.75, 'cm'),
                        legend.title = element_text(size = 14),
                        legend.text = element_text(size = 14)) + ylab("") +
  geom_vline(xintercept = quantile(arrivals_unscaled$EPPO_reporting_service, prob = 0.25), linetype = 2) +
  geom_vline(xintercept = quantile(arrivals_unscaled$EPPO_reporting_service, prob = 0.75), linetype = 2)


plot_int_with_surv
#surv_plot <- plot(surv_int_plot, plot = FALSE)[[1]] 



# to reverse scaling and centering:
# surv_scale <- 2 * sd(log(arrivals_unscaled$EPPO_reporting_service), na.rm = TRUE)
# surv_center <- mean(log(arrivals_unscaled$EPPO_reporting_service), na.rm = TRUE)
# # unscale the axes and labels for plotting
# surv_plot$data$EPPO_reporting_service <- exp(surv_plot$data$EPPO_reporting_service * surv_scale + surv_center)
# surv_plot$data$effect1__ <- exp(surv_plot$data$effect1__ * surv_scale + surv_center)
# library(scales) # avoid scientific notation of numbers
# surv_plot <- surv_plot + scale_x_continuous("Importer plant health surveillance effort\n(EPPO Reporting Service entries 1967-2021)", 
#                                               trans = "log2",
#                                               breaks = c(1, 5, 10, 50,  100, 500, 1000),
#                                               labels = comma) +
#   scale_y_continuous("")
# surv_plot <- surv_plot + scale_fill_manual("survival structures",
#                                              values = alpha(c("0" = "steelblue",
#                                                               "1" = "orange",
#                                                               "2" = "red4"), 0.5)) +
#   scale_color_manual("survival structures",
#                      values = c("0" = "steelblue",
#                                 "1" = "orange",
#                                 "2" = "red4"))
# surv_plot <- surv_plot + theme(panel.grid.major = element_blank(),
#                                panel.grid.minor = element_blank(),
#                                panel.background = element_blank(),
#                                axis.line = element_line(colour = "black"),
#                                axis.title.y = element_blank(),
#                                axis.title.x = element_text(size = 16),
#                                axis.text = element_text(size = 14),#
#                                legend.position = c(0.35,0.75),
#                                legend.key.size = unit(0.75, 'cm'),
#                                legend.title = element_text(size = 14),
#                                legend.text = element_text(size = 14))
# # show differences in number of arrivals
# #sum_arrivals <- table(arrivals_unscaled$EPPO_reporting_service, arrivals_unscaled$arrival)
# #arrival_pressure <- data.frame(EPPO_reporting_service = rownames(sum_arrivals),
# #                               nonarrival = sum_arrivals[,1],
# #                               arrival = sum_arrivals[,2])
# 
# surv_CH_OO_interaction <- surv_plot + ylim(c(0,1))#+ geom_point(data = arrivals_unscaled, # add the observed data points, if needed
#                                     #               aes(x = EPPO_reporting_service, y = arrival),
#                                      #              inherit.aes = FALSE)
# 
# surv_CH_OO_interaction

ggsave(filename = "surveillance_CHOO_interaction.tiff",
       height = 6, width = 5,
       dpi = 600,
       units = "in")


ggarrange(plot_temp_range_main,
          plot_int_with_trade,
          plot_int_with_surv,
          labels = "auto",
          nrow = 1,
          align = "v",
          font.label = list(size = 18), common.legend = TRUE,
          legend = "right")
ggsave(paste0(minRE_spp, "/interaction_plot.png"),
       height = 6, width = 18, units = "in",
       dpi = 500)

