rm(list = ls())


# load packages



library(brms)
library(rslurm)

# load processed nursery data

arrivals <- readRDS("arrivals.rds")
dredge_arrival <- readRDS("dredge_arrival.rds")
Phytophthora_covariance <- readRDS("PhytophthoraA.rds")

# check for missing models
missing <- setdiff(1:nrow(dredge_arrival), as.integer(gsub("fit|\\.rds", "", list.files(dir(pattern = "_rslurm_"), patter = "fit"))))
dredge_arrival <- dredge_arrival[missing,]
# function to fit model
slurm_brms <- function(model_formula, index){
  if(index %in% "fit0001"){
    prior <-  
      prior(student_t(3, 0, 20), "sd")
  }
  else {
    prior = c(
      prior(normal(0, 10), "b"),
      prior(normal(0, 50), "Intercept"),
      prior(student_t(3, 0, 20), "sd"))
  }
  
  fit <-
    brm(bf(model_formula),
        data = arrivals, family = zero_inflated_binomial() ,
        data2 = list(Phytophthora_covariance = Phytophthora_covariance),
        prior = prior,
        control = list(adapt_delta = 0.999,
                       max_treedepth = 15),
        iter = 4000,
        chains = 3
    )
  fit <- add_criterion(fit, criterion = c("loo", "loo_R2", "bayes_R2"))
  
  saveRDS(fit, file = paste0(index, ".rds"))
}

sjob <- slurm_apply(f = slurm_brms,
                    params = dredge_arrival,
                    jobname = 'trade_hort_minREspp0',
                    nodes = nrow(dredge_arrival),
                    cpus_per_node = 1,
                    submit = TRUE,
                    global_objects = c('Phytophthora_covariance', 'arrivals'),
                    slurm_options = list(time = '23:59:00',
                                         #mem = 20000,
                                         partition = 'short-serial',
                                         error = '%a.err'))