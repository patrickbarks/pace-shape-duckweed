

##### 1. Required libraries
library(readr)
library(dplyr)
library(tibble)
library(broom)
library(tidyr)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ape)
library(geosphere)
library(raster)
library(mgcv)
library(loo)




##### 2. Load RStan and set relevent options

# In the R code below, sections that require the RStan library are commented
#  out. Our analyses can be replicated without installing RStan, by using model
#  output that we have already saved to .csv files. If users wish to re-run the
#  RStan models, uncomment the relevent sections (and ensure that RStan is
#  installed correctly).

# For RStan installation instructions, see:
# (Windows) https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Windows
# (Mac/Linux) https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Mac-or-Linux

# # load RStan and set relevent options
# library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())




##### 3. Ensure working directory is set to '/pace-shape-duckweed/'
# setwd('~/.../pace-shape-duckweed/')




##### 4. Read and organize data files
# common garden demographic data
dat <- read_csv('dat/data_demographic.csv') %>% 
  filter(!site %in% c('dbn', 'han', 'skf')) %>%   # omitted sites
  filter(discard != TRUE)                # omit fronds that could not be tracked

# frond areas from supplementary experiment with all available strains
dat_areas_supp <- read_csv('dat/data_areas_supplementary.csv') %>% 
  filter(site != 'skf')           # omit site 'skf' (L. minor)

# site-specific environmental data
dat_sites <- read_csv("dat/data_sites.csv") %>% 
  filter(site != 'skf') %>%       # omit site 'skf' (L. minor)
  mutate(nutrient_pc1 = princomp(data.frame(log(tdn), log(tdp)))$scores[,1])




##### 5. Check final sample sizes by site and site-by-block
group_by(dat, strain) %>% summarize(n()) %>% as.data.frame()
group_by(dat, strain, block) %>% summarize(n()) %>% as.data.frame()




##### 6. Assign site-specific colours for plotting
replic_col <- data.frame(
  site = dat_sites$site[dat_sites$replicated == TRUE],
  col = brewer.pal(n = 5, name = "Accent"),
  stringsAsFactors = F
)

# join replicated colours to dat
dat <- left_join(dat, replic_col, by = 'site') %>% 
  mutate(col = ifelse(is.na(col), '#ffffff', col)) # non-replic sites in white




##### 7. Extract subset containing reproduction data only
# number of offspring produced each day by each frond
dat_raw_repro <- subset(dat, select = Jun_01_2014:Aug_31_2014)




##### 8. Determine dates of first and last reproduction
# get column name (i.e. date) associated with first (min) or last (max) repro
GetDateFirstRepro <- function (x) names(dat_raw_repro)[min(which(x > 0))]
GetDateLastRepro <- function (x) names(dat_raw_repro)[max(which(x > 0))]
dat$date_first_repro <- apply(dat_raw_repro, 1, GetDateFirstRepro)
dat$date_last_repro <- apply(dat_raw_repro, 1, GetDateLastRepro)

# convert dates to R's date format
dat$date_birth <- as.Date(dat$date_birth)
dat$date_first_repro <- as.Date(dat$date_first_repro, format = "%b_%d_%Y")
dat$date_last_repro <- as.Date(dat$date_last_repro, format = "%b_%d_%Y")




##### 9. Calculate frond-level demographic statistics
# 'ffr' is from first reproduction, 'latency' is latency to reproduce
dat$lifespan <- as.numeric(dat$date_last_repro - dat$date_birth + 1)
dat$lifespan_ffr <- as.numeric(dat$date_last_repro - dat$date_first_repro + 1)
dat$latency <- as.numeric(dat$date_first_repro - dat$date_birth + 1)
dat$total_offspring <- rowSums(dat_raw_repro, na.rm = T)
dat$total_offspring[which(dat$uncertain_repro == TRUE)] <- NA
dat$fecund_mean <- dat$total_offspring / dat$lifespan
dat$fecund_mean_ffr <- dat$total_offspring / dat$lifespan_ffr




##### 10. Convert demographic data to flat (tidy) format
# create variable for daily survival and fecundity from age 1 to death
dat_tidy <- dat %>% 
  dplyr::select(-area, -discard, -tray, -col) %>% 
  gather(date, fecund, ends_with('2014')) %>%
  arrange(id) %>% 
  mutate(date = as.Date(date, format = "%b_%d_%Y")) %>% 
  filter(date >= date_birth & date <= date_last_repro) %>%
  mutate(age = as.numeric(date - date_birth + 1),
         died = ifelse(date == date_last_repro, 1, 0),
         surv = ifelse(date == date_last_repro, 0, 1)) %>% 
  dplyr::select(-contains('date'))




##### 11. Cohort-level parameters by site, strain, and strain-by-block
# lifespan_group is cohort-specific pace from birth
# lifespan_ffr_group is cohort-specific pace from age of first repro
cohort_par_site <- dat %>% 
  group_by(site) %>% 
  summarize(lifespan_group = mean(lifespan, na.rm = T),
            lifespan_ffr_group = mean(lifespan_ffr, na.rm = T))

cohort_par_strain <- dat %>% 
  group_by(strain) %>% 
  summarize(lifespan_group = mean(lifespan, na.rm = T),
            lifespan_ffr_group = mean(lifespan_ffr, na.rm = T))

cohort_par_block <- dat %>% 
  group_by(strain, block) %>% 
  summarize(lifespan_group = mean(lifespan, na.rm = T),
            lifespan_ffr_group = mean(lifespan_ffr, na.rm = T)) %>% 
  ungroup()




##### 12. Mean standardize fecundity and pace-standardize age, by site, strain, etc.
# from birth
dat_tidy_std_site <- dat_tidy %>%
  left_join(cohort_par_site, by = 'site') %>%
  mutate(age_std = age / lifespan_group,
         fecund_std = fecund / fecund_mean)

dat_tidy_std_strain <- dat_tidy %>%
  left_join(cohort_par_strain, by = 'strain') %>%
  mutate(age_std = age / lifespan_group,
         fecund_std = fecund / fecund_mean)

dat_tidy_std_block <- dat_tidy %>%
  left_join(cohort_par_block, by = c('strain', 'block')) %>%
  mutate(age_std = age / lifespan_group,
         fecund_std = fecund / fecund_mean)

# from age of first reproduction
dat_tidy_std_ffr_site <- dat_tidy %>%
  group_by(id) %>% 
  mutate(first_repro = which(fecund > 0)[1]) %>% 
  filter(age >= first_repro) %>% 
  ungroup() %>% 
  left_join(dplyr::select(cohort_par_site, site, lifespan_ffr_group),
            by = 'site') %>% 
  mutate(age_std = age / lifespan_ffr_group,
         fecund_std = fecund / fecund_mean_ffr)

dat_tidy_std_ffr_strain <- dat_tidy %>%
  group_by(id) %>% 
  mutate(first_repro = which(fecund > 0)[1]) %>% 
  filter(age >= first_repro) %>% 
  ungroup() %>% 
  left_join(dplyr::select(cohort_par_strain, strain, lifespan_ffr_group),
            by = 'strain') %>% 
  mutate(age_std = age / lifespan_ffr_group,
         fecund_std = fecund / fecund_mean_ffr)

dat_tidy_std_ffr_block <- dat_tidy %>%
  group_by(id) %>% 
  mutate(first_repro = which(fecund > 0)[1]) %>% 
  filter(age >= first_repro) %>% 
  ungroup() %>% 
  left_join(dplyr::select(cohort_par_block, strain, block, lifespan_ffr_group),
            by = c('strain', 'block')) %>% 
  mutate(age_std = age / lifespan_ffr_group,
         fecund_std = fecund / fecund_mean_ffr)




##### 13. Empirical survivorship and hazard by age and strain

# function to get proportional survivorship
# x is age, t is vector of ages at death, n is total starting sample size
Survivorship <- function (x, t, n) { length(which(t >= x)) / n }

# wrapper function to get proportional survivorship by strain
SurvivorshipWrap <- function(age_vec, t) {
  n <- length(t)                                      # starting sample size
  surv <- sapply(age_vec, Survivorship, t = t, n = n) # prop. surivorship at age
  return(data.frame(age = age_vec, surv = surv))
}

Hazard <- function(t) {
  n <- length(t)                                                 # starting sample size
  max_age <- max(t) - 1                                          # max age - 1
  n_last_repro <- tabulate(t, nbins = max_age)                   # number of deaths by age
  n_surv_cum <-  c(n, n - cumsum(n_last_repro[1:(max_age - 1)])) # cumulat. survivorship
  prob_surv <- (n_surv_cum - n_last_repro) / n_surv_cum          # probability survival
  haz <- -log(prob_surv)
  return(data.frame(age = 1:max_age, n_last_repro, n_surv_cum, prob_surv, haz))
}

# survivorship by absolute age
surv_strain <- dat %>% 
  group_by(strain) %>%
  do(SurvivorshipWrap(age_vec = 1:40, t = .$lifespan)) %>% 
  ungroup() %>% 
  mutate(site = substr(strain, 1, 3))

# survivorship by standardized age
surv_strain_std <- dat %>% 
  group_by(strain) %>%
  mutate(lifespan_std = lifespan / mean(lifespan)) %>% 
  do(SurvivorshipWrap(age_vec = seq(0, 1.8, 0.02), t = .$lifespan_std)) %>% 
  ungroup() %>% 
  mutate(site = substr(strain, 1, 3))

# hazard by absolute age
hazard_strain <- dat %>%
  group_by(strain) %>%
  do(Hazard(t = .$lifespan)) %>% 
  ungroup()




##### 14. Non-parametric hazard by age and strain
GetHazardCurve <- function(data, k) {
  
  # mean lifespan
  lifespan_group <- mean(filter(data, died == 1)$age)
  
  # fit models for absolute and standardized age
  mod <- gam(surv ~ s(age, k = k), family = 'binomial', data = data)
  mod_std <- gam(surv ~ s(age_std, k = k), family = 'binomial', data = data)
  
  # vector of ages for prediction
  newage <- seq(0, max(data$age), length = 100)
  newagestd <- seq(0, max(data$age_std), length = 100)
  
  # predict absolute hazard and SE
  pred_y <- predict(mod, newdata = data.frame(age = newage), type = 'response')
  pred_y_se <- predict(mod, newdata = data.frame(age = newage), type = 'response', se = T)$se.fit
  
  # predict standardized hazard and SE
  pred_y_std <- predict(mod_std, newdata = data.frame(age_std = newagestd), type = 'response')
  pred_y_std_se <- predict(mod_std, newdata = data.frame(age_std = newagestd), type = 'response', se = T)$se.fit
  
  # organize in df and return
  df <- data.frame(
    newage,
    newagestd,
    haz = -log(pred_y),
    haz_upper = -log(pred_y + pred_y_se),
    haz_lower = -log(pred_y - pred_y_se),
    haz_std = (-log(pred_y_std)) * lifespan_group,
    haz_std_upper = (-log(pred_y_std + pred_y_std_se)) * lifespan_group,
    haz_std_lower = (-log(pred_y_std - pred_y_std_se)) * lifespan_group
  )
  
  return(df)
}

hazard_spline <- dat_tidy_std_strain %>%
  group_by(strain) %>%
  do(GetHazardCurve(., k = 3)) %>% 
  ungroup()




##### 15. Parametric mortality models

# read model data already written to file, obtained via code below
#  (commented-out code below takes ~10 minutes to run, and requires
#  installation of RStan)
mortality_param <- read_csv('analysis/mortality_parametric_strain.csv')

# # function to fit parametric mortality models
# ParametricMortality <- function(t) {
# 
#   # prepare data
#   t_new <- seq(0, max(t), length = 100) # sequence of ages for prediction
#   N <- length(t)
#   N_pred <- length(t_new)
# 
#   # stan parameter init functions
#   inits_expo <- function() { list(l = runif(1, 0, 1))}
#   inits_weib <- function() { list(l = runif(1, 0, 1), b = runif(1, 0, 1))}
#   inits_gomp <- function() { list(l = runif(1, 0, 0.1), g = runif(1, 0, 1))}
#   inits_logi <- function() { list(l = runif(1, 0, 1), g = runif(1, 0, 1), s = runif(1, 0, 1))}
# 
#   # fit stan models
#   FitStan <- function(model, init) {
#     fit <- sampling(
#       model,
#       data = list(t = t, N = N, N_pred = N_pred, t_new = t_new),
#       init = init,
#       warmup = 1500,
#       iter = 3000,
#       thin = 2,
#       chains = 2,
#       control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 12)
#     )
#   }
# 
#   fit_expo <- FitStan(mod_expo, inits_expo)
#   fit_weib <- FitStan(mod_weib, inits_weib)
#   fit_gomp <- FitStan(mod_gomp, inits_gomp)
#   fit_logi <- FitStan(mod_logi, inits_logi)
# 
#   # check for divergent transitions
#   Diverg <- function(fit) {
#     args <- get_sampler_params(fit, inc_warmup = F)
#     diverg <- do.call(rbind, args)[,5]
#     return(n_diverg = length(which(diverg > 0)))
#   }
# 
#   diverg_expo <- Diverg(fit_expo)
#   diverg_weib <- Diverg(fit_weib)
#   diverg_gomp <- Diverg(fit_gomp)
#   diverg_logi <- Diverg(fit_logi)
# 
#   # rerun if any divergent transitions
#   if (any(c(diverg_expo, diverg_weib, diverg_gomp, diverg_logi)) > 0) {
#     fit_expo <- FitStan(mod_expo, inits_expo)
#     fit_weib <- FitStan(mod_weib, inits_weib)
#     fit_gomp <- FitStan(mod_gomp, inits_gomp)
#     fit_logi <- FitStan(mod_logi, inits_logi)
# 
#     diverg_expo <- Diverg(fit_expo)
#     diverg_weib <- Diverg(fit_weib)
#     diverg_gomp <- Diverg(fit_gomp)
#     diverg_logi <- Diverg(fit_logi)
#   }
# 
#   # extract posterior samples of survivorship and hazard
#   surv_expo <- rstan::extract(fit_expo, pars = "surv_pred")$surv_pred
#   surv_weib <- rstan::extract(fit_weib, pars = "surv_pred")$surv_pred
#   surv_gomp <- rstan::extract(fit_gomp, pars = "surv_pred")$surv_pred
#   surv_logi <- rstan::extract(fit_logi, pars = "surv_pred")$surv_pred
# 
#   haz_expo <- rstan::extract(fit_expo, pars = "haz_pred")$haz_pred
#   haz_weib <- rstan::extract(fit_weib, pars = "haz_pred")$haz_pred
#   haz_gomp <- rstan::extract(fit_gomp, pars = "haz_pred")$haz_pred
#   haz_logi <- rstan::extract(fit_logi, pars = "haz_pred")$haz_pred
# 
#   # calculate waic
#   waic_expo <- waic(extract_log_lik(fit_expo))$waic
#   waic_weib <- waic(extract_log_lik(fit_weib))$waic
#   waic_gomp <- waic(extract_log_lik(fit_gomp))$waic
#   waic_logi <- waic(extract_log_lik(fit_logi))$waic
# 
#   # summary table for each model type
#   df_summary <- data.frame(
#     model = c('Exponential', 'Weibull', 'Gompertz', 'Logistic'),
#     waic = c(waic_expo, waic_weib, waic_gomp, waic_logi),
#     diverg = c(diverg_expo, diverg_weib, diverg_gomp, diverg_logi),
#     stringsAsFactors = F
#   )
# 
#   # get posterior median and 95% CI of survivorship and hazard
#   GetQuantiles <- function(model, surv, haz, t) {
#     surv_low <- apply(surv, 2, function(x) quantile(x, probs = 0.025))
#     surv_med <- apply(surv, 2, function(x) quantile(x, probs = 0.500))
#     surv_upp <- apply(surv, 2, function(x) quantile(x, probs = 0.975))
# 
#     haz_low <- apply(haz, 2, function(x) quantile(x, probs = 0.025))
#     haz_med <- apply(haz, 2, function(x) quantile(x, probs = 0.500))
#     haz_upp <- apply(haz, 2, function(x) quantile(x, probs = 0.975))
# 
#     cbind.data.frame(surv_low, surv_med, surv_upp, haz_low, haz_med, haz_upp) %>%
#       mutate(t = t, model = model)
#   }
# 
#   posterior_quantiles <- rbind.data.frame(
#     GetQuantiles('Exponential', surv_expo, haz_expo, t_new),
#     GetQuantiles('Weibull', surv_weib, haz_weib, t_new),
#     GetQuantiles('Gompertz', surv_gomp, haz_gomp, t_new),
#     GetQuantiles('Logistic', surv_logi, haz_logi, t_new)
#   )
# 
#   # join posterior_quantiles to df_summary and return
#   posterior_quantiles %>%
#     left_join(df_summary, by = 'model') %>%
#     mutate(model = factor(model, levels = c('Exponential',
#                                             'Weibull',
#                                             'Gompertz',
#                                             'Logistic')))
# }
# 
# # pre-compile stan models
# mod_expo <- stan_model('stan/surv-exponential.stan')
# mod_weib <- stan_model('stan/surv-weibull.stan')
# mod_gomp <- stan_model('stan/surv-gompertz.stan')
# mod_logi <- stan_model('stan/surv-logistic.stan')
# 
# # apply ParametricMortality function by strain
# mortality_param <- dat %>%
#   group_by(strain) %>%
#   do(ParametricMortality(t = .$lifespan))

# check for divergent transitions (sign of poor model convergence)
mortality_param %>%
  group_by(strain, model) %>%
  summarize(diverg = unique(diverg)) %>%
  as.data.frame()

# identify best mortality model for each strain
mortality_param_best <- mortality_param %>%
  group_by(strain) %>%
  summarize(best_model = model[which.min(waic)]) %>%
  as.data.frame()

# summary table of waic values (Table S2)
mortality_param_summary <- mortality_param %>%
  dplyr::select(strain, model, waic) %>% 
  unique() %>% 
  spread(model, waic) %>%
  left_join(mortality_param_best, by = 'strain') %>% 
  subset(select = c(1, 2, 5, 3, 4, 6)) %>%
  as.data.frame()

# # write results to file
# write.csv(mortality_param, 'analysis/mortality_parametric_strain.csv', row.names = F)
# write.csv(mortality_param_summary, 'analysis/mortality_parametric_summary.csv', row.names = F)




##### 16. Bootstrap shape_mortality
BootShape <- function(data, n_rep) {
  t <- data$lifespan
  Resample <- replicate(n_rep, sample(t, length(t), replace = T), simplify = F)
  shape_mort_boot <- sapply(Resample, function(y) 1 - (sd(y) / mean(y)))
  return(data.frame(shape_mort_boot))
}

set.seed(08); shape_site <- group_by(dat, site) %>%
  do(BootShape(., n_rep = 5e4)) %>% ungroup()

set.seed(80); shape_strain <- group_by(dat, strain) %>%
  do(BootShape(., n_rep = 5e4)) %>% ungroup()

set.seed(88); shape_block <- group_by(dat, strain, block) %>%
  do(BootShape(., n_rep = 5e4)) %>% ungroup()

# add colors to shape_strain df
shape_strain <- mutate(shape_strain, site = substr(strain, 1, 3)) %>% 
  left_join(replic_col, by = 'site') %>% 
  mutate(col = ifelse(is.na(col), '#ffffff', col))




##### 17. Calculate shape_fecundity
ShapeFecundity <- function(data) {
  mod <- glm(fecund ~ age, data = data)
  mod_std <- glm(fecund_std ~ age_std, data = data)
  out <- data.frame(shape_fec_nonstd = coef(mod)[2],
                    shape_fec = coef(mod_std)[2])
  return(out)
}

fecund_slope_site <- dat_tidy_std_ffr_site %>% 
  filter(uncertain_repro == FALSE) %>% group_by(id) %>%
  do(ShapeFecundity(data = .)) %>%
  ungroup() %>% rename(shape_fec_site = shape_fec)

fecund_slope_strain <- dat_tidy_std_ffr_strain %>% 
  filter(uncertain_repro == FALSE) %>% group_by(id) %>%
  do(ShapeFecundity(data = .)) %>%
  ungroup() %>% rename(shape_fec_strain = shape_fec)

fecund_slope_block <- dat_tidy_std_ffr_block %>% 
  filter(uncertain_repro == FALSE) %>% group_by(id) %>%
  do(ShapeFecundity(data = .)) %>%
  ungroup() %>% rename(shape_fec_block = shape_fec)

# join to dat
dat_join <- dat %>% 
  dplyr::select(-(Jun_01_2014:Aug_31_2014)) %>% 
  left_join(dplyr::select(fecund_slope_site, id, shape_fec_site), by = 'id') %>% 
  left_join(dplyr::select(fecund_slope_strain, id, shape_fec_strain), by = 'id') %>% 
  left_join(dplyr::select(fecund_slope_block, id, shape_fec_block), by = 'id')

# print mean shape_fecundity by strain, from highest to lowest
dat_join %>% 
  group_by(strain) %>% 
  summarize(shape_fecund = mean(shape_fec_strain, na.rm = T)) %>% 
  arrange(desc(shape_fecund))




##### 18. Compare strain-specific life history trait values across blocks
BlockCompare <- function(data, var_focal, label, shape_mort = FALSE) {
# for shape_mort, SE is SD of bootstrap values
# for all other traits, SE calculated as normal
  
  df <- data %>% 
    rename_(y = var_focal) %>% 
    group_by(strain, block) %>% 
    summarise(y_mean = mean(y, na.rm = TRUE),
              y_se = ifelse(shape_mort == TRUE,
                            sd(y, na.rm = TRUE),
                            sd(y, na.rm = TRUE) / sqrt(length(y[!is.na(y)])))) %>% 
    setDT() %>%
    dcast(strain ~ block,value.var = c('y_mean', 'y_se')) %>% 
    mutate(plot_min = min(c(y_mean_a - y_se_a, y_mean_b - y_se_b)),
           plot_max = max(c(y_mean_a + y_se_a, y_mean_b + y_se_b)),
           name = label)
  
  return(df)
}

trait_strain_block <- rbind.data.frame(
  BlockCompare(dat_join, 'lifespan', 'lifespan'),
  BlockCompare(shape_block, 'shape_mort_boot', 'shape_mort', TRUE),
  BlockCompare(dat_join, 'shape_fec_block', 'shape_fecund'),
  BlockCompare(dat_join, 'total_offspring', 'total_offspring'),
  BlockCompare(dat_join, 'area', 'area1'),
  BlockCompare(dat_areas_supp, 'area', 'area2')
)




##### 19. Life history trait variance components (block, strain, site)
# (Table 1)

# read model data already written to file, from code below (commented-out code
#  below takes ~5 minutes to run, and requires installation of RStan
var_strain_df <- read_csv('analysis/var_strain_summary.csv')
var_site_df <- read_csv('analysis/var_site_summary.csv')

# # Function to obtain variance components
# StanVarComp <- function(data, y, J, K, label, shape_mort = FALSE) {
# 
#   # organize data
#   data <- data %>%
#     rename_(y = y, J = J, K = K) %>%
#     filter(!is.na(y)) %>%
#     mutate(K = paste0(J, K))
# 
#   if (shape_mort == TRUE) {  # if trait is shape_mortality
#     data <- data %>%
#       group_by(J, K) %>%
#       summarize(y_mean = mean(y), y_sd = sd(y)) %>%
#       ungroup()
# 
#     dat_stan <- list(
#       N = nrow(data),
#       N_J = length(unique(data$J)),
#       N_K = length(unique(data$K)),
#       y = data$y_mean,
#       J = as.numeric(as.factor(data$J)),
#       K = as.numeric(as.factor(data$K)),
#       sigma_within = mean(data$y_sd)
#     )
# 
#     model <- stan_varcomp_shape
# 
#   } else {                   # if trait is not shape_mortality
#     dat_stan <- list(
#       N = nrow(data),
#       N_J = length(unique(data$J)),
#       N_K = length(unique(data$K)),
#       y = data$y,
#       J = as.numeric(as.factor(data$J)),
#       K = as.numeric(as.factor(data$K))
#     )
# 
#     model <- stan_varcomp
#   }
# 
#   # fit model
#   fit <- sampling(
#     model,
#     data = dat_stan,
#     pars = c('eta_J', 'eta_K'),  # don't store fitted values for eta
#     include = FALSE,             # don't store fitted values for eta
#     warmup = 2000,
#     iter = 4000,
#     thin = 2,
#     chains = 2,
#     control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 12)
#   )
# 
#   # extract posterior samples of variance parameters
#   var_within <- rstan::extract(fit, 'var_within')$var_within
#   var_J <- rstan::extract(fit, 'var_J')$var_J
#   var_K <- rstan::extract(fit, 'var_K')$var_K
#   var_total <- rstan::extract(fit, 'var_total')$var_total
# 
#   # calculate variance components (also called intraclass correlation)
#   icc_within <- var_within / var_total * 100
#   icc_J <- var_J / var_total * 100
#   icc_K <- var_K / var_total * 100
#   icc_ratio <- icc_J / icc_K
# 
#   return(data.frame(icc_within, icc_J, icc_K, icc_ratio, label))
# }
# 
# # pre-compile stan models
# stan_varcomp <- stan_model('stan/varcomp-nested.stan')
# stan_varcomp_shape <- stan_model('stan/varcomp-nested-shape-mort.stan')
# 
# # strain/block
# var_strain_life <- dat_join %>%
#   StanVarComp('lifespan', 'strain', 'block', 'lifespan')
# 
# var_strain_shape_mort <- shape_block %>%
#   StanVarComp('shape_mort_boot', 'strain', 'block', 'shape_mort', TRUE)
# 
# var_strain_shape_fecund <- dat_join %>%
#   StanVarComp('shape_fec_block', 'strain', 'block', 'shape_fecund')
# 
# var_strain_total_offspr <- dat_join %>%
#   StanVarComp('total_offspring', 'strain', 'block', 'total_offspring')
# 
# var_strain_area1 <- dat_join %>%
#   StanVarComp('area', 'strain', 'block', 'area1')
# 
# var_strain_area2 <- dat_areas_supp %>%
#   StanVarComp('area', 'strain', 'block', 'area2')
# 
# # site/strain
# var_site_life <- filter(dat_join, site %in% replic_col$site) %>%
#   StanVarComp('lifespan', 'site', 'strain', 'lifespan')
# 
# var_site_shape_mort <- filter(shape_strain, site %in% replic_col$site) %>%
#   StanVarComp('shape_mort_boot', 'site', 'strain', 'shape_mort', TRUE)
# 
# var_site_shape_fecund <- filter(dat_join, site %in% replic_col$site) %>%
#   StanVarComp('shape_fec_strain', 'site', 'strain', 'shape_fecund')
# 
# var_site_total_offspr <- filter(dat_join, site %in% replic_col$site) %>%
#   StanVarComp('total_offspring', 'site', 'strain', 'total_offspring')
# 
# var_site_area1 <- filter(dat_join, site %in% replic_col$site) %>%
#   StanVarComp('area', 'site', 'strain', 'area1')
# 
# var_site_area2 <- dat_areas_supp %>%
#   StanVarComp('area', 'site', 'strain', 'area2')
# 
# # summarize strain/block
# var_strain_df <- rbind.data.frame(
#   var_strain_life,
#   var_strain_shape_mort,
#   var_strain_shape_fecund,
#   var_strain_total_offspr,
#   var_strain_area1,
#   var_strain_area2
# ) %>% group_by(label) %>%
#   summarize(strain_med = quantile(icc_J, 0.500),
#             strain_low = quantile(icc_J, 0.025),
#             strain_upp = quantile(icc_J, 0.975),
#             block_med = quantile(icc_K, 0.500),
#             block_low = quantile(icc_K, 0.025),
#             block_upp = quantile(icc_K, 0.975),
#             ratio_med = quantile(icc_ratio, 0.500),
#             ratio_low = quantile(icc_ratio, 0.025),
#             ratio_upp = quantile(icc_ratio, 0.975))
# 
# # summarize site/strain
# var_site_df <- rbind.data.frame(
#   var_site_life,
#   var_site_shape_mort,
#   var_site_shape_fecund,
#   var_site_total_offspr,
#   var_site_area1,
#   var_site_area2
# ) %>% group_by(label) %>%
#   summarize(site_med = quantile(icc_J, 0.500),
#             site_low = quantile(icc_J, 0.025),
#             site_upp = quantile(icc_J, 0.975),
#             strain_med = quantile(icc_K, 0.500),
#             strain_low = quantile(icc_K, 0.025),
#             strain_upp = quantile(icc_K, 0.975),
#             ratio_med = quantile(icc_ratio, 0.500),
#             ratio_low = quantile(icc_ratio, 0.025),
#             ratio_upp = quantile(icc_ratio, 0.975))
# 
# # write to file
# write.csv(var_strain_df, 'analysis/var_strain_summary.csv', row.names = F)
# write.csv(var_site_df, 'analysis/var_site_summary.csv', row.names = F)




##### 20. Correlations among life history traits, within and among strains

## among-strain correlations
AmongTest <- function(data, var_x, var_y) {
  dat_sub <- data %>% 
    dplyr::select(strain, var_x , var_y) %>% 
    setNames(c('strain', 'x' ,'y')) %>% 
    filter(!is.na(x) & !is.na(y))
  
  mod <- lme(y ~ scale(x), random = ~ 1|strain, data = dat_sub,
             method = 'ML', control = lmeControl(opt = 'optim'))
  
  beta_a <- summary(mod)$tTable[2,'Value']
  t_a <- summary(mod)$tTable[2,'t-value']
  DF_a <- summary(mod)$tTable[2,'DF']
  Pval_a <- summary(mod)$tTable[2,'p-value']
  sig_a <- ifelse(Pval_a < 0.05/10, 'sig', NA)
  
  return(data.frame(beta_a, t_a, DF_a, Pval_a, sig_a, stringsAsFactors = F))
}

# focal variables
vars_foc <- c('lifespan',
              'shape_mort',
              'shape_fecund',
              'total_offspring',
              'area')

# strain-specific means of shape_mortality
strain_means_shape <- shape_strain %>% 
  group_by(strain) %>% 
  summarize(g_shape_mort = mean(shape_mort_boot))

# strain-specific means of all other traits
strain_means_full <- dat_join %>% 
  group_by(strain) %>% 
  summarize(g_lifespan = mean(lifespan, na.rm = T),
            g_shape_fecund = mean(shape_fec_strain, na.rm = T),
            g_total_offspring = mean(total_offspring, na.rm = T),
            g_area = mean(area, na.rm = T)) %>%
  ungroup() %>% 
  left_join(strain_means_shape, by = 'strain') # join shape_mort data

# join strain-specific means to frond-specific data
# 'g_' prefix represents group-level (i.e. cohort-level) traits
dat_among <- dat_join %>% 
  rename(shape_fecund = shape_fec_strain) %>% 
  dplyr::select(id, strain, lifespan, shape_fecund, total_offspring, area) %>% 
  left_join(strain_means_full, by = 'strain')

# create every combination of response and predictor
var_among_comb <- expand.grid(vars_foc, vars_foc, stringsAsFactors = F) %>% 
  dplyr::select(x_var = Var2, y_var = Var1) %>% 
  filter(x_var != y_var) %>% 
  filter(y_var != 'shape_mort') %>% 
  mutate(x_var = paste0('g_', x_var))

# estimate among-strain correlation between all pairs of traits
trait_corr_among <- var_among_comb %>%
  group_by(y_var, x_var) %>% 
  do(AmongTest(data = dat_among, var_x = .$x_var, var_y = .$y_var)) %>% 
  ungroup() %>% 
  mutate(x_var = gsub('g_', '', x_var))

## within-strain correlation
WithinTest <- function(data, var_x, var_y) {
  ScaleFn <- function(x, y) {
    data.frame(x = scale(x), y = scale(y))
  }
  
  dat_sub <- data %>% 
    dplyr::select(id, strain, var_x, var_y) %>% 
    setNames(c('id', 'strain', 'x' ,'y')) %>% 
    filter(!is.na(x) & !is.na(y)) %>%
    group_by(strain) %>%
    do(ScaleFn(x = .$x, y = .$y)) %>%
    ungroup()
  
  mod <- lme(y ~ x, random = ~ x|strain, data = dat_sub,
             method = 'ML', control = lmeControl(opt = 'optim'))
  
  beta_w <- summary(mod)$tTable[2,'Value']
  t_w <- summary(mod)$tTable[2,'t-value']
  DF_w <- summary(mod)$tTable[2,'DF']
  Pval_w <- summary(mod)$tTable[2,'p-value']
  sig_w <- ifelse(Pval_w < 0.05/6, 'sig', NA)
  
  return(data.frame(beta_w, t_w, DF_w, Pval_w, sig_w, stringsAsFactors = F))
}

# create every combination of response and predictor
var_within_comb <- expand.grid(vars_foc[-2], vars_foc[-2], stringsAsFactors = F) %>% 
  dplyr::select(x_var = Var2, y_var = Var1) %>% 
  filter(x_var != y_var)

# organize data for within-strain correlations
dat_within <- dat_join %>% 
  rename(shape_fecund = shape_fec_strain) %>% 
  dplyr::select(id, strain, lifespan, shape_fecund, total_offspring, area)

# estimate within-strain correlation between all pairs of traits
trait_corr_within <- var_within_comb %>%
  group_by(y_var, x_var) %>% 
  do(WithinTest(data = dat_within, var_x = .$x_var, var_y = .$y_var)) %>% 
  ungroup()

# combine within- and among-strain results (Table S3)
trait_corr_full <- trait_corr_within %>% 
  right_join(trait_corr_among, by = c('y_var', 'x_var')) %>% 
  mutate(x_var = factor(x_var, levels = vars_foc)) %>%
  mutate(y_var = factor(y_var, levels = vars_foc[-2])) %>%
  as.data.frame() %>% 
  slice(order(y_var, x_var))

# write to file
# write.csv(trait_corr_full, 'analysis/life_history_relationships.csv', row.names = F)




##### 21. Relationship between site-level traits and site-level environmental
# characteristics

# read model data already written to file, obtained via commented-out code below
#  (which takes ~5 minutes to run, and requires installation of RStan)
df_pred <- read_csv('analysis/trait_vs_env_pred.csv')
df_alpha <- read_csv('analysis/trait_vs_env_alpha.csv')
beta_summary <- read_csv('analysis/trait_vs_env_beta.csv')

# ## function to fit stan models
# StanFn <- function(data, var, X, seed, shape_mort = FALSE) {
# 
#   # remove NAs
#   which_na <- which(is.na(data[,var]))
#   if(length(which_na)) data <- data[-which_na,]
# 
#   # arrange data for stan
#   dat_stan <- list(
#     n = nrow(data),
#     k = ncol(X),
#     nsites = length(unique(data$site)),
#     site = as.numeric(as.factor(data$site)),
#     y = unlist(data[,var]),
#     X = X
#   )
# 
#   # if trait is shape_mort
#   if(shape_mort == TRUE) {
#     dat_stan$sigma_y = data$y_sd
#     model <- stan_multilevel_shape
#   } else {
#     model <- stan_multilevel 
#   }
# 
#   # fit stan model
#   sampling(
#     model,
#     data = dat_stan,
#     par = c('alpha', 'beta', 'beta_std'),
#     warmup = 2000,
#     iter = 4000,
#     thin = 2,
#     chains = 2,
#     control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 12)
#   )
# }
# 
# ## site-specific environmental characteristics
# # 22 sites in common garden
# X <- filter(dat_sites, dropped == FALSE) %>%
#   model.matrix(~ log(conductivity) + nutrient_pc1 + log(degree_days), data = .)
# 
# # 24 sites for supplementary experiment on frond size
# X_full <- dat_sites %>%
#   model.matrix(~ log(conductivity) + nutrient_pc1 + log(degree_days), data = .)
# 
# ## site specific mean and se in shape_mort
# shape_site_summary <- shape_site %>%
#   group_by(site) %>%
#   summarize(y_mean = mean(shape_mort_boot), y_sd = sd(shape_mort_boot))
# 
# # pre-compile stan models
# stan_multilevel <- stan_model('stan/multilevel.stan')
# stan_multilevel_shape <- stan_model('stan/multilevel-shape-mort.stan')
# 
# # fit models
# fit_life <- StanFn(dat_join, 'lifespan', X, 123450)
# fit_shapem <- StanFn(shape_site_summary, 'y_mean', X, 234501, TRUE)
# fit_shapef <- StanFn(dat_join, 'shape_fec_site', X, 345012)
# fit_to <- StanFn(dat_join, 'total_offspring', X, 450123)
# fit_area1 <- StanFn(dat_join, 'area', X, 501234)
# fit_area2 <- StanFn(dat_areas_supp, 'area', X_full, 012345)
# 
# # extract posterior samples for parameters of interest
# beta_life <- rstan::extract(fit_life, pars = "beta")$beta
# beta_shapef <- rstan::extract(fit_shapef, pars = "beta")$beta
# beta_shapem <- rstan::extract(fit_shapem, pars = "beta")$beta
# beta_to <- rstan::extract(fit_to, pars = "beta")$beta
# beta_area1 <- rstan::extract(fit_area1, pars = "beta")$beta
# beta_area2 <- rstan::extract(fit_area2, pars = "beta")$beta
# 
# beta_std_life <- rstan::extract(fit_life, pars = "beta_std")$beta_std
# beta_std_shapef <- rstan::extract(fit_shapef, pars = "beta_std")$beta_std
# beta_std_shapem <- rstan::extract(fit_shapem, pars = "beta_std")$beta_std
# beta_std_to <- rstan::extract(fit_to, pars = "beta_std")$beta_std
# beta_std_area1 <- rstan::extract(fit_area1, pars = "beta_std")$beta_std
# beta_std_area2 <- rstan::extract(fit_area2, pars = "beta_std")$beta_std
# 
# alpha_life <- rstan::extract(fit_life, pars = 'alpha', permuted = T)
# alpha_shapef <- rstan::extract(fit_shapef, pars = 'alpha', permuted = T)
# alpha_shapem <- rstan::extract(fit_shapem, pars = 'alpha', permuted = T)
# alpha_to <- rstan::extract(fit_to, pars = 'alpha', permuted = T)
# alpha_area1 <- rstan::extract(fit_area1, pars = 'alpha', permuted = T)
# alpha_area2 <- rstan::extract(fit_area2, pars = 'alpha', permuted = T)
# 
# ## posterior summaries for standardized regression coefficients
# # fns to get posterior probability mass > 0, and 95% CI
# ProbG0 <- function(y) length(which(y < 0)) / length(y) * 100
# Cred95 <- function(z) quantile(z, probs = c(0.5, 0.025, 0.975))
# 
# # wrapper fn for ProbG0 and Cred95
# BetaSummary <- function(beta, label) {
#   t(apply(beta, 2, Cred95)) %>%
#     round(3) %>%
#     as.data.frame() %>%
#     setNames(c('med', 'low95', 'upp95')) %>%
#     mutate(pg0 = round(apply(beta, 2, ProbG0), 3)) %>%
#     mutate(y = label) %>%
#     mutate(x = c('conduct', 'nutrient', 'degree_days'))
# }
# 
# # get posterior probability mass > 0, and 95% CI (Table S5)
# beta_summary <- BetaSummary(beta_std_life, 'lifespan') %>%
#   rbind.data.frame(BetaSummary(beta_std_shapem, 'shape_mort')) %>%
#   rbind.data.frame(BetaSummary(beta_std_shapef, 'shape_fecund')) %>%
#   rbind.data.frame(BetaSummary(beta_std_to, 'total_offspring')) %>%
#   rbind.data.frame(BetaSummary(beta_std_area1, 'area1')) %>%
#   rbind.data.frame(BetaSummary(beta_std_area2, 'area2'))
# 
# ## regular sequence of predictor values for prediction lines
# # conductivity varies, others held at mean
# new_X1 <- cbind(
#   b0 = 1,
#   b1 = seq(min(log(dat_sites$conductivity)), max(log(dat_sites$conductivity)), length = 100),
#   b2 = mean(dat_sites$nutrient_pc1),
#   b3 = mean(log(dat_sites$degree_days))
# )
# 
# # nutrient pc1 varies, others held at mean
# new_X2 <- cbind(
#   b0 = 1,
#   b1 = mean(log(dat_sites$conductivity)),
#   b2 = seq(min(dat_sites$nutrient_pc1), max(dat_sites$nutrient_pc1), length = 100),
#   b3 = mean(log(dat_sites$degree_days))
# )
# 
# # degree-days varies, others held at mean
# new_X3 <- cbind(
#   b0 = 1,
#   b1 = mean(log(dat_sites$conductivity)),
#   b2 = mean(dat_sites$nutrient_pc1),
#   b3 = seq(min(log(dat_sites$degree_days)), max(log(dat_sites$degree_days)), length = 100)
# )
# 
# ## predicted values for site-level models
# pred_life <- tibble(
#   p1 = list(new_X1 %*% t(beta_life)), # conductivity
#   p2 = list(new_X2 %*% t(beta_life)), # nutrient pc1
#   p3 = list(new_X3 %*% t(beta_life))  # degree-days
# )
# 
# pred_shapem <- tibble(
#   p1 = list(new_X1 %*% t(beta_shapem)),
#   p2 = list(new_X2 %*% t(beta_shapem)),
#   p3 = list(new_X3 %*% t(beta_shapem))
# )
# 
# pred_shapef <- tibble(
#   p1 = list(new_X1 %*% t(beta_shapef)),
#   p2 = list(new_X2 %*% t(beta_shapef)),
#   p3 = list(new_X3 %*% t(beta_shapef))
# )
# 
# pred_to <- tibble(
#   p1 = list(new_X1 %*% t(beta_to)),
#   p2 = list(new_X2 %*% t(beta_to)),
#   p3 = list(new_X3 %*% t(beta_to))
# )
# 
# pred_area2 <- tibble(
#   p1 = list(new_X1 %*% t(beta_area2)),
#   p2 = list(new_X2 %*% t(beta_area2)),
#   p3 = list(new_X3 %*% t(beta_area2))
# )
# 
# ## posterior median prediction lines and 95% credible intervals
# # x1 is conductivity, x2 is nutrient pc1, and x3 is degree-days
# PredQuantiles <- function(x, label) {
#   tibble(
#     low1 = apply(x$p1[[1]], 1, quantile, probs = 0.025),
#     med1 = apply(x$p1[[1]], 1, quantile, probs = 0.500),
#     upp1 = apply(x$p1[[1]], 1, quantile, probs = 0.975),
#     low2 = apply(x$p2[[1]], 1, quantile, probs = 0.025),
#     med2 = apply(x$p2[[1]], 1, quantile, probs = 0.500),
#     upp2 = apply(x$p2[[1]], 1, quantile, probs = 0.975),
#     low3 = apply(x$p3[[1]], 1, quantile, probs = 0.025),
#     med3 = apply(x$p3[[1]], 1, quantile, probs = 0.500),
#     upp3 = apply(x$p3[[1]], 1, quantile, probs = 0.975)
#   ) %>% setNames(paste(label, names(.), sep = '_'))
# }
# 
# df_pred <- data.frame(x1 = new_X1[,2], x2 = new_X2[,3], x3 = new_X3[,4]) %>%
#   cbind.data.frame(PredQuantiles(pred_life, 'life')) %>%
#   cbind.data.frame(PredQuantiles(pred_shapem, 'shapem')) %>%
#   cbind.data.frame(PredQuantiles(pred_shapef, 'shapef')) %>%
#   cbind.data.frame(PredQuantiles(pred_to, 'to')) %>%
#   cbind.data.frame(PredQuantiles(pred_area2, 'area2'))
# 
# ## posterior median site-specific intercepts and 95% credible intervals
# AlphaQuantiles <- function(x, label) {
#   tibble(
#     low = apply(x, 2, quantile, probs = 0.025),
#     med = apply(x, 2, quantile, probs = 0.500),
#     upp = apply(x, 2, quantile, probs = 0.975)
#   ) %>% setNames(paste(label, names(.), sep = '_'))
# }
# 
# # main common garden experiment (Nsite = 22)
# df_alpha_main <- dat_sites %>%
#   filter(dropped == FALSE) %>%
#   dplyr::select(site) %>%
#   cbind.data.frame(AlphaQuantiles(alpha_life$alpha, 'life')) %>%
#   cbind.data.frame(AlphaQuantiles(alpha_shapem$alpha, 'shapem')) %>%
#   cbind.data.frame(AlphaQuantiles(alpha_shapef$alpha, 'shapef')) %>%
#   cbind.data.frame(AlphaQuantiles(alpha_to$alpha, 'to'))
# 
# # supplementary experiment with all available strains (Nsite = 24)
# df_alpha_supp <- dat_sites %>%
#   dplyr::select(site) %>%
#   cbind.data.frame(X_full[,2:4]) %>%
#   as_tibble() %>%
#   setNames(c('site', 'conduct', 'nutrient', 'degdays')) %>%
#   cbind.data.frame(AlphaQuantiles(alpha_area2$alpha, 'area2'))
# 
# # join
# df_alpha <- df_alpha_supp %>%
#   left_join(df_alpha_main, by = 'site')
# 
# # write to file
# write.csv(beta_summary, 'analysis/trait_vs_env_beta.csv', row.names = F)
# write.csv(df_pred, 'analysis/trait_vs_env_pred.csv', row.names = F)
# write.csv(df_alpha, 'analysis/trait_vs_env_alpha.csv', row.names = F)




##### 22. Spatial autocorrelation in traits or environmental characteristics
# (Table S4)

# check for shared climate stations
dat_sites %>% 
  arrange(climate_station) %>% 
  as.data.frame()

# compile median life history traits by site
meds_main <- dat_join %>% 
  group_by(site) %>% 
  summarize(med_lifespan = median(lifespan, na.rm = T),
            med_shapef = median(shape_fec_site, na.rm = T),
            med_total_offspring = median(total_offspring, na.rm = T),
            med_area1 = median(area, na.rm = T))

meds_shape_mort <- shape_site %>% 
  group_by(site) %>% 
  summarize(med_shapem = median(shape_mort_boot, na.rm = T))

meds_areas_supp <- dat_areas_supp %>% 
  group_by(site) %>% 
  summarize(med_area2 = median(area, na.rm = T))

# combine site-specific medians into single df
meds_full <- dat_sites %>% 
  left_join(meds_main, by = 'site') %>% 
  left_join(meds_shape_mort, by = 'site') %>% 
  left_join(meds_areas_supp, by = 'site')

# distance matrix
dists <- matrix(ncol = nrow(meds_full), nrow = nrow(meds_full))

for(x in 1:nrow(meds_full)) {
  for(y in 1:nrow(meds_full)) {
    dists[x,y]<- distCosine(c(meds_full$lon[x], meds_full$lat[x]),
                            c(meds_full$lon[y], meds_full$lat[y]))
  }
}

dists_inv <- 1/dists  # invert distance matrix
diag(dists_inv) <- 0  # set diagonal to zero

# calculate Moran's I and relevent stats for each variable of interest
MoranFn <- function(value, dists_inv) {
  Moran.I(value, dists_inv, na.rm = T) %>% 
    as.data.frame() %>%
    mutate(zI = (observed - expected) / sd)
}

moran_summary <- meds_full %>% 
  dplyr::select(site, conductivity, nutrient_pc1, degree_days, starts_with('med')) %>% 
  mutate(degree_days = log(degree_days),
         conductivity = log(conductivity)) %>% 
  gather(var, value, conductivity:med_area2) %>% 
  group_by(var) %>% 
  do(MoranFn(.$value, dists_inv)) %>% 
  ungroup() %>% 
  # rearange rows and columns
  dplyr::select(var, observed, expected, sd, zI, p.value) %>% 
  slice(c(1, 9, 2, 5, 7, 6, 8, 3, 4))

# # write to file
# write.csv(moran_summary, 'analysis/moran_summary.csv', row.names = F)






### -------------------------------------------------------------------------- #
### Figure 1 ----------------------------------------------------------------- #
### survivorship, mortality, and fecundity trajectories by strain ------------ #
### -------------------------------------------------------------------------- #

# ggplot theme
tt1 <- theme_bw() +
  theme(axis.title = element_text(size = 13),
        plot.title = element_text(size = 13, hjust = 0.5),
        text = element_text(size = 12),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.3, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .3, 0, 0, unit = 'cm')))

# function for labelling axes
LogLab <- function(x) ifelse(x %in% c(0.1, 1), as.character(x), '')

# generate all six plot panels
p1_1 <- ggplot(surv_strain, aes(age, surv, group = strain)) +
  geom_step(size = 0.75, alpha = 0.3, direction = 'hv') +
  coord_cartesian(xlim = c(0, 36), ylim = c(0.02, 1)) +
  scale_y_log10(breaks = c((1:10)/100, (1:10)/10), labels = LogLab) +
  xlab(NULL) +
  ylab('Survivorship') +
  ggtitle('Non-standardized') +
  annotate('text', x = 37, y = Inf, label = '(a)', hjust = 1, vjust = 1.6, size = 4.5) +
  tt1

p1_2 <- ggplot(surv_strain_std, aes(age, surv, group = strain)) +
  geom_step(size = 0.75, alpha = 0.3, direction = 'hv') +
  coord_cartesian(xlim = c(0, 1.8), ylim = c(0.02, 1)) +
  scale_x_continuous(labels = function(x) formatC(x, width = 1)) +
  scale_y_log10(breaks = c((1:10)/100, (1:10)/10), labels = LogLab) +
  coord_cartesian(xlim = c(0, 1.8), ylim = c(0.02, 1)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle('Standardized') +
  annotate('text', x = 1.85, y = Inf, label = '(b)', hjust = 1, vjust = 1.6, size = 4.5) +
  tt1

p1_3 <- ggplot(hazard_spline) +
  geom_ribbon(aes(x = newage, ymin = haz_lower, ymax = haz_upper, group = strain), fill = 'grey40', alpha = 0.1)+
  geom_line(aes(x = newage, y = haz, group = strain), col = 'darkred', size = 0.6) +
  coord_cartesian(xlim = c(0, 36), ylim = c(0, 1.3)) +
  scale_y_continuous(labels = function(x) formatC(x, width = 1)) +
  xlab(NULL) +
  ylab('Mortality') +
  annotate('text', x = 37, y = Inf, label = '(c)', hjust = 1, vjust = 1.6, size = 4.5) +
  tt1

p1_4 <- ggplot(hazard_spline) +
  geom_ribbon(aes(x = newagestd, ymin = haz_std_lower, ymax = haz_std_upper, group = strain), fill = 'grey40', alpha = 0.1)+
  geom_line(aes(x = newagestd, y = haz_std, group = strain), col = 'darkred', size = 0.6) +
  scale_x_continuous(labels = function(x) formatC(x, width = 1)) +
  scale_y_continuous(breaks = c(0, 10, 20)) +
  coord_cartesian(xlim = c(0, 1.8), ylim = c(0, 25)) +
  xlab(NULL) +
  ylab(NULL) +
  annotate('text', x = 1.85, y = Inf, label = '(d)', hjust = 1, vjust = 1.6, size = 4.5) +
  tt1

p1_5 <- ggplot(filter(dat_tidy_std_strain, uncertain_repro == FALSE), aes(age, fecund, group = strain)) +
  geom_smooth(method = 'loess', se = T, size = 0, fill = 'grey40', alpha = 0.1) +
  geom_smooth(method = 'loess', se = F, size = 0.6, col = 'darkred', fill = 1) +
  coord_cartesian(xlim = c(0, 36), ylim = c(0, 1.05)) +
  scale_y_continuous(breaks = seq(0, 1, 0.5), labels = function(x) formatC(x, width = 1)) +
  xlab('Age (days)') +
  ylab('Fecundity') +
  annotate('text', x = 37, y = Inf, label = '(e)', hjust = 1, vjust = 1.6, size = 4.5) +
  tt1

p1_6 <- ggplot(filter(dat_tidy_std_strain, uncertain_repro == FALSE), aes(age_std, fecund_std, group = strain)) +
  geom_smooth(method = 'loess', se = T, size = 0, fill = 'grey40', alpha = 0.1) +
  geom_smooth(method = 'loess', se = F, size = 0.6, col = 'darkred', fill = 1) +
  coord_cartesian(xlim = c(0, 1.8), ylim = c(0, 1.75)) +
  scale_x_continuous(labels = function(x) formatC(x, width = 1)) +
  scale_y_continuous(labels = function(x) formatC(x, width = 1)) +
  xlab('Age (life expectancies)') +
  ylab(NULL) +
  annotate('text', x = 1.85, y = Inf, label = '(f)', hjust = 1, vjust = 1.6, size = 4.5) +
  tt1

# convert panels to ggplot grobs
g1_1 <- ggplotGrob(p1_1)
g1_2 <- ggplotGrob(p1_2)
g1_3 <- ggplotGrob(p1_3)
g1_4 <- ggplotGrob(p1_4)
g1_5 <- ggplotGrob(p1_5)
g1_6 <- ggplotGrob(p1_6)

# standardize panel dimensions
g1_3$widths <- g1_1$widths
g1_5$widths <- g1_1$widths
g1_2$widths <- g1_6$widths
g1_4$widths <- g1_6$widths

# generate full figure
fig_1 <- arrangeGrob(g1_1, g1_2, g1_3, g1_4, g1_5, g1_6, nrow = 3,
                     heights = c(1.09, 0.98, 1.09), widths = c(1.062, 1))

# # print to Mac OSX device (change 'quartz' to 'window' if running Windows)
# dev.off()
# quartz(height = 7, width = 7)
# grid.arrange(fig_1)

# write to file
ggsave('img/Fig_1.png', fig_1, height = 7, width = 7, units = 'in', dpi = 300)




### -------------------------------------------------------------------------- #
### Figure 2 ----------------------------------------------------------------- #
### correlation between strain-specific trait values in block a vs b --------- #
### -------------------------------------------------------------------------- #

# ensure panel x and y limit are equal, then add 10% extra space at top
panel_lims <- trait_strain_block %>% group_by(strain, name) %>% 
  summarize(plot_min = unique(plot_min),
            plot_max = unique(plot_max)) %>% 
  mutate(buffer = plot_max + (plot_max - plot_min) * 0.1)

# ggplot theme
tt2 <- theme_bw() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.3, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .3, 0, 0, unit = 'cm')))

# generate all six panels in plot
p2_1 <- ggplot(filter(trait_strain_block, name == 'lifespan')) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
  geom_linerange(aes(x = y_mean_a, ymin = y_mean_b - y_se_b, ymax = y_mean_b + y_se_b), col = 'grey60') +
  geom_errorbarh(aes(y = y_mean_b, x = y_mean_a, xmin = y_mean_a - y_se_a, xmax = y_mean_a + y_se_a), col = 'grey60', height = 0) +
  annotate('text', x = -Inf, y = Inf, label = 'Lifespan', hjust = -0.13, vjust = 1.6, size = 4) +
  annotate('text', x = Inf, y = -Inf, label = '(a)', hjust = 1.25, vjust = -0.6, size = 4) +
  scale_x_continuous(breaks = seq(17, 25, 2)) +
  scale_y_continuous(breaks = seq(17, 25, 2)) +
  geom_point(aes(y_mean_a, y_mean_b), size = 1.2) +
  geom_blank(data = filter(panel_lims, name == 'lifespan'), aes(plot_min, plot_min)) +
  geom_blank(data = filter(panel_lims, name == 'lifespan'), aes(plot_max, buffer)) +
  xlab(NULL) + ylab(NULL) +
  tt2

p2_2 <- ggplot(filter(trait_strain_block, name == 'shape_mort')) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
  geom_linerange(aes(x = y_mean_a, ymin = y_mean_b - y_se_b, ymax = y_mean_b + y_se_b), col = 'grey60') +
  geom_errorbarh(aes(y = y_mean_b, x = y_mean_a, xmin = y_mean_a - y_se_a, xmax = y_mean_a + y_se_a), col = 'grey60', height = 0) +
  annotate('text', x = -Inf, y = Inf, label = 'shape[mortality]', hjust = -0.08, vjust = 1.37, size = 4, parse = T) +
  annotate('text', x = Inf, y = -Inf, label = '(b)', hjust = 1.25, vjust = -0.6, size = 4) +
  geom_point(aes(y_mean_a, y_mean_b), size = 1.2) +
  geom_blank(data = filter(panel_lims, name == 'shape_mort'), aes(plot_min, plot_min)) +
  geom_blank(data = filter(panel_lims, name == 'shape_mort'), aes(plot_max, buffer)) +
  xlab(NULL) + ylab(NULL) +
  tt2

p2_3 <- ggplot(filter(trait_strain_block, name == 'shape_fecund')) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
  geom_linerange(aes(x = y_mean_a, ymin = y_mean_b - y_se_b, ymax = y_mean_b + y_se_b), col = 'grey60') +
  geom_errorbarh(aes(y = y_mean_b, x = y_mean_a, xmin = y_mean_a - y_se_a, xmax = y_mean_a + y_se_a), col = 'grey60', height = 0) +
  annotate('text', x = -Inf, y = Inf, label = 'shape[fecundity]', hjust = -0.08, vjust = 1.37, size = 4, parse = T) +
  annotate('text', x = Inf, y = -Inf, label = '(c)', hjust = 1.25, vjust = -0.6, size = 4) +
  scale_x_continuous(breaks = c(-0.8, -0.4, 0, 0.4)) +
  scale_y_continuous(breaks = c(-0.8, -0.4, 0, 0.4)) +
  geom_point(aes(y_mean_a, y_mean_b), size = 1.2) +
  geom_blank(data = filter(panel_lims, name == 'shape_fecund'), aes(plot_min, plot_min)) +
  geom_blank(data = filter(panel_lims, name == 'shape_fecund'), aes(plot_max, buffer)) +
  xlab(NULL) + ylab(NULL) +
  tt2

p2_4 <- ggplot(filter(trait_strain_block, name == 'total_offspring')) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
  geom_linerange(aes(x = y_mean_a, ymin = y_mean_b - y_se_b, ymax = y_mean_b + y_se_b), col = 'grey60') +
  geom_errorbarh(aes(y = y_mean_b, x = y_mean_a, xmin = y_mean_a - y_se_a, xmax = y_mean_a + y_se_a), col = 'grey60', height = 0) +
  annotate('text', x = -Inf, y = Inf, label = 'Cumulative fecundity', hjust = -0.05, vjust = 1.6, size = 4) +
  annotate('text', x = Inf, y = -Inf, label = '(d)', hjust = 1.25, vjust = -0.6, size = 4) +
  scale_x_continuous(breaks = seq(9, 15, 2)) +
  scale_y_continuous(breaks = seq(9, 15, 2)) +
  geom_point(aes(y_mean_a, y_mean_b), size = 1.2) +
  geom_blank(data = filter(panel_lims, name == 'total_offspring'), aes(plot_min, plot_min)) +
  geom_blank(data = filter(panel_lims, name == 'total_offspring'), aes(plot_max, buffer)) +
  xlab(NULL) + ylab('Block B') +
  tt2 + theme(axis.title.y = element_text(hjust = 1.42))

p2_5 <- ggplot(filter(trait_strain_block, name == 'area1')) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
  geom_linerange(aes(x = y_mean_a, ymin = y_mean_b - y_se_b, ymax = y_mean_b + y_se_b), col = 'grey60') +
  geom_errorbarh(aes(y = y_mean_b, x = y_mean_a, xmin = y_mean_a - y_se_a, xmax = y_mean_a + y_se_a), col = 'grey60', height = 0) +
  annotate('text', x = -Inf, y = Inf, label = 'Surface area', hjust = -0.08, vjust = 1.6, size = 4, parse = F) +
  annotate('text', x = Inf, y = -Inf, label = '(e)', hjust = 1.25, vjust = -0.6, size = 4) +
  scale_x_continuous(breaks = seq(4, 8, 1)) +
  scale_y_continuous(breaks = seq(4, 8, 1)) +
  geom_point(aes(y_mean_a, y_mean_b), size = 1.2) +
  geom_blank(data = filter(panel_lims, name == 'area1'), aes(plot_min, plot_min)) +
  geom_blank(data = filter(panel_lims, name == 'area1'), aes(plot_max, buffer)) +
  xlab('Block A') + ylab(NULL) +
  tt2

p2_6 <- ggplot(filter(trait_strain_block, name == 'area2')) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
  geom_linerange(aes(x = y_mean_a, ymin = y_mean_b - y_se_b, ymax = y_mean_b + y_se_b), col = 'grey60') +
  geom_errorbarh(aes(y = y_mean_b, x = y_mean_a, xmin = y_mean_a - y_se_a, xmax = y_mean_a + y_se_a), col = 'grey60', height = 0) +
  annotate('text', x = -Inf, y = Inf, label = 'Surface area (all\navailable strains)', hjust = -0.06, vjust = 1.2, size = 4, parse = F) +
  annotate('text', x = Inf, y = -Inf, label = '(f)', hjust = 1.3, vjust = -0.6, size = 4) +
  scale_x_continuous(breaks = seq(3, 7, 1)) +
  scale_y_continuous(breaks = seq(3, 7, 1)) +
  geom_point(aes(y_mean_a, y_mean_b), size = 1.2) +
  geom_blank(data = filter(panel_lims, name == 'area2'), aes(plot_min, plot_min)) +
  geom_blank(data = filter(panel_lims, name == 'area2'), aes(plot_max, buffer)) +
  xlab(NULL) + ylab(NULL) +
  tt2

# convert panels to ggplot grobs
g2_1 <- ggplotGrob(p2_1)
g2_2 <- ggplotGrob(p2_2)
g2_3 <- ggplotGrob(p2_3)
g2_4 <- ggplotGrob(p2_4)
g2_5 <- ggplotGrob(p2_5)
g2_6 <- ggplotGrob(p2_6)

# standardize panel dimensions
g2_1$widths <- g2_4$widths
g2_4$heights <- g2_5$heights
g2_6$heights <- g2_5$heights

# generate full figure
fig_2 <- arrangeGrob(
  g2_1, g2_2, g2_3, g2_4, g2_5, g2_6,
  nrow = 2, heights = c(1, 1.12), widths = c(1.1, 1, 1)
)

# # print to Mac OSX device (change 'quartz' to 'window' if running Windows)
# dev.off()
# quartz(height = 4.5, width = 7.5)
# grid.arrange(fig_2)

# write to file
ggsave('img/Fig_2.png', fig_2, height = 4.5, width = 7.5, units = 'in', dpi = 300)




### -------------------------------------------------------------------------- #
### Figure 3 ----------------------------------------------------------------- #
### site-specific life history traits vs environmental characteristics ------- #
### -------------------------------------------------------------------------- #

# ggplot theme
tt3 <- theme_bw() +
  theme(axis.title = element_text(size = 12),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.15, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .15, 0, .1, unit = 'cm')),
        axis.ticks = element_line(size = 0.5))

## generate all panels
# pace
p3_01 <- ggplot() +
  geom_ribbon(inherit.aes = F, data = df_pred, aes(x = x1, ymin = life_low1, ymax = life_upp1), fill = 'grey82') +
  geom_line(inherit.aes = F, data = df_pred, aes(x1, life_med1), size = 0.8) +
  geom_linerange(inherit.aes = F, data = df_alpha, aes(x = conduct, ymin = life_low, ymax = life_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(inherit.aes = F, data = df_alpha, aes(conduct, life_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab(NULL) + ylab('Lifespan') +
  scale_y_continuous(limits = c(17.3, 24)) +
  tt3 + theme(axis.text.x = element_blank())

p3_02 <- ggplot() +
  geom_ribbon(inherit.aes = F, data = df_pred, aes(x = x2, ymin = life_low2, ymax = life_upp2), fill = 'grey82') +
  geom_line(inherit.aes = F, data = df_pred, aes(x2, life_med2), size = 0.8) +
  geom_linerange(inherit.aes = F, data = df_alpha, aes(x = nutrient, ymin = life_low, ymax = life_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(inherit.aes = F, data = df_alpha, aes(nutrient, life_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab(NULL) + ylab(NULL) +
  scale_y_continuous(limits = c(17.3, 24)) +
  tt3 + theme(axis.text.y = element_blank()) + theme(axis.text.x = element_blank())

p3_03 <- ggplot() +
  geom_ribbon(inherit.aes = F, data = df_pred, aes(x = x3, ymin = life_low3, ymax = life_upp3), fill = 'grey82') +
  geom_line(inherit.aes = F, data = df_pred, aes(x3, life_med3), size = 0.8) +
  geom_linerange(inherit.aes = F, data = df_alpha, aes(x = degdays, ymin = life_low, ymax = life_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(inherit.aes = F, data = df_alpha, aes(degdays, life_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab(NULL) + ylab(NULL) +
  scale_y_continuous(limits = c(17.3, 24)) +
  tt3 + theme(axis.text.y = element_blank()) + theme(axis.text.x = element_blank())

# shape mortality
p3_04 <- ggplot() +
  geom_ribbon(inherit.aes = F, data = df_pred, aes(x = x1, ymin = shapem_low1, ymax = shapem_upp1), fill = 'grey82') +
  geom_line(inherit.aes = F, data = df_pred, aes(x1, shapem_med1), size = 0.8) +
  geom_linerange(inherit.aes = F, data = df_alpha, aes(x = conduct, ymin = shapem_low, ymax = shapem_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(inherit.aes = F, data = df_alpha, aes(conduct, shapem_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab(NULL) + ylab(expression(shape[mortality])) +
  tt3 + theme(axis.text.x = element_blank())

p3_05 <- ggplot() +
  geom_ribbon(inherit.aes = F, data = df_pred, aes(x = x2, ymin = shapem_low2, ymax = shapem_upp2), fill = 'grey82') +
  geom_line(inherit.aes = F, data = df_pred, aes(x2, shapem_med2), size = 0.8) +
  geom_linerange(inherit.aes = F, data = df_alpha, aes(x = nutrient, ymin = shapem_low, ymax = shapem_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(inherit.aes = F, data = df_alpha, aes(nutrient, shapem_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab(NULL) + ylab(NULL) +
  tt3 + theme(axis.text.y = element_blank()) + theme(axis.text.x = element_blank())

p3_06 <- ggplot() +
  geom_ribbon(inherit.aes = F, data = df_pred, aes(x = x3, ymin = shapem_low3, ymax = shapem_upp3), fill = 'grey82') +
  geom_line(inherit.aes = F, data = df_pred, aes(x3, shapem_med3), size = 0.8) +
  geom_linerange(inherit.aes = F, data = df_alpha, aes(x = degdays, ymin = shapem_low, ymax = shapem_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(inherit.aes = F, data = df_alpha, aes(degdays, shapem_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab(NULL) + ylab(NULL) +
  tt3 + theme(axis.text.y = element_blank()) + theme(axis.text.x = element_blank())

# shape fecundity
p3_07 <- ggplot() +
  geom_ribbon(inherit.aes = F, data = df_pred, aes(x = x1, ymin = shapef_low1, ymax = shapef_upp1), fill = 'grey82') +
  geom_line(inherit.aes = F, data = df_pred, aes(x1, shapef_med1), size = 0.8) +
  geom_linerange(inherit.aes = F, data = df_alpha, aes(x = conduct, ymin = shapef_low, ymax = shapef_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(inherit.aes = F, data = df_alpha, aes(conduct, shapef_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab(NULL) + ylab(expression(shape[fecundity])) +
  scale_y_continuous(breaks = c(-0.8, -0.5, -0.2)) +
  tt3 + theme(axis.text.x = element_blank())

p3_08 <- ggplot() +
  geom_ribbon(inherit.aes = F, data = df_pred, aes(x = x2, ymin = shapef_low2, ymax = shapef_upp2), fill = 'grey82') +
  geom_line(inherit.aes = F, data = df_pred, aes(x2, shapef_med2), size = 0.8) +
  geom_linerange(inherit.aes = F, data = df_alpha, aes(x = nutrient, ymin = shapef_low, ymax = shapef_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(inherit.aes = F, data = df_alpha, aes(nutrient, shapef_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab(NULL) + ylab(NULL) +
  scale_y_continuous(breaks = c(-0.8, -0.5, -0.2)) +
  tt3 + theme(axis.text.y = element_blank()) + theme(axis.text.x = element_blank())

p3_09 <- ggplot() +
  geom_ribbon(inherit.aes = F, data = df_pred, aes(x = x3, ymin = shapef_low3, ymax = shapef_upp3), fill = 'grey82') +
  geom_line(inherit.aes = F, data = df_pred, aes(x3, shapef_med3), size = 0.8) +
  geom_linerange(inherit.aes = F, data = df_alpha, aes(x = degdays, ymin = shapef_low, ymax = shapef_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(inherit.aes = F, data = df_alpha, aes(degdays, shapef_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab(NULL) + ylab(NULL) +
  scale_y_continuous(breaks = c(-0.8, -0.5, -0.2)) +
  tt3 + theme(axis.text.y = element_blank()) + theme(axis.text.x = element_blank())

# cumulative fecundity
p3_10 <- ggplot() +
  geom_ribbon(data = df_pred, aes(x = x1, ymin = to_low1, ymax = to_upp1), fill = 'grey82') +
  geom_line(data = df_pred, aes(x1, to_med1), size = 0.8) +
  geom_linerange(data = df_alpha, aes(x = conduct, ymin = to_low, ymax = to_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(data = df_alpha, aes(conduct, to_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab(NULL) + ylab('Cumulative\nfecundity') +
  scale_y_continuous(limits = c(9, 16.1), breaks = seq(10, 16, 2)) +
  tt3 + theme(axis.text.x = element_blank())

p3_11 <- ggplot() +
  geom_ribbon(data = df_pred, aes(x = x2, ymin = to_low2, ymax = to_upp2), fill = 'grey82') +
  geom_line(data = df_pred, aes(x2, to_med2), size = 0.8) +
  geom_linerange(data = df_alpha, aes(x = nutrient, ymin = to_low, ymax = to_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(data = df_alpha, aes(nutrient, to_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab(NULL) + ylab(NULL) +
  scale_y_continuous(limits = c(9, 16.1), breaks = seq(10, 16, 2)) +
  tt3 + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank())

p3_12 <- ggplot() +
  geom_ribbon(data = df_pred, aes(x = x3, ymin = to_low3, ymax = to_upp3), fill = 'grey82') +
  geom_line(data = df_pred, aes(x3, to_med3), size = 0.8) +
  geom_linerange(data = df_alpha, aes(x = degdays, ymin = to_low, ymax = to_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(data = df_alpha, aes(degdays, to_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab(NULL) + ylab(NULL) +
  scale_y_continuous(limits = c(9, 16.1), breaks = seq(10, 16, 2)) +
  tt3 + theme(axis.text.x = element_blank()) + theme(axis.text.y = element_blank())

# frond surface area (supplementary study with all available strains)
p3_13 <- ggplot() +
  geom_ribbon(data = df_pred, aes(x = x1, ymin = area2_low1, ymax = area2_upp1), fill = 'grey82') +
  geom_line(data = df_pred, aes(x1, area2_med1), size = 0.8) +
  geom_linerange(data = df_alpha, aes(x = conduct, ymin = area2_low, ymax = area2_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(data = df_alpha, aes(conduct, area2_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab('Conductivity') + ylab('Frond area\n(all strains)') +
  scale_y_continuous(limits = c(3.5, 6.5), breaks = 4:6) +
  tt3

p3_14 <- ggplot() +
  geom_ribbon(data = df_pred, aes(x = x2, ymin = area2_low2, ymax = area2_upp2), fill = 'grey82') +
  geom_line(data = df_pred, aes(x2, area2_med2), size = 0.8) +
  geom_linerange(data = df_alpha, aes(x = nutrient, ymin = area2_low, ymax = area2_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(data = df_alpha, aes(nutrient, area2_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab('Nutrients') + ylab(NULL) +
  scale_y_continuous(limits = c(3.5, 6.5), breaks = 4:6) +
  tt3 + theme(axis.text.y = element_blank())

p3_15 <- ggplot() +
  geom_ribbon(data = df_pred, aes(x = x3, ymin = area2_low3, ymax = area2_upp3), fill = 'grey82') +
  geom_line(data = df_pred, aes(x3, area2_med3), size = 0.8) +
  geom_linerange(data = df_alpha, aes(x = degdays, ymin = area2_low, ymax = area2_upp), size = 0.5, col = 'blue', alpha = 0.7) +
  geom_point(data = df_alpha, aes(degdays, area2_med), size = 1.3, col = 'darkblue', shape = 16) +
  xlab('Degree-days') + ylab(NULL) +
  scale_y_continuous(limits = c(3.5, 6.5), breaks = 4:6) +
  tt3 + theme(axis.text.y = element_blank())

# convert panels to ggplot grobs
g3_01 <- ggplotGrob(p3_01)
g3_02 <- ggplotGrob(p3_02)
g3_03 <- ggplotGrob(p3_03)
g3_04 <- ggplotGrob(p3_04)
g3_05 <- ggplotGrob(p3_05)
g3_06 <- ggplotGrob(p3_06)
g3_07 <- ggplotGrob(p3_07)
g3_08 <- ggplotGrob(p3_08)
g3_09 <- ggplotGrob(p3_09)
g3_10 <- ggplotGrob(p3_10)
g3_11 <- ggplotGrob(p3_11)
g3_12 <- ggplotGrob(p3_12)
g3_13 <- ggplotGrob(p3_13)
g3_14 <- ggplotGrob(p3_14)
g3_15 <- ggplotGrob(p3_15)

# standardize panel dimensions
g3_01$widths <- g3_07$widths
g3_04$widths <- g3_07$widths
g3_10$widths <- g3_07$widths
g3_13$widths <- g3_07$widths
g3_14$heights <- g3_13$heights
g3_15$heights <- g3_13$heights

# generate full figure
fig_3 <- arrangeGrob(
  g3_01, g3_02, g3_03,
  g3_04, g3_05, g3_06,
  g3_07, g3_08, g3_09,
  g3_10, g3_11, g3_12,
  g3_13, g3_14, g3_15,
  nrow = 5, widths = c(1.46, 1, 1), heights = c(1, 1, 1, 1, 1.27)
)

# # print to Mac OSX device (change 'quartz' to 'window' if running Windows)
# dev.off()
# quartz(height = 6.5, width = 4.5)
# grid.arrange(fig_3)

# write to file
ggsave('img/Fig_3.png', fig_3, height = 7, width = 4.5, units = 'in', dpi = 300)




### -------------------------------------------------------------------------- #
### Figure S1 ---------------------------------------------------------------- #
### map of study sites ------------------------------------------------------- #
### -------------------------------------------------------------------------- #

# get shapefiles for Alberta and North America
alberta <- raster::getData('GADM', country = 'CAN', level = 1, path = 'dat/') %>% 
  raster::subset(NAME_1 == 'Alberta')

north_america <- borders('world',
                         regions = c('canada', 'usa', 'mexico'),
                         colour = 'gray70', fill = 'gray70')

# lat/lon for Alberta cities to plot
cities <- data.frame(
  name = c('Edmonton', 'Calgary', 'Fort McMurray'),
  lat = c(53.53, 51.05, 56.73),
  lon = c(-113.5, -114.07, -111.38)
)

# ggplot theme
tts1 <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = 'grey70', size = 0.3),
  axis.ticks = element_line(color = 'grey70', size = 0.3),
  axis.title.x = element_text(size = 16, margin = margin(.5, 0, 0, 0, unit = 'cm')),
  axis.title.y = element_text(size = 16, margin = margin(0, .4, 0, 0, unit = 'cm')),
  axis.text = element_text(size = 12, colour = "grey30"),
  plot.margin = unit(c(1, 1, 6, 4), "mm")
)

# generate all panels
ps1_a <- ggplot(alberta) +
  north_america +
  coord_quickmap(xlim = c(-172.1945, -51.65366)) +
  geom_path(aes(long, lat), size = 0.6) +
  scale_x_continuous(labels = seq(175, 50, -25), breaks = seq(-175, -50, 25), expand = c(0.05, 0.05)) +
  xlab("Longitude (°W)") + ylab("Latitude (°N)") +
  tts1

ps1_b <- ggplot(alberta, aes(x = long, y = lat)) +
  geom_polygon(fill = "grey95") +
  geom_path(size = 0.8) +
  coord_map() +
  geom_point(aes(x = lon, y = lat), data = cities, size = 1.9, shape = 1, colour = "black") +
  geom_point(aes(x = lon, y = lat), data = dat_sites, size = 4, shape = 17, colour = "black") +
  geom_point(aes(x = lon, y = lat), data = dat_sites, size = 2.8, shape = 17, colour = "grey40") +
  xlab("Longitude (°W)") + ylab("Latitude (°N)") +
  scale_x_continuous(labels = c("120", "115", "110"), breaks = c(-120, -115, -110), expand = c(0.05, 0.05)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  annotate("text", label = cities$name, x = cities$lon + c(-0.3, -0.2, -0.6), y = cities$lat + c(-0.25, 0.31, 0.25), size = 4.2, colour = "black") +
  annotate("text", label = "Alberta", x = -114.7, y = 59.45, size = 6, colour = 'black') +
  tts1

# convert panels to ggplot grobs
gs1_a <- ggplotGrob(ps1_a)
gs1_b <- ggplotGrob(ps1_b)

# generate full figure
fig_s1 <- arrangeGrob(gs1_a, gs1_b, nrow = 1, widths = c(0.65, 0.35))

# # print to Mac OSX device (change 'quartz' to 'window' if running Windows)
# quartz(height = 8, width = 14)
# grid.arrange(fig_s1)

# write to file
ggsave('img/Fig_S1.png', fig_s1, height = 8, width = 14, units = 'in')




### -------------------------------------------------------------------------- #
### Figure S2 ---------------------------------------------------------------- #
### parametric mortality models by strain ------------------------------------ #
### -------------------------------------------------------------------------- #

# ggplot theme
tts2 <- theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.3, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .4, 0, 0, unit = 'cm')),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 10))

# subset data to best model by strain for plotting
mortality_param_plot <- mortality_param %>%
  left_join(mortality_param_best, by = 'strain') %>% 
  filter(model == best_model)

# function to format strain labels
FormatName <- function(x) {
  paste0(toupper(substring(x, 1, 1)), substring(x, 2, 5)) %>% 
    gsub('\\.', ' ', .)
}

strain_labels <- tibble(strain = sort(unique(dat$strain))) %>% 
  left_join(mortality_param_best, by = 'strain') %>% 
  mutate(lab = FormatName(strain)) %>% 
  mutate(ast = ifelse(best_model == 'Logistic', '', '***'))

# generate figure
fig_s2 <- ggplot(mortality_param_plot) +
  geom_ribbon(aes(x = t, ymin = haz_low, ymax = haz_upp), fill = '#bdc9e1', alpha = 0.7) +
  geom_line(aes(t, haz_med), size = 0.5, col = 'darkblue') +
  geom_point(data = filter(hazard_strain, n_surv_cum > 0), aes(age, haz, size = n_surv_cum), shape = 1, alpha = 0.7) +
  geom_text(data = strain_labels, aes(label = lab, 0, Inf), hjust = 0, vjust = 1.5, size = 3.6) +
  geom_text(data = strain_labels, aes(label = ast, 0, Inf), hjust = 0, vjust = 2.7, size = 4.5, fontface = 'bold') +
  coord_cartesian(xlim = c(0, 31), c(0, 1.4)) +
  scale_y_continuous(breaks = c(0, 1)) +
  scale_size(range = c(0.2, 2), name = 'N') +
  facet_wrap(~ strain, nrow = 7) +
  xlab('Age (days)') + ylab('Mortality') +
  tts2

# # print to Mac OSX device (change 'quartz' to 'window' if running Windows)
# dev.off()
# quartz(height = 6.5, width = 6)
# print(fig_s2)

# write to file
ggsave('img/Fig_S2.png', fig_s2, height = 6.5, width = 6, units = 'in')




### -------------------------------------------------------------------------- #
### Figure S3 ---------------------------------------------------------------- #
### boxplots of life history traits by strain -------------------------------- #
### -------------------------------------------------------------------------- #

## organize data for plotting
# exclude one outlier with very high shape_fecundity (noted in caption)
dat_join_boxplot <- filter(dat_join, shape_fec_strain < 3)

# subsample 1k bootstrap values per strain, rather than plotting full 50k
set.seed(532)
shape_strain_boxplot <- shape_strain %>% 
  group_by(strain) %>% 
  do(slice(., sample(1:n(), 1000))) %>% 
  ungroup()

## function to order strains based on median/mean of given variable
RankStrains <- function(data, var) {
  data %>% 
    rename_(y = var) %>% 
    group_by(strain) %>%
    summarize(med = median(y, na.rm = T), mean = mean(y, na.rm = T)) %>% 
    arrange(med, mean) %>%
    dplyr::select(strain) %>%
    unlist() %>% as.character()
}

## function to format site names for y-axis
FormatSite <- function(x) {
  paste0(toupper(substring(x,1,1)), substring(x, 2, 3))
}


## ggplot theme
tts3 <- theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 7.5),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.2, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .2, 0, 0, unit = 'cm')))

## generate all panels
ps3_1 <- dat_join %>% 
  mutate(strain = factor(strain, RankStrains(., 'lifespan'))) %>% 
  ggplot(aes(strain, lifespan)) +
  geom_boxplot(aes(fill = col), size = 0.3, outlier.size = 0.7) +
  scale_x_discrete(labels = FormatSite) +
  scale_fill_identity() +
  coord_flip() + 
  xlab(NULL) + ylab('Lifespan (days)') +
  annotate('text', x = Inf, y = -Inf, label = '(a)', hjust = -0.3, vjust = 1.5, size = 4.6) +
  tts3

ps3_2 <- shape_strain_boxplot %>% 
  mutate(strain = factor(strain, RankStrains(., 'shape_mort_boot'))) %>% 
  ggplot(aes(strain, shape_mort_boot)) +
  geom_boxplot(aes(fill = col), size = 0.3, outlier.size = 0.7) +
  scale_x_discrete(labels = FormatSite) +
  scale_fill_identity() +
  coord_flip() +
  xlab(NULL) + ylab(expression(shape[mortality])) +
  annotate('text', x = Inf, y = -Inf, label = '(b)', hjust = -0.3, vjust = 1.5, size = 4.6) +
  tts3

ps3_3 <- dat_join_boxplot %>% 
  mutate(strain = factor(strain, RankStrains(., 'shape_fec_strain'))) %>% 
  ggplot(aes(strain, shape_fec_strain)) +
  geom_boxplot(aes(fill = col), size = 0.3, outlier.size = 0.7) +
  scale_x_discrete(labels = FormatSite) +
  scale_y_continuous(breaks = seq(-2, 1, 1)) +
  scale_fill_identity() +
  coord_flip() +
  xlab(NULL) + ylab(expression(shape[fecundity])) +
  annotate('text', x = Inf, y = -Inf, label = '(c)', hjust = -0.3, vjust = 1.5, size = 4.6) +
  tts3

ps3_4 <- dat_join %>% 
  filter(!is.na(total_offspring)) %>% 
  mutate(strain = factor(strain, RankStrains(., 'total_offspring'))) %>% 
  ggplot(aes(strain, total_offspring)) +
  geom_boxplot(aes(fill = col), size = 0.3, outlier.size = 0.7) +
  scale_x_discrete(labels = FormatSite) +
  scale_y_continuous(breaks = seq(5, 20, 5)) +
  scale_fill_identity() +
  coord_flip() +
  xlab(NULL) + ylab('Cumulative fecundity') +
  annotate('text', x = Inf, y = -Inf, label = '(d)', hjust = -0.3, vjust = 1.5, size = 4.6) +
  tts3

ps3_5 <- dat_join %>% 
  filter(!is.na(area)) %>%
  mutate(strain = factor(strain, RankStrains(., 'area'))) %>% 
  ggplot(aes(strain, area)) +
  geom_boxplot(aes(fill = col), size = 0.3, outlier.size = 0.7) +
  scale_x_discrete(labels = FormatSite) +
  scale_y_continuous(breaks = c(4, 6, 8)) +
  scale_fill_identity() +
  coord_flip() +
  xlab(NULL) + ylab('Surface area (mm²)') +
  annotate('text', x = Inf, y = -Inf, label = '(e)', hjust = -0.3, vjust = 1.5, size = 4.6) +
  tts3

# convert panels to ggplot grobs
gs3_1 <- ggplotGrob(ps3_1)
gs3_2 <- ggplotGrob(ps3_2)
gs3_3 <- ggplotGrob(ps3_3)
gs3_4 <- ggplotGrob(ps3_4)
gs3_5 <- ggplotGrob(ps3_5)

# standardize panel dimensions
gs3_1$heights <- gs3_2$heights
gs3_3$heights <- gs3_2$heights
gs3_4$heights <- gs3_2$heights
gs3_5$heights <- gs3_2$heights

# generate full figure
fig_s3 <- arrangeGrob(gs3_1, gs3_2, gs3_3, gs3_4, gs3_5, nrow = 2)

# # print to Mac OSX device (change 'quartz' to 'window' if running Windows)
# dev.off()
# quartz(height = 6.5, width = 9)
# grid.arrange(fig_s3)

# write to file
ggsave('img/Fig_S3.png', fig_s3, height = 6.5, width = 9, units = 'in')




### -------------------------------------------------------------------------- #
### Figure S4 ---------------------------------------------------------------- #
### variance components (site/strain/block) by life history trait ------------ #
### -------------------------------------------------------------------------- #

# prepare data for plotting
var_strain_plot <- var_strain_df %>% 
  mutate(label = factor(label, levels = rev(label)))

var_site_plot <- var_site_df %>% 
  mutate(label = factor(label, levels = rev(label)))

# ggplot theme
tts4 <- theme_bw() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 13),
        plot.title = element_text(size = 13),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.3, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .3, 0, 0, unit = 'cm')))

## generate all panels
ps4_a <- ggplot(var_strain_plot) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbarh(aes(x = ratio_med, xmin = ratio_low, xmax = ratio_upp, y = label), height = 0) +
  geom_point(aes(ratio_med, label)) +
  annotate('text', x = Inf, y = Inf, label = '(a)', hjust = 1.5, vjust = 2, size = 4.6) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c('0.01', '0.1', '1', '10', '100')) +
  xlab(expression(over(italic(var[strain]), italic(var[block])))) + ylab(NULL) +
  coord_cartesian(xlim = c(0.07, 100)) +
  scale_y_discrete(labels = c('Surface area (all\navailable strains)',
                              'Surface area',
                              'Cumulative fecundity',
                              expression(shape[fecundity]),
                              expression(shape[mortality]),
                              'Lifespan')) +
  tts4

ps4_b <- ggplot(var_site_plot) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbarh(aes(x = ratio_med, xmin = ratio_low, xmax = ratio_upp, y = label), height = 0) +
  geom_point(aes(ratio_med, label)) +
  annotate('text', x = Inf, y = Inf, label = '(b)', hjust = 1.5, vjust = 2, size = 4.6) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c('0.01', '0.1', '1', '10', '100')) +
  xlab(expression(over(italic(var[site]), italic(var[strain])))) + ylab(NULL) +
  coord_cartesian(xlim = c(0.07, 100)) +
  scale_y_discrete(labels = NULL) +
  tts4

# generate full figure
fig_s4 <- arrangeGrob(ps4_a, ps4_b, nrow = 1, widths = c(1.35, 1))

# # print to Mac OSX device (change 'quartz' to 'window' if running Windows)
# dev.off()
# quartz(height = 5, width = 10)
# grid.arrange(fig_s4)

# write to file
ggsave('img/Fig_S4.png', fig_s4, height = 5, width = 10, units = 'in')




### -------------------------------------------------------------------------- #
### Figure S5 ---------------------------------------------------------------- #
### among-strain correlations between all pairs of life history traits ------- #
### -------------------------------------------------------------------------- #

# ggplot theme
tts5 <- theme_bw() +
  theme(axis.title = element_text(size = 13),
        plot.title = element_text(size = 13, hjust = 0.5),
        text = element_text(size = 12),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.3, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .3, 0, 0, unit = 'cm')))

## organize data
# strain-specific mean and se of all traits except shape_mort
df_means_main <- dat_join %>%
  group_by(strain) %>%
  summarize(
    mean_lifespan = mean(lifespan, na.rm = T),
    se_lifespan = sd(lifespan, na.rm = T) / sqrt(length(lifespan[!is.na(lifespan)])),
    mean_shape_fecund = mean(shape_fec_strain, na.rm = T),
    se_shape_fecund = sd(shape_fec_strain, na.rm = T) / sqrt(length(shape_fec_strain[!is.na(shape_fec_strain)])),
    mean_total_offspring = mean(total_offspring, na.rm = T),
    se_total_offspring = sd(total_offspring, na.rm = T) / sqrt(length(total_offspring[!is.na(total_offspring)])),
    mean_area = mean(area, na.rm = T),
    se_area = sd(area, na.rm = T) / sqrt(length(area[!is.na(area)]))
  )

# strain-specific mean and se of shape_mort
df_means_shape_mort <- shape_strain %>%
  group_by(strain) %>% 
  summarize(mean_shape_mort = mean(shape_mort_boot),
            se_shape_mort = sd(shape_mort_boot))

# join strain-specific mean and se for all traits
df_means <- df_means_main %>%
  full_join(df_means_shape_mort, by = 'strain')

## list of axis breaks for labels
df_breaks <- list(lifespan = c(19, 21, 23),
                  shape_mort = c(0.75, 0.82, 0.89),
                  shape_fecund = c(-0.6, -0.2, 0.2),
                  total_offspring = c(10, 12, 14),
                  area = c(4, 5, 6, 7))

## function to produce trait correlation plots
TraitCorrelations <- function(data, var_x, var_y, pos, lab_x = NULL, lab_y = NULL) {
  # organize data
  data_sub <- data[,c('strain',
                        paste0('mean_', var_x),
                        paste0('se_', var_x),
                        paste0('mean_', var_y),
                        paste0('se_', var_y))] %>% 
    setNames(c('strain', 'x_mean', 'x_se', 'y_mean', 'y_se'))
  
  # axis limits
  lim_x <- range(c(data_sub$x_mean - data_sub$x_se, data_sub$x_mean + data_sub$x_se))
  lim_y <- range(c(data_sub$y_mean - data_sub$y_se, data_sub$y_mean + data_sub$y_se))
  
  # base plot
  p <- ggplot(data_sub) +
    geom_blank() +
    xlab(lab_x) +
    ylab(lab_y) +
    coord_cartesian(xlim = lim_x, ylim = lim_y) +
    scale_x_continuous(limits = lim_x, breaks = df_breaks[[var_x]]) +
    scale_y_continuous(limits = lim_y, breaks = df_breaks[[var_y]]) + 
    tts5
  
  # if x != y, add data to plot (panels on diagonal are blank)
  if(var_x != var_y) {
    p <- p +
      geom_linerange(data = data_sub, aes(x = x_mean, ymin = y_mean - y_se, ymax = y_mean + y_se), col = 'grey50', alpha = 0.5) +
      geom_errorbarh(data = data_sub, aes(x = x_mean, y = y_mean, xmin = x_mean - x_se, xmax = x_mean + x_se), col = 'grey50', alpha = 0.5) +
      geom_point(data = data_sub, aes(x = x_mean, y = y_mean), size = 1.2) +
      geom_smooth(data = data_sub, aes(x = x_mean, y = y_mean), method = 'lm', col = 'darkblue', fill = 'darkblue', alpha = 0.16)
  }
  
  # adjust axis text depending on panel position
  # 'b' = bottom-most row, 'l' = leftmost column, 'm' = middle
  if (pos == 'l') p <- p + theme(axis.text.x = element_blank())
  if (pos == 'b') p <- p + theme(axis.text.y = element_blank())
  if (pos == 'm') p <- p + theme(axis.text = element_blank())
  
  return(list(var_x = var_x, var_y = var_y, p = p))
}

## generate all panels
trait_corr11 <- TraitCorrelations(df_means, 'lifespan', 'lifespan', pos = 'l', lab_y = 'Lifespan')
trait_corr12 <- TraitCorrelations(df_means, 'lifespan', 'shape_mort', pos = 'l', lab_y = expression(Shape[mortality]))
trait_corr13 <- TraitCorrelations(df_means, 'lifespan', 'shape_fecund', pos = 'l', lab_y = expression(Shape[fecundity]))
trait_corr14 <- TraitCorrelations(df_means, 'lifespan', 'total_offspring', pos = 'l', lab_y = 'Cumul. fecund.')
trait_corr15 <- TraitCorrelations(df_means, 'lifespan', 'area', pos = 'bl', lab_x = 'Lifespan', lab_y = 'Frond area')

trait_corr21 <- TraitCorrelations(df_means, 'shape_mort', 'lifespan', pos = 'm')
trait_corr22 <- TraitCorrelations(df_means, 'shape_mort', 'shape_mort', pos = 'm')
trait_corr23 <- TraitCorrelations(df_means, 'shape_mort', 'shape_fecund', pos = 'm')
trait_corr24 <- TraitCorrelations(df_means, 'shape_mort', 'total_offspring', pos = 'm')
trait_corr25 <- TraitCorrelations(df_means, 'shape_mort', 'area', lab_x = expression(Shape[mortality]), pos = 'b')

trait_corr31 <- TraitCorrelations(df_means, 'shape_fecund', 'lifespan', pos = 'm')
trait_corr32 <- TraitCorrelations(df_means, 'shape_fecund', 'shape_mort', pos = 'm')
trait_corr33 <- TraitCorrelations(df_means, 'shape_fecund', 'shape_fecund', pos = 'm')
trait_corr34 <- TraitCorrelations(df_means, 'shape_fecund', 'total_offspring', pos = 'm')
trait_corr35 <- TraitCorrelations(df_means, 'shape_fecund', 'area', lab_x = expression(Shape[fecundity]), pos = 'b')

trait_corr41 <- TraitCorrelations(df_means, 'total_offspring', 'lifespan', pos = 'm')
trait_corr42 <- TraitCorrelations(df_means, 'total_offspring', 'shape_mort', pos = 'm')
trait_corr43 <- TraitCorrelations(df_means, 'total_offspring', 'shape_fecund', pos = 'm')
trait_corr44 <- TraitCorrelations(df_means, 'total_offspring', 'total_offspring', pos = 'm')
trait_corr45 <- TraitCorrelations(df_means, 'total_offspring', 'area', lab_x = 'Cumul. fecund.', pos = 'b')

trait_corr51 <- TraitCorrelations(df_means, 'area', 'lifespan', pos = 'm')
trait_corr52 <- TraitCorrelations(df_means, 'area', 'shape_mort', pos = 'm')
trait_corr53 <- TraitCorrelations(df_means, 'area', 'shape_fecund', pos = 'm')
trait_corr54 <- TraitCorrelations(df_means, 'area', 'total_offspring', pos = 'm')
trait_corr55 <- TraitCorrelations(df_means, 'area', 'area', lab_x = 'Frond area', pos = 'b')

# convert panels to ggplot grobs
gs5_11 <- ggplotGrob(trait_corr11$p)
gs5_12 <- ggplotGrob(trait_corr12$p)
gs5_13 <- ggplotGrob(trait_corr13$p)
gs5_14 <- ggplotGrob(trait_corr14$p)
gs5_15 <- ggplotGrob(trait_corr15$p)

gs5_21 <- ggplotGrob(trait_corr21$p)
gs5_22 <- ggplotGrob(trait_corr22$p)
gs5_23 <- ggplotGrob(trait_corr23$p)
gs5_24 <- ggplotGrob(trait_corr24$p)
gs5_25 <- ggplotGrob(trait_corr25$p)

gs5_31 <- ggplotGrob(trait_corr31$p)
gs5_32 <- ggplotGrob(trait_corr32$p)
gs5_33 <- ggplotGrob(trait_corr33$p)
gs5_34 <- ggplotGrob(trait_corr34$p)
gs5_35 <- ggplotGrob(trait_corr35$p)

gs5_41 <- ggplotGrob(trait_corr41$p)
gs5_42 <- ggplotGrob(trait_corr42$p)
gs5_43 <- ggplotGrob(trait_corr43$p)
gs5_44 <- ggplotGrob(trait_corr44$p)
gs5_45 <- ggplotGrob(trait_corr45$p)

gs5_51 <- ggplotGrob(trait_corr51$p)
gs5_52 <- ggplotGrob(trait_corr52$p)
gs5_53 <- ggplotGrob(trait_corr53$p)
gs5_54 <- ggplotGrob(trait_corr54$p)
gs5_55 <- ggplotGrob(trait_corr55$p)

# standardize panel dimensions
gs5_11$widths <- gs5_13$widths
gs5_12$widths <- gs5_13$widths
gs5_13$widths <- gs5_13$widths
gs5_14$widths <- gs5_13$widths
gs5_15$widths <- gs5_13$widths
gs5_15$heights <- gs5_25$heights
gs5_25$heights <- gs5_25$heights
gs5_35$heights <- gs5_25$heights
gs5_45$heights <- gs5_25$heights
gs5_55$heights <- gs5_25$heights

# arrange full figure
fig_s5 <- arrangeGrob(gs5_11, gs5_21, gs5_31, gs5_41, gs5_51,
                      gs5_12, gs5_22, gs5_32, gs5_42, gs5_52,
                      gs5_13, gs5_23, gs5_33, gs5_43, gs5_53,
                      gs5_14, gs5_24, gs5_34, gs5_44, gs5_54,
                      gs5_15, gs5_25, gs5_35, gs5_45, gs5_55,
                      nrow = 5,
                      heights = c(1, 1, 1, 1, 1.35),
                      widths = c(1.3, 1, 1, 1, 1))

# # print to Mac OSX device (change 'quartz' to 'window' if running Windows)
# dev.off()
# quartz(height = 6.75, width = 8)
# grid.arrange(fig_s5)

# write to file
ggsave('img/Fig_S5.png', fig_s5, height = 6.75, width = 8, units = 'in')




### -------------------------------------------------------------------------- #
### Figure A1 ---------------------------------------------------------------- #
### illustrate calculation of shape_fecundity, part 1 ------------------------ #
### -------------------------------------------------------------------------- #

# ggplot theme
tta1 <- theme_bw() +
  theme(axis.title = element_text(size = 16),
        text = element_text(size = 15),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.3, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .3, 0, 0, unit = 'cm')))

FormatName <- function(x) {
  paste0(toupper(substring(x, 1, 1)), substring(x, 2, 5)) %>% 
    gsub('\\.', ' ', .)
}

# organize data for plotting
plot_shape_fecund <- filter(dat_tidy_std_ffr_strain, uncertain_repro == FALSE) %>%
  filter(id %in% c('C0334', 'C0481', 'C0009', 'C1140')) %>%
  mutate(strain = FormatName(strain)) %>%
  mutate(id_lab = paste0(id, ' (', strain, ')')) %>%
  mutate(age = age - first_repro,
         age_std = age_std - (first_repro/lifespan_ffr_group))

plot_labs <- fecund_slope_strain %>%
  filter(id %in% plot_shape_fecund$id) %>%
  mutate(fecund_slope = paste0('β = ', round(shape_fec_nonstd, 3))) %>%
  mutate(shape_fec_strain = paste0('β = ', round(shape_fec_strain, 2))) %>%
  left_join(dplyr::select(plot_shape_fecund, id, id_lab), by = 'id')

# generate all panels
pa1_a <- ggplot(plot_shape_fecund, aes(age, fecund)) +
  geom_point() +
  geom_smooth(method = 'lm', se = T, col = 'black') +
  geom_hline(aes(yintercept = fecund_mean_ffr), linetype = 2) +
  geom_vline(aes(xintercept = lifespan_ffr_group), linetype = 2) +
  scale_y_continuous(breaks = c(0, 1, 2)) +
  facet_wrap(~ id_lab, nrow = 1) +
  geom_text(data = plot_labs, inherit.aes = F, aes(x = 0, y = -Inf, label = fecund_slope), hjust = 0, vjust = -1, size = 4.5) +
  xlab('Age (days)') + ylab('Fecundity') +
  tta1

pa1_b <- ggplot(plot_shape_fecund, aes(age_std, fecund_std)) +
  geom_point() +
  geom_smooth(method = 'lm', se = T, col = 'black') +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_vline(aes(xintercept = 1), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 1.2, 0.2)) +
  facet_wrap(~ id_lab, nrow = 1) +
  geom_text(data = plot_labs, inherit.aes = F, aes(x = 0, y = -Inf, label = shape_fec_strain), hjust = 0, vjust = -1, size = 4.5) +
  xlab('Standardized age (life expectancies)') + ylab('Standardized fecundity') +
  tta1

# convert panels to ggplot grobs
ga1_a <- ggplotGrob(pa1_a)
ga1_b <- ggplotGrob(pa1_b)

# standardize panel dimensions
ga1_b$widths <- ga1_a$widths

# arrange full figure
fig_a1 <- arrangeGrob(ga1_a, ga1_b, nrow = 2)

# # print to Mac OSX device (change 'quartz' to 'window' if running Windows)
# dev.off()
# quartz(height = 7, width = 13)
# grid.arrange(fig_a1)

# write to file
ggsave('img/Fig_A1.png', fig_a1, height = 7, width = 13, units = 'in')




### -------------------------------------------------------------------------- #
### Figure A2 ---------------------------------------------------------------- #
### illustrate effect of pace-standardizing age in calculation of shape_fecund #
### -------------------------------------------------------------------------- #

# generate two hypothetical fecundity trajectories
a <- 100
b1 <- -14
b2 <- -7
x1 <- 0:5
x2 <- 0:10
y1 <- a + b1 * x1
y2 <- a + b2 * x2

spp_A <- data.frame(x1, y1)
spp_B <- data.frame(x2, y2)

# ggplot theme
tta2 <- theme_bw() +
  theme(axis.title = element_text(size = 17),
        plot.title = element_text(size = 13, hjust = 0.5),
        text = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.3, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .3, 0, 0, unit = 'cm')))

# generate all panels
pa2_a<- ggplot() +
  geom_line(data = spp_A, aes(x1, y1), size = 1.5) +
  geom_line(data = spp_B, aes(x2, y2), size = 1.5) +
  geom_hline(aes(yintercept = 65), linetype = 2) +
  coord_cartesian(xlim = c(0, 10.5), ylim = c(0, 100)) +
  scale_x_continuous(breaks = c(0, 5, 10)) +
  annotate('text', x = 4.5, y = 25, label = 'Species A', size = 6) +
  annotate('text', x = 9.5, y = 25, label = 'Species B', size = 6) +
  annotate('text', x = Inf, y = Inf, label = '(a)', hjust = 1.4, vjust = 1.75, size = 7.5) +
  xlab('Age (absolute time, e.g. days)') + ylab('Fecundity') +
  tta2

pa2_b <- ggplot() +
  geom_line(data = spp_A, aes(x1, y1), size = 1.5) +
  geom_hline(aes(yintercept = 65), linetype = 2) +
  coord_cartesian(xlim = c(0, 5.5), ylim = c(0, 100)) +
  scale_x_continuous(breaks = c(0, 5), labels = c(0, 1)) +
  annotate('text', x = 5, y = 25, label = 'Species A', size = 6) +
  annotate('text', x = 5, y = 20, label = 'Species B', size = 6) +
  annotate('text', x = Inf, y = Inf, label = '(b)', hjust = 1.4, vjust = 1.75, size = 7.5) +
  xlab('Age (life expectancies)') + ylab('Fecundity') +
  tta2

# arrange full figure
fig_a2 <- arrangeGrob(pa2_a, pa2_b, nrow = 1)

# # print to Mac OSX device (change 'quartz' to 'window' if running Windows)
# dev.off()
# quartz(height = 6, width = 12)
# grid.arrange(fig_a2)

# write to file
ggsave('img/Fig_A2.png', fig_a2, height = 6, width = 12, units = 'in')

