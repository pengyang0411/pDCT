##-----------------------------
## Run LDT
## Author: Peng Yang
##-----------------------------
rm(list=ls(all=TRUE))
args <- commandArgs(trailingOnly = TRUE)
i = as.numeric(args[1])

# library(randomForest)
library(Metrics)
library(tsbart)
# library(mmrm)
library(ggplot2)
library(coda)
library(MASS)
library(aod)

source('source/source.R')


# sim = c(1:1000)[id]

##-----------------------
## Simulation scenarios
##-----------------------

n <- 200

sigma_onsite = c(1, 2)
d_grid <- c(0.5, 1, 3)
sigma_offsite = c(2, 3, 5)

n_p = 5
n_T = 7
ntree = 200
nsamp = 2000

df_grid <- expand.grid(
  d_grid = d_grid,
  sigma_onsite = sigma_onsite,
  sigma_offsite = sigma_offsite,
  sim = 1:1000
)
df_grid <- df_grid[df_grid$sigma_offsite > df_grid$sigma_onsite, ]

Methods <- c('True',
             'f(t) + X',
             'f(t) + X + S', 'f(t) + X + S + Z',
             'Digital Twins')

n.Methods <- length(Methods)

## For each of the grid, we running models
d_bias        <- df_grid[i, 1]
sigma_onsite  <- df_grid[i, 2]
sigma_offsite <- df_grid[i, 3]
sim           <- df_grid[i, 4]

##-----------------------------------
## Generate the data baseline
##-----------------------------------

##------------------
## Define functions
##------------------
# fx <- function(t) 4*(-t^2 + 3*t)
fx <- function(t){
  return(3 + 10/(1 + exp(- 5*(t - 0.05))) - 7.37)
}

Delta <- function(t){
  return(0 * t)
}

Gamma_bias <- function(t){
  
  # return(d_bias*t^2 + t)
  return(0 * t)
  
}

fx_a0 <- function(t, n_p) return(replicate(n_p, rep(2.5, length(t))))
fx_b0 <- function(t, n_p) return(replicate(n_p, rep(1, length(t))))
fx_a1 <- function(t, n_p) return(replicate(n_p, rep(2.5, length(t))))
fx_b1 <- function(t, n_p) return(replicate(n_p, rep(1, length(t))))

data <- sim_data(n = n,
                 n_T = n_T, n_p = n_p,
                 sig_1 = sigma_onsite, 
                 sig_2 = sigma_offsite, 
                 d_s = d_bias,
                 if_null = T, if_plot = F)

##------------------------------
## run models on full cohort
##------------------------------
res_all <- c()
## The proposed model
if(!dir.exists('Figures/Baseline/DT/')) dir.create('Figures/Baseline/DT/')
res_LDT <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                   if_covariates = F, if_digital_twins = T,
                   if_model_diagnostic = T, 
                   plot_path = paste0('Figures/Baseline/DT/sim_', sim))
res_LDT <- rbind(res_LDT$post_ATE,
                 res_LDT$Prob)

## Y = f(t) + X 
if(!dir.exists('Figures/Baseline/f_X/')) dir.create('Figures/Baseline/f_X/')
res_LDT_1 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = F,
                     if_covariates = F, if_digital_twins = F,
                     if_model_diagnostic = T, 
                     plot_path = paste0('Figures/Baseline/f_X/sim_', sim))
res_LDT_1 <- rbind(res_LDT_1$post_ATE,
                   res_LDT_1$Prob)

## Y = f(t) + X + S
if(!dir.exists('Figures/Baseline/f_X_S/')) dir.create('Figures/Baseline/f_X_S/')
res_LDT_2 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = T,
                     if_covariates = F, if_digital_twins = F,
                     if_model_diagnostic = T, 
                     plot_path = paste0('Figures/Baseline/f_X_S/sim_', sim))
res_LDT_2 <- rbind(res_LDT_2$post_ATE,
                   res_LDT_2$Prob)

## Y = f(t) + X + S + Z
if(!dir.exists('Figures/Baseline/f_X_S_Z/')) dir.create('Figures/Baseline/f_X_S_Z/')
res_LDT_3 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = T,
                     if_covariates = T, if_digital_twins = F,
                     if_model_diagnostic = T, 
                     plot_path = paste0('Figures/Baseline/f_X_S_Z/sim_', sim))
res_LDT_3 <- rbind(res_LDT_3$post_ATE,
                   res_LDT_3$Prob)


## True ATE at each time
res_True <- c()
for(t in 1:n_T){
  Indx <- seq(t, n_T * n, by = n_T)
  res_True <- c(res_True,
                mean(data$Y1[Indx] - data$Y0[Indx]))
}


meta.data <- c()
for(j in (1:n.Methods)){
  for(mm in c('ATE', 'Prob')){
    meta.data <- rbind(meta.data,  c('Baseline', 'Full', sim, d_bias, sigma_onsite, sigma_offsite, Methods[j], mm))
  }
}


res_tmp <- cbind(meta.data, 
                 rbind(res_LDT_1, res_LDT_2, res_LDT_3, res_LDT))

res_all <- rbind(res_all, res_tmp)

save(res_all, file = paste0('Res/Sim_ID_', i, '.RData'))

##-----------------------------------
## Generate the data Linear
##-----------------------------------

fx_a0 <- function(t, n_p) return(replicate(n_p, rep(0, length(t))))
fx_b0 <- function(t, n_p) return(replicate(n_p, rep(1, length(t))))
fx_a1 <- function(t, n_p) return(replicate(n_p, rep(0, length(t))))
fx_b1 <- function(t, n_p) return(replicate(n_p, rep(1, length(t))))

Gamma_bias <- function(t){
  
  return(d_bias*t^2 + t)
  # return(0)
  
}

data <- sim_data(n = n,
                 n_T = n_T, n_p = n_p,
                 sig_1 = sigma_onsite, 
                 sig_2 = sigma_offsite, 
                 d_s = d_bias,
                 if_null = T, if_plot = F)

##--------------------
## run models
##--------------------

## The proposed model
if(!dir.exists('Figures/Linear/DT/')) dir.create('Figures/Linear/DT/')
res_LDT <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                   if_covariates = F, if_digital_twins = T,
                   if_model_diagnostic = T, 
                   plot_path = paste0('Figures/Linear/DT/sim_', sim))
res_LDT <- rbind(res_LDT$post_ATE,
                 res_LDT$Prob)

## Y = f(t) + X 
if(!dir.exists('Figures/Linear/f_X/')) dir.create('Figures/Linear/f_X/')
res_LDT_1 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = F,
                     if_covariates = F, if_digital_twins = F,
                     if_model_diagnostic = T, 
                     plot_path = paste0('Figures/Linear/f_X/sim_', sim))
res_LDT_1 <- rbind(res_LDT_1$post_ATE,
                   res_LDT_1$Prob)

## Y = f(t) + X + S
if(!dir.exists('Figures/Linear/f_X_S/')) dir.create('Figures/Linear/f_X_S/')
res_LDT_2 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = T,
                     if_covariates = F, if_digital_twins = F,
                     if_model_diagnostic = T, 
                     plot_path = paste0('Figures/Linear/f_X_S/sim_', sim))
res_LDT_2 <- rbind(res_LDT_2$post_ATE,
                   res_LDT_2$Prob)

## Y = f(t) + X + S + Z
if(!dir.exists('Figures/Linear/f_X_S_Z/')) dir.create('Figures/Linear/f_X_S_Z/')
res_LDT_3 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = T,
                     if_covariates = T, if_digital_twins = F,
                     if_model_diagnostic = T, 
                     plot_path = paste0('Figures/Linear/f_X_S_Z/sim_', sim))
res_LDT_3 <- rbind(res_LDT_3$post_ATE,
                   res_LDT_3$Prob)

## True ATE at each time
res_True <- c()
for(t in 1:n_T){
  Indx <- seq(t, n_T * n, by = n_T)
  res_True <- c(res_True,
                mean(data$Y1[Indx] - data$Y0[Indx]))
}

meta.data <- c()
for(j in (1:n.Methods)){
  for(mm in c('ATE', 'Prob')){
    meta.data <- rbind(meta.data,  c('Linear', 'Full', sim, d_bias, sigma_onsite, sigma_offsite, Methods[j], mm))
  }
}

res_tmp <- cbind(meta.data, 
                 rbind(res_LDT_1, res_LDT_2, res_LDT_3, res_LDT))

res_all <- rbind(res_all, res_tmp)

save(res_all, file = paste0('Res/Sim_ID_', i, '.RData'))


##-----------------------------------
## Generate the data Linearly 
## Increased ATE
##-----------------------------------

##------------------
## Define function
##------------------
Delta <- function(t){
  return(t)
}

fx_a0 <- function(t, n_p) return(replicate(n_p, rep(2.5, length(t))))
fx_b0 <- function(t, n_p) return(replicate(n_p, rep(1, length(t))))
fx_a1 <- function(t, n_p) return(replicate(n_p, rep(2.5, length(t))))
fx_b1 <- function(t, n_p) return(replicate(n_p, rep(1, length(t))))

data <- sim_data(n = n,
                 n_T = n_T, n_p = n_p,
                 sig_1 = sigma_onsite, 
                 sig_2 = sigma_offsite, 
                 d_s = d_bias,
                 if_null = F, if_plot = F)


##--------------------
## run models
##--------------------

## The proposed model
if(!dir.exists('Figures/Linear_Increased_ATE/DT/')) dir.create('Figures/Linear_Increased_ATE/DT/')
res_LDT <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                   if_covariates = F, if_digital_twins = T,
                   if_model_diagnostic = T, 
                   plot_path = paste0('Figures/Linear_Increased_ATE/DT/sim_', sim))
res_LDT <- rbind(res_LDT$post_ATE,
                 res_LDT$Prob)

## Y = f(t) + X 
if(!dir.exists('Figures/Linear_Increased_ATE/f_X/')) dir.create('Figures/Linear_Increased_ATE/f_X/')
res_LDT_1 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = F,
                     if_covariates = F, if_digital_twins = F,
                     if_model_diagnostic = T, 
                     plot_path = paste0('Figures/Linear_Increased_ATE/f_X/sim_', sim))
res_LDT_1 <- rbind(res_LDT_1$post_ATE,
                   res_LDT_1$Prob)

## Y = f(t) + X + S
if(!dir.exists('Figures/Linear_Increased_ATE/f_X_S/')) dir.create('Figures/Linear_Increased_ATE/f_X_S/')
res_LDT_2 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = T,
                     if_covariates = F, if_digital_twins = F,
                     if_model_diagnostic = T, 
                     plot_path = paste0('Figures/Linear_Increased_ATE/f_X_S/sim_', sim))
res_LDT_2 <- rbind(res_LDT_2$post_ATE,
                   res_LDT_2$Prob)

## Y = f(t) + X + S + Z
if(!dir.exists('Figures/Linear_Increased_ATE/f_X_S_Z/')) dir.create('Figures/Linear_Increased_ATE/f_X_S_Z/')
res_LDT_3 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = T,
                     if_covariates = T, if_digital_twins = F,
                     if_model_diagnostic = T, 
                     plot_path = paste0('Figures/Linear_Increased_ATE/f_X_S_Z/sim_', sim))
res_LDT_3 <- rbind(res_LDT_3$post_ATE,
                   res_LDT_3$Prob)



## True ATE at each time
res_True <- c()
for(t in 1:n_T){
  Indx <- seq(t, n_T * n, by = n_T)
  res_True <- c(res_True,
                mean(data$Y1[Indx] - data$Y0[Indx]))
}

meta.data <- c()
for(j in (1:n.Methods)){
  for(mm in c('ATE', 'Prob')){
    meta.data <- rbind(meta.data,  c('L_Increased_ATE', 'Full', sim, d_bias, sigma_onsite, sigma_offsite, Methods[j], mm))
  }
}

res_tmp <- cbind(meta.data, 
                 rbind(res_LDT_1, res_LDT_2, res_LDT_3, res_LDT,
                       res_LDT_onsite_all, res_LDT_onsite_1_all))


res_all <- rbind(res_all, res_tmp)

save(res_all, file = paste0('Res/Sim_ID_', i, '.RData'))

##-----------------------------------
## Generate the data none-linearly
## Increased ATE
##-----------------------------------

##------------------
## Define function
##------------------

Delta <- function(t){
  return(1.5*log(t + 1))
}

data <- sim_data(n = n,
                 n_T = n_T, n_p = n_p,
                 sig_1 = sigma_onsite, 
                 sig_2 = sigma_offsite, 
                 d_s = d_bias,
                 if_null = F, if_plot = F)


##--------------------
## run models
##--------------------

## The proposed model
if(!dir.exists('Figures/nonLinear_Increased_ATE/DT/')) dir.create('Figures/nonLinear_Increased_ATE/DT/')
res_LDT <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                   if_covariates = F, if_digital_twins = T,
                   if_model_diagnostic = T, 
                   plot_path = paste0('Figures/nonLinear_Increased_ATE/DT/sim_', sim))
res_LDT <- rbind(res_LDT$post_ATE,
                 res_LDT$Prob)

## Y = f(t) + X 
if(!dir.exists('Figures/nonLinear_Increased_ATE/f_X/')) dir.create('Figures/nonLinear_Increased_ATE/f_X/')
res_LDT_1 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = F,
                     if_covariates = F, if_digital_twins = F,
                     if_model_diagnostic = T, 
                     plot_path = paste0('Figures/nonLinear_Increased_ATE/f_X/sim_', sim))
res_LDT_1 <- rbind(res_LDT_1$post_ATE,
                   res_LDT_1$Prob)

## Y = f(t) + X + S
if(!dir.exists('Figures/nonLinear_Increased_ATE/f_X_S/')) dir.create('Figures/nonLinear_Increased_ATE/f_X_S/')
res_LDT_2 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = T,
                     if_covariates = F, if_digital_twins = F,
                     if_model_diagnostic = T, 
                     plot_path = paste0('Figures/nonLinear_Increased_ATE/f_X_S/sim_', sim))
res_LDT_2 <- rbind(res_LDT_2$post_ATE,
                   res_LDT_2$Prob)

## Y = f(t) + X + S + Z
if(!dir.exists('Figures/nonLinear_Increased_ATE/f_X_S_Z/')) dir.create('Figures/nonLinear_Increased_ATE/f_X_S_Z/')
res_LDT_3 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = T,
                     if_covariates = T, if_digital_twins = F,
                     if_model_diagnostic = T, 
                     plot_path = paste0('Figures/nonLinear_Increased_ATE/f_X_S_Z/sim_', sim))
res_LDT_3 <- rbind(res_LDT_3$post_ATE,
                   res_LDT_3$Prob)



## True ATE at each time
res_True <- c()
for(t in 1:n_T){
  Indx <- seq(t, n_T * n, by = n_T)
  res_True <- c(res_True,
                mean(data$Y1[Indx] - data$Y0[Indx]))
}


meta.data <- c()
for(j in (1:n.Methods)){
  for(mm in c('ATE', 'Prob')){
    meta.data <- rbind(meta.data,  c('NL_Increased_ATE', 'Full', sim, d_bias, sigma_onsite, sigma_offsite, Methods[j], mm))
  }
}

res_tmp <- cbind(meta.data, 
                 rbind(res_LDT_1, res_LDT_2, res_LDT_3, res_LDT))


res_all <- rbind(res_all, res_tmp)

save(res_all, file = paste0('Res/Sim_ID_', i, '.RData'))

##-----------------------------------
## Generate the data Time-varying
## coefficient
##-----------------------------------

##------------------
## Define function
##------------------

fx_a0 <- function(t, n_p) return(replicate(n_p, seq(0, 2.5, length.out = n_T)))
fx_b0 <- function(t, n_p) return(replicate(n_p, seq(0, 1, length.out = n_T)))
fx_a1 <- function(t, n_p) return(replicate(n_p, seq(0, 2.5, length.out = n_T)))
fx_b1 <- function(t, n_p) return(replicate(n_p, seq(0, 1, length.out = n_T)))

data <- sim_data(n = n,
                 n_T = n_T, n_p = n_p,
                 sig_1 = sigma_onsite, 
                 sig_2 = sigma_offsite, 
                 d_s = d_bias,
                 if_null = F, if_plot = F)

##--------------------
## run models
##--------------------

## The proposed model
if(!dir.exists('Figures/Time_Varying_Coefficient/DT/')) dir.create('Figures/Time_Varying_Coefficient/DT/')
res_LDT <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                   if_covariates = F, if_digital_twins = T,
                   if_model_diagnostic = T, 
                   plot_path = paste0('Figures/Time_Varying_Coefficient/DT/sim_', sim))
res_LDT <- rbind(res_LDT$post_ATE,
                 res_LDT$Prob)

## Y = f(t) + X 
if(!dir.exists('Figures/Time_Varying_Coefficient/f_X/')) dir.create('Figures/Time_Varying_Coefficient/f_X/')
res_LDT_1 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = F,
                     if_covariates = F, if_digital_twins = F, 
                     plot_path = paste0('Figures/Time_Varying_Coefficient/f_X/sim_', sim))
res_LDT_1 <- rbind(res_LDT_1$post_ATE,
                   res_LDT_1$Prob)

## Y = f(t) + X + S
if(!dir.exists('Figures/Time_Varying_Coefficient/f_X_S/')) dir.create('Figures/Time_Varying_Coefficient/f_X_S/')
res_LDT_2 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = T,
                     if_covariates = F, if_digital_twins = F, 
                     plot_path = paste0('Figures/Time_Varying_Coefficient/f_X_S/sim_', sim))
res_LDT_2 <- rbind(res_LDT_2$post_ATE,
                   res_LDT_2$Prob)

## Y = f(t) + X + S + Z
if(!dir.exists('Figures/Time_Varying_Coefficient/f_X_S_Z/')) dir.create('Figures/Time_Varying_Coefficient/f_X_S_Z/')
res_LDT_3 <- run_LDT(data = data, nsamp = nsamp, if_plot = F, 
                     if_site = T,
                     if_covariates = T, if_digital_twins = F, 
                     plot_path = paste0('Figures/Time_Varying_Coefficient/f_X_S_Z/sim_', sim))
res_LDT_3 <- rbind(res_LDT_3$post_ATE,
                   res_LDT_3$Prob)


## True ATE at each time
res_True <- c()
for(t in 1:n_T){
  Indx <- seq(t, n_T * n, by = n_T)
  res_True <- c(res_True,
                mean(data$Y1[Indx] - data$Y0[Indx]))
}


meta.data <- c()
for(j in (1:n.Methods)){
  for(mm in c('ATE', 'Prob')){
    meta.data <- rbind(meta.data,  c('Time_varying', 'Full', sim, d_bias, sigma_onsite, sigma_offsite, Methods[j], mm))
  }
}

res_tmp <- cbind(meta.data, 
                 rbind(res_LDT_1, res_LDT_2, res_LDT_3, res_LDT))

res_all <- rbind(res_all, res_tmp)

save(res_all, file = paste0('Res/Sim_ID_', i, '.RData'))




