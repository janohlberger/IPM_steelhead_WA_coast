##=================================================================##
##                                                                 ##
## Load coastal steelhead population data and fit hierarchical IPM ##
##                                                                 ##
##=================================================================##

##========================================================## packages
pkg<-c("here","dplyr","tidyverse","rstan","readr","readxl","tibble", "dataRetrieval","posterior","ggsidekick","RColorBrewer","ggplot2", "officer","MuMIn","ncdf4","reshape2","pracma","relaimpo","visreg", "Hmisc","bayesdfa","MARSS","faraway","gtools","gridExtra","gsl","rcartocolor","bayesplot","rstanarm","distributional")
if(length(setdiff(pkg,rownames(installed.packages())))>0){install.packages(setdiff(pkg,rownames(installed.packages())),dependencies=T)}
invisible(lapply(pkg,library,character.only=T))
## devtools::install_github("ebuhle/salmonIPM@dev",auth_token=PAT)
library(salmonIPM)

##========================================================## settings
theme_set(theme_sleek())
options(rstudio.help.showDataPreview=FALSE)

##===================================================## rstan options
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

##================================================## output directory
out_dir<-paste0(home,"/output/")
if(!file.exists(out_dir)) dir.create(file.path(out_dir))
setwd(file.path(out_dir))

##======================================================## read data
fish_dat<-read.csv("IPM_fish_dat_all.csv")
covar_dat<-read.csv("IPM_covar_dat_all.csv")

##=================================================================##
##=====================================## Integrated Population Model
##=================================================================##
covar_effects<-TRUE ## TRUE or FALSE
##-----------------------------------## stock-recruitment function
SR_mod<-"Ricker" 
##-------------------------------------------------------## priors
mu_p_vec<-c(0.01,0.02,0.47,0.44,0.05,0.01)
priors<- list(
   mu_alpha ~ normal(1.5,0.5), ## hyper mean/SD of log productivity
   mu_p ~ dirichlet(mu_p_vec*100), ## mean maiden age distribution
   # tau ~ gnormal(1,0.85,30), ## default: gnormal(1,0.85,30))
   mu_SS ~ beta(1.5,3) ## mean kelt survival rate
)
##-------------------------------------------## add covariate data
if(covar_effects) { ## no NAs in covariates
   fish_dat<-fish_dat %>%
      left_join(pinks,by='year') %>%
      left_join(sst,by='year') %>%
      left_join(sst_cst,by='year') %>%
      left_join(npgo,by='year') %>%
      # left_join(seals,by='year') %>% ## ONLY FOR TESTS (FEW DATA)
      na.omit()
   par_models<-list(s_SS~SST+pinks,R~NPGO_2+SST_4+pinks_4)
} else {
   par_models<-NULL
}
##-----------------------------------------------## years included
dat_years<-sort(unique(fish_dat$year))
nY<-length(dat_years)
##---------------------------------## save fish data used in model
if(covar_effects) { 
   write.csv(fish_dat,"IPM_fish_dat_with_covars.csv",row.names=F)
}else{ 
   write.csv(fish_dat,"IPM_fish_dat_without_covars.csv",row.names=F)
}
##----------------------------------------------------## fit model
IPM_fit <- salmonIPM(
   life_cycle="SSiter",
   pool_pops=TRUE,
   SR_fun=SR_mod, 
   par_models=par_models, 
   fish_data=fish_dat,
   prior=priors,
   pars=c(stan_pars("IPM_SSiter_pp")),
   chains=3, 
   iter=1000, ## tests: 1000 | final: 2000
   warmup=500, ## tests: 500 | final: 1000
   control=list(
      adapt_delta=0.95,  ## tests: 0.95 | final: 0.98
      max_treedepth=10 ## tests: 10 | final: 12 
   )
)
##-------------------------------------------------## time elapsed
print(paste(round(max(rowSums(get_elapsed_time(IPM_fit))/60)),"min"))
##-------------------------------------------------## save results
if(covar_effects) { 
   saveRDS(IPM_fit,"IPM_fit_with_covars.Rdata") 
}else{ 
   saveRDS(IPM_fit,"IPM_fit_without_covars.Rdata") 
}

##=================================================================##
##=================================================================##
##=================================================================##