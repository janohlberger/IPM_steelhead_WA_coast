##=================================================================##
##                                                                 ##
##    Fit hierarchical IPM for WA coastal steelhead populations    ##
##                                                                 ##
##=================================================================##

##========================================================## packages
pkg<-c("here","dplyr","tidyverse","rstan","readr","readxl","tibble", "dataRetrieval","posterior","ggsidekick","RColorBrewer","ggplot2", "officer","MuMIn","ncdf4","reshape2","pracma","relaimpo","visreg", "Hmisc","bayesdfa","MARSS","faraway","gtools","gridExtra","gsl","rcartocolor","bayesplot","rstanarm","distributional","salmonIPM")
if(length(setdiff(pkg,rownames(installed.packages())))>0){install.packages(setdiff(pkg,rownames(installed.packages())),dependencies=T)}
invisible(lapply(pkg,library,character.only=T))
home<-here::here()

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

##=======================================================## read data
fish_dat<-read.csv("IPM_fish_dat_all.csv")
covar_dat<-read.csv("IPM_covar_dat_all.csv")

##=========================================================## fit IPM
covar_effects<-TRUE ## TRUE or FALSE
year_effect<-TRUE ## TRUE or FALSE
SR_mod<-"Ricker" ## SR function (Ricker/BH)
##----------------------------------------------------------## priors
priors<- list(
   ## hyper mean and SD of log productivity
   mu_alpha ~ normal(1.5,0.5), 
   ## mean maiden spawner age distribution
   mu_p ~ dirichlet(c(1,2,47,44,5,1)), 
   ## mean kelt survival rate across years
   mu_SS ~ beta(1.5,3) 
)
##----------------------------------------------## add covariate data
if(covar_effects) { ## no NAs in covariates
   fish_dat<-fish_dat %>%
      left_join(covar_dat,by='year') %>%
      na.omit()
   par_models<-list(s_SS~NPGO+SST+pinks,
                    R~NPGO_2+SST_4+pinks_4+av_CFS_min_1)
   if(year_effect) { 
      par_models<-list(s_SS~year+NPGO+SST+pinks,
                       R~year+NPGO_2+SST_4+pinks_4+av_CFS_min_1)
   }
} else {
   if(year_effect) { 
      par_models<-list(s_SS~year,R~year)
   }else{
      par_models<-NULL
   }
}
##--------------------------------------------------## years included
dat_years<-sort(unique(fish_dat$year))
nY<-length(dat_years)
##-----------------------------------------## save data used in model
if(covar_effects) { 
   if(year_effect) { 
      write.csv(fish_dat,"IPM_fish_dat_with_covars_with_year.csv",
                row.names=F)
   }else{
      write.csv(fish_dat,"IPM_fish_dat_with_covars.csv",
                row.names=F)
   }
}else{ 
   if(year_effect) { 
      write.csv(fish_dat,"IPM_fish_dat_without_covars_with_year.csv",
                row.names=F)
   }else{
      write.csv(fish_dat,"IPM_fish_dat_without_covars.csv",
                row.names=F)
   }
}
##-------------------------------------------------------## fit model
IPM_fit<-salmonIPM(
   life_cycle="SSiter",
   pool_pops=TRUE,
   SR_fun=SR_mod, 
   par_models=par_models, 
   fish_data=fish_dat,
   prior=priors,
   pars=c(stan_pars("IPM_SSiter_pp")),
   chains=3, 
   iter=2000,
   warmup=1000, 
   control=list(adapt_delta=0.95,max_treedepth=10)
) ## adapt_delta/max_treedepth: tests=0.95/10 | final=0.99/12
##----------------------------------------------------## time elapsed
print(paste(round(max(rowSums(get_elapsed_time(IPM_fit))/60)),"min"))
##----------------------------------------------------## save results
if(covar_effects) { 
   if(year_effect) { 
      saveRDS(IPM_fit,"IPM_fit_with_covars_with_year.Rdata")
   }else{
      saveRDS(IPM_fit,"IPM_fit_with_covars.Rdata")
   }
}else{ 
   if(year_effect) { 
      saveRDS(IPM_fit,"IPM_fit_without_covars_with_year.Rdata") 
   }else{
      saveRDS(IPM_fit,"IPM_fit_without_covars.Rdata") 
   }
}

##=================================================================##
##=================================================================##
##=================================================================##