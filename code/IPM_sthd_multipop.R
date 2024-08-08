##=================================================================##
##                                                                 ##
##                 Fit integrated population model                 ##
##                                                                 ##
##=================================================================##
pacman::p_load(here,tidyverse,rstan,gtools,salmonIPM)

##===================================================## rstan options
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

##================================================## output directory
out_dir<-paste0(here(),"/output/")
if(!file.exists(out_dir)) dir.create(file.path(out_dir))
setwd(file.path(out_dir))

##========================================================## settings
covar_effects<-TRUE ## TRUE/FALSE
year_effect<-FALSE ## TRUE/FALSE
biastest<-FALSE ## TRUE/FALSE

##============================================================## data
covar_dat<-read.csv(paste0(here(),"/data/IPM_covar_dat_selected.csv"))

if(biastest){
   fish_dat<-read.csv(paste0(here(),"/data/IPM_fish_dat_biastest.csv"))%>%
      dplyr::select(-F_rate_NA) ## column used for plotting only
}else{
   fish_dat<-read.csv(paste0(here(),"/data/IPM_fish_dat_all.csv")) %>%
      dplyr::select(-F_rate_NA) ## column used for plotting only
}

## no NAs allowed in covariate data
if(covar_effects) { 
   fish_dat<-fish_dat %>%
      left_join(covar_dat,by='year') %>%
      na.omit() 
}

dat_years<-sort(unique(fish_dat$year)) ## years included
nY<-length(dat_years) ## number of years

##======================================## regressions to be included
if(covar_effects) {
   par_models<-list(s_SS~SST+pinks+av_CFS_min,R~NPGO_2+SST_4+pinks_4)
   if(year_effect) { 
      par_models<-list(s_SS~year+SST+pinks+av_CFS_min,R~year+NPGO_2+SST_4+pinks_4)
   }
} else {
   if(year_effect) { 
      par_models<-list(s_SS~year,R~year)
   }else{
      par_models<-NULL
   }
}

##=========================================================## fit IPM
IPM_fit<-salmonIPM(
   life_cycle="SSiter",
   pool_pops=TRUE,
   SR_fun="Ricker", 
   par_models=par_models, 
   fish_data=fish_dat,
   prior=list(
      mu_alpha~normal(1.5,0.5),
      mu_p~dirichlet(c(1,2,47,44,5,1)),
      mu_SS~beta(1.5,3)),
   pars=c(stan_pars("IPM_SSiter_pp")),
   chains=3, 
   iter=2000,
   warmup=1000, 
   control=list(adapt_delta=0.98,max_treedepth=12)
)

##=====================================## save model data and results
if(covar_effects) { 
   if(year_effect) { 
      write.csv(fish_dat,"IPM_fish_dat_with_covars_with_year.csv",row.names=F)
      saveRDS(IPM_fit,"IPM_fit_with_covars_with_year.Rdata")
   }else{
      if(biastest){
         write.csv(fish_dat,"IPM_fish_dat_with_covars.csv",row.names=F)
         saveRDS(IPM_fit,"IPM_fit_with_covars_biastest.Rdata")
      }else{
         write.csv(fish_dat,"IPM_fish_dat_with_covars.csv",row.names=F)
         saveRDS(IPM_fit,"IPM_fit_with_covars.Rdata")  
      }
   }
}else{ 
   if(year_effect) { 
      write.csv(fish_dat,"IPM_fish_dat_without_covars_with_year.csv",row.names=F)
      saveRDS(IPM_fit,"IPM_fit_without_covars_with_year.Rdata") 
   }else{
      write.csv(fish_dat,"IPM_fish_dat_without_covars.csv",row.names=F)
      saveRDS(IPM_fit,"IPM_fit_without_covars.Rdata") 
   }
}

##=================================================================##