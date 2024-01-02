##=================================================================##
##                                                                 ##
##          Compare models with vs. without year effect            ##
##                                                                 ##
##=================================================================##

##========================================================## packages
pkg<-c("here","tidyverse","cowplot","ggsidekick","gridExtra","RColorBrewer","salmonIPM","viridis","posterior")
if(length(setdiff(pkg,rownames(installed.packages())))>0){install.packages(setdiff(pkg,rownames(installed.packages())),dependencies=T)}
invisible(lapply(pkg,library,character.only=T))
theme_set(theme_sleek())
home<-here::here()

##================================================## output directory
out_dir<-paste0(home,"/output/")
if(!file.exists(out_dir)) dir.create(file.path(out_dir))
setwd(file.path(out_dir))

##========================================================## IPM fits

##-------------------------------## without covariates no year effect
IPM_fit_without_covars<-readRDS(paste0(home,"/output/IPM_fit_without_covars.Rdata"))

##----------------------------------## with covariates no year effect
IPM_fit_with_covars<-readRDS(paste0(home,"/output/IPM_fit_with_covars.Rdata"))

##-----------------------------## without covariates with year effect
IPM_fit_without_covars_with_year<-readRDS(paste0(home,"/output/IPM_fit_without_covars_with_year.Rdata"))

##--------------------------------## with covariates with year effect
IPM_fit_with_covars_with_year<-readRDS(paste0(home,"/output/IPM_fit_with_covars_with_year.Rdata"))

##=======================================================## fish data

##-------------------------------------------------## with covariates 
fish_dat_with_covars<-read.csv(paste0(home,"/output/IPM_fish_dat_with_covars.csv"))
dat_years_with_covars<-sort(unique(fish_dat_with_covars$year))

##----------------------------------------------## without covariates 
fish_dat_without_covars<-read.csv(paste0(home,"/output/IPM_fish_dat_without_covars.csv"))
dat_years_without_covars<-sort(unique(fish_dat_without_covars$year))

##================================================## model comparison
probs<-c(0.05,0.25,0.5,0.75,0.95)

##---------------------------------------## log recruitment anomalies

eta_R_post1<-extract1(IPM_fit_without_covars,"eta_year_R")
eta_R_1<-t(apply(eta_R_post1,2,function(x)quantile(x,prob=probs)))%>%
   data.frame() %>% dplyr::select(median_eta_base=X50.) %>%
   add_column(year=dat_years_without_covars) %>%
   mutate(year=as.numeric(year))

eta_R_post2<-extract1(IPM_fit_without_covars_with_year,"eta_year_R")
eta_R_2<-t(apply(eta_R_post2,2,function(x)quantile(x,prob=probs)))%>%
   data.frame() %>% dplyr::select(median_eta_year=X50.) %>%
   add_column(year=dat_years_without_covars) %>%
   mutate(year=as.numeric(year))

eta_R_post3<-extract1(IPM_fit_with_covars,"eta_year_R")
eta_R_3<-t(apply(eta_R_post3,2,function(x)quantile(x,prob=probs)))%>%
   data.frame() %>% dplyr::select(median_eta_covars=X50.) %>%
   add_column(year=dat_years_with_covars) %>%
   mutate(year=as.numeric(year))

eta_R_post4<-extract1(IPM_fit_with_covars_with_year,"eta_year_R")
eta_R_4<-t(apply(eta_R_post4,2,function(x)quantile(x,prob=probs)))%>%
   data.frame() %>% dplyr::select(median_eta_covars_year=X50.) %>%
   add_column(year=dat_years_with_covars) %>%
   mutate(year=as.numeric(year))

eta_df<-eta_R_1 %>%
   left_join(eta_R_2) %>%
   left_join(eta_R_3) %>%
   left_join(eta_R_4) %>%
   drop_na() %>%
   data.frame() %>%
   dplyr::select(-year)

##-----------------------------------## logit kelt survival anomalies

eta_SS_post1<-extract1(IPM_fit_without_covars,"eta_year_SS")
eta_SS_1<-t(apply(eta_SS_post1,2,function(x)quantile(x,prob=probs)))%>%
   data.frame() %>% dplyr::select(median_eta_base=X50.) %>%
   add_column(year=dat_years_without_covars) %>%
   mutate(year=as.numeric(year))

eta_SS_post2<-extract1(IPM_fit_without_covars_with_year,"eta_year_SS")
eta_SS_2<-t(apply(eta_SS_post2,2,function(x)quantile(x,prob=probs)))%>%
   data.frame() %>% dplyr::select(median_eta_year=X50.) %>%
   add_column(year=dat_years_without_covars) %>%
   mutate(year=as.numeric(year))

eta_SS_post3<-extract1(IPM_fit_with_covars,"eta_year_SS")
eta_SS_3<-t(apply(eta_SS_post3,2,function(x)quantile(x,prob=probs)))%>%
   data.frame() %>% dplyr::select(median_eta_covars=X50.) %>%
   add_column(year=dat_years_with_covars) %>%
   mutate(year=as.numeric(year))

eta_SS_post4<-extract1(IPM_fit_with_covars_with_year,"eta_year_SS")
eta_SS_4<-t(apply(eta_SS_post4,2,function(x)quantile(x,prob=probs)))%>%
   data.frame() %>% dplyr::select(median_eta_covars_year=X50.) %>%
   add_column(year=dat_years_with_covars) %>%
   mutate(year=as.numeric(year))

eta_df<-eta_SS_1 %>%
   left_join(eta_SS_2) %>%
   left_join(eta_SS_3) %>%
   left_join(eta_SS_4) %>%
   drop_na() %>%
   data.frame() %>%
   dplyr::select(-year)

##=================================================================##