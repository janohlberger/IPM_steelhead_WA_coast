##=================================================================##
##                                                                 ##
##          Compare models with vs. without year effect            ##
##                                                                 ##
##=================================================================##

##========================================================## packages
pkg<-c("here","tidyverse","salmonIPM") ## ,"posterior"
if(length(setdiff(pkg,rownames(installed.packages())))>0){install.packages(setdiff(pkg,rownames(installed.packages())),dependencies=T)}
invisible(lapply(pkg,library,character.only=T))
theme_set(theme_sleek())
home<-here::here()

##======================================## IPM fits with year effects

IPM_fit_without_covars_with_year<-readRDS(paste0(home,"/output/IPM_fit_without_covars_with_year.Rdata"))

IPM_fit_with_covars_with_year<-readRDS(paste0(home,"/output/IPM_fit_with_covars_with_year.Rdata"))

##=====================## compare year effects on recruitment anomaly

beta_R_yr_qs<-extract1(IPM_fit_without_covars_with_year,"beta_R")%>%
   apply(.,2,function(x) quantile(x,prob=c(0.05,0.5,0.95))) %>%
   data.frame() %>% mutate(across(everything(),\(x) round(x,4))) %>%
   rename_with(~c("Year")) 
## strong and significant year effect

beta_R_cov_yr_qs<-extract1(IPM_fit_with_covars_with_year,"beta_R")%>%
   apply(.,2,function(x) quantile(x,prob=c(0.05,0.5,0.95))) %>%
   data.frame() %>% mutate(across(everything(),\(x) round(x,4))) %>%
   rename_with(~c("Year","NPGO","SST","Pinks","Flow")) 
## weaker and non-significant year effect

##===================## compare year effects on kelt survival anomaly

beta_SS_yr_qs<-extract1(IPM_fit_without_covars_with_year,"beta_SS")%>%
   apply(.,2,function(x) quantile(x,prob=c(0.05,0.5,0.95))) %>%
   data.frame() %>% mutate(across(everything(),\(x) round(x,4))) %>%
   rename_with(~c("Year")) 
## strong and significant year effect

beta_SS_cov_yr_qs<-extract1(IPM_fit_with_covars_with_year,"beta_SS")%>%
   apply(.,2,function(x) quantile(x,prob=c(0.05,0.5,0.95))) %>%
   data.frame() %>% mutate(across(everything(),\(x) round(x,4))) %>%
   rename_with(~c("Year","NPGO","SST","Pinks")) 
## weaker but still significant year effect

##=================================================================##