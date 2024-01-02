##=================================================================##
##                                                                 ##
##                   Prepare supplementary tables                  ##
##                                                                 ##
##=================================================================##

##========================================================## packages
pkg<-c("here","tidyverse","salmonIPM","posterior")
if(length(setdiff(pkg,rownames(installed.packages())))>0){install.packages(setdiff(pkg,rownames(installed.packages())),dependencies=T)}
invisible(lapply(pkg,library,character.only=T))
home<-here::here()

##================================================## output directory
out_dir<-paste0(home,"/output/")
if(!file.exists(out_dir)) dir.create(file.path(out_dir))
setwd(file.path(out_dir))

##========================================## load data and model fit
covar_dat<-read.csv("IPM_covar_dat_all.csv")
IPM_fit<-readRDS(paste0(home,"/output/IPM_fit_with_covars.Rdata"))
fish_dat<-read.csv(paste0(home,"/output/IPM_fish_dat_with_covars.csv"))
pops<-unique(fish_dat$pop)

##=================================================## hyperparameters
cis<-c(0.05,0.5,0.95) 

para_names<-c("mu_alpha","sigma_alpha","mu_Rmax","sigma_Rmax", "rho_alphaRmax","rho_R","sigma_year_R","sigma_R","mu_SS","rho_SS", "sigma_year_SS","sigma_SS","tau") 

post<-as_draws_rvars(IPM_fit) %>%
   mutate_variables(
      mu_alpha=exp(as.vector(mu_alpha)),
      sigma_alpha=exp(as.vector(sigma_alpha)),
      mu_Rmax=exp(as.vector(mu_Rmax)),
      sigma_Rmax=exp(as.vector(sigma_Rmax))
   ) %>%
   purrr::keep(names(.) %in% para_names) %>%
   as_draws_matrix() %>%
   as.data.frame() 

quants<-data.frame(t(apply(post,2,function(x) quantile(x,prob=cis))))%>%
   mutate_all(round,2) %>%
   rename_with(~c("lower","median","upper"))

write.table(quants,"IPM_table_S1.csv",row.names=T,sep=",")

##============================## population productivity and capacity

alphas_post<-data.frame(extract1(IPM_fit,"alpha"))
names(alphas_post)<-pops
alpha_qs<-apply(alphas_post,2,function(x) quantile(x,prob=cis))

Rmax_post<-data.frame(extract1(IPM_fit,"Rmax"))
names(Rmax_post)<-pops
Rmax_qs<-apply(Rmax_post,2,function(x) quantile(x,prob=cis))

t1<-data.frame(t(alpha_qs)) %>%
   dplyr::select(lower=X5.,median=X50.,upper=X95.) %>%
   mutate_all(round,2)

t2<-data.frame(t(Rmax_qs)) %>%
   dplyr::select(lower=X5.,median=X50.,upper=X95.) %>%
   mutate_all(round,2)

tab<-t1 %>% cbind(t2)

write.table(tab,"IPM_table_S2.csv",row.names=T,sep=",")

##=================================================================##
##=================================================================##
##=================================================================##