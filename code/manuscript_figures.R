##=================================================================##
##                                                                 ##
##                    Plot manuscript figures                      ##
##                                                                 ##
##=================================================================##

##========================================================## packages
pkg<-c("here","tidyverse","cowplot","ggsidekick","gridExtra","salmonIPM","posterior")
if(length(setdiff(pkg,rownames(installed.packages())))>0){install.packages(setdiff(pkg,rownames(installed.packages())),dependencies=T)}
invisible(lapply(pkg,library,character.only=T))
theme_set(theme_sleek())
home<-here::here()

##================================================## output directory
out_dir<-paste0(home,"/output/")
if(!file.exists(out_dir)) dir.create(file.path(out_dir))
setwd(file.path(out_dir))

##=================================================================##
##========================================## load data and model fit
##=================================================================##
# covar_dat<-read.csv(paste0(home,"/data/IPM_covar_dat_all.csv"))
covar_dat<-read.csv(paste0(home,"/data/IPM_covar_dat_selected.csv"))

##========================================================## IPM fits
##----------------------------------------------## without covariates
file1<-"IPM_fit_without_covars.Rdata"
IPM_fit_without_covars<-readRDS(paste0(home,"/output/",file1))
##-------------------------------------------------## with covariates
file2<-"IPM_fit_with_covars.Rdata"
IPM_fit_with_covars<-readRDS(paste0(home,"/output/",file2))

##=======================================================## fish data
##----------------------------------------------## without covariates
file3<-"IPM_fish_dat_without_covars.csv"
fish_dat_without_covars<-read.csv(paste0(home,"/output/",file3))
dat_years_without_covars<-sort(unique(fish_dat_without_covars$year))
##-------------------------------------------------## with covariates
file4<-"IPM_fish_dat_with_covars.csv"
fish_dat_with_covars<-read.csv(paste0(home,"/output/",file4))
dat_years_with_covars<-sort(unique(fish_dat_with_covars$year))

##=====================================## fit and data for main plots
IPM_fit_plot<-IPM_fit_with_covars
fish_dat<-fish_dat_with_covars
pops<-unique(fish_dat$pop)
nP<-length(pops)
N<-dim(fish_dat)[1]

##===================================================## plot settings
colors<-rev(c("#A86260","#A5B1C4","#3B4D6B","#89A18D","#747260","#B3872D","#774C2C")) 
colors<-colorRampPalette(colors)(nP)
##-------------------------------------------------------## functions
summary_CI90<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.05),na.rm=T),ymax=quantile(x,prob=c(0.95),na.rm=T))) }
summary_CI50<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.25),na.rm=T),ymax=quantile(x,prob=c(0.75),na.rm=T))) }
##-------------------------------------------------------------## CIs
probs<-c(0.05,0.25,0.5,0.75,0.95) ## median with 50% and 95% CIs

##=======================================================## posterior
df_post<-as.data.frame(IPM_fit_plot) ## full posterior
nS<-dim(df_post)[1] ## number of samples
df_out<-data.frame(summary(IPM_fit_plot)$summary) ## summary 

##===============================================## observation error
tau_post<-extract1(IPM_fit_plot,"tau")
tau_quants<-quantile(tau_post,prob=probs)

##=================================## variance explained by covariates
## calculate variance explained (R2) by covariates by dividing sigma2 for process errors with covariates by the sigma2 without covariates
## except that the model is fit to slightly different time series due to the required lag of covariates >> 'approximate variance explained'?

##------------------------------------## shared recruitment anomalies
sigma_R_wo_cov<-extract1(IPM_fit_without_covars,"sigma_year_R")
sigma_R_w_cov<-extract1(IPM_fit_with_covars,"sigma_year_R")
sigma_R_median_ratio<-median(sigma_R_w_cov^2/sigma_R_wo_cov^2) 
sigma_R_ratio_median<-median(sigma_R_w_cov)^2/median(sigma_R_wo_cov)^2

##----------------------------------## shared kelt survival anomalies
sigma_SS_wo_cov<-extract1(IPM_fit_without_covars,"sigma_year_SS")
sigma_SS_w_cov<-extract1(IPM_fit_with_covars,"sigma_year_SS")
sigma_SS_median_ratio<-median(sigma_SS_w_cov^2/sigma_SS_wo_cov^2) 
sigma_SS_ratio_median<-median(sigma_SS_w_cov)^2/median(sigma_SS_wo_cov)^2

##=================================================================##
##=================## Figure 2 - spawners, recruitment, kelt survival
##=================================================================##

##==============================================## estimated spawners
S_post<-extract1(IPM_fit_plot,"S")

S_qs<-t(apply(S_post,2,function(x) quantile(x,prob=probs)))
S_est_qs<-S_qs %>% data.frame() %>%
   add_column(fish_dat%>%dplyr::select(pop,year,S_obs)%>%data_frame()) 

## posterior predictive distribution of spawners
quants<-c(0.05,0.5,0.95)
S_draws<-as_draws_rvars(IPM_fit_plot) %>%
   mutate_variables(S_ppd=rvar_rng(rlnorm,length(S),log(S),tau))
S_ppd<-summarise_draws(S_draws$S_ppd,~quantile(.x,probs=quants))

S_all<-S_est_qs %>% add_column(S_ppd)

S_all<-S_all %>% mutate(TF=ifelse(between(S_obs,`5%`,`95%`),T,F)) 

pTRUE<-S_all %>% filter(TF==TRUE) %>% nrow() / nrow(S_all) * 100

p1 <- S_all %>%
   ggplot(aes(x=year,y=X50.,color=pop,fill=pop)) +
   ## estimated states
   geom_line(aes(y=X50.),lwd=0.5,alpha=1) + 
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.4) +
   ## posterior predictive distribution
   geom_ribbon(aes(ymin=`5%`,ymax=`95%`),color=NA,alpha=0.2,lwd=0.1) +
   ## observations
   geom_point(aes(y=S_obs),size=0.75,alpha=1) +   
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   labs(x="Spawning year",y="Spawner abundance") + 
   scale_y_continuous(limits=c(0,NA)) +
   facet_wrap(vars(pop),ncol=nP,scales="free_y") + 
   theme_sleek() +
   theme(legend.position="none") +
   theme(
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15)
   ) +
   NULL

##==============================================## estimated recruits
R_post<-extract1(IPM_fit_plot,"R")
R_qs<-t(apply(R_post,2,function(x) quantile(x,prob=probs))) 
R_est_qs<-R_qs %>% 
   data.frame() %>%
   add_column(fish_dat%>%dplyr::select(pop,year)%>%data.frame())

p2 <- R_est_qs %>%
   ggplot(aes(x=year,y=X50.,color=pop,fill=pop)) +
   geom_line(aes(y=X50.),lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.4) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   labs(x="Spawning year",y="Recruitment") + 
   scale_y_continuous(limits=c(0,NA)) +
   facet_wrap(vars(pop),ncol=nP,scales="free_y") + 
   theme_sleek() +
   theme(legend.position="none") +
   theme(
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15)
   ) +
   NULL

##=============================================## kelt survival rates
mu_surv_post<-extract1(IPM_fit_plot,"mu_SS")
mu_surv_qs<-t(quantile(mu_surv_post,prob=probs)) ## quantiles
surv_post<-extract1(IPM_fit_plot,"s_SS") ## survival rates by population

##--------------------------## median survival by population and year
surv_med<-apply(surv_post,2,median) %>% 
   data.frame() %>%
   rename(surv='.') %>%
   add_column(fish_dat%>%dplyr::select(pop,year)%>%data.frame())

##---------------------------------## quantiles by population and year
surv_qs<-t(apply(surv_post,2,function(x) quantile(x,prob=probs))) %>%
   data.frame() %>%
   add_column(fish_dat%>%dplyr::select(pop,year)%>%data.frame())

##---------------------------------## mean survival rate across years
surv_mean_of_annual_medians<-surv_med %>%
   group_by(pop) %>%
   dplyr::summarize(surv=mean(surv)) %>%
   data.frame()

##---------------## mean of medians across past five years and rivers
surv_med_recent_mean<-surv_med %>%
   filter(year>(max(surv_med$year-5))) %>%
   group_by(pop) %>%
   dplyr::summarize(recent_mean=mean(surv))

##------------------------------------## plot kelt survival over time
p3 <- surv_qs %>%
   ggplot(aes(x=year,y=X50.,color=pop,fill=pop)) +
   geom_line(aes(y=X50.),lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.4) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   labs(x="Year",y="Kelt survival rate") + 
   facet_wrap(vars(pop),ncol=nP,scales="free_y") + 
   theme_sleek() +
   theme(legend.position="none") +
   theme(
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15),
      plot.margin=unit(c(0.4,0.4,0.4,1.2),"lines"), ## t-r-b-l
      panel.spacing=unit(1,"lines")
   ) +
   NULL

##===================================================## combined plot
pp<-grid.arrange(p1,p2,p3)
ggsave("IPM-sthd-spawner-recruit-kelt.pdf",pp,width=1+2*nP,height=7.5)

##=================================================================##
##===========================================## recruitment anomalies
##=================================================================##

##==============================================## without covariates
eta_R_post<-extract1(IPM_fit_without_covars,"eta_year_R")
eta_R_qsr<-t(apply(eta_R_post,2,function(x) quantile(x,prob=probs)))
eta_R_qs<-data.frame(eta_R_qsr) %>% add_column(year=dat_years_without_covars)

sigma_R_post<-extract1(IPM_fit_without_covars,"sigma_R")
sigma_R_without_covars<-quantile(sigma_R_post,prob=probs)
sigma_year_R_post<-extract1(IPM_fit_without_covars,"sigma_year_R")
sigma_year_R_without_covars<-quantile(sigma_year_R_post,prob=probs)

##----------------------------------------## estimate trend over time
etas_R<-data.frame(median=apply(eta_R_post,2,median),sd=apply(eta_R_post,2,sd)) %>% add_column(year=dat_years_without_covars)
newdata<-data.frame(year=dat_years_without_covars)
lm_fit<-lm(etas_R$median~etas_R$year,weights=1/(etas_R$sd^2))
pred<-predict(lm_fit,newdata,interval="confidence",level=0.95)
pred<-data.frame(pred) %>% add_column(year=etas_R$year) 

##------------------------------------------------------------## plot
p2a <- eta_R_qs %>%
   ggplot(aes(x=year,y=X50.)) +
   geom_line(lwd=0.4,alpha=1) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.2) +
   labs(x="Spawning year",y="Recruitment anomaly") + 
   scale_y_continuous(limits=c(-1.1,1.1),expand=c(0,0)) +
   #geom_hline(yintercept=0,lwd=0.2,linetype="dashed") +
   geom_line(data=pred,aes(y=fit),lwd=0.2,linetype="dashed") +
   theme_sleek() +
   theme(axis.text.x=element_text(size=12),
         axis.title.x=element_text(size=15),
         axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=15)) +
   NULL

##=================================================## with covariates
eta_R_post<-extract1(IPM_fit_with_covars,"eta_year_R")
eta_R_qsr<-t(apply(eta_R_post,2,function(x) quantile(x,prob=probs)))
eta_R_qs<-data.frame(eta_R_qsr) %>% add_column(year=dat_years_with_covars)

sigma_R_post<-extract1(IPM_fit_with_covars,"sigma_R")
sigma_R_with_covars<-quantile(sigma_R_post,prob=probs)
sigma_year_R_post<-extract1(IPM_fit_with_covars,"sigma_year_R")
sigma_year_R_with_covars<-quantile(sigma_year_R_post,prob=probs)

##----------------------------------------## estimate trend over time
etas_R<-data.frame(median=apply(eta_R_post,2,median),sd=apply(eta_R_post,2,sd)) %>% add_column(year=dat_years_with_covars)
newdata<-data.frame(year=dat_years_with_covars)
lm_fit<-lm(etas_R$median~etas_R$year,weights=1/(etas_R$sd^2))
pred<-predict(lm_fit,newdata,interval="confidence",level=0.95)
pred<-data.frame(pred) %>% add_column(year=etas_R$year) 

##------------------------------------------------------------## plot
p2c<-eta_R_qs %>%
   ggplot(aes(x=year,y=X50.)) +
   geom_line(lwd=0.4,alpha=1) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),fill="goldenrod1",color=NA,alpha=0.3,lwd=0.1) +
   labs(x="Spawning year",y="Recruitment anomaly") + 
   scale_y_continuous(limits=c(-1.1,1.1),expand=c(0,0)) +
   # geom_hline(yintercept=0,lwd=0.2,linetype="dashed") +
   geom_line(data=pred,aes(y=fit),lwd=0.2,linetype="dashed") +
   theme_sleek() +
   theme(axis.text.x=element_text(size=12),
         axis.title.x=element_text(size=15),
         axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=15)) +
   NULL

##===============================================## covariate effects
beta_R_post<-data.frame(extract1(IPM_fit_plot,"beta_R"))
names(beta_R_post)<-c("NPGO","SST","Pinks")#,"Flow")
beta_R_med<-apply(beta_R_post,2,median) ## median effect sizes

##---------------------------------------------------## effects plots
p2b<-beta_R_post %>% 
   pivot_longer(col=everything(),names_to="name",values_to="value")%>%
   ggplot(aes(x=name,y=value)) +
   # geom_violin(lwd=0.1,col="white") +
   stat_summary(fun.data=summary_CI90,size=0.25,color="goldenrod1") +
   stat_summary(fun.data=summary_CI50,size=1.0,color="goldenrod1") +  
   geom_hline(yintercept=0,linetype="dashed",linewidth=0.2) +
   scale_y_continuous(limits=c(-0.3,0.3),expand=c(0,0)) +
   labs(x="",y="Effect size") + 
   theme_sleek() +
   theme(axis.title.x=element_blank(),
         axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1),
         axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=15)) +
   NULL
                            
##---------------------------------## probability above/below zero
my_pnorm<-function(x,output){ return(pnorm(0,mean=x[1],sd=x[2])) }
beta_df<-data.frame(apply(beta_R_post,2,median),apply(beta_R_post,2,sd))
prob_below_zero<-apply(beta_df,1,my_pnorm)
prob_above_zero <- 1-prob_below_zero

##=================================================================##
##=========================================## kelt survival anomalies
##=================================================================##

##==============================================## without covariates
eta_SS_post<-extract1(IPM_fit_without_covars,"eta_year_SS")
eta_SS_qs<-t(apply(eta_SS_post,2,function(x) quantile(x,prob=probs)))
eta_SS_qs<-eta_SS_qs %>% data.frame() %>% add_column(year=dat_years_without_covars)

sigma_SS_post<-extract1(IPM_fit_without_covars,"sigma_SS")
sigma_SS_without_covars<-quantile(sigma_SS_post,prob=probs)
sigma_year_SS_post<-extract1(IPM_fit_without_covars,"sigma_year_SS")
sigma_year_SS_without_covars<-quantile(sigma_year_SS_post,prob=probs)

##----------------------------------------## estimate trend over time
etas_SS<-data.frame(median=apply(eta_SS_post,2,median),sd=apply(eta_SS_post,2,sd)) %>% add_column(year=dat_years_without_covars)
newdata<-data.frame(year=dat_years_without_covars)
lm_fit<-lm(etas_SS$median~etas_SS$year,weights=1/(etas_SS$sd^2))
pred<-predict(lm_fit,newdata,interval="confidence",level=0.95)
pred<-data.frame(pred) %>% add_column(year=etas_SS$year) 

##------------------------------------------------------------## plot
p1a <- eta_SS_qs %>%
   ggplot(aes(x=year,y=X50.)) +
   geom_line(lwd=0.4,alpha=1) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.2) +
   labs(x="Year",y="Kelt survival anomaly") + 
   scale_y_continuous(limits=c(-1.1,1.1),expand=c(0,0)) +
   # geom_hline(yintercept=0,lwd=0.2,linetype="dashed") +
   geom_line(data=pred,aes(y=fit),lwd=0.2,linetype="dashed") +
   theme_sleek() +
   theme(axis.text.x=element_text(size=12),
         axis.title.x=element_text(size=15),
         axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=15)) +
   NULL

##=================================================## with covariates
eta_SS_post<-extract1(IPM_fit_with_covars,"eta_year_SS")
eta_SS_qs<-t(apply(eta_SS_post,2,function(x) quantile(x,prob=probs)))
eta_SS_qs<-eta_SS_qs %>% data.frame() %>% add_column(year=dat_years_with_covars)

sigma_SS_post<-extract1(IPM_fit_with_covars,"sigma_SS")
sigma_SS_with_covars<-quantile(sigma_SS_post,prob=probs)
sigma_year_SS_post<-extract1(IPM_fit_with_covars,"sigma_year_SS")
sigma_year_SS_with_covars<-quantile(sigma_year_SS_post,prob=probs)

##----------------------------------------## estimate trend over time
etas_SS<-data.frame(median=apply(eta_SS_post,2,median),sd=apply(eta_SS_post,2,sd)) %>% add_column(year=dat_years_with_covars)
newdata<-data.frame(year=dat_years_with_covars)
lm_fit<-lm(etas_SS$median~etas_SS$year,weights=1/(etas_SS$sd^2))
pred<-predict(lm_fit,newdata,interval="confidence",level=0.95)
pred<-data.frame(pred) %>% add_column(year=etas_SS$year) 

##------------------------------------------------------------## plot
p1c <- eta_SS_qs %>%
   ggplot(aes(x=year,y=X50.)) +
   geom_line(lwd=0.4,alpha=1) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),fill="goldenrod1",color=NA,alpha=0.3,lwd=0.2) +
   labs(x="Year",y="Kelt survival anomaly") + 
   scale_y_continuous(limits=c(-1.1,1.1),expand=c(0,0)) +
   # geom_hline(yintercept=0,lwd=0.2,linetype="dashed") +
   geom_line(data=pred,aes(y=fit),lwd=0.2,linetype="dashed") +
   theme_sleek() +
   theme(axis.text.x=element_text(size=12),
         axis.title.x=element_text(size=15),
         axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=15)) +
   NULL

##===============================================## covariate effects
cov_eff_post<-data.frame(extract1(IPM_fit_with_covars,"beta_SS"))
# cov_eff_post$NPGO<-NA; names(cov_eff_post)<-c("SST","Pinks","NPGO")
names(cov_eff_post)<-c("SST","Pinks","Flow")
cov_eff_med<-apply(cov_eff_post,2,median) ## median effect sizes

##---------------------------------------------------## effects plots
p1b<-cov_eff_post %>% 
   pivot_longer(col=everything(),names_to="name",values_to="value")%>%
   # mutate(name=case_when(
   #    name=="Flow"~"Flow",
   #    name=="SST"~"SST",
   #    name=="Pinks"~"Pinks"
   #    )) %>%
   ggplot(aes(x=name,y=value)) +
   stat_summary(fun.data=summary_CI90,size=0.25,col="goldenrod1") +
   stat_summary(fun.data=summary_CI50,size=1.0,col="goldenrod1") +  
   geom_hline(yintercept=0,linetype="dashed",linewidth=0.2) +
   scale_y_continuous(limits=c(-0.3,0.3),expand=c(0,0)) +
   labs(x="",y="Effect size") + 
   theme_sleek() +
   theme(axis.title.x=element_blank(),
         axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=1),
         axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=15)) +
   NULL

##------------------------------------## probability above/below zero
my_pnorm<-function(x,output){ return(pnorm(0,mean=x[1],sd=x[2])) }
beta_df<-data.frame(apply(cov_eff_post,2,median),apply(cov_eff_post,2,sd))
prob_below_zero<-apply(beta_df,1,my_pnorm)
prob_above_zero <- 1-prob_below_zero

##=================================================================##
##==========## Figure 3 - recruitment and kelt survival anomaly plots
##=================================================================##
title<-paste0("          Models without environmental drivers     ",
             "              Environmental effects                   ",
             "    Models with environmental drivers")

fig<-ggdraw() +
   theme(plot.margin=margin(20,0,0,0)) +
   draw_label(title,x=0.5,y=1,hjust=0.5,vjust=0,size=15,fontface="bold")+
   draw_plot(p1a,x=0.0,y=0,width=0.4,height=0.5,scale=0.95) +
   draw_plot(p1b,x=0.4,y=0.02,width=0.2,height=0.48,scale=0.95) +
   draw_plot(p1c,x=0.6,y=0,width=0.4,height=0.5,scale=0.95) +
   draw_plot(p2a,x=0.0,y=0.5,width=0.4,height=0.5,scale=0.95) +
   draw_plot(p2b,x=0.4,y=0.52,width=0.2,height=0.48,scale=0.95) +
   draw_plot(p2c,x=0.6,y=0.5,width=0.4,height=0.5,scale=0.95) +
   draw_plot_label(label=c("A","B","C","D","E","F"),
                   x=c(0,0.4,0.6,0,0.4,0.6),
                   y=c(1,1,1,0.5,0.5,0.5),
                   size=15)

save_plot("IPM-sthd-covariate-effects-recruitment-and-kelt-survival.pdf",fig,ncol=3,nrow=2,base_height=4.5,base_width=4.5)

##=================================================================##
##==============## Figure 4 - anomalies as functions of SST and pinks
##=================================================================##

##========================## kelt survival anomaly without covariates
eta_SS_post<-extract1(IPM_fit_without_covars,"eta_year_SS") 

etas_SS_med<-data.frame(anomaly=apply(eta_SS_post,2,median)) %>% 
   add_column(year=dat_years_without_covars) %>%
   mutate(color=factor(ifelse(anomaly>0,"positive","negative")))

pdat_SS<-etas_SS_med %>% 
   slice(-n()) %>%
   left_join(covar_dat%>%dplyr::select(year,SST,pinks),by="year") %>%
   filter(!is.na(pinks)) %>%
   filter(!is.na(SST)) %>%
   add_column(name="Kelt survival anomaly") %>%
   dplyr::select(name,anomaly,year,color,SST,Pinks=pinks)

##==========================## recruitment anomaly without covariates
eta_R_post<-extract1(IPM_fit_without_covars,"eta_year_R") 

etas_R_med<-data.frame(anomaly=apply(eta_R_post,2,median)) %>% 
   add_column(year=dat_years_without_covars) %>%
   mutate(color=factor(ifelse(anomaly>0,"positive","negative")))

pdat_R<-etas_R_med %>%
   left_join(covar_dat%>%dplyr::select(year,SST_4,pinks_4),by="year")%>%
   filter(!is.na(pinks_4)) %>%
   filter(!is.na(SST_4)) %>%
   add_column(name="Recuitment anomaly") %>%
   dplyr::select(name,anomaly,year,color,SST=SST_4,Pinks=pinks_4)

##===================================================## combined plot
# pdat_R<-pdat_R %>% filter(year %in% pdat_SS$year)
# pdat_SS<-pdat_SS %>% filter(year %in% pdat_R$year)

pdat<-rbind(pdat_R,pdat_SS) 

xmin<-min(pdat$SST)-0.2;xmax<-max(pdat$SST)+0.2
ymin<-min(pdat$Pinks)-80;ymax<-max(pdat$Pinks)+80

xmin<-8.6;xmax<-12.0
ymin<-110;ymax<-875

p<-pdat%>% 
   ggplot(aes(x=SST,y=Pinks,size=abs(anomaly),color=color))+
   #geom_hline(yintercept=mean(pdat$Pinks),linetype="dashed",linewidth=0.2)+
   #geom_vline(xintercept=mean(pdat$SST),linetype="dashed",linewidth=0.2)+
   annotate("text",x=xmin,y=ymin,
            size=2.5,color="black",hjust=0,vjust=0,lineheight=0.9, 
            label="cold and low competition\nocean environment") +
   annotate("text",x=xmax,y=ymax,
            size=2.5,color="black",hjust=1,vjust=1,lineheight=0.9,
            label="warm and high competition\nocean environment") +
   geom_point(alpha=0.75)+
   geom_point(shape=1,fill=NA,color="black") + 
   scale_radius(range=c(0.5,9))+
   scale_color_manual(values=c("red","blue"),
                      labels=c("negative","positive")) +
   scale_x_continuous(limits=c(xmin,xmax),breaks=seq(8,15,0.5))+
   scale_y_continuous(limits=c(ymin,ymax),breaks=seq(0,800,100))+
   theme_classic()+
   labs(x="SST (°C)",y="Pink salmon abundance (millions)")+
   labs(size="anomaly",color="")+
   theme(strip.background=element_blank(),
         strip.text=element_text(size=14),
         axis.line=element_line(linewidth=0.1),
         axis.text=element_text(size=12),
         axis.title=element_text(size=14),
         panel.border=element_rect(fill=NA,linewidth=1),
         legend.key.size=unit(0.5,'cm'),
         legend.title=element_text(size=10),
         legend.text=element_text(size=10))+
   guides(color=guide_legend(order=1,override.aes=list(size=5)),
          size=guide_legend(order=2)) +
   facet_wrap(vars(factor(name,c("Recuitment anomaly","Kelt survival anomaly"))),ncol=2,scales="free") +
   NULL

ggsave("IPM-sthd-anomalies-vs-covariates.pdf",p,width=8.5,height=4)
       
##=================================================================##
##=================================================================##
##=================================================================##