##=================================================================##
##                                                                 ##
## Plots for hierarchical IPM for winter steelhead on the WA coast ##
##                                                                 ##
##=================================================================##

##========================================================## packages
pkg<-c("here","tidyverse","ggplot2","cowplot","ggsidekick","gridExtra","RColorBrewer","salmonIPM")
if(length(setdiff(pkg,rownames(installed.packages())))>0){install.packages(setdiff(pkg,rownames(installed.packages())),dependencies=T)}
invisible(lapply(pkg,library,character.only=T))

##========================================================## settings
theme_set(theme_sleek())
options(rstudio.help.showDataPreview=FALSE)

##=======================================================## functions
home<-here::here()
fxn<-list.files(paste0(home,"/functions"))
invisible(sapply(FUN=source,paste0(home,"/functions/",fxn)))

##================================================## output directory
out_dir<-paste0(home,"/output/")
if(!file.exists(out_dir)) dir.create(file.path(out_dir))
setwd(file.path(out_dir))

##=================================================================##
##========================================## load data and model fit
##=================================================================##
covar_effects<-TRUE ## plot model with or without covariate effects

##=========================================================## IPM fit
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
##======================================================## main plots
if(covar_effects){
   IPM_fit<-IPM_fit_with_covars
   fish_dat<-fish_dat_with_covars
}else{
   IPM_fit<-IPM_fit_without_covars
   fish_dat<-fish_dat_without_covars
}
pops<-unique(fish_dat$pop)
nP<-length(pops)
N<-dim(fish_dat)[1]
   
##=======================================================## posterior
df_post<-as.data.frame(IPM_fit) ## full posterior
nS<-dim(df_post)[1] ## number of samples
df_out<-data.frame(summary(IPM_fit)$summary) ## summary 

##===================================================## plot settings
# colors<-colorRampPalette(brewer.pal(name="Set1",n=9))(nP)
##----------------------------------------------------------## colors
colors<-rev(c("#A86260","#A5B1C4","#3B4D6B","#89A18D","#747260","#B3872D","#774C2C")) 
colors <- colorRampPalette(colors)(nP)
##-------------------------------------------------------## functions
summary_CI95<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.025),na.rm=T),ymax=quantile(x,prob=c(0.975),na.rm=T))) }
summary_CI90<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.05),na.rm=T),ymax=quantile(x,prob=c(0.95),na.rm=T))) }
summary_CI50<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.25),na.rm=T),ymax=quantile(x,prob=c(0.75),na.rm=T))) }
##-------------------------------------------------------------## CIs
probs<-c(0.05,0.25,0.5,0.75,0.95) ## median with 50% and 90% CIs

##=================================================================##
##=======================## observation error, productivity, capacity
##=================================================================##

##===============================================## observation error
tau_post<-extract1(IPM_fit,"tau")
tau_median<-median(tau_post)

##============================================## productivity (alpha)
alphas_post<-data.frame(extract1(IPM_fit,"alpha"))
names(alphas_post)<-pops
alpha_qs<-apply(alphas_post,2,function(x) quantile(x,prob=probs))
alpha_medians<-round(apply(alphas_post,2,median),1)

mu_alpha_post<-data.frame(exp(extract1(IPM_fit,"mu_alpha")))
names(mu_alpha_post)<-" mu_alpha "
mu_alpha_qs<-t(apply(mu_alpha_post,2,function(x) quantile(x,prob=probs)))

##=================================================## capacity (Rmax)
Rmax_post<-data.frame(extract1(IPM_fit,"Rmax"))
names(Rmax_post)<-pops
Rmax_qs<-apply(Rmax_post,2,function(x) quantile(x,prob=probs))
Rmax_medians<-round(apply(Rmax_post,2,median))

mu_Rmax_post<-data.frame(exp(extract1(IPM_fit,"mu_Rmax")))
names(mu_Rmax_post)<-" mu_Rmax"
mu_Rmax_qs<-t(apply(mu_Rmax_post,2,function(x) quantile(x,prob=probs)))

##=================================================================##
##============================## spawners, recruitment, kelt survival
##=================================================================##

##==============================================## estimated spawners
S_post<-extract1(IPM_fit,"S")
S_qs<-t(apply(S_post,2,function(x) quantile(x,prob=probs)))
S_est_qs<-S_qs %>% data.frame() %>%
   add_column(fish_dat%>%dplyr::select(pop,year,S_obs)%>%data.frame())

p1 <- S_est_qs %>%
   # filter(pop=='Humptulips') %>%
   ggplot(aes(x=year,y=X50.,color=pop,fill=pop)) +
   geom_line(aes(y=X50.),lwd=0.5,alpha=1) +
   geom_point(aes(y=S_obs),size=0.6,alpha=1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.4) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.2) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   labs(x="Brood year",y="Spawner abundance") + 
   scale_y_continuous(limits=c(0,NA)) +
   facet_wrap(vars(pop),ncol=nP,scales="free_y") + ## ,scales="free_x"
   theme_sleek() +
   theme(legend.position="none") +
   theme(
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15)
   ) +
   NULL

##==============================================## estimated recruits
R_post<-extract1(IPM_fit,"R")
R_qs<-t(apply(R_post,2,function(x) quantile(x,prob=probs))) 
R_est_qs<-R_qs %>% 
   data.frame() %>%
   add_column(fish_dat%>%dplyr::select(pop,year)%>%data.frame())

p2 <- R_est_qs %>%
   ggplot(aes(x=year,y=X50.,color=pop,fill=pop)) +
   geom_line(aes(y=X50.),lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.4) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.2) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   labs(x="Brood year",y="Recruitment") + 
   scale_y_continuous(limits=c(0,NA)) +
   facet_wrap(vars(pop),ncol=nP,scales="free_y") +  # scales="free_x"
   theme_sleek() +
   theme(legend.position="none") +
   theme(
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15)
   ) +
   NULL

##============================================## recruits per spawner
RpS_post<-extract1(IPM_fit,"R")/extract1(IPM_fit,"S")
RpS_qs<-t(apply(RpS_post,2,function(x) quantile(x,prob=probs))) %>%
   data.frame() %>%
   rename_all(list(~stringr::str_replace_all(.,'X','RpS_'))) %>%
   add_column(fish_dat%>%dplyr::select(pop,year)%>%data.frame())

##------------------------------------## mean age structure by cohort
## p: cohort (maiden) age distributions
p_table<-df_out[grep("^(?=.*p)(?!.*_)",rownames(df_out),perl=TRUE),]
p_table<-p_table[!grepl("alpha",rownames(p_table)),]
m_names<-colnames(fish_dat)[grepl("M_obs",colnames(fish_dat))]
m_ages<-gsub("n_age","pM",gsub("M_obs","",m_names))
nm_age<-length(m_ages)
p_est<-matrix(p_table$mean,ncol=nm_age,nrow=N,byrow=TRUE) %>% 
   data.frame() %>%
   rename_with(~m_ages,all_of(names(.))) %>%
   add_column(fish_dat%>%dplyr::select(pop,year))

##----------## kelt survival rate and recruits/spawner at equilibrium
## s_SS: true kelt survival by out-migration year
s_SS_post<-extract1(IPM_fit,"s_SS")
surv_med<-data.frame(surv=apply(s_SS_post,2,median)) %>% 
   add_column(fish_dat%>%dplyr::select(pop,year)%>%data.frame()) %>%
   dplyr::select(pop,year,surv) %>%
   full_join(p_est,by=c("pop","year")) %>%
   group_by(pop) %>% ## calculate kelt survival rate X years later
   mutate( ## RpS by brood year so kelt survival needs to be lagged
      surv3=lead(surv,3),
      surv4=lead(surv,4), 
      surv5=lead(surv,5), 
      surv6=lead(surv,6), 
      surv7=lead(surv,7), 
      surv8=lead(surv,8), 
   ) %>% ## now apply kelt survival for dominant maiden ages
   ungroup() %>% ## spawn year same as kelt out-migration year
   mutate(RpS_M3=1/(1+surv3+surv3*surv4+surv3*surv4*surv5)) %>%
   mutate(RpS_M4=1/(1+surv4+surv4*surv5+surv4*surv5*surv6)) %>%
   mutate(RpS_M5=1/(1+surv5+surv5*surv6+surv5*surv6*surv7)) %>%
   mutate(RpS_M6=1/(1+surv6+surv6*surv7+surv6*surv7*surv8)) %>%
   mutate(RpS_eq=RpS_M3*pM3+RpS_M4*pM4+RpS_M5*pM5+RpS_M6*pM6) %>%
   dplyr::select(pop,year,RpS_eq) %>%
   data.frame()

##---------------------------------## recruits/spawner at equilibrium
RpS_plot<-RpS_qs %>% 
   left_join(surv_med,by=c("pop","year"))

if(covar_effects) { 
   RpS_plot<-RpS_plot 
}else{ 
   RpS_plot<-RpS_plot %>% 
      group_by(pop) %>%
      slice(1:(n()-5)) %>%
      data.frame()
} 
## drop years with incomplete recruitment when no covariates in model
## i.e., only recruitment for which dominant ages have been observed

##------------------------------## plot time-varying recruits/spawner
p3 <- RpS_plot  %>%
   ggplot(aes(x=year,y=RpS_50.,color=pop,fill=pop)) +
   ## estimated recruits per spawner median, 50%, and 90% CIs
   geom_line(aes(y=RpS_50.),lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=RpS_25.,ymax=RpS_75.),color=NA,alpha=0.4) +
   geom_ribbon(aes(ymin=RpS_5.,ymax=RpS_95.),color=NA,alpha=0.2) +
   ## estimated recruits per spawner at equilibrium
   geom_line(aes(y=RpS_eq),lwd=0.5,alpha=0.5) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   labs(x="Brood year",y="Recruits/spawner") + 
   scale_y_continuous(limits=c(0,NA)) + ## ,expand=c(0,0)
   facet_wrap(vars(pop),ncol=nP,scales="free_y") + #scales="free_y") +
   theme_sleek() +
   theme(legend.position="none") +
   theme(
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15),
      plot.margin=unit(c(0.4,0.4,0.4,1.8),"lines"), ## t-r-b-l
      panel.spacing=unit(1.4,"lines")
   ) +
   NULL

##=============================================## kelt survival rates
mu_surv_post<-extract1(IPM_fit,"mu_SS")
mu_surv_qs<-t(quantile(mu_surv_post,prob=probs)) ## quantiles
surv_post<-extract1(IPM_fit,"s_SS") ## survival rates by population

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
p4 <- surv_qs %>%
   ggplot(aes(x=year,y=X50.,color=pop,fill=pop)) +
   geom_line(aes(y=X50.),lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.4) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.2) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   # scale_x_continuous(limits=c(1978,2022)) +
   labs(x="Year",y="Kelt survival rate") + 
   facet_wrap(vars(pop),ncol=nP,scales="free_y") + 
   theme_sleek() +
   theme(legend.position="none") +
   theme(
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
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
pp<-grid.arrange(p1,p2,p4)
ggsave("IPM-sthd-spawner-recruit-kelt.pdf",pp,width=1+2*nP,height=7.5)

##=================================================================##
##=========================================## kelt survival anomalies
##=================================================================##

##==============================================## without covariates
eta_SS_post<-extract1(IPM_fit_without_covars,"eta_year_SS")
eta_SS_qs<-t(apply(eta_SS_post,2,function(x) quantile(x,prob=probs)))
eta_SS_qs<-eta_SS_qs %>% data.frame() %>% add_column(year=dat_years_without_covars)

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
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.1) +
   labs(x="Year",y="Kelt survival anomaly") + 
   scale_y_continuous(limits=c(-1,1)) +
   ## add weighted linear fit based on median recruitment residuals
   geom_line(data=pred,aes(y=fit),lwd=0.2,linetype="dashed") +
   #geom_line(data=pred,aes(y=lwr),lwd=0.2,color="darkgray",alpha=1) +
   #geom_line(data=pred,aes(y=upr),lwd=0.2,color="darkgray",alpha=1) +
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
   geom_ribbon(aes(ymin=X25.,ymax=X75.),fill="goldenrod1",color="darkgray",alpha=0.3,lwd=0.2) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),fill="goldenrod1",color="darkgray",alpha=0.1,lwd=0.2) +
   labs(x="Year",y="Kelt survival anomaly") + 
   scale_y_continuous(limits=c(-1,1)) +
   ## add weighted linear fit based on median recruitment residuals
   geom_line(data=pred,aes(y=fit),lwd=0.2,linetype="dashed") +
   #geom_line(data=pred,aes(y=lwr),lwd=0.2,color="darkgray",alpha=1) +
   #geom_line(data=pred,aes(y=upr),lwd=0.2,color="darkgray",alpha=1) +
   theme_sleek() +
   theme(axis.text.x=element_text(size=12),
         axis.title.x=element_text(size=15),
         axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=15)) +
   NULL

##===============================================## covariate effects
cov_eff_post<-data.frame(extract1(IPM_fit_with_covars,"beta_SS"))
apply(cov_eff_post,2,median)
names(cov_eff_post)<-c("NPGO","SST","Pinks")
##---------------------------------------------------## effects plots
p1b <- cov_eff_post %>% 
   pivot_longer(col=everything(),names_to="name",values_to="value") %>%
   ggplot(aes(x=name,y=value)) +
   # geom_violin(lwd=0.1,col="gray") +
   stat_summary(fun.data=summary_CI90,size=0.25,col="goldenrod1") +
   # stat_summary(fun.data=summary_CI95,size=0.25,col="goldenrod1") +
   stat_summary(fun.data=summary_CI50,size=1.0,col="goldenrod1") +  
   geom_hline(yintercept=0,linetype="dashed",size=0.2) +
   scale_y_continuous(limits=c(-0.28,0.11)) +
   labs(x="",y="Effect size") + 
   theme_sleek() +
   theme(axis.title.x=element_blank(),
         axis.text.x=element_text(size=12),
         axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=15)) +
   NULL
##------------------------------------## probability above/below zero
my_pnorm<-function(x,output){ return(pnorm(0,mean=x[1],sd=x[2])) }
beta_df<-data.frame(apply(cov_eff_post,2,median),
                    apply(cov_eff_post,2,sd))
prob_below_zero<-apply(beta_df,1,my_pnorm)
prob_above_zero <- 1-prob_below_zero

##=================================================================##
##===========================================## recruitment anomalies
##=================================================================##

##==============================================## without covariates
eta_R_post<-extract1(IPM_fit_without_covars,"eta_year_R")
eta_R_qsr<-t(apply(eta_R_post,2,function(x) quantile(x,prob=probs)))
eta_R_qs<-data.frame(eta_R_qsr) %>% add_column(year=dat_years_without_covars)

##----------------------------------------## estimate trend over time
etas_R<-data.frame(median=apply(eta_R_post,2,median),sd=apply(eta_R_post,2,sd)) %>% add_column(year=dat_years_without_covars)
newdata<-data.frame(year=dat_years_without_covars)
lm_fit<-lm(etas_R$median~etas_R$year,weights=1/(etas_R$sd^2))
pred<-predict(lm_fit,newdata,interval="confidence",level=0.95)
pred<-data.frame(pred) %>% add_column(year=etas_R$year) 

##------------------------------------------------------------## plot
p2a <- eta_R_qs %>%
   ## plot median recruitment residuals and 95% CIs
   ggplot(aes(x=year,y=X50.)) +
   geom_line(lwd=0.4,alpha=1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.1) +
   labs(x="Brood year",y="Recruitment anomaly") + 
   scale_y_continuous(limits=c(-1,0.7)) +
   ## add weighted linear fit based on median recruitment residuals
   geom_line(data=pred,aes(y=fit),lwd=0.2,linetype="dashed") +
   #geom_line(data=pred,aes(y=lwr),lwd=0.2,color="darkgray",alpha=1) +
   #geom_line(data=pred,aes(y=upr),lwd=0.2,color="darkgray",alpha=1) +
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

##----------------------------------------## estimate trend over time
etas_R<-data.frame(median=apply(eta_R_post,2,median),sd=apply(eta_R_post,2,sd)) %>% add_column(year=dat_years_with_covars)
newdata<-data.frame(year=dat_years_with_covars)
lm_fit<-lm(etas_R$median~etas_R$year,weights=1/(etas_R$sd^2))
pred<-predict(lm_fit,newdata,interval="confidence",level=0.95)
pred<-data.frame(pred) %>% add_column(year=etas_R$year) 

##------------------------------------------------------------## plot
p2c <- eta_R_qs %>%
   ## plot median recruitment residuals and 95% CIs
   ggplot(aes(x=year,y=X50.)) +
   geom_line(lwd=0.4,alpha=1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),fill="goldenrod1",color="darkgray",alpha=0.3,lwd=0.1) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),fill="goldenrod1",color="darkgray",alpha=0.1,lwd=0.1) +
   labs(x="Brood year",y="Recruitment anomaly") + 
   scale_y_continuous(limits=c(-1,0.7)) +
   ## add weighted linear fit based on median recruitment residuals
   geom_line(data=pred,aes(y=fit),lwd=0.2,linetype="dashed") +
   #geom_line(data=pred,aes(y=lwr),lwd=0.2,color="darkgray",alpha=1) +
   #geom_line(data=pred,aes(y=upr),lwd=0.2,color="darkgray",alpha=1) +
   theme_sleek() +
   theme(axis.text.x=element_text(size=12),
         axis.title.x=element_text(size=15),
         axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=15)) +
   NULL

##===============================================## covariate effects
beta_R_post<-data.frame(extract1(IPM_fit,"beta_R")) ## posterior
apply(beta_R_post,2,median)
names(beta_R_post)<-c("NPGO","SST","Pinks")
## new names
# cov1<-c("Pink salmon effect estimate")
# cov2<-c(expression("Summer SST"[arc]*" (°C) effect estimate"))
# cov3<-c("North Pacific Gyre Oscillation")
# names(beta_R_post)<-c(cov1,cov2,cov3)

##---------------------------------------------------## effects plots
p2b <- beta_R_post %>% 
   pivot_longer(col=everything(),names_to="name",values_to="value") %>%
   ggplot(aes(x=name,y=value)) +
   geom_violin(lwd=0.1,col="white") +
   stat_summary(fun.data=summary_CI90,size=0.25,color="goldenrod1") +
   # stat_summary(fun.data=summary_CI95,size=0.25,color="goldenrod1") +
   stat_summary(fun.data=summary_CI50,size=1.0,color="goldenrod1") +  
   geom_hline(yintercept=0,linetype="dashed",size=0.2) +
   labs(x="",y="Effect size") + 
   theme_sleek() +
   theme(axis.title.x=element_blank(),
         axis.text.x=element_text(size=12),
         axis.text.y=element_text(size=12),
         axis.title.y=element_text(size=15)) +
   NULL
                            
##---------------------------------## probability above/below zero
my_pnorm<-function(x,output){ return(pnorm(0,mean=x[1],sd=x[2])) }
beta_df<-data.frame(apply(beta_R_post,2,median),
                    apply(beta_R_post,2,sd))
prob_below_zero<-apply(beta_df,1,my_pnorm)
prob_above_zero <- 1-prob_below_zero

##=================================================================##
##==========## combined recruitment and kelt survival anomaly figures
##=================================================================##
title<-paste0("      Models without covariates                    ",
             "     Covariate effects                             ",
             " Models with covariates")

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

save_plot("IPM-sthd-covariate-effects-recruitment-and-kelt-survival.pdf",fig,ncol=3,nrow=2,base_height=4,base_width=4)

##=================================================================##
##=================================================================##
##=================================================================##