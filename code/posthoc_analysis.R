##=================================================================##
##                                                                 ##
##                      Run posthoc analysis                       ##
##                                                                 ##
##=================================================================##

##========================================================## packages
pkg<-c("here","tidyverse","MuMIn","relaimpo","visreg","salmonIPM")
if(length(setdiff(pkg,rownames(installed.packages())))>0){install.packages(setdiff(pkg,rownames(installed.packages())),dependencies=T)}
invisible(lapply(pkg,library,character.only=T))
home<-here::here()

##================================================## output directory
out_dir<-paste0(home,"/output/")
if(!file.exists(out_dir)) dir.create(file.path(out_dir))
setwd(file.path(out_dir))

##=================================================================##
##=========================## load IPM fit, covariates, and fish data
##=================================================================##
IPM_fit<-readRDS(paste0(home,"/output/IPM_fit_without_covars.Rdata"))
cov_dat<-read.csv(paste0(home,"/data/IPM_covar_dat_all.csv"))
fish_dat<-read.csv(paste0(home,"/data/IPM_fish_dat_all.csv"))
dat_years<-sort(unique(fish_dat$year))
pops<-unique(fish_dat$pop)
nP<-length(pops)
N<-dim(fish_dat)[1]

##=================================================================##
##===========================## linear model of recruitment anomalies
##=================================================================##
probs<-c(0.05,0.25,0.5,0.75,0.95) ## median with 50% and 90% CIs

##===========================================## recruitment anomalies
eta_R_post<-extract1(IPM_fit,"eta_year_R")
etas_R<-data.frame(median=apply(eta_R_post,2,median),sd=apply(eta_R_post,2,sd)) %>% add_column(year=dat_years)

##==============================## merge residuals and covariate data
df<-data.frame(year=seq(min(dat_years-4),max(dat_years),1)) %>%
   left_join(etas_R %>% dplyr::select(year,median,sd)) %>%
   rename(residuals=median) %>%
   left_join(cov_dat,by='year') 

##------------------------------------------------## scale covariates
df1<-df %>% dplyr::select(year,residuals,sd) 
df2<-df %>% dplyr::select(-year,-residuals,-sd) 
df_means<-sapply(df2,function(x) mean(x,na.rm=T))
df_sds<-sapply(df2,function(x) sd(x,na.rm=T))
df<-data.frame(cbind(df1,scale(df2))) %>% na.omit()
##----------------------------------------------------## correlations
out<-data.frame(cor(as.matrix(df),use="pairwise.complete.obs"))
res<-data.frame(round(out[-c(1:3),1:2],2)) %>% dplyr::select(-year)

##=================================================## model selection
options(na.action="na.fail")
##------------------------------------------------## year effect only
mod_lm_year<-lm(residuals~year,data=df)
data.frame(summary(mod_lm_year)$coefficients)[2,]
##------------------------------------------------------## model form
form<-formula(residuals~pinks_2+pinks_3+pinks_4+NPGO_2+NPGO_3+NPGO_4+SST_2+SST_3+SST_4+SST_cst_2+SST_cst_3+SST_cst_4+av_CFS_min+av_CFS_max+av_CFS_min_1+av_CFS_max_1) 
##---------------------------------------------------## model fitting
mod_lm_full<-lm(form,data=df)
mod_select<-dredge(mod_lm_full,trace=F,rank="AICc")
aic_table<-data.frame(mod_select)[1:10,-1]
mod_delta2<-aic_table[aic_table$delta<2,] ## models with delta_AIC<2
index<-which(mod_delta2$df==min(mod_delta2$df)) 
mod_sel<-get.models(mod_select,subset=1)[[1]] ## lowest AICc
summary(mod_sel) 
ncovar<-dim(summary(mod_sel)$coefficients)[1]-1
par(mfcol=c(1,ncovar),mar=c(4,4,1,1));visreg(mod_sel)

##==================================================## selected model
mod_rec<-mod_sel
##---------------------------------------------------## model results
residuals<-residuals(mod_rec,type="response")
fitted<-fitted(mod_rec)
out_mod<-summary(mod_rec) 
as.numeric(pacf(residuals(mod_rec),lag=9,plot=F)$acf) 
car::vif(mod_rec) ## VIF
##----------------------------------## pairwise covariate correlations
mod_terms<-as.character(names(mod_rec[[1]])[-1])
nterms<-length(mod_terms)
test_data<-df[,colnames(df) %in% mod_terms]
test<-apply(test_data,2,function(x) as.numeric(as.character(x)))
round(rcorr(test,type="pearson")$r,2)

##=============================================## variable importance
relimp<-calc.relimp(mod_rec,type="lmg") ## print(relimp)
v1<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"pinks")]))
v2<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"SST")]))
v3<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"NPGO")]))
v4<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"CFS")]))
relimp_df<-data.frame(pinks=v1,sst=v2,npgo=v3,flow=v4)
relimp_df<-round(relimp_df*100,2) ## % of response variance explained
round(sum(relimp_df),2) ## total % of response variance explained

write.csv(relimp_df,"percent_recruit_var_expl_posthoc.csv")

##=================================================## partial effects
pdf("IPM-sthd-recruitment-anamaly-covariates-effects.pdf",width=4*nterms,height=4)
layout(matrix(c(1:length(mod_terms)),nrow=1,byrow=T))
par(mar=c(4,4,1,1),mgp=c(2.5,0.5,0),tcl=-0.3,cex.lab=1.5,cex.axis=1.2)
xlab<-c("NPGO","Pink salmon abundance (millions)","Summer SST (°C)","Summer low flow (cf/s)")
xtrans<-function(x) { x*df_sds[names(df_sds)==covar]+df_means[names(df_means)==covar] }
for(i in 1:length(mod_terms)) {
   covar<-mod_terms[i]
   visreg(mod_rec,xvar=covar,xlab=xlab[i],partial=T,ylab="Partial effect on recruitment anomaly",scale="response",xtrans=xtrans,points.par=list(cex=1,pch=16,col=1),fill.par=list(col=alpha(1,0.1)),line.par=list(col=1))
}
dev.off()

##=================================================================##
##====================## linear model of kelt survival rate anomalies
##=================================================================##

##===========================================## kelt survival anomaly
eta_SS_post<-extract1(IPM_fit,"eta_year_SS")
etas_SS<-data.frame(median=apply(eta_SS_post,2,median),sd=apply(eta_SS_post,2,sd)) %>% add_column(year=dat_years)

##=======================================## merge with covariate data
## leads/lags not used for kelt survival
df<-data.frame(year=seq(min(dat_years-4),max(dat_years),1)) %>%
   left_join(etas_SS %>% dplyr::select(year,median,sd)) %>%
   rename(survival=median) %>%
   left_join(cov_dat,by='year') %>%
   dplyr::select(-contains("_1")) %>%
   dplyr::select(-contains("_2")) %>%
   dplyr::select(-contains("_3")) %>%
   dplyr::select(-contains("_4")) 
##------------------------------------------------## scale covariates
df1<-df %>% dplyr::select(year,survival,sd) 
df2<-df %>% dplyr::select(-year,-survival,-sd) 
df_means<-sapply(df2,function(x) mean(x,na.rm=T))
df_sds<-sapply(df2,function(x) sd(x,na.rm=T))
df<-data.frame(cbind(df1,scale(df2))) %>% na.omit()
##----------------------------------------------------## correlations
out<-data.frame(cor(as.matrix(df),use="pairwise.complete.obs"))
res<-data.frame(round(out[-c(1:3),1:2],2)) %>% dplyr::select(-year)

##=================================================## model selection
options(na.action="na.fail")
##------------------------------------------------## year effect only
mod_lm_year<-lm(survival~year,data=df)
data.frame(summary(mod_lm_year)$coefficients)[2,]
##--------------------------------------------## model form (no lags)
form<-formula(survival~pinks+NPGO+SST+SST_cst+av_CFS_min+av_CFS_max)
##---------------------------------------------------## model fitting
mod_lm_full<-lm(form,data=df)
mod_select<-dredge(mod_lm_full,trace=F,rank="AICc")
aic_table<-data.frame(mod_select)[1:10,-1]
mod_delta2<-aic_table[aic_table$delta<2,] ## models with delta_AIC<2
index<-which(mod_delta2$df==min(mod_delta2$df))
mod_sel<-get.models(mod_select,subset=1)[[1]] ## lowest AICc
summary(mod_sel)
ncovar<-dim(summary(mod_sel)$coefficients)[1]-1
par(mfcol=c(1,ncovar),mar=c(4,4,1,1));visreg(mod_sel)

##==================================================## selected model
mod_surv<-mod_sel
##---------------------------------------------------## model results
resid<-residuals(mod_surv,type="response")
fitted<-fitted(mod_surv)
out_mod<-summary(mod_surv) 
as.numeric(pacf(residuals(mod_surv),lag=9,plot=F)$acf) 
mod_terms<-as.character(names(mod_surv[[1]])[-1])
nterms<-length(mod_terms)
if(nterms>1) car::vif(mod_surv) ## VIF
##----------------------------------## pairwise covariate correlations
if(nterms>1) {
   test_data<-df[,colnames(df) %in% mod_terms]
   test<-apply(test_data,2,function(x) as.numeric(as.character(x)))
   round(rcorr(test,type="pearson")$r,4)
}

##=============================================## variable importance
relimp<-calc.relimp(mod_surv,type="lmg") ## print(relimp)
v1<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"pinks")]))
v2<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"SST")]))
v3<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"NPGO")]))
relimp_df<-data.frame(pinks=v1,sst=v2,npgo=v3)
relimp_df<-round(relimp_df*100,2) ## % of response variance explained
round(sum(relimp_df),2) ## total % of response variance explained

write.csv(relimp_df,"percent_keltsurv_var_expl_posthoc.csv")

##=================================================## partial effects
pdf("IPM-sthd-kelt-survival-anomaly-covariate-effects.pdf",width=4*nterms,height=4)
layout(matrix(c(1:length(mod_terms)),nrow=1,byrow=T))
par(mar=c(4,4,1,1),mgp=c(2.5,0.5,0),tcl=-0.3,cex.axis=1.2,cex.lab=1.5)
xlabels<-data.frame(term=mod_terms) %>%
   mutate(name=case_when(
      grepl("pink",term) ~ "Pink salmon abundance (millions)", 
      grepl("SST",term) ~ "Summer SST (°C)", 
      grepl("NPGO",term) ~ "NPGO", 
   ))
xtrans<-function(x) { x*df_sds[names(df_sds)==covar]+df_means[names(df_means)==covar] }
for(i in 1:length(mod_terms)) {
   covar<-mod_terms[i]
   visreg(mod_surv,xvar=covar,xlab=xlabels$name[i],partial=T,ylab="Partial effect on kelt survival anomaly",scale="response",xtrans=xtrans,points.par=list(cex=1,pch=16,col=1),fill.par=list(col=alpha(1,0.1)),line.par=list(col=1))
} 
dev.off()

##=================================================================##
##=================================================================##
##=================================================================##