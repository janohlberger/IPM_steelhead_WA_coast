##=================================================================##
##                                                                 ##
## Load coastal steelhead population data and fit hierarchical IPM ##
##                                                                 ##
##=================================================================##

##========================================================## packages
pkg<-c("here","dplyr","tidyverse","rstan","readr","readxl","tibble", "dataRetrieval","posterior","ggsidekick","RColorBrewer","ggplot2", "officer","MuMIn","ncdf4","reshape2","pracma","relaimpo","visreg", "Hmisc","bayesdfa","MARSS","faraway","gtools","gridExtra","gsl","rcartocolor","bayesplot","rstanarm","distributional") ##,"TMB"
if(length(setdiff(pkg,rownames(installed.packages())))>0){install.packages(setdiff(pkg,rownames(installed.packages())),dependencies=T)}
invisible(lapply(pkg,library,character.only=T))
## devtools::install_github("ebuhle/salmonIPM@dev",auth_token=PAT)
library(salmonIPM)
theme_set(theme_sleek())
options(rstudio.help.showDataPreview=FALSE)

##===================================================## rstan options
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

##=======================================================## functions
home<-here::here()
fxn<-list.files(paste0(home,"/R/functions"))
invisible(sapply(FUN=source,paste0(home,"/R/functions/",fxn)))

##=================================================================##
##=================================================================##
##============================================================## data
##=================================================================##
##=================================================================##

##=================================================================##
##=======================================================## fish data
##=================================================================##
spatial_extent<-"northcoast" ## 'OP',"coastal", "northcoast"
cdata_type<-"rate" ## 'count' or 'rate'
##=======================================================## load data
file_dir<-age_dir<-paste0(home,"/R/data/RiverFiles")
##--------------------------------------------------------------## OP
if(spatial_extent=="OP"){
pops<-c("Hoh","Queets","Quillayute","Quinault") ## OP (alphabetical)
# skm<-c(10431,12144,19571,14723) ## stream km (OP sthd petition Tab4)
# areas<-c(299,445,629,325) ## square miles watershed
areas<-c(4.38,6.95,9.47,6.61) ## km^2 habitat (Ohlberger et al. 2017)
areaD<-data.frame(cbind(pop=pops,area=areas))
# write.csv(areaD,"IPM_habitat_areas.csv",row.names=F)
}
##---------------------------------------------------------## coastal
if(spatial_extent=="coastal"){
pops<-c("Chehalis","Hoh","Humptulips","Queets","Quillayute","Quinault","Willapa","North","Palix","Nemah","Naselle","Bear")
areaD<-NULL
}
##---------------------------------------------------------## coastal
if(spatial_extent=="northcoast"){
pops<-c("Chehalis","Hoh","Humptulips","Queets","Quillayute","Quinault")
areaD<-NULL
}

##-----------------------------------------------## prepare fish data
pops<-sort(pops)
nP<-length(pops)
for(i in 1:nP) {
   pop_dat<-data_prep(
      system=pops[i],
      covars=FALSE,
      file_dir=file_dir,
      age_dir=age_dir)
   if(i==1){fish_dat<-pop_dat}else{fish_dat<-bind_rows(fish_dat,pop_dat)}
}
if(spatial_extent=="OPswWA"){
fish_dat$S_obs[fish_dat$pop=="Grays River" & fish_dat$year==1999]<-1000
}
##-----------------------------## only years with observed escapement
fish_dat<-fish_dat %>% filter(S_obs>0) 
# fish_dat<-fish_dat[!is.na(fish_dat$year),]
all_years<-sort(unique(fish_dat$year))
##------------------------------------------## add area by population
if(spatial_extent %in% c("OP")){
fish_dat<-fish_dat %>% 
   left_join(data.frame(cbind(pop=pops,area=areas))) %>%
   mutate_at('area',as.numeric) %>%
   dplyr::select(-A) %>%
   rename(A=area)
}

##=====================## limit number of initial years without catch
##-------------------------## get first year catch data by population
frst_yrs<-fish_dat %>%
   filter(!is.na(F_rate)) %>% 
   group_by(pop) %>% 
   slice_min(year) %>% 
   dplyr::select(pop,min_year=year)
##---------------## drop data prior to first catch year minus min age
fish_dat<-fish_dat %>% 
   left_join(frst_yrs) %>%
   filter(year>min_year-3) %>% ## min_age=2
   dplyr::select(-min_year)
## now early year NAs in F_rate/B_take_obs can be set to zero because post-removal spawner abundance in years 1:min_age is drawn from prior
##----------------------------## set NAs in F_rate/B_take_obs to zero
if(cdata_type=="count") { 
   fish_dat<-fish_dat %>% 
      mutate(F_rate=0) %>%
      mutate(B_take_obs=replace_na(B_take_obs,0))
   }
if(cdata_type=="rate") {
   fish_dat<-fish_dat %>% 
      mutate(B_take_obs=0) %>%
      mutate(F_rate=replace_na(F_rate,0))
}

##==============================================## correct age format
age_cols<-names(fish_dat)[grepl("n_age",names(fish_dat))]
##-------------------## drop age columns with only zeros if necessary
ind<-which(colSums(fish_dat[,names(fish_dat)%in%age_cols],na.rm=T)==0)
drop<-which(names(fish_dat) %in% names(ind))
fish_dat<-fish_dat %>% dplyr::select(-all_of(drop))
##-------------------------------------## names of all all age groups
m_names<-sort(names(fish_dat)[grepl("M_obs",names(fish_dat))])
m_ages<-as.numeric(gsub("n_age","",gsub("M_obs","",m_names)))
r_names<-sort(names(fish_dat)[grepl("K_obs",names(fish_dat))])
r_ages<-as.numeric(gsub("n_age","",gsub("K_obs","",r_names)))
##---------------------------## drop youngest maiden age if necessary
if(min(m_ages)+1==min(r_ages)){ 
   fish_dat<-fish_dat 
}else{ 
   fish_dat<-fish_dat %>% 
      dplyr::select(-m_names[1]) 
}
##-------------------## create repeat spawner plus-group if necessary
index<-which(r_ages>max(m_ages)) 
if(length(index)>0){
   r_plus_group<-r_names[index[1]]
   r_older_ages<-r_names[index[-1]]
   fish_dat[names(fish_dat)==r_plus_group]<-rowSums(fish_dat[names(fish_dat) %in% c(r_plus_group,r_older_ages)],na.rm=T)
   fish_dat<-fish_dat %>% dplyr::select(-all_of(r_older_ages))
}
nages<-length(colnames(fish_dat)[grepl("M_obs",colnames(fish_dat))])
##--------------------------------------------------------## reorder
m_names<-sort(names(fish_dat)[grepl("M_obs",names(fish_dat))])
r_names<-sort(names(fish_dat)[grepl("K_obs",names(fish_dat))])
##-----------------------------------## convert NA in ages into zeros
index<-which(names(fish_dat) %in% c(m_names,r_names)) ## "B_take_obs"
fish_dat<-fish_dat %>% mutate_at(index, ~replace_na(.,0))
##---------------------------------------------## columns in fish_dat
fish_dat<-fish_dat %>% dplyr::select(pop,A,year,S_obs,all_of(m_names),all_of(r_names),n_W_obs,n_H_obs,B_take_obs,fit_p_HOS,F_rate)

##==============================================## set sub-directory
setwd(file.path(paste0(home,"/R/output/")))

##=================================================================##
##==================================================## covariate data
##=================================================================##

##================================================## pink salmon data
file_path<-paste0(home,"/R/data/pink_salmon_total_abundance.xlsx")

pink_salmon<-data.frame(read_excel(file_path))
pinks<-pink_salmon %>%
   dplyr::select(year=Year,pinks=Total)
add<-data.frame(year=seq(2016,2021,1),pinks=c(432,513,701,639,315,730))
pinks<-data.frame(rbind(pinks,add)) ## recent years (2021 preliminary)

pinks<-pinks %>%
   mutate(
      pinks_1=lead(pinks,1),
      pinks_2=lead(pinks,2),
      pinks_3=lead(pinks,3),
      pinks_4=lead(pinks,4)
   ) %>%
   mutate(across(-1,round,2))

##===================================================## get NPGO data
npgo<-read_table("http://www.o3d.org/npgo/npgo.php",skip=29,col_names=F,comment="#") %>%
   filter(!is.na(X2)) %>%
   dplyr::rename(Year=X1,Month=X2,NPGO=X3) %>%
   mutate(year=as.numeric(Year)) %>%
   group_by(year) %>%
   add_tally() %>%
   filter(!n<12) %>% ## use only complete years
   # filter(Month>3 & Month<9) %>% 
   group_by(year) %>%
   dplyr::summarise(NPGO=mean(NPGO)) %>%
   data.frame()

npgo<-npgo %>% 
   mutate(
      NPGO_1=lead(NPGO,1),
      NPGO_2=lead(NPGO,2),
      NPGO_3=lead(NPGO,3),
      NPGO_4=lead(NPGO,4),
   ) %>%
   mutate(across(-1,round,4))

##==================================## get ERSST data for coastal SST
data.dir<-paste0(home,"/R/data/")
thisYr<-as.numeric(format(Sys.Date(),"%Y"))
sst_dat<-get_ersst_data(years=c(1960,thisYr),data.dir=data.dir ,ncfilename="sst.mnmean.nc",latrange=c(46,49),lonrange=c(-124,-127))

sst_cst<-sst_dat %>% ## select 'meanSST' or 'gmeanSST' or 'resid':
   dplyr::select(year,month,meanSST) %>%
   data.frame() %>%
   mutate(month=paste0("m_",month)) %>%
   pivot_wider(names_from=month,values_from=meanSST) %>%
   rowwise() %>%
   mutate(SST_cst=(m_06+m_07+m_08)/3) %>%
   dplyr::select(year,SST_cst) %>%
   data.frame()

# sst_cst<-sst_dat %>% ## annual only
#    dplyr::select(year,meanSST) %>%
#    group_by(year) %>%
#    add_tally()%>%
#    filter(!n<12)%>% #use only complete years
#    dplyr::summarize(SST_cst=mean(meanSST))

sst_cst<-sst_cst %>%
   mutate(
      SST_cst_1=lead(SST_cst,1),
      SST_cst_2=lead(SST_cst,2),
      SST_cst_3=lead(SST_cst,3),
      SST_cst_4=lead(SST_cst,4),
   ) %>%
   mutate(across(-1,round,2))

##==================================## get ERSST data for summer SST
## SSTs are representative of the thermal conditions experienced by steelhead in the ocean, because steelhead are surface oriented and remain in the upper 20m with periodic dives to 40-60m (Burgner et al 1992, Walker et al. 2000) 

# sst_dat<-read.csv(paste0(home,"/R/data/ersstArc.raw.csv"))
# sst_dat<-sst_dat %>%
#    filter(year>1960) %>%
#    dplyr::select(year,month,sstarc) %>%
#    data.frame() %>%
#    mutate(month=paste0("m_",month)) %>%
#    pivot_wider(names_from=month,values_from=sstarc) %>%
#    rowwise() %>%
#    mutate(sst=(m_6+m_7+m_8)/3) %>% ## use summer (June/July/August)
#    dplyr::select(year,sst) %>%
#    data.frame()

sst_dat<-read.csv(paste0(home,"/R/data/ersst_sthd.csv"))

sst<-sst_dat %>%
   rename(SST=sst) %>%
   mutate(
      SST_1=lead(SST,1),
      SST_2=lead(SST,2),
      SST_3=lead(SST,3),
      SST_4=lead(SST,4),
   ) %>%
   mutate(across(-1,round,4))

##================================================## harbor seal data
## ONLY FOR TESTING _ MUCH FEWER YEARS OF DATA
file_path<-paste0(home,"/R/data/harbor_seal_abundances_updated.xlsx")
seal_data<-data.frame(read_excel(file_path)) %>%
   rename(year=Year) %>%
   rename(seals=Abundance) %>%
   dplyr::select(year,seals)

seals<-data.frame(year=all_years) %>%
   left_join(seal_data) %>%
   mutate(
      seals_2=lead(seals,2),
      seals_3=lead(seals,3),
      seals_4=lead(seals,4),
   ) %>%
   mutate(across(-1,round,2))

##----------------------------------------------------## correlations
df_dat<-pinks %>%
   left_join(npgo,by='year') %>%
   left_join(sst,by='year') %>%
   left_join(sst_cst,by='year') %>%
   dplyr::select(-c(matches("[[:digit:]]"))) %>%
   na.omit() ## no NAs in covariates

res<-data.frame(cor(as.matrix(df_dat),use="pairwise.complete.obs"))
res[-1,-1]

##=================================================================##
##=================================================================##
##======================================================## data plots
##=================================================================##
##=================================================================##

##=================================================================##
##==================================================## plot fish data
##=================================================================##
##---------------------------------------------------## color palette
# pic_path<-paste0(home,"/R/other/SteelheadPicture.jpg")
# cols_pic<-create_palette(image_path=pic_path,number_of_colors=30, type_of_variable="categorical")
# print(cols_pic)
colors<-rev(c("#A86260","#A5B1C4","#3B4D6B","#89A18D","#747260","#B3872D","#774C2C")) 

##================================================## spawners/runsize
colors <- colorRampPalette(colors)(nP)

##-------------------------------------## spawner abundances by river
pd1 <- fish_dat %>% 
   dplyr::select(pop,year,S_obs) %>%
   ggplot(aes(x=year,y=S_obs,color=pop)) +
   geom_line(aes(y=S_obs),lwd=0.5,alpha=1) +
   geom_point(aes(y=S_obs),size=1,alpha=1) +
   scale_color_manual(values=colors) +
   # scale_fill_manual(values=colors) +
   labs(x="Year",y="Spawner abundance") + 
   theme_sleek() +
   scale_y_continuous(limits=c(0,NA)) +
   theme(
      legend.position="none",
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      axis.title=element_text(size=12),
      axis.text=element_text(size=10)
   ) +
   facet_wrap(vars(pop),ncol=nP,scales="free_y") +
   NULL
ggsave("coastal-steelhead-spawners.pdf",pd1,width=3+2*nP,height=3)

##-----------------------------------------------## run size by river
pd2 <- fish_dat %>%
   dplyr::select(pop,year,S_obs,F_rate) %>%
   filter(!is.na(S_obs)) %>%
   mutate(runsize=round(S_obs/(1-F_rate))) %>%
   filter(!is.na(runsize)) %>%
   mutate(runsize=ifelse(runsize==S_obs,NA,runsize)) %>%
   ggplot(aes(x=year,y=runsize,color=pop)) +
   geom_line(aes(y=runsize),lwd=1,alpha=1) +
   scale_color_manual(values=colors) +
   labs(x="Year",y="Run size") +
   theme_sleek() +
   theme(
      legend.position=c(0.9,0.9),
      legend.title=element_blank(),
      axis.title=element_text(size=12),
      axis.text=element_text(size=10)
   ) +
   NULL
ggsave("coastal-steelhead-runsizes-oneplot.pdf",pd2,width=6,height=4)

pd3 <- fish_dat %>% 
   mutate(F_rate=ifelse(F_rate==0,NA,F_rate)) %>%
   ggplot(aes(x=year,y=F_rate,color=pop)) +
   geom_line(aes(y=F_rate),lwd=0.5,alpha=1) +
   geom_point(aes(y=F_rate),size=1,alpha=1) +
   scale_color_manual(values=colors) +
   # scale_fill_manual(values=colors) +
   labs(x="Year",y="Harvest rate") + 
   theme_sleek() +
   scale_y_continuous(limits=c(0,0.7)) +
   theme(
      legend.position="none",
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      axis.title=element_text(size=12),
      axis.text=element_text(size=10)
   ) +
   facet_wrap(vars(pop),ncol=nP,scales="free_y") +
   NULL
ggsave("coastal-steelhead-harvest-rates.pdf",pd3,width=3+2*nP,height=3)

pd4 <- fish_dat %>% 
   dplyr::select(pop,year,S_obs,F_rate) %>%
   filter(!is.na(S_obs)) %>%
   mutate(runsize=round(S_obs/(1-F_rate))) %>%
   filter(!is.na(runsize)) %>%
   mutate(runsize=ifelse(runsize==S_obs,NA,runsize)) %>%
   ggplot(aes(x=year,y=runsize,color=pop)) +
   #geom_line(aes(y=runsize),lwd=0.7,col="black") +
   geom_line(aes(y=runsize),lwd=0.6,alpha=0.8) +
   #geom_line(aes(y=S_obs),lwd=0.7,color="black",linetype="dashed") +
   geom_line(aes(y=S_obs),lwd=0.6,linetype="dashed") +
   #geom_point(aes(y=runsize),size=1,alpha=1) +
   scale_color_manual(values=colors) +
   # scale_fill_manual(values=colors) +
   labs(x="Year",y="Abundance") + 
   theme_sleek() +
   scale_y_continuous(limits=c(0,NA)) +
   theme(
      legend.position="none",
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      axis.title=element_text(size=12),
      axis.text=element_text(size=10)
   ) +
   facet_wrap(vars(pop),ncol=nP,scales="free_y") +
   NULL
ggsave("coastal-steelhead-runsizes.pdf",pd4,width=3+2*nP,height=3)

##---------------------------------------------## run size cumulative
fish_dat_wide<-fish_dat %>% 
   mutate(runsize=round(S_obs/(1-F_rate))) %>%
   dplyr::select(pop,year,runsize) %>%
   pivot_wider(names_from=pop,values_from=runsize) %>%
   na.omit() %>%
   data.frame()

pdf("coastal-steelhead-cumulative-runsize.pdf",width=6,height=4)
par(mar=c(4,4,1,1),mgp=c(2,0.5,0),tck=-0.02,cex.lab=1.2,cex.axis=0.8, xaxs="i",yaxs="i")
ymax<-max(rowSums(fish_dat_wide))
yrs<-unique(fish_dat_wide$year)
plot(NA,NA,xlim=c(min(yrs),max(yrs)),ylim=c(0,ymax),xlab="Year",ylab="Run size")
poly.x<-c(yrs,rev(yrs))
for(i in 1:nP){
   if(i==1) {
      run_added<-as.vector(fish_dat_wide[,i+1])
      poly.y<-c(run_added,rep(0,length(yrs)))
   }
   if(i!=1) {
      run_prev<-run_added
      run_added<-run_prev+fish_dat_wide[,i+1]
      poly.y<-c(run_added,rev(run_prev))
   }
   polygon(poly.x,poly.y,lwd=0.5,col=colors[i],border=colors[i])
}
box()
legend("topright",rev(names(fish_dat_wide)[-1]),pch=15,pt.cex=1,cex=0.6,adj=0,bty="n",col=rev(colors),xpd=T,inset=c(0,0),y.intersp=1.1)
dev.off()

##==================================================## age proportions
ages<-colnames(fish_dat)[grepl("age",colnames(fish_dat))]

##---------------------## only plot years with at least XX age samples
min_age_samples<-10
age_dat<-fish_dat
age_dat[which(rowSums(age_dat[,ages])<min_age_samples),ages]<-NA

##--------------------------## age proportions for maiden and repeats
age_prop <- data.frame(cbind(age_dat %>% dplyr::select(pop,year),prop.table(age_dat[,names(age_dat) %in% ages] %>% replace(is.na(.),0) %>% as.matrix(),1)))

# plot_pops<-c("Quillayute")
# npops<-length(plot_pops)
# age_prop<-age_prop %>% filter(pop %in% plot_pops)

##-------------------## only keep age groups with prop>0.001 for plot
av_prop<-round(colMeans(age_prop[,-c(1,2)],na.rm=T),4)
drop<-names(av_prop)[av_prop<0.001] ## at least 0.1% of samples
age_prop<-age_prop %>% dplyr::select(-all_of(drop))

##------------------## drop years with certain age props equal to one
# drop <- age_prop %>% filter(if_any(.cols=everything(),.fns=~.==1))
drop <- age_prop %>% filter_all(any_vars(.==1))
age_prop <- anti_join(age_prop, drop, by = c("pop","year"))

##--------------------------------------------------------## reformat
age_data <- age_prop  %>%
   pivot_longer(3:last_col(),names_to="age",values_to="prop") %>%
   mutate(MR=substr(age,7,7)) %>%
   mutate(MR=gsub("K","repeats",MR)) %>%
   mutate(MR=gsub("M","maiden",MR)) %>%
   mutate(age=substr(age,6,6)) 

nA<-length(unique(age_data$age))
colors <- colorRampPalette(brewer.pal(name="Set1",n=9))(nA+1)[-1]

##----------------------------## plot age proportions as points/lines
## proportions across all spawners
pd5<-age_data %>%
   ggplot(aes(x=year,y=prop,color=age)) +
   geom_line(lwd=0.2,alpha=1) +
   geom_point(pch=16,size=1,alpha=1) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   labs(x="Year",y="Proportion")+ 
   theme_sleek() +
   theme(
      axis.title=element_text(size=12),
      axis.text=element_text(size=6),
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)
   ) +
   facet_grid(MR~pop,scales="free_y") +
   NULL
ggsave("coastal-steelhead-age-data.pdf",pd5,width=4+2*nP,height=4)

##=================================================================##
##============================================## plot covariate data
##=================================================================##
pc1<-ggplot(pinks,aes(x=year,y=pinks))+geom_line()+theme_sleek()
pc2<-ggplot(npgo,aes(x=year,y=NPGO))+geom_line()+theme_sleek()
pc3<-ggplot(sst_cst,aes(x=year,y=SST_cst))+geom_line()+theme_sleek()
pc4<-ggplot(sst,aes(x=year,y=SST))+geom_line()+theme_sleek()
# pc5<-ggplot(seals,aes(x=year,y=seals_OP))+geom_line()+theme_sleek()
# pc6<-ggplot(seals,aes(x=year,y=seals_CE))+geom_line()+theme_sleek()
# gg<-grid.arrange(pc1,pc2,pc3,pc4,pc5,pc6,ncol=2)
gg<-grid.arrange(pc1,pc2,pc3,pc4,ncol=2)
ggsave("covariate-data.pdf",gg,width=6,height=6)

##=================================================================##
##=================================================================##
##=====================================## Integrated Population Model
##=================================================================##
##=================================================================##

##=================================================================##
##=========================================================## fit IPM
##=================================================================##
SR_mod<-"Ricker" ## "Ricker" or "BH"
covar_effects<-TRUE
##----------------------------------------------------------## priors
mu_p_vec<-c(0.01,0.02,0.47,0.44,0.05,0.01)
priors<- list(
   mu_alpha ~ normal(1.5,0.5), ## hyper mean / SD of log productivity
   mu_p ~ dirichlet(mu_p_vec*100), ## mean maiden age distribution
   # tau ~ gnormal(1,0.85,30), ## default: gnormal(1,0.85,30))
   mu_SS ~ beta(1.5,3) ## mean kelt survival rate
)
##----------------------------------------------## add covariate data
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
##--------------------------------------------------## years included
dat_years<-sort(unique(fish_dat$year))
nY<-length(dat_years)
##------------------------------------## save fish data used in model
write.csv(fish_dat,"IPM_fish_dat.csv",row.names=F)
##-------------------------------------------## read previous IPM fit
# IPM_fit<-readRDS(paste0(home,"/R/output/IPM_fit.Rdata"))
# fish_dat<-read.csv(paste0(home,"/R/output/IPM_fish_dat.csv"))

##=======================================================## fit model
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
## stan_pars("IPM_SSiter_pp") ## parameters monitored
##-----------------------------------------------------## time elapsed
## ~1h for test run; ~3h for a final run
print(paste(round(max(rowSums(get_elapsed_time(IPM_fit))/60),2),"min"))
##-----------------------------------------------------## save results
saveRDS(IPM_fit,"IPM_fit.Rdata")

##=================================================================##
##=========================================================## results
##=================================================================##
df_post<-as.data.frame(IPM_fit) ## posterior of all parameters
nS<-dim(df_post)[1] ## number of samples
posterior<-as.array(IPM_fit) ## get full posterior 
monitored<-dimnames(posterior)$parameters
# monitored[grepl("_R",monitored)] 
df_out<-data.frame(summary(IPM_fit)$summary) ## summary all parameters
draws<-as_draws_rvars(IPM_fit) ## all draws as multidimensional array
# probs<-c(0.025,0.25,0.5,0.75,0.975) ## median with 50% and 95% CIs
probs<-c(0.05,0.25,0.5,0.75,0.95) ## median with 50% and 90% CIs

##=====================================================## plot colors
colors<-colorRampPalette(brewer.pal(name="Set1",n=9))(nP)
# colors <- c("goldenrod1","chocolate2","firebrick2", "darkorchid3", "midnightblue","lightblue","slategray")

##==================================================## plot functions
summary_CI95<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.025),na.rm=T),ymax=quantile(x,prob=c(0.975),na.rm=T))) }
summary_CI90<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.05),na.rm=T),ymax=quantile(x,prob=c(0.95),na.rm=T))) }
summary_CI50<-function(x) { return(data.frame(y=median(x,na.rm=T),ymin=quantile(x,prob=c(0.25),na.rm=T),ymax=quantile(x,prob=c(0.75),na.rm=T))) }

##=====================================## plot diagnostics and priors
plot_diagnostics(IPM_fit,fish_dat)
plot_priors_vs_posteriors(IPM_fit,fish_dat,pdf=TRUE)

tau_post<-extract1(IPM_fit,"tau")
tau_median<-median(tau_post)

##=================================================================##
##=========================================## population productivity
##=================================================================##

##-------------------------------------------## alpha hyper parameter
mu_alpha_post<-data.frame(exp(extract1(IPM_fit,"mu_alpha")))
names(mu_alpha_post)<-" mu_alpha "
mu_alpha_qs<-t(apply(mu_alpha_post,2,function(x) quantile(x,prob=probs)))

##-------------------------------------## population parameters alpha 
alphas_post<-data.frame(extract1(IPM_fit,"alpha"))
names(alphas_post)<-pops
alpha_medians<-apply(alphas_post,2,median)
alpha_means<-apply(alphas_post,2,mean)

alphas_qs<-t(apply(alphas_post,2,function(x) quantile(x,prob=probs))) %>% data.frame() %>% add_column(pop=pops)

p01 <- alphas_post %>% 
   pivot_longer(col=everything(),names_to="name",values_to="value") %>%
   ggplot(aes(x=name,y=value,color=name)) +
   #geom_violin(lwd=0.1) +
   stat_summary(fun.data=summary_CI90,size=0.25,shape=20) +
   # stat_summary(fun.data=summary_CI95,size=0.25,shape=20) +
   stat_summary(fun.data=summary_CI50,size=0.75,shape=20) +  
   scale_y_continuous(limits=c(0,5)) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   geom_hline(yintercept=1,linetype="dashed",size=0.2) +
   labs(x="",y="Popuation productivity") + 
   theme_sleek() +
   theme(
      legend.position="none",
      axis.title.x=element_blank(),
      plot.margin=unit(c(1,1,1,1),"lines"),
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)
      ) +
   NULL
ggsave("IPM-sthd-productivity-all-pops.pdf",p01,width=3,height=3)

p02 <- mu_alpha_post %>%
   pivot_longer(col=everything(),names_to="name",values_to="value") %>%
   ggplot(aes(x=name,y=value)) +
   #geom_violin(lwd=0.1) +
   stat_summary(fun.data=summary_CI90,size=0.25,shape=20) +
   # stat_summary(fun.data=summary_CI95,size=0.25,shape=20) +
   stat_summary(fun.data=summary_CI50,size=0.75,shape=20) +  
   scale_y_continuous(limits=c(0,5)) +
   scale_color_manual(values=c(colors,"black")) +
   scale_fill_manual(values=c(colors,"black")) +
   geom_hline(yintercept=1,linetype="dashed",size=0.2) +
   labs(x=paste0("hyper-\nparameter"),y="") + 
   theme_sleek() +
   theme(
      legend.position="none",
      plot.margin=unit(c(1,1,0.25,-0.5),"lines"), ## 3rd mar: 1.9 | 0.3
      axis.text.x=element_blank(),
      axis.title.x=element_text(angle=90,vjust=0.5,hjust=1),
      axis.ticks.x=element_blank()
   ) +
   NULL
ggsave("IPM-sthd-productivity-hyperparameter.pdf",p02,width=2,height=3)

p03<-grid.arrange(p01,p02,nrow=1,widths=c(3,1))
ggsave("IPM-sthd-productivity-combined.pdf",p03,width=5,height=3)

##=================================================================##
##=============================================## population capacity
##=================================================================##
per_area<-length(unique(fish_dat$A))>1 ## capacity per area or not ?
   
##--------------------------------------------## Rmax hyper parameter
mu_Rmax_post<-data.frame(exp(extract1(IPM_fit,"mu_Rmax")))
names(mu_Rmax_post)<-" mu_Rmax"
mu_Rmax_qs<-t(apply(mu_Rmax_post,2,function(x) quantile(x,prob=probs)))

##----------------------------------------## population-specific Rmax
Rmax_post<-data.frame(extract1(IPM_fit,"Rmax"))
names(Rmax_post)<-pops

# for(i in 1:nP){ Rmax_post[,i]<-Rmax_post[,i]*as.numeric(areaD$area[areaD$pop==pops[i]]) } ## convert to total capacity 

Rmax_medians<-apply(Rmax_post,2,median)
Rmax_means<-apply(Rmax_post,2,mean)

Rmax_qs<-apply(Rmax_post,2,function(x) quantile(x,prob=probs))

if(per_area) { 
   ylab<-expression(paste("Popuation capacity per km"^"2")) 
}else{ 
   ylab<-"Popuation capacity" 
} 

p04 <- Rmax_post %>% 
   pivot_longer(col=everything(),names_to="pop",values_to="Rmax") %>%
   ggplot(aes(x=pop,y=Rmax,color=pop)) +
   # geom_violin(lwd=0.1) +
   stat_summary(fun.data=summary_CI90,size=0.2,shape=20) +
   # stat_summary(fun.data=summary_CI95,size=0.2,shape=20) +
   stat_summary(fun.data=summary_CI50,size=0.7,shape=20) +  
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   geom_hline(yintercept=0,linetype="dashed",size=0.2) +
   labs(x="",y=ylab) + 
   theme_sleek() +
   theme(
      legend.position="none",
      axis.title.x=element_blank(),
      plot.margin=unit(c(1,1,1,1),"lines"),
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)
   ) +
   NULL
ggsave("IPM-sthd-capacity-all-pops.pdf",p04,width=3,height=3)

p05 <- data.frame(mu_Rmax_post) %>% 
   ggplot(aes(x=1,y=X.mu_Rmax)) +
   # geom_violin(lwd=0.1) +
   stat_summary(fun.data=summary_CI90,size=0.2,shape=20) +
   # stat_summary(fun.data=summary_CI95,size=0.2,shape=20) +
   stat_summary(fun.data=summary_CI50,size=0.7,shape=20) +  
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   geom_hline(yintercept=0,linetype="dashed",size=0.2) +
   # labs(x=expression(mu[" Rmax"]),y="") + 
   labs(x=paste0("hyper-\nparameter"),y="") + 
   theme_sleek() +
   theme(
      legend.position="none",
      plot.margin=unit(c(1,1,0.25,-0.5),"lines"), ## 3rd mar: 1.9 | 0.3
      axis.text.x=element_blank(),
      axis.title.x=element_text(angle=90,vjust=0.5,hjust=1),
      axis.ticks.x=element_blank()
   ) +
   NULL
ggsave("IPM-sthd-capacity-hyperparameter.pdf",p05,width=3,height=3)

p05b<-grid.arrange(p04,p05,nrow=1,widths=c(3,1))
ggsave("IPM-sthd-capacity-combined.pdf",p05b,width=5,height=3)

if(per_area) {
   p06<-Rmax_post %>% 
      pivot_longer(col=everything(),names_to="pop",values_to="Rmax")%>%
      mutate(system_Rmax=case_when(
         pop==pops[1]~Rmax*as.numeric(areaD$area[areaD$pop==pops[1]]),
         pop==pops[2]~Rmax*as.numeric(areaD$area[areaD$pop==pops[2]]),
         pop==pops[3]~Rmax*as.numeric(areaD$area[areaD$pop==pops[3]]),
         pop==pops[4]~Rmax*as.numeric(areaD$area[areaD$pop==pops[4]])
         #pop==pops[5]~Rmax*as.numeric(areaD$area[areaD$pop==pops[5]]),
         #pop==pops[6]~Rmax*as.numeric(areaD$area[areaD$pop==pops[6]]),
         #pop==pops[7]~Rmax*as.numeric(areaD$area[areaD$pop==pops[7]])
      )) %>%
      ggplot(aes(x=pop,y=system_Rmax,color=pop)) +
      # geom_violin(lwd=0.1) +
      stat_summary(fun.data=summary_CI90,size=0.2,shape=20) +
      # stat_summary(fun.data=summary_CI95,size=0.2,shape=20) +
      stat_summary(fun.data=summary_CI50,size=0.7,shape=20) +  
      scale_color_manual(values=colors) +
      scale_fill_manual(values=colors) +
      geom_hline(yintercept=0,linetype="dashed",size=0.2) +
      labs(x="",y="Popuation capacity") + 
      theme_sleek() +
      theme(
         legend.position="none",
         axis.title.x=element_blank(),
         plot.margin=unit(c(1,1,1,1),"lines"),
         axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)
      ) +
      NULL
   ggsave("IPM-sthd-capacity-full-system.pdf",p06,width=3,height=3)
   
   p06b<-grid.arrange(p04,p06,nrow=1)
   ggsave("IPM-sthd-capacity-w-wo-per-area.pdf",p06b,width=5,height=3)
   
}

##=================================================================##
##==============================================## plot age structure
##=================================================================##
mu_p_post<-data.frame(extract1(IPM_fit,"mu_p"))
mu_p_med<-round(data.frame(t(apply(mu_p_post,2,median))),2)
mu_p_med
mu_p_vec ## vector of proportions used for prior
## if prior is mu_p_vec*1000, posterior is basically equal to prior

mu_p_qs<-t(apply(mu_p_post,2,function(x) quantile(x,prob=probs))) %>%
   data.frame() %>%
   mutate(across(where(is.numeric),round,3))

sigma_p_pop_post<-data.frame(extract1(IPM_fit,"sigma_pop_p"))#n_age-1
sigma_p_pop_med<-data.frame(t(apply(sigma_p_pop_post,2,median)))

sigma_p_post<-data.frame(extract1(IPM_fit,"sigma_p"))
sigma_p_med<-round(data.frame(t(apply(sigma_p_post,2,median))),3)

##-------------------------------------------## age structure by year
## q: true maiden and kelt age distributions
q_table<-df_out[grep("^(?=.*q)(?!.*_)",rownames(df_out),perl=TRUE),]
##--------------------------------------------## get mean proportions
N<-dim(fish_dat)[1]
m_names<-colnames(fish_dat)[grepl("M_obs",colnames(fish_dat))]
r_names<-colnames(fish_dat)[grepl("K_obs",colnames(fish_dat))]
all_ages<-gsub("n_age","",gsub("_obs","",c(m_names,r_names)))
n_age<-length(all_ages)
   
q_est<-matrix(q_table$mean,ncol=n_age,nrow=N,byrow=TRUE) %>% 
   data.frame() %>%
   rename_with(~all_ages,all_of(names(.))) %>%
   add_column(fish_dat%>%dplyr::select(pop,year)) %>%
   pivot_longer(1:(last_col()-2),names_to="age",values_to="prop") %>%
   data.frame()
##---------------------------------------## get age and maiden/repeat
age_est<-q_est %>%
   mutate(MR=substr(age,2,2)) %>%
   mutate(MR=gsub("K","repeats",MR)) %>%
   mutate(MR=gsub("M","maiden",MR)) %>%
   mutate(age=substr(age,1,1))
##-----------------------------------------## add full table with CIs
age_df<-data.frame(cbind(age_est,q_table)) %>% 
   dplyr::select(-prop) %>%
   filter(age>2)
##-----------------------------------------------------## plot colors
nA<-length(unique(age_df$age))
cols <- colorRampPalette(brewer.pal(name="Set1",n=9))(nA+1)[-1]
# cols <- c("goldenrod1","chocolate2","firebrick2", "darkorchid3", "midnightblue","lightblue")
##----------------------------------## plot estimated age proportions
p07<-age_df %>%
   ggplot(aes(x=year,y=mean,color=age,fill=age)) +
   geom_line(lwd=0.2,alpha=1) +
   geom_point(pch=16,size=1,alpha=1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=X2.5.,ymax=X97.5.),color=NA,alpha=0.1) +
   scale_color_manual(values=cols) +
   scale_fill_manual(values=cols) +
   labs(x="Year",y="Proportion")+ 
   theme_sleek() +
   theme(
      axis.title=element_text(size=12),
      axis.text=element_text(size=6),
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)
   ) +
   facet_grid(MR~pop,scales="free_y") + 
   NULL
ggsave("IPM-sthd-age-proportions.pdf",p07,width=2+2*nP,height=4)

##=================================================================##
##===========================================## spawners and recruits
##=================================================================##

##================================================## recruit/spawners
## iteroparous steelhead: recruits/spawner at equilibrium is < 1 (assuming kelt survival rate is larger than 0)
## only the expected life-time recruits/spawner is 1 at equilibrium (includes recruits from repeat spawning events)
## life-time RpS at equilibrium = maiden RpS - kelt survival rate (given that kelt survival rate is time-invariant and age-independent)

##------------------------------------------## get posterior for S/R
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

## drop years with incomplete recruitment when no covariates in model
## i.e., only recruitment for which dominant ages have been observed

if(covar_effects) { 
   RpS_plot<-RpS_plot 
}else{ 
   RpS_plot<-RpS_plot %>% 
      group_by(pop) %>%
      slice(1:(n()-5)) %>%
      data.frame()
} 

##------------------------------## plot time-varying recruits/spawner
p <- RpS_plot  %>%
   ggplot(aes(x=year,y=RpS_50.,color=pop,fill=pop)) +
   # geom_hline(yintercept=1.0,lwd=0.1,col="gray",linetype="dashed") +
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
ggsave("IPM-sthd-recruits-per-spawner-w-RpSeq.pdf",p,width=1+2*nP,height=2.5)

##------------------------------## plot time-varying recruits/spawner
p <- RpS_plot %>% 
   # left_join(RpSeq,by=c("pop","year")) %>%
   ggplot(aes(x=year,y=RpS_50.,color=pop,fill=pop)) +
   ## lines for time-invariant kelt survival rates
   geom_hline(yintercept=1.0,lwd=0.2,col="gray",alpha=0.1) +
   geom_hline(yintercept=0.9,lwd=0.2,col="gray",alpha=0.2) +
   geom_hline(yintercept=0.8,lwd=0.2,col="gray",alpha=0.3) +
   geom_hline(yintercept=0.7,lwd=0.2,col="gray",alpha=0.4) +
   geom_hline(yintercept=0.6,lwd=0.2,col="gray",alpha=0.5) +
   geom_hline(yintercept=0.5,lwd=0.2,col="gray",alpha=0.6) +
   geom_hline(yintercept=0.4,lwd=0.2,col="gray",alpha=0.7) +
   geom_hline(yintercept=0.3,lwd=0.2,col="gray",alpha=0.8) +
   geom_hline(yintercept=0.2,lwd=0.2,col="gray",alpha=0.9) +
   geom_hline(yintercept=0.1,lwd=0.2,col="gray",alpha=1.0) +
   ## estimated recruits per spawner median, 50%, and 90% CIs
   geom_line(aes(y=RpS_50.),lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=RpS_25.,ymax=RpS_75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=RpS_5.,ymax=RpS_95.),color=NA,alpha=0.1) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   labs(x="Brood year",y="Recruits/spawner") + 
   scale_y_continuous(limits=c(0,NA),expand=c(0,0)) +
   facet_wrap(vars(pop),ncol=nP) + #,scales="free_y") +
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
ggsave("IPM-sthd-recruits-per-spawner.pdf",p,width=1+2*nP,height=2.5)

##-------------------------## plot time-varying log(recruits/spawner)
p <- RpS_plot %>%
   mutate(across(!year & !pop, log)) %>%
   ggplot(aes(x=year,y=RpS_50.,color=pop,fill=pop)) +
   geom_line(aes(y=RpS_50.),lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=RpS_25.,ymax=RpS_75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=RpS_5.,ymax=RpS_95.),color=NA,alpha=0.1) +
   ## add iteroparous recruits/spawner at equilibrium/replacement ??? 
   geom_line(aes(y=RpS_eq),lwd=0.5,alpha=0.5) +
   #geom_line(aes(y=RpSeq_50.),lwd=0.25,alpha=0.5) +
   #geom_ribbon(aes(ymin=RpSeq_25.,ymax=RpSeq_75.),color=NA,alpha=0.1)+
   #geom_ribbon(aes(ymin=RpSeq_5.,ymax=RpSeq_95.),color=NA,alpha=0.05)+
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   labs(x="Year",y="log(recruits/spawner)") + 
   # scale_y_continuous(limits=c(0,NA)) +
   facet_wrap(vars(pop),ncol=nP) + #,scales="free_y") +
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
ggsave("IPM-sthd-recruits-per-spawner-log.pdf",p,width=1+2*nP,height=2.5)

##----------------------------------------------## estimated spawners
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
ggsave("IPM-sthd-spawners.pdf",p1,width=1+2*nP,height=2.5)

##----------------------## plot recruits/spawner vs spawner abundance 
p1b <- RpS_qs %>% 
   left_join(S_est_qs) %>%
   ggplot(aes(x=X50.,y=RpS_50.,color=pop,fill=pop)) +
   geom_point(lwd=0.5,alpha=1) + ## aes(x=X50.,y=RpS_50.)
scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   labs(x="Year",y="Recruits/spawner") + 
   scale_x_continuous(limits=c(0,NA)) +
   scale_y_continuous(limits=c(0,NA)) +
   facet_wrap(vars(pop),ncol=nP,scales="free_x") +
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
ggsave("IPM-sthd-RpS-vs-spawner-abundance.pdf",p1b,width=1+2*nP,height=2.5)

##----------------------------------------------## estimated recruits
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
ggsave("IPM-sthd-recruits.pdf",p2,width=1+2*nP,height=2.5)


##================================================## reference points

##--------------------------------------------------## time-invariant
df_MSY<-ref_points_itero_multi_pop(IPM_fit,SR_mod,fish_dat,per_area,areaD)

Smsy_qs<-df_MSY %>%
   data.frame() %>%
   group_by(pop) %>%
   summarise(q2.5=quantile(S_MSY,0.025),
             q5=quantile(S_MSY,0.05),
             q25=quantile(S_MSY,0.25),
             q50=quantile(S_MSY,0.50),
             q75=quantile(S_MSY,0.75),
             q95=quantile(S_MSY,0.95),
             q97.5=quantile(S_MSY,0.975)) 

MSY_medians<-df_MSY %>%
   data.frame() %>%
   group_by(pop) %>%
   summarise(S_MSY=median(S_MSY),
             U_MSY=median(U_MSY),
             R_MSY=median(R_MSY),
             S_zero=median(S_zero),
             R_zero=median(R_zero))

##----------------------------------------------------## time-varying
df_MSY_tv<-ref_points_itero_multi_pop_time_var(IPM_fit,SR_mod,fish_dat,per_area,areaD)

Smsy_tv_qs<-df_MSY_tv %>%
   data.frame() %>%
   group_by(pop,year) %>%
   summarise(q2.5=quantile(S_MSY,0.025),
             q5=quantile(S_MSY,0.05),
             q25=quantile(S_MSY,0.25),
             q50=quantile(S_MSY,0.50),
             q75=quantile(S_MSY,0.75),
             q95=quantile(S_MSY,0.95),
             q97.5=quantile(S_MSY,0.975)) %>%
   ungroup() %>%
   mutate_at('year',as.numeric)
   
MSY_tv_medians<-df_MSY_tv %>%
   data.frame() %>%
   group_by(pop,year) %>%
   summarise(S_MSY=median(S_MSY),
             U_MSY=median(U_MSY),
             R_MSY=median(R_MSY),
             S_zero=median(S_zero),
             R_zero=median(R_zero)) %>%
   ungroup() %>%
   mutate_at('year',as.numeric)

##------------------------------------------## plot time-varying Smsy
p <- Smsy_tv_qs %>%
   ggplot(aes(x=year,y=q50,color=pop,fill=pop)) +
   geom_line(aes(y=q50),lwd=0.5,alpha=1) +
   #geom_hline(yintercept=as.numeric(surv_med_recent_mean),lwd=0.5) +
   geom_ribbon(aes(ymin=q25,ymax=q75),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=q5,ymax=q95),color=NA,alpha=0.1) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   #scale_x_continuous(limits=c(1978,2022)) +
   #scale_y_continuous(limits=c(0,0.65),breaks=seq(0,0.6,0.1),expand=c(0,0)) +
   labs(x="Year",y="S_msy") + 
   facet_wrap(vars(pop),ncol=nP,scales="free_y") + 
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
ggsave("IPM-sthd-Smsy-time-varying.pdf",p,width=2+2*nP,height=2.9)

##================================================## SR relationships
Rmax_post<-data.frame(extract1(IPM_fit,"Rmax"))
names(Rmax_post)<-pops
if(per_area){ 
   for(i in 1:nP){ 
      Rmax_post[,i]<-Rmax_post[,i]*as.numeric(areaD$area[areaD$pop==pops[i]]) 
   } 
}
Rmax_medians<-apply(Rmax_post,2,median)
Rmax_means<-apply(Rmax_post,2,mean)
Rmax_qs<-apply(Rmax_post,2,function(x) quantile(x,prob=probs))

alpha_Rmax_pops<-data.frame(alpha=alpha_medians,R_max=Rmax_medians)

# p3<-plot_instr( code={ ## instead of next line for ppt slides
pdf("IPM-sthd-spawner-recruit-curves.pdf",width=3*nP,height=3.5)
par(mar=c(1,1,3,1),oma=c(4,4,0,0),mgp=c(2,0.5,0))
par(tck=-0.02,cex.lab=1.2,xaxs="i",yaxs="i")
layout(matrix(1:nP,ncol=nP,nrow=1,byrow=T))
for(j in 1:nP){
   ## observed fish data
   pop_dat<-fish_dat[fish_dat$pop==pops[j],]
   ## estimated spawners and recruits with uncertainty
   S_est<-S_est_qs[S_est_qs$pop==pops[j],]
   R_est<-R_est_qs[R_est_qs$pop==pops[j],]
   # xylim<-c(0,max(1.5*max(pop_dat$S_obs,na.rm=T),Rmax_qs[4,j]))
   # xylim<-c(0,round(1.1*max(rbind(S_est[,1:5],R_est[,1:5]),na.rm=T)))
   # plot(NA,NA,xlab="",ylab="",xlim=xylim,ylim=xylim)
   xlim<-c(0,round(1.1*max(S_est[,1:5],na.rm=T)))
   ylim<-c(0,round(1.1*max(R_est[,1:5],na.rm=T)))
   plot(NA,NA,xlab="",ylab="",xlim=xlim,ylim=ylim)
   ## random draws from the posterior
   nSp<-100
   index<-sample(seq(nS),nSp,replace=F)
   for(i in 1:nSp){ curve(SR(SR_fun=SR_mod,alpha=alphas_post[index[i],j],Rmax=Rmax_post[index[i],j],S=x),add=TRUE,col="gray",lwd=0.1) }
   ## median stock-recruit relationship
   curve(SR(SR_fun=SR_mod,alpha=alpha_medians[j],Rmax=Rmax_medians[j],S=x), lwd=2,add=TRUE)
   abline(0,1,lwd=1,lty=3)
   ## add vertical lines for escapement at MSY estimates
   smsy<-MSY_medians$S_MSY[MSY_medians$pop==pops[j]]
   abline(v=smsy,lwd=1,col="goldenrod")
   ## estimated spawners and recruits with uncertainty
   S_est<-S_est_qs[S_est_qs$pop==pops[j],]
   R_est<-R_est_qs[R_est_qs$pop==pops[j],]
   points(S_est[,3],R_est[,3],pch=16,cex=1.2)
   text(S_est[,3],R_est[,3],labels=paste0(S_est$year,"\n "),cex=0.5,pos=4,offset=0.2)
   segments(S_est[,1],R_est[,3],S_est[,5],R_est[,3],lwd=0.1)
   segments(S_est[,2],R_est[,3],S_est[,4],R_est[,3],lwd=1)
   segments(S_est[,3],R_est[,1],S_est[,3],R_est[,5],lwd=0.1)
   segments(S_est[,3],R_est[,2],S_est[,3],R_est[,4],lwd=1)
   box()
   mtext(paste(pops[j]),side=3,line=0.5,cex=1,font=2,outer=F)
}
mtext("Spawners",side=1,line=2,cex=1.2,font=1,outer=T)
mtext("Recruits",side=2,line=2,cex=1.2,font=1,outer=T)
dev.off()
# }) ## instead of previous line for ppt slides

# pp <- read_pptx("../template.pptx") %>%
#   ph_with(value=p3,ph_location(left=0.5,top=1,width=12,height=2.4))%>%
#   print(target="IPM_sthd_model_SR_plots.pptx")

##----------------## median spawners and recruits with year indicated
SR_est_qs<-S_est_qs %>% 
   rename_all(list(~stringr::str_replace_all(.,'X','S_'))) %>%
   full_join(R_est_qs,by=c('pop','year')) %>%
   rename_all(list(~stringr::str_replace_all(.,'X','R_'))) %>%
   data.frame()

p3 <- SR_est_qs %>%
   ggplot(aes(x=S_50.,y=R_50.,color=year,fill=year)) +
   geom_point(size=2,alpha=1) +
   labs(x="Spawners",y="Recruitment") + 
   scale_x_continuous(limits=c(0,NA),expand=c(0,1000)) +
   scale_y_continuous(limits=c(0,NA),expand=c(0,1000)) +
   facet_wrap(vars(pop),ncol=nP,scales="free") +  # scales="free_x"
   theme_sleek() +
   # theme(legend.position="none") +
   theme(
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15)
   ) +
   NULL
ggsave("IPM-sthd-spawner-recruit-curve.pdf",p3,width=1+3*nP,height=3.5)

##=================================================================##
##=============================================## kelt survival rates
##=================================================================##
mu_surv_post<-extract1(IPM_fit,"mu_SS")
mu_surv_qs<-t(quantile(mu_surv_post,prob=probs)) ## quantiles
surv_post<-extract1(IPM_fit,"s_SS") ## survival rates by population

##---------------------## median survival rate by population and year
surv_med<-apply(surv_post,2,median) %>% 
   data.frame() %>%
   rename(surv='.') %>%
   add_column(fish_dat%>%dplyr::select(pop,year)%>%data.frame())

surv_annual_means_of_medians<-surv_med %>%
   pivot_wider(names_from=pop,values_from=surv) %>%
   arrange(year) %>% 
   mutate(mean_across_pops=rowMeans(across(!year))) %>%
   dplyr::select(year,mean_across_pops) %>%
   data.frame()

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
surv_recent_mean<-surv_med %>%
   filter(year>(max(surv_med$year-5))) 

surv_med_recent_mean<-surv_recent_mean %>%
   group_by(pop) %>%
   dplyr::summarize(recent_mean=mean(surv))

## select specific populations (those with single-population models)
# keep_pops<-c("Quillayute","Queets","Quinault","Humptulips","Chehalis")
# surv_qs<-surv_qs %>% filter(pop %in% keep_pops)

##------------------------------------## plot kelt survival over time
p4 <- surv_qs %>%
   ggplot(aes(x=year,y=X50.,color=pop,fill=pop)) +
   geom_line(aes(y=X50.),lwd=0.5,alpha=1) +
   #geom_hline(yintercept=as.numeric(surv_med_recent_mean),lwd=0.5) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.4) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.2) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   scale_x_continuous(limits=c(1978,2022)) +
   #scale_y_continuous(limits=c(0,0.65),breaks=seq(0,0.6,0.1),expand=c(0,0)) +
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
ggsave("IPM-sthd-kelt-survival-rates.pdf",p4,width=1+2*nP,height=2.5)

##===========================================## kelt survival anomaly
## 'eta_year_SS' are logit kelt survival anomalies
eta_SS_post<-extract1(IPM_fit,"eta_year_SS")
eta_SS_qs<-t(apply(eta_SS_post,2,function(x) quantile(x,prob=probs)))
eta_SS_qs<-eta_SS_qs %>%
   data.frame() %>%
   add_column(year=dat_years) # %>%
   # slice(-1) %>% ## first year
   # slice(-n()) ## last year

## estimate trend over time
etas_SS<-data.frame(median=apply(eta_SS_post,2,median),sd=apply(eta_SS_post,2,sd)) %>% add_column(year=dat_years)

newdata<-data.frame(year=dat_years)
lm_fit<-lm(etas_SS$median~etas_SS$year,weights=1/(etas_SS$sd^2))

pred<-predict(lm_fit,newdata,interval="confidence",level=0.95)
pred<-data.frame(pred) %>% add_column(year=etas_SS$year) 

p5 <- eta_SS_qs %>%
   ggplot(aes(x=year,y=X50.)) +
   geom_line(lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.1) +
   labs(x="Year",y="Kelt survival anomaly") + 
   ## add weighted linear fit based on median recruitment residuals
   geom_line(data=pred,aes(y=fit),lwd=0.2,linetype="dashed") +
   #geom_line(data=pred,aes(y=lwr),lwd=0.2,color="darkgray",alpha=1) +
   #geom_line(data=pred,aes(y=upr),lwd=0.2,color="darkgray",alpha=1) +
   theme_sleek() +
   NULL
ggsave("IPM-sthd-kelt-survival-anomaly.pdf",p5,width=4,height=3)

##===============================================## covariate effects
if(covar_effects) {
   cov_eff_post<-data.frame(extract1(IPM_fit,"beta_SS")) ## posterior
   apply(cov_eff_post,2,median)
   names(cov_eff_post)<-c("SST","Pinks")
   #cov1<-c("Summer SST (°C)")
   #cov2<-c("Pink salmon")
   ##------------------------------------------------## effects plots
   p6 <- cov_eff_post %>% 
      pivot_longer(col=everything(),names_to="name",values_to="value") %>%
      ggplot(aes(x=name,y=value)) +
      geom_violin(lwd=0.1,col="white") +
      stat_summary(fun.data=summary_CI95,size=0.25) +
      stat_summary(fun.data=summary_CI50,size=1.0) +  
      geom_hline(yintercept=0,linetype="dashed",size=0.2) +
      labs(x="",y="Effect size") + 
      theme_sleek() +
      theme(axis.title.x=element_blank()) +
      NULL
   ggsave("IPM-sthd-kelt-survival-cov-effect.pdf",p6,width=3,height=3)
   ##---------------------------------## probability above/below zero
   my_pnorm<-function(x,output){ return(pnorm(0,mean=x[1],sd=x[2])) }
   beta_df<-data.frame(apply(cov_eff_post,2,median),
                       apply(cov_eff_post,2,sd))
   prob_below_zero<-apply(beta_df,1,my_pnorm)
   prob_above_zero <- 1-prob_below_zero
   
}

##-------------------## similar plots for single-population model fits
# select_pops<-c("Quinault","Quillayute","Queets","Humptulips","Chehalis")
# populations<-pops[pops %in% select_pops]
# nP<-length(populations)
# for(i in 1:nP) {
#    ## get fish data for system
#    file_dir<-paste0(home,"/R/data/River_Files")
#    age_dir<-paste0(home,"/R/data/age_data_for_IPM")
#    fish_dat<-data_prep(system=populations[i],covars=F,file_dir=file_dir,age_dir=age_dir,type=cdata_type) %>% filter(year>=1980)
#    ## get IPM fit for system
#    system_dir<-paste0(home,"/R/output/",populations[i])
#    # setwd(file.path(system_dir))
#    fit<-readRDS(paste0(home,"/R/output/",populations[i],"/IPM_fit.Rdata"))
#    ## plot kelt survival rate
#    plot_kelt_survival_rate(fit,fish_dat,trend=F,cov=F,pdf=T)
# }

##=================================================================##
##=============================================## recruitment process
##=================================================================##
## auto-correlated residuals represent environmental stochasticity in survival occurring after density dependence, most likely during marine residence (in the multi-population context the auto-correlated errors represent shared environmental fluctuations occurring after population-level density dependence)
##==========================================## process error hyper-SD
## 'sigma_year_R' = hyper-SD of brood year log productivity anomalies
sigma_R_post<-extract1(IPM_fit,"sigma_year_R")
sigma_R_qs<-quantile(sigma_R_post,prob=probs)
sigma_R_est<-data.frame(value=sigma_R_post) %>% add_column(name="sigma_R")

##=================================================## AR1 coefficient
## 'rho_R' = AR(1) coef for log productivity anomalies
rho_R_post<-extract1(IPM_fit,"rho_R")
rho_R_qs<-quantile(rho_R_post,prob=probs)
rho_R_est<-data.frame(value=rho_R_post) %>% add_column(name="rho_R")

##-----------------------------------------------------------## plots
p7a<-data.frame(rbind(sigma_R_est,rho_R_est)) %>%
   ggplot(aes(x=name,y=value)) +
   # geom_violin(lwd=0.1,col="darkgray") +
   stat_summary(fun.data=summary_CI95,size=0.25) +
   stat_summary(fun.data=summary_CI50,size=0.5) +  
   scale_y_continuous(limits=c(-0.2,0.6)) +
   geom_hline(yintercept=0,linetype="dashed",size=0.2) +
   labs(x="",y="Effect size") + 
   theme_sleek() +
   theme(axis.title.x=element_blank()) +
   NULL
ggsave("IPM-sthd-rho-AR1-and-sigma-hyperSD.pdf",p7a,width=3,height=3)

##===========================================## recruitment anomalies
## 'eta_year_R' = log brood year productivity anomalies
eta_R_post<-extract1(IPM_fit,"eta_year_R")
eta_R_qsr<-t(apply(eta_R_post,2,function(x) quantile(x,prob=probs)))

eta_R_qs<-data.frame(eta_R_qsr) %>% add_column(year=dat_years)

## estimate trend over time
etas_R<-data.frame(median=apply(eta_R_post,2,median),sd=apply(eta_R_post,2,sd)) %>% add_column(year=dat_years)

newdata<-data.frame(year=dat_years)
lm_fit<-lm(etas_R$median~etas_R$year,weights=1/(etas_R$sd^2))

pred<-predict(lm_fit,newdata,interval="confidence",level=0.95)
pred<-data.frame(pred) %>% add_column(year=etas_R$year) 

p7 <- eta_R_qs %>%
   ## plot median recruitment residuals and 95% CIs
   ggplot(aes(x=year,y=X50.)) +
   geom_line(lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.1) +
   labs(x="Brood year",y="Recruitment anomaly") + 
   ## add weighted linear fit based on median recruitment residuals
   geom_line(data=pred,aes(y=fit),lwd=0.2,linetype="dashed") +
   #geom_line(data=pred,aes(y=lwr),lwd=0.2,color="darkgray",alpha=1) +
   #geom_line(data=pred,aes(y=upr),lwd=0.2,color="darkgray",alpha=1) +
   theme_sleek() +
   NULL
ggsave("IPM-sthd-recruitment-anomaly.pdf",p7,width=4,height=3)

write.csv(etas_R,"IPM_etas_R.csv",row.names=F)

##---------------------------## % deviation from expected recruitment
p_dev_exp_rec<-data.frame(100*(exp(eta_R_qsr)-1))
p_dev_exp_rec$year=eta_R_qs$year

## estimate trend over time
etas_exp_R<-data.frame(median=p_dev_exp_rec$X50.) %>% 
   add_column(year=dat_years)

newdata<-data.frame(year=dat_years)
lm_fit<-lm(etas_exp_R$median~etas_exp_R$year)

pred<-predict(lm_fit,newdata,interval="confidence",level=0.95)
pred<-data.frame(pred) %>% add_column(year=etas_exp_R$year) 

p7b <- p_dev_exp_rec %>%
   ## plot median recruitment residuals and 95% CIs
   ggplot(aes(x=year,y=X50.)) +
   geom_line(lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.1) +
   labs(x="Brood year",y="Recruitment deviation (%)") + 
   ## add weighted linear fit based on median recruitment residuals
   geom_line(data=pred,aes(y=fit),lwd=0.2,linetype="dashed") +
   #geom_line(data=pred,aes(y=lwr),lwd=0.2,color="darkgray",alpha=1) +
   #geom_line(data=pred,aes(y=upr),lwd=0.2,color="darkgray",alpha=1) +
   theme_sleek() +
   NULL
ggsave("IPM-sthd-recruitment-deviation-percent.pdf",p7b,width=4,height=3)

##=======================## recruitment anomaly vs proportion repeats
propR<-age_est %>%
   group_by(pop,year,MR) %>%
   dplyr::summarize(prop=sum(prop)) %>%
   ungroup() %>%
   filter(!MR=="maiden") %>%
   dplyr::select(-MR) %>%
   rename(prop_repeats=prop) %>%
   data.frame()

etas_vs_propR<-etas_R %>% ## recruitment anomalies by by brood year
   slice(1:(n()-6)) %>% ## dominant ages missing for last six years
   left_join(propR,by="year") %>%
   group_by(pop) %>%
   slice(-1) %>%
   ungroup() %>%
   dplyr::select(year,pop,median,sd,prop_repeats) %>%
   data.frame()

p7c<-etas_vs_propR %>%
   ggplot(aes(x=prop_repeats,y=median,color=pop,fill=pop)) +
   geom_point(size=2) +
   geom_smooth(method="lm",se=FALSE,formula=y~x) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   #scale_x_continuous(limits=c(0,max(etas_vs_propR$prop_repeats))) +
   #scale_y_continuous(limits=c(-0.8,0.8)) +
   labs(x="Proportion repeat spawners",y="Recruitment anomaly") + 
   facet_wrap(vars(pop),ncol=nP,scales="free_x") + 
   theme_sleek() +
   theme(
      legend.position="none",
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15)
   ) +
   NULL
ggsave("IPM-sthd-recruitment-anomaly-vs-proportion-repeats.pdf",p7c, width=2+2*nP,height=3)

##========================## recruitment anomaly vs mean spawners age
meanage<-age_est %>%
   mutate_at('age',as.numeric) %>%
   group_by(pop,year) %>%
   dplyr::summarize(mean_age=mean(rep(age,round(100*prop)))) %>%
   ungroup() %>%
   data.frame()

etas_vs_meanage<-etas_R %>% ## recruitment anomalies by by brood year
   slice(1:(n()-6)) %>% ## dominant ages missing for last six years
   left_join(meanage,by="year") %>%
   group_by(pop) %>%
   slice(-1) %>%
   ungroup() %>%
   dplyr::select(year,pop,median,sd,mean_age) %>%
   data.frame()

p7d<-etas_vs_meanage %>%
   ggplot(aes(x=mean_age,y=median,color=pop,fill=pop)) +
   geom_point(size=2) +
   geom_smooth(method="lm",se=FALSE,formula=y~x) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   scale_y_continuous(limits=c(-0.8,0.8)) +
   labs(x="Mean spawner age",y="Recruitment anomaly") + 
   facet_wrap(vars(pop),ncol=nP,scales="free_x") + 
   theme_sleek() +
   theme(
      legend.position="none",
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15)
   ) +
   NULL
ggsave("IPM-sthd-recruitment-anomaly-vs-proportion-repeats.pdf",p7c, width=2+2*nP,height=3)


##===============================================## covariate effects
## 'beta_R' regression coefficients for log recruitment
if(covar_effects) {
   beta_R_post<-data.frame(extract1(IPM_fit,"beta_R")) ## posterior
   apply(beta_R_post,2,median)
   names(beta_R_post)<-c("NPGO","SST","Pinks")
   ## new names
   # cov1<-c("Pink salmon effect estimate")
   # cov2<-c(expression("Summer SST"[arc]*" (°C) effect estimate"))
   # cov3<-c("North Pacific Gyre Oscillation")
   ##------------------------------------------------## effects plots
   p8 <- beta_R_post %>% 
      pivot_longer(col=everything(),names_to="name",values_to="value") %>%
      ggplot(aes(x=name,y=value)) +
      geom_violin(lwd=0.1,col="white") +
      stat_summary(fun.data=summary_CI95,size=0.25) +
      stat_summary(fun.data=summary_CI50,size=1.0) +  
      geom_hline(yintercept=0,linetype="dashed",size=0.2) +
      labs(x="",y="Effect size") + 
      theme_sleek() +
      theme(axis.title.x=element_blank()) +
      NULL
   ggsave("IPM-sthd-recruit-anomaly-cov-effect.pdf",p8,width=3,height=3)
   ##---------------------------------## probability above/below zero
   my_pnorm<-function(x,output){ return(pnorm(0,mean=x[1],sd=x[2])) }
   beta_df<-data.frame(apply(beta_R_post,2,median),
                       apply(beta_R_post,2,sd))
   prob_below_zero<-apply(beta_df,1,my_pnorm)
   prob_above_zero <- 1-prob_below_zero

}

##=================================================================##
##==========================## time-varying productivity and capacity
##=================================================================##

##====================================================## productivity

##-----------------------------------------------------------## alpha
mu_alpha_tv<-array(NA,dim=c(nY,nS),dimnames=list(dat_years,seq(nS)))
for(y in 1:nY) { mu_alpha_tv[y,]<-as.vector(exp(eta_R_post[,y]+log(mu_alpha_post)))$' mu_alpha ' }

mu_alpha_tv_qs<-apply(mu_alpha_tv,1,function(x)quantile(x,prob=probs))
mu_alpha_tv_qs<-data.frame(t(data.frame(mu_alpha_tv_qs)))
max_y<-max(mu_alpha_tv_qs)
mu_alpha_tv_qs$year=dat_years

##--------------------## plot time-varying productivity by population
p9a <- mu_alpha_tv_qs %>%
   ggplot(aes(x=year,y=X50.)) +
   geom_line(lwd=0.5,alpha=1) +
   geom_hline(yintercept=1,color="gray",size=0.1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.1) +
   # scale_x_continuous(limits=c(1978,2022)) +
   scale_y_continuous(limits=c(0,max_y)) +
   labs(x="Year",y="Population productivity") + 
   theme_sleek() +
   theme(
      legend.position="none",
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15)
   ) +
   NULL
ggsave("IPM-sthd-productivity-time-varying-hyperparameter.pdf",p9a,width=5,height=4)

##-----------------------------------------------## population alphas
alpha_tv<-array(NA,dim=c(nP,nY,nS),dimnames=list(pops,dat_years,seq(nS)))

for(p in 1:nP) {
   for(y in 1:nY) {   
      alpha_tv[p,y,]<-exp(eta_R_post[,y]+log(alphas_post[,p]))
   }
}

alpha_tv_qs<-apply(alpha_tv,c(1,2),function(x) quantile(x,prob=probs))
alpha_tv_qs<-data.frame(t(data.frame(alpha_tv_qs)))
max_y<-max(alpha_tv_qs)
alpha_tv_qs$population<-as.vector(rep(pops,nY))
alpha_tv_qs$year=as.vector(rep(dat_years,each=nP))

##--------------------## plot time-varying productivity by population
p9b <- alpha_tv_qs %>%
   ggplot(aes(x=year,y=X50.,color=population,fill=population)) +
   geom_line(lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.1) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   scale_x_continuous(limits=c(1978,2022)) +
   scale_y_continuous(limits=c(0,max_y)) +
   labs(x="Year",y="Population productivity") + 
   facet_wrap(vars(population),ncol=nP,scales="free_x") + 
   theme_sleek() +
   theme(
      legend.position="none",
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15)
   ) +
   NULL
ggsave("IPM-sthd-productivity-time-varying.pdf",p9b,width=15,height=3)

##========================================## iteroparous productivity
## ITEROPAROUS PRODUCTIVITY CORRECTED WITH KELT SURVIVAL RATE
## alpha_itero = intrinsic productivity of iteroparous population 
## mutate(alpha=alpha*(1+surv_av/(1-surv_av))) %>% 

##-----------------------------------------------------------## alpha
mu_alpha_itero_tv<-array(NA,dim=c(nY,nS),dimnames=list(dat_years,seq(nS)))
for(y in 1:nY) { 
   mu_alpha_itero_tv[y,]<-as.vector(exp(eta_R_post[,y]+log(mu_alpha_post)*(1+mu_surv_post/(1-mu_surv_post))))$' mu_alpha ' 
   }

mu_alpha_itero_tv_qs<-apply(mu_alpha_itero_tv,1,function(x)quantile(x,prob=probs))

mu_alpha_itero_tv_qs<-data.frame(t(data.frame(mu_alpha_itero_tv_qs)))
max_y<-max(mu_alpha_itero_tv_qs)
mu_alpha_itero_tv_qs$year=dat_years

##--------------------## plot time-varying productivity by population
p9c <- mu_alpha_itero_tv_qs %>%
   ggplot(aes(x=year,y=X50.)) +
   geom_line(lwd=0.5,alpha=1) +
   geom_hline(yintercept=1,color="gray",size=0.1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.1) +
   # scale_x_continuous(limits=c(1978,2022)) +
   scale_y_continuous(limits=c(0,max_y)) +
   labs(x="Year",y="Population productivity") + 
   theme_sleek() +
   theme(
      legend.position="none",
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15)
   ) +
   NULL
ggsave("IPM-sthd-productivity-iteroparous-time-varying-hyperparameter.pdf",p9c,width=5,height=4)

## NEED THE SAME FOR POPULATION-SPECIFIC ALPHA ESTIMATES BUT REQUIRES SURVIVIAL RATES BY YEAR AND POPULATION FROM POSTERIOR TO BE FOLDED IN AS 3-D ARRAY

##========================================================## capacity

##-----------------------------------------------------------## Rmax
mu_Rmax_tv<-array(NA,dim=c(nY,nS),dimnames=list(dat_years,seq(nS)))
for(y in 1:nY) { mu_Rmax_tv[y,]<-as.vector(exp(eta_R_post[,y]+log(mu_Rmax_post)))$' mu_Rmax' }

mu_Rmax_tv_qs<-apply(mu_Rmax_tv,1,function(x)quantile(x,prob=probs))
mu_Rmax_tv_qs<-data.frame(t(data.frame(mu_Rmax_tv_qs)))
max_y<-max(mu_Rmax_tv_qs)
mu_Rmax_tv_qs$year=dat_years

##--------------------------------------## plot time-varying capacity
p10a <- mu_Rmax_tv_qs %>%
   ggplot(aes(x=year,y=X50.)) +
   geom_line(lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.1) +
   labs(x="Year",y="Population capacity") + 
   theme_sleek() +
   theme(
      legend.position="none",
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15)
   ) +
   NULL
ggsave("IPM-sthd-capacity-time-varying-hyperparameter.pdf",p10a,width=5,height=4)

##------------------------------------------------## population Rmax
Rmax_tv<-array(NA,dim=c(nP,nY,nS),dimnames=list(pops,dat_years,seq(nS)))
for(p in 1:nP) {
   for(y in 1:nY) {   
      Rmax_tv[p,y,]<-exp(eta_R_post[,y]+log(Rmax_post[,p]))
   }
}

Rmax_tv_qs<-apply(Rmax_tv,c(1,2),function(x) quantile(x,prob=probs))
Rmax_tv_qs<-data.frame(t(data.frame(Rmax_tv_qs)))
max_y<-max(Rmax_tv_qs)
Rmax_tv_qs$population<-as.vector(rep(pops,nY))
Rmax_tv_qs$year=as.vector(rep(dat_years,each=nP))

##------------------------## plot time-varying capacity by population
p10 <- Rmax_tv_qs %>%
   ggplot(aes(x=year,y=X50.,color=population,fill=population)) +
   geom_line(lwd=0.5,alpha=1) +
   geom_ribbon(aes(ymin=X25.,ymax=X75.),color=NA,alpha=0.3) +
   geom_ribbon(aes(ymin=X5.,ymax=X95.),color=NA,alpha=0.1) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   scale_x_continuous(limits=c(1978,2022)) +
   scale_y_continuous(limits=c(0,max_y)) +
   labs(x="Year",y="Population capacity") + 
   facet_wrap(vars(population),ncol=nP,scales="free_x") + 
   theme_sleek() +
   theme(
      legend.position="none",
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      panel.grid.minor=element_blank(),
      panel.grid.major.y=element_blank(),
      strip.background=element_rect(fill=NA),
      strip.text=element_text(margin=margin(b=3,t=3)),
      legend.box.margin=margin(0,-10,0,-15)
   ) +
   NULL
ggsave("IPM-sthd-capacity-time-varying.pdf",p10,width=15,height=3)

##=================================================================##
##==========================================## make PowerPoint slides
##=================================================================##
## using an empty template with the correct widescreen dimensions
pp <- read_pptx("../template.pptx") %>%
   ph_with(value=p01,ph_location(left=1,top=2,width=4,height=3))%>%
   add_slide()%>%
   ph_with(value=p02,ph_location(left=1,top=2,width=1.5,height=3))%>%
   add_slide()%>%
   ph_with(value=p04,ph_location(left=1,top=2,width=4,height=3))%>%
   add_slide()%>%
   ph_with(value=p05,ph_location(left=1,top=2,width=1.5,height=3))%>%
   add_slide()%>%
   ph_with(value=p1,ph_location(left=1,top=2,width=9,height=2))%>%
   #add_slide()%>%
   ph_with(value=p2,ph_location(left=1,top=4,width=9,height=2))%>%
   #add_slide()%>%
   #ph_with(value=p3,ph_location(left=1,top=5,width=12,height=2.2))%>%
   add_slide()%>%
   ph_with(value=p4,ph_location(left=1,top=2,width=9,height=2.4))%>%
   # add_slide()%>%
   ph_with(value=p5,ph_location(left=1,top=4,width=2,height=2))%>%
   #add_slide()%>%
   #ph_with(value=p6,ph_location(left=1,top=2,width=2.5,height=3))%>%
   add_slide()%>%
   ph_with(value=p7,ph_location(left=1,top=2,width=4,height=3))%>%
   #add_slide()%>%
   #ph_with(value=p8,ph_location(left=1,top=2,width=3,height=3))%>%
   print(target="IPM_sthd_model_plots.pptx")

# pp <- read_pptx("../template.pptx") %>%
#    ph_with(value=p3,ph_location(left=1,top=5,width=12,height=2.2))%>%
#    print(target="IPM_sthd_model_SR_plots.pptx")

##=================================================================##
##=================================================================##
##===============================================## post-hoc analyses
##=================================================================##
##=================================================================##
## continue with post-hoc analyses if covariates not included in IPM
if(covar_effects) {
   stop("covariates already included in the model") 
}else{ 
   
##=======================================================## comments
## 1) Multiple linear regression model of recruitment anomaly 
## 2) Multiple linear regression model of kelt survival rate anomaly
## >> relevant when no covariate effects were included in the model
## 3) Dynamic Factor Analysis of kelt survival rates by population

##=================================================================##
##===========================## linear model of recruitment anomalies
##=================================================================##
## lags: dominant FW age 2 so test ocean conditions in years 2,3,4 ##

##=======================================## proportion repeat spawners
propR<-age_est %>%
   group_by(year,MR) %>% ## averaged across populations
   dplyr::summarize(prop=sum(prop)) %>%
   ungroup() %>%
   filter(!MR=="maiden") %>%
   dplyr::select(-MR) %>%
   rename(prop_repeats=prop) %>%
   slice(-1) %>%
   data.frame()   
   
##==============================## merge residuals and covariate data
df<-data.frame(year=seq(min(all_years-4),max(all_years),1)) %>%
   left_join(etas_R %>% dplyr::select(year,median,sd)) %>%
   rename(residuals=median) %>%
   left_join(propR,by="year") %>%
   left_join(pinks,by='year') %>%
   left_join(npgo,by='year') %>%
   left_join(sst,by='year') %>%
   left_join(sst_cst,by='year') #%>%
   #na.omit()
##---------------------------------------------------## add seal data
# df<-df %>%
#    left_join(seals %>% dplyr::select(year,seals_2),by='year') %>% 
#    # dplyr::select(year,residuals,seals_2) %>%
#    na.omit()
##------------------------------------------------## scale covariates
df1<-df %>% dplyr::select(year,residuals,sd) 
df2<-df %>% dplyr::select(-year,-residuals,-sd) 
df_means<-sapply(df2,function(x) mean(x,na.rm=T))
df_sds<-sapply(df2,function(x) sd(x,na.rm=T))
df<-data.frame(cbind(df1,scale(df2))) %>% na.omit()
##----------------------------------------------------## correlations
out<-data.frame(cor(as.matrix(df),use="pairwise.complete.obs"))
out[,1:2]
##--------------------------------------------------## even-odd dummy
# is.even<-function(x){ x%%2 ==0 }
# df$evenodd<-factor(is.even(df$year),labels=c("odd","even"))

##=================================================## model selection
options(na.action="na.fail")
##----------------------------------------## highest correlation lags
form<-formula(residuals~pinks_4+NPGO_2+SST_3+SST_cst_1) 
##--------------------------------------------------------## all lags
## models w/o NPGO not as good (but depends on spatial scale)
form<-formula(residuals~pinks_2+pinks_3+pinks_4+NPGO_2+NPGO_3+NPGO_4+SST_2+SST_3+SST_4+SST_cst_2+SST_cst_3+SST_cst_4) 
##---------------------------------------------------## model fitting
# mod_lm_full<-lm(form,data=df,weights=1/(sd^2))
mod_lm_full<-lm(form,data=df)
mod_select<-dredge(mod_lm_full,trace=F,rank="AICc")
aic_table<-data.frame(mod_select)[1:10,-1];aic_table
mod_delta2<-aic_table[aic_table$delta<2,] ## models with delta_AIC<2
index<-which(mod_delta2$df==min(mod_delta2$df)) 
mod_sel<-get.models(mod_select,subset=index)[[1]] ## use simple model
# mod_sel<-get.models(mod_select,subset=1)[[1]] ## use lowest AICc
summary(mod_sel)
ncovar<-dim(summary(mod_sel)$coefficients)[1]-1
par(mfcol=c(1,ncovar),mar=c(4,4,1,1));visreg(mod_sel)

##==================================================## selected model
mod<-mod_sel
##---------------------------------------------------## model results
residuals<-residuals(mod,type="response")
fitted<-fitted(mod)
out_mod<-summary(mod) 
as.numeric(pacf(residuals(mod),lag=9,plot=F)$acf) ## autocorrelation
car::vif(mod) ## VIF
##----------------------------------## pairwise covariate correlations
mod_terms<-as.character(names(mod[[1]])[-1])
nterms<-length(mod_terms)
test_data<-df[,colnames(df) %in% mod_terms]
test<-apply(test_data,2,function(x) as.numeric(as.character(x)))
round(rcorr(test,type="pearson")$r,2) ## SST_acr~pinks corr: 0.14

##=============================================## variable importance
relimp<-calc.relimp(mod,type="lmg") ## print(relimp)
v1<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"pinks")]))
v2<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"SST")]))
v3<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"NPGO")]))
relimp_df<-data.frame(pinks=v1,sst=v2,npgo=v3)
round(relimp_df*100,2) ## % of response variance explained by variable
round(sum(relimp_df),2) ## total % of response variance explained

##=================================================## partial effects
# pp1<-plot_instr( code={ ## instead of next line for ppt slides
pdf("IPM-sthd-recruitment-anamaly-covariates-effects.pdf",width=4*nterms,height=4)
layout(matrix(c(1:length(mod_terms)),nrow=1,byrow=T))
par(mar=c(4,4,1,1),mgp=c(2.5,0.5,0),tcl=-0.3,cex.lab=1.2)
xlab<-c("NPGO","Pink salmon abundance (millions)","Summer SST (°C)")
xtrans<-function(x) { x*df_sds[names(df_sds)==covar]+df_means[names(df_means)==covar] }
for(i in 1:length(mod_terms)) {
   covar<-mod_terms[i]
   visreg(mod,xvar=covar,xlab=xlab[i],partial=T,ylab="Partial effect on recruitment anomaly",scale="response",xtrans=xtrans,points.par=list(cex=1,pch=16,col=1),fill.par=list(col=alpha(1,0.1)),line.par=list(col=1))
}
dev.off()
# }) ## instead of previous line for ppt slides

##===============================================## model predictions
newD<-list(pinks_4=df$pinks_4,SST_4=df$SST_4,NPGO_2=df$NPGO_2)
# newD<-list(pinks_4=df$pinks_4,SST_3=df$SST_3,NPGO_2=df$NPGO_2)
# newD<-list(pinks_4=df$pinks_4,SST_4=df$SST_4,NPGO_2=df$NPGO_2)
# newD<-list(pinks_4=df$pinks_4,SST_coast_2=df$SST_coast_2,seals_2=df$seals_2,NPGO_2=df$NPGO_2)
predicted<-predict(mod,newdata=newD,se.fit=T,interval="confidence",level=0.95,type="response")
pp2<-data.frame(predicted$fit) %>% 
   add_column(year=df$year) %>%
   ggplot(aes(x=year,y=fit)) +
   geom_line(aes(y=fit),lwd=0.5,color="gray") +
   geom_ribbon(aes(ymin=lwr,ymax=upr),color="gray",alpha=0.1) +
   geom_line(aes(y=residuals),data=df,lwd=0.25) +
   geom_point(aes(y=residuals),data=df,size=1) +
   labs(x="Year",y="Recruitment anomaly") + 
   theme_sleek() +
   NULL
ggsave("IPM-sthd-recruitment-anomaly-predicted-and-observed.pdf",pp2,width=5,height=4)

##=================================================================##
##====================## linear model of kelt survival rate anomalies
##=================================================================##

##==================## merge kelt survival anomaly and covariate data
df<-data.frame(year=seq(min(all_years-4),max(all_years),1)) %>%
   left_join(etas_SS %>% dplyr::select(year,median,sd)) %>%
   rename(survival=median) %>%
   left_join(pinks,by='year') %>%
   left_join(npgo,by='year') %>%
   left_join(sst,by='year') %>%
   left_join(sst_cst,by='year') %>%
   dplyr::select(-contains("_1")) %>%
   dplyr::select(-contains("_2")) %>%
   dplyr::select(-contains("_3")) %>%
   dplyr::select(-contains("_4")) 
##---------------------------------------------------## add seal data
# df<-df %>%
#    left_join(seals %>% dplyr::select(year,seals_2),by='year') %>% 
#    # dplyr::select(year,survival,seals_2) %>%
#    na.omit()
##------------------------------------------------## scale covariates
df1<-df %>% dplyr::select(year,survival,sd) 
df2<-df %>% dplyr::select(-year,-survival,-sd) 
df_means<-sapply(df2,function(x) mean(x,na.rm=T))
df_sds<-sapply(df2,function(x) sd(x,na.rm=T))
df<-data.frame(cbind(df1,scale(df2))) %>% na.omit()
##----------------------------------------------------## correlations
out<-data.frame(cor(as.matrix(df),use="pairwise.complete.obs"))
out[,1:2]
##--------------------------------------------------## even-odd dummy
# is.even<-function(x){ x%%2 ==0 }
# df$evenodd<-factor(is.even(df$year),labels=c("odd","even"))

##=================================================## model selection
options(na.action="na.fail")
##---------------------------------------------------------## no lags
form<-formula(survival~pinks+NPGO+SST+SST_cst) 
##---------------------------------------------------## model fitting
# mod_lm_full<-lm(form,data=df,weights=1/(sd^2))
mod_lm_full<-lm(form,data=df)
mod_select<-dredge(mod_lm_full,trace=F,rank="AICc")
aic_table<-data.frame(mod_select)[1:10,-1];aic_table
mod_delta2<-aic_table[aic_table$delta<2,] ## models with delta_AIC<2
index<-which(mod_delta2$df==min(mod_delta2$df))
mod_sel<-get.models(mod_select,subset=index)[[1]] ## use simple model
# mod_sel<-get.models(mod_select,subset=1)[[1]] ## use lowest AICc
summary(mod_sel)
ncovar<-dim(summary(mod_sel)$coefficients)[1]-1
par(mfcol=c(1,ncovar),mar=c(4,4,1,1));visreg(mod_sel)

##==================================================## selected model
mod<-mod_sel
##---------------------------------------------------## model results
resid<-residuals(mod,type="response")
fitted<-fitted(mod)
out_mod<-summary(mod) 
as.numeric(pacf(residuals(mod),lag=9,plot=F)$acf) ## autocorrelation
mod_terms<-as.character(names(mod[[1]])[-1])
nterms<-length(mod_terms)
if(nterms>1) car::vif(mod) ## VIF
##----------------------------------## pairwise covariate correlations
if(nterms>1) {
   test_data<-df[,colnames(df) %in% mod_terms]
   test<-apply(test_data,2,function(x) as.numeric(as.character(x)))
   round(rcorr(test,type="pearson")$r,4)
}

##=============================================## variable importance
relimp<-calc.relimp(mod,type="lmg") ## print(relimp)
v1<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"pinks")]))
v2<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"SST")]))
# v3<-as.numeric(sum(relimp$lmg[str_detect(names(relimp$lmg),"NPGO")]))
relimp_df<-data.frame(pinks=v1,sst=v2)#,npgo=v3)
round(relimp_df*100,2) ## % of response variance explained by variable
round(sum(relimp_df),3) ## total % of response variance explained

##=================================================## partial effects
## pp3<-plot_instr( code={ ## instead of next line for ppt slides
pdf("IPM-sthd-kelt-survival-anomaly-covariate-effects.pdf",width=4*nterms,height=4)
layout(matrix(c(1:length(mod_terms)),nrow=1,byrow=T))
par(mar=c(4,4,1,1),mgp=c(2.5,0.5,0),tcl=-0.3,cex.lab=1.2)
#xlab<-c("NPGO","Pink salmon abundance (millions)","Summer SST (°C)")
xlabels<-data.frame(term=mod_terms) %>%
   mutate(name=case_when(
      grepl("pink",term) ~ "Pink salmon abundance (millions)", 
      grepl("SST",term) ~ "Summer SST (°C)", 
      grepl("NPGO",term) ~ "NPGO", 
   ))
xtrans<-function(x) { x*df_sds[names(df_sds)==covar]+df_means[names(df_means)==covar] }
for(i in 1:length(mod_terms)) {
   covar<-mod_terms[i]
   visreg(mod,xvar=covar,xlab=xlabels$name[i],partial=T,ylab="Partial effect on kelt survival anomaly",scale="response",xtrans=xtrans,points.par=list(cex=1,pch=16,col=1),fill.par=list(col=alpha(1,0.1)),line.par=list(col=1))
} 
dev.off()
## }) ## instead of previous line for ppt slides

##===============================================## model predictions
newD<-list(SST=df$SST,pinks=df$pinks)
predicted<-predict(mod,newdata=newD,se.fit=T,interval="confidence",level=0.95,type="response")
pp4<-data.frame(predicted$fit) %>% 
   add_column(year=df$year) %>%
   ggplot(aes(x=year,y=fit)) +
   geom_line(aes(y=fit),lwd=0.5,color="gray") +
   geom_ribbon(aes(ymin=lwr,ymax=upr),color="gray",alpha=0.1) +
   geom_line(aes(y=survival),data=df,lwd=0.25) +
   geom_point(aes(y=survival),data=df,size=1) +
   labs(x="Year",y="Kelt survival anomaly") + 
   theme_sleek() +
   NULL
ggsave("IPM-sthd-kelt-survival-anomaly-predicted-and-observed.pdf",pp4,width=5,height=4)

##=================================================================##
##==========================================## make PowerPoint slides
##=================================================================##
## using an empty template with the correct widescreen dimensions
pp <- read_pptx("../template.pptx") %>%
   ph_with(value=pp1,ph_location(left=1,top=2,width=10,height=4))%>%
   add_slide()%>%
   ph_with(value=pp2,ph_location(left=1,top=2,width=5,height=4))%>%
   add_slide()%>%
   ph_with(value=pp3,ph_location(left=1,top=2,width=10,height=4))%>%
   add_slide()%>%
   ph_with(value=pp4,ph_location(left=1,top=2,width=5,height=4))%>%
   print(target="IPM_sthd_model_posthoc_analyses.pptx")

##=================================================================##
##======================================## DFA of kelt survival rates
##=================================================================##
## 'bayesdfa' options:
## extremes: estimate_nu=TRUE
## AR(1) process: estimate_trend_ar=TRUE
## moving average process: estimate_trend_ma=TRUE 
## smooth spline model: trend_model="spline",n_knots=xx 
## Gaussian process: trend_model="gp",n_knots=xx 

##============================================## survival time series
## gaussian/lognormal/gamma/binomial/poisson/nbinom2
mvts<-surv_med %>%
   mutate(logit_surv=logit(surv)) %>%
   dplyr::select(-surv) %>%
   pivot_wider(names_from=year,values_from=logit_surv,names_sort=T) %>%
   dplyr::select(-pop) %>%  ## time (years) across columns
   data.frame()
years<-as.numeric(gsub("X","",colnames(mvts)))

##-------------------------------## correlations logit survival rates
res<-cor(as.matrix(t(mvts)),use="pairwise.complete.obs")
min(res) ## all pairwise correlations in logit survival rates > 0.5

##---------------------------------------## plot logit survival rates
p<-surv_med %>% 
   mutate(logit_surv=logit(surv)) %>%
   ggplot(aes(x=year,y=logit_surv,group=pop,color=pop)) +
   geom_line(size=1) +
   scale_color_manual(values=colors) +
   labs(x="Year",y="Median survival rate (logit scale)") + 
   theme_sleek() +
   theme(
      legend.position=c(0.85,0.8),
      legend.title=element_blank()
   ) +
   NULL
ggsave("IPM-sthd-kelt-survival-logit-scale.pdf",p,width=5,height=4)

##===========================================## covariate time series
## lag-1 covariate time series (out-migration year)
cov_ts<- #data.frame(year=c(min(years)-1,years)) %>%
   data.frame(year=years) %>%
   left_join(npgo,by='year')  %>%
   left_join(pinks,by='year') %>% 
   left_join(sst_arc,by='year') %>% 
   left_join(sst_cst,by='year') %>% 
   mutate( ## all lag-1 means previous, i.e. out-migration year
      NPGO=scale(NPGO),
      SST_cst=scale(SST_cst),
      SST_arc=scale(SST_arc),
      pinks=scale(pinks)
      # NPGO=scale(lag(NPGO,1)),
      # SST_cst=scale(lag(SST_cst,1)),
      # SST_arc=scale(lag(SST_arc,1)),
      # pinks=scale(lag(pinks,1))
   ) %>%
   na.omit()
covars<-c("NPGO","pinks","SST_arc","SST_cst")
nC<-length(covars)
##----------------------------------------------------## correlations
data.frame(cor(as.matrix(cov_ts[,-1]),use="pairwise.complete.obs"))
## only include one climate covariate in each model

##================================================## fit Bayesian DFA
# fit_1<-fit_dfa(y=mvts,num_trends=1,scale="zscore",iter=900,chains=4)
# is_converged(fit_1,threshold=1.05)
# rot_1<-rotate_trends(fit_1)
# plot_trends(rot_1) + theme_bw()
# plot_loadings(rot_1) + theme_bw()
# # plot_fitted(fit_1) + theme_bw()
# loo_1<-loo(fit_1)
# loo_1$estimates

##==================================================## fit DFA in TMB
source(paste0(home,"/R/other/DFA_TMB.R")) ## Tim's code
mvts<-t(as.matrix(apply(mvts,1,scale)))
res<-list()
## model without covariates
res[[1]]<-runDFA(mvts,NumStates=1,ErrStruc="DUE",EstCovar=F)
## models with covariates
covs<-list("NPGO","pinks","SST_arc","SST_cst",c("pinks","NPGO"),c("pinks","SST_arc"),c("pinks","SST_cst"))
nC<-length(covs)
for(i in 1:nC) {
   cov_ts_mod<-cov_ts[names(cov_ts) %in% as.vector(covs[[i]])]
   covts<-t(as.matrix(cov_ts_mod))
   dfa<-runDFA(mvts,NumStates=1,ErrStruc="DUE",EstCovar=T,Covars=covts)
   res[[1+i]]<-dfa
}
nM<-length(res)
##-------------------------------------------## model selection table
tbl<-data.frame(num=seq(1:nM),AIC=sapply(res,function(x) x[["AIC"]]))
tbl$deltaAIC<-round(tbl$AIC-min(tbl$AIC),3) ## delta-AIC
wt<-exp(-0.5*tbl$deltaAIC) ## Akaike weights
tbl$wt<-round(wt/sum(wt),3) ## add to table
tbl<-tbl[order(tbl$AIC),] ## sort by AIC
tbl$cov<-NA
for(i in 1:nC) tbl$cov[1+i]<-paste0(as.vector(covs[[i]]),collapse=" ")
##-----------------------------------------------## select best model
sel_dfa<-res[[tbl[1,"num"]]]
##-----------------------------------------------## covariate effects
effects<-data.frame(sel_dfa$Estimates$D)
colnames(effects)<-tbl$cov[tbl[1,"num"]]
effects$pop<-pops
cov_effs<-effects %>%
   pivot_longer(cols=-pop,names_to="covariate",values_to="effect") %>%
   data.frame()
##------------------------------------------## plot covariate effects
# cols<-c("dodgerblue3","chocolate2","forestgreen","goldenrod1")
p<-cov_effs %>% 
   ggplot(aes(x=effect,y=pop,group=covariate,color=covariate)) +
   geom_point(size=2) +
   geom_vline(xintercept=0) +
   # scale_color_manual(values=cols) +
   # scale_fill_manual(values=cols) +
   labs(y="",x="Covariate effect") + 
   theme_sleek() +
   theme(legend.position="none") +
   facet_wrap(vars(covariate),ncol=nC) +
   NULL
ggsave("DFA-covariate-effects.pdf",p,width=5,height=4)
##---------------------------------------------## trends and loadings
trend<-ilogit(as.vector(sel_dfa$Estimates$u))
loadings<-as.vector(sel_dfa$Estimates$Z)
##--------------------------------------------## plot trends/loadings
pdf("DFA-trends-loadings.pdf",width=8,height=4)
par(mar=c(4,4,1,1),mgp=c(2,0.5,0),tcl=-0.3)
layout(matrix(c(1:2),ncol=2,byrow=T),widths=c(0.6,0.4))
plot(years,trend,type='l',lwd=2,xlab="Year",ylab="Trend",cex.lab=1.2, ylim=c(min(0,trend),max(0,trend)),xlim=c(min(years),max(years)))
abline(h=0,lwd=0.5,lty=2,col=1)
par(mar=c(6,4,1,1))
max.load<-max(abs(loadings))
plot(seq(nP),loadings,type='h',lwd=3,xaxt="n",xlab="",ylab="Loading",ylim=c(-max.load*1.2,max.load*1.2),col=colors,cex.lab=1.2) 
axis(side=1,at=seq(nP),labels=pops,las=3,cex.axis=1.2)
abline(h=0,lwd=1,lty=1,col=1)
dev.off()

} ## end if statement on covariate effects
##=================================================================##
##=================================================================##
##=================================================================##