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

##=======================================================## functions
home<-here::here()
fxn<-list.files(paste0(home,"/functions"))
invisible(sapply(FUN=source,paste0(home,"/functions/",fxn)))

##================================================## output directory
our_dir<-paste0(home,"/output/")
if(!file.exists(our_dir)) dir.create(file.path(our_dir))
setwd(file.path(our_dir))

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
file_dir<-age_dir<-paste0(home,"/data/RiverFiles")
##--------------------------------------------------------------## OP
if(spatial_extent=="OP"){
pops<-c("Hoh","Queets","Quillayute","Quinault") ## OP (alphabetical)
# skm<-c(10431,12144,19571,14723) ## stream km (OP sthd petition Tab4)
# areas<-c(299,445,629,325) ## square miles watershed
areas<-c(4.38,6.95,9.47,6.61) ## km^2 habitat (Ohlberger et al. 2017)
areaD<-data.frame(cbind(pop=pops,area=areas))
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

write.csv(areaD,"IPM_habitat_areas.csv",row.names=F)

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

##==================================================## save fish data
write.csv(fish_dat,"IPM_fish_dat_all.csv",row.names=F)

##=================================================================##
##==================================================## covariate data
##=================================================================##

##================================================## pink salmon data
file_path<-paste0(home,"/data/pink_salmon_total_abundance.xlsx")

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
sst_cst<-read.csv(paste0(home,"/data/SST_coast_JJA.csv"))

sst_cst<-sst_cst %>%
   mutate(
      SST_cst_1=lead(SST_cst,1),
      SST_cst_2=lead(SST_cst,2),
      SST_cst_3=lead(SST_cst,3),
      SST_cst_4=lead(SST_cst,4),
   ) %>%
   mutate(across(-1,round,2))

##==================================## get ERSST data for summer SST
sst_dat<-read.csv(paste0(home,"/data/ersst_sthd.csv"))

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
file_path<-paste0(home,"/data/harbor_seal_abundances_updated.xlsx")
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

##----------------------------------------------------## combine data
covar_dat<-pinks %>%
   left_join(npgo,by='year') %>%
   left_join(sst,by='year') %>%
   left_join(sst_cst,by='year')

##=================================================================##
write.csv(covar_dat,"IPM_covar_dat_all.csv",row.names=F)

##----------------------------------------------------## correlations
df_dat<-covar_dat %>%
   dplyr::select(-c(matches("[[:digit:]]"))) %>%
   na.omit() ## no NAs in covariates

res<-data.frame(cor(as.matrix(df_dat),use="pairwise.complete.obs"))
res[-1,-1]

##=================================================================##
##=================================================================##
##=====================================## Integrated Population Model
##=================================================================##
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