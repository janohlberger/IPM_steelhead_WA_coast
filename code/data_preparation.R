##=================================================================##
##                                                                 ##
## Load coastal steelhead population data and fit hierarchical IPM ##
##                                                                 ##
##=================================================================##

##========================================================## packages
pkg<-c("here","tidyverse","readr","readxl","dataRetrieval","lubridate")
if(length(setdiff(pkg,rownames(installed.packages())))>0){install.packages(setdiff(pkg,rownames(installed.packages())),dependencies=T)}
invisible(lapply(pkg,library,character.only=T))
home<-here::here()
source(paste0(home,"/code/data_prep.R"))

##================================================## output directory
out_dir<-paste0(home,"/output/")
if(!file.exists(out_dir)) dir.create(file.path(out_dir))
setwd(file.path(out_dir))

##=================================================================##
##=======================================================## fish data
##=================================================================##
file_dir<-age_dir<-paste0(home,"/data/fish_data")
pops<-c("Chehalis","Hoh","Humptulips","Queets","Quillayute","Quinault")
areas<-c(1965,242,224,332,658,253)
write.csv(areas,"IPM_habitat_areas.csv",row.names=F)
##-----------------------------------------------## prepare fish data
pops<-sort(pops)
nP<-length(pops)
fish_dat<-NULL
for(i in 1:nP) {
   pop_dat<-data_prep(
      system=pops[i],
      file_dir=file_dir,
      age_dir=age_dir)
   fish_dat<-bind_rows(fish_dat,pop_dat)
}

##-----------------------------## only years with observed escapement
fish_dat<-fish_dat %>% filter(S_obs>0) 
all_years<-sort(unique(fish_dat$year))
##------------------------------------------## add area by population
fish_dat<-fish_dat %>% 
   full_join(data.frame(cbind(pop=pops,area=areas))) %>%
   mutate_at('area',as.numeric) %>%
   dplyr::select(-A) %>%
   rename(A=area)

##=====================## limit number of initial years without catch
##-------------------------## get first year catch data by population
frst_yrs<-fish_dat %>%
   filter(!is.na(F_rate)) %>% 
   group_by(pop) %>% 
   slice_min(year) %>% 
   dplyr::select(pop,min_year=year)
##---------------## drop data prior to first catch year minus min age
## so NAs in F_rate during early years can be set to zero because the post-removal spawner abundance in years 1:min_age is drawn from prior
fish_dat<-fish_dat %>% 
   full_join(frst_yrs) %>%
   filter(year>min_year-3) %>% ## min_age=2
   dplyr::select(-min_year)

##----------------------------## set NAs in F_rate/B_take_obs to zero
fish_dat<-fish_dat %>% 
   mutate(B_take_obs=0) %>%
   mutate(F_rate=replace_na(F_rate,0))
   
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
index<-which(names(fish_dat) %in% c(m_names,r_names)) 
fish_dat<-fish_dat %>% mutate_at(index, ~replace_na(.,0))
##---------------------------------------------## columns in fish_dat
fish_dat<-fish_dat %>% dplyr::select(pop,A,year,S_obs,all_of(m_names),all_of(r_names),n_W_obs,n_H_obs,B_take_obs,fit_p_HOS,F_rate)

##==================================================## save fish data
write.csv(fish_dat,"IPM_fish_dat_all.csv",row.names=F)

##=================================================================##
##==================================================## covariate data
##=================================================================##

##==========================================## pink salmon abundances
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

##============================================================## NPGO
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

##======================================## ERSST data for coastal SST
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
## Only for test runs - much fewer years of data
file_path<-paste0(home,"/data/harbor_seal_abundances.xlsx")
seal_data<-data.frame(read_excel(file_path)) %>%
   rename(year=Year) %>%
   rename(seals=Abundance) %>%
   dplyr::select(year,seals)

seals<-data.frame(year=all_years) %>%
   full_join(seal_data) %>%
   mutate(
      seals_2=lead(seals,2),
      seals_3=lead(seals,3),
      seals_4=lead(seals,4),
   ) %>%
   mutate(across(-1,round,2))

##=====================================================## river flows
flow_data<-read.csv(paste0(home,"/data/flow_data.csv"))

flows<-flow_data %>%
   mutate(av_CFS_min_1=lead(av_CFS_min,1),
          av_CFS_max_1=lead(av_CFS_max,1)
          ) %>%
   mutate(across(-1,round,4))

##===========================================## combine data and save
covar_dat<-pinks %>%
   full_join(npgo,by='year') %>%
   full_join(sst,by='year') %>%
   full_join(sst_cst,by='year') %>%
   full_join(flows,by='year') %>%
   filter(year>min(fish_dat$year))

write.csv(covar_dat,"data/IPM_covar_dat_all.csv",row.names=F)

##===================================================## include flows
## can't be included in post-hoc test of shared residuals
## can be included to model recruitment in the IPM though
##----------------------------------------------------## combine data
# covar_dat<-flows %>%
#    full_join(pinks,by='year') %>%
#    full_join(npgo,by='year') %>%
#    full_join(sst,by='year') %>%
#    full_join(sst_cst,by='year') %>%
#    filter(year>min(fish_dat$year))
# 
# write.csv(covar_dat,"IPM_covar_dat_with_flows.csv",row.names=F)

##=================================================================##
##=================================================================##
##=================================================================##