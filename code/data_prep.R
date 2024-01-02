#' Function to prepare fish data as input to IPM
#' 
#' @author Jan Ohlberger
#' @param system river system when using a single population
#' @param file_dir directory for escapement and catch data
#' @param age_dir directory for age data 
#' @import dataRetrieval dplyr tidyverse readr tibble
#' @return fish data for use in IPM
#'
data_prep<-function(system,file_dir,age_dir){
   ##================================================## load fish data
   syst<-substr(system,1,4)
   LCR_systems<-c("Mill Creek","Elochaman","Grays River")
   Will_systems<-c("Willapa","North","Palix","Nemah","Naselle","Bear")
   ##---------------------------------------------------## escapements
   if(system %in% Will_systems) {
      esc_file<-paste0(file_dir,"/Will_Esc_by_System.csv")
      esc_data<-read_csv(esc_file,col_types="n") %>%
         dplyr::select(year,escapement=all_of(system))
   } else{
      esc_file<-paste0(file_dir,"/",syst,"_Esc.csv")
      esc_data<-read_csv(esc_file,col_types="n")  
   }
   ##---------------------------------------## combine with catch data
   ## works with catch as total count or as rate (e.g. 'impact_rate')
   ## or split into 'impact_rate_tribs' and 'impact_rate_main'
   if(system %in% LCR_systems) {
      catch_file<-paste0(file_dir,"/LCR_Impact_Rates.csv")
      catch_data<-read_csv(catch_file,col_types="n")
      rate_names<-names(catch_data)[-1]
      nrates<-length(rate_names)
      data<-merge(esc_data,catch_data,by="year",all=TRUE)
      if(nrates==1) {
         data<-data %>% 
            rowwise() %>% 
            mutate(ocean_rec=round(escapement/(1-impact_rate))) %>%
            mutate(catch=round(ocean_rec-escapement))
         }
      if(nrates==2) { 
         data<-data %>% 
            rowwise() %>% 
            mutate(river_rec=round(escapement/(1-impact_rate_tribs))) %>% 
            mutate(ocean_rec=round(river_rec/(1-impact_rate_main))) %>%
            mutate(catch=round(ocean_rec-escapement)) 
         }
      } else {
         if(system %in% Will_systems) {
            catch_file<-paste0(file_dir,"/Will_Catch_by_System.csv")
            catch_data<-read_csv(catch_file,col_types="n") %>%
               dplyr::select(year,catch=all_of(system))
            data<-merge(esc_data,catch_data,by="year",all=TRUE)
         } else {
         catch_file<-paste0(file_dir,"/",syst,"_Catch.csv")
         catch_data<-read_csv(catch_file,col_types="n") %>%
            filter(if_any(!year,~!is.na(.))) %>%
            rowwise()%>% 
            mutate(catch=sum(c_across(-year),na.rm=T))
         data<-merge(esc_data,catch_data,by="year",all=TRUE)
         }
      }
   data<-data %>% 
      data.frame() %>%
      replace_na(list(catch=0)) %>%
      rowwise() %>%
      mutate(runsize=catch+escapement) %>%
      dplyr::select(year,escapement,catch,runsize)
   ##------------------------------------------------------## age data
   if(system %in% LCR_systems) {
      age_dat1<-read_csv(paste0(age_dir,"/Aber_Age.csv"),col_types="n")
      age_dat2<-read_csv(paste0(age_dir,"/Lower_CR.csv"),col_types="n")
      age_data<-bind_rows(age_dat2,age_dat1) 
      if("type" %in% names(age_data)) {
         age_data<-age_data %>%
            filter(!type=="drop") %>% ## drop age data when indicated
            dplyr::select(-type) ## remove 'type'
         }
      age_data<-age_data %>% 
         replace(is.na(.),0) %>%
         group_by(year) %>% 
         summarise_all(sum) %>%
         dplyr::select(-'...1')
   } else {
      if(system %in% Will_systems) {
         age_data<-data.frame(cbind(year=NA))
      } else {
      age_file<-paste0(age_dir,"/",syst,"_Age.csv")
      age_data<-read_csv(age_file,col_types="n")
      if("type"%in%names(age_data)){
         age_data<-age_data %>%
            filter(!type=="drop") %>% ## drop age data when indicated
            dplyr::select(-type) ## remove 'type'
      }
      }
   }
   ##----------------------------------------------## combine all data
   dat<-merge(data,age_data,by="year",all=TRUE)
   dat<-dat %>%
      mutate(across(everything(),round)) %>%
      mutate(F_rate=round(catch/runsize,4))%>%
      mutate(B_take_obs=catch) %>%
      mutate(Year=year) %>%
      add_column(pop=system)
   ##==================================## spawner age counts and total
   dat$S_obs <- dat$escapement
   dat <- dat %>% 
      rename_with(.cols = contains("age"), .fn = ~ paste0("n_", ., "_obs"))
   ##==================================## rename columns for new model
   dat<-dat %>%
      rename_with(.cols=starts_with("n_age"),.fn=~gsub("_r_","K_",.))%>%
      rename_with(.cols=starts_with("n_age"),.fn=~gsub("_m_","M_",.))
   ##============================================## other data columns
   dat$A <- 1 ## spawning habitat size
   dat$n_W_obs <- 0 ## frequency of natural-origin spawners
   dat$n_H_obs <- 0 ## frequency of hatchery-origin spawners
   dat$fit_p_HOS <- FALSE ## estimate p_HOS (TRUE/FALSE)
   ##======================================## fish data for stan model
   fish_dat<-dat %>% as.data.frame() %>% 
      dplyr::select(pop,A,year,S_obs,contains("_obs"),fit_p_HOS,F_rate)
   ##===================================================## return data
   return(fish_dat)
}
