##=================================================================##
##                                                                 ##
## Create a map of WA coast with steelhead rivers and river names  ##
##                                                                 ##
##=================================================================##
pkg<-c("here","dplyr","tidyr","ggplot2","maps","mapdata","maptools", "rgdal","raster","rgeos","tigris","RColorBrewer","ggplotify","ggsidekick","patchwork","cowplot","grid","marmap")
if(length(setdiff(pkg,rownames(installed.packages())))>0){install.packages(setdiff(pkg,rownames(installed.packages())),dependencies=T)}
invisible(lapply(pkg,library,character.only=T))

##========================================================## settings
setwd(file.path(paste0(here(),"/R/data/map")))

##=================================================================##
##=======================================================## fish data
##=================================================================##
data<-read.csv("../../output/IPM_fish_dat_without_covars.csv")
pops<-unique(data$pop)
nP<-length(pops)

##----------------------------------------------------------## colors
cols<-c("#A86260","#A5B1C4","#3B4D6B","#89A18D","#747260","#B3872D","#774C2C") 
colors<-colorRampPalette(rev(cols))(nP)
colors<-rep("deepskyblue4",nP)

##==============================================## spawner abundances
gg_esc<-data %>% 
   dplyr::select(pop,year,S_obs) %>%
   ggplot(aes(x=year,y=S_obs,color=pop)) +
   geom_line(aes(y=S_obs),lwd=0.5,alpha=1) +
   geom_point(aes(y=S_obs),size=1.2,alpha=1) +
   scale_color_manual(values=colors) +
   # scale_fill_manual(values=colors) +
   labs(x="Year",y="Spawner abundance") + 
   theme_sleek() +
   scale_y_continuous(limits=c(0,NA)) +
   theme(
      legend.position="none",
      plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), ## t-r-b-l
      strip.text=element_text(size=15),
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      axis.title=element_text(size=15),
      axis.text=element_text(size=12)
   ) +
   facet_wrap(vars(pop),ncol=nP) + #,scales="free_y") +
   NULL

##==================================================## age proportions
ages<-colnames(data)[grepl("age",colnames(data))]

##---------------------## only plot years with at least XX age samples
min_age_samples<-10
age_dat<-data
age_dat[which(rowSums(age_dat[,ages])<min_age_samples),ages]<-NA

##--------------------------## age proportions for maiden and repeats
age_prop <- data.frame(cbind(age_dat %>% dplyr::select(pop,year),prop.table(age_dat[,names(age_dat) %in% ages] %>% replace(is.na(.),0) %>% as.matrix(),1)))

##-------------------## only keep age groups with prop>0.001 for plot
av_prop<-round(colMeans(age_prop[,-c(1,2)],na.rm=T),4)
drop<-names(av_prop)[av_prop<0.001] ## at least 0.1% of samples
age_prop<-age_prop %>% dplyr::select(-all_of(drop))

##------------------## drop years with certain age props equal to one
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
gg_age<-age_data %>%
   ggplot(aes(x=year,y=prop,color=age)) +
   geom_line(lwd=0.25,alpha=1) +
   geom_point(pch=16,size=1.5,alpha=1) +
   scale_color_manual(values=colors) +
   scale_fill_manual(values=colors) +
   labs(x="Year",y="Age proportions")+ 
   theme_sleek() +
   theme(
      plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), ## t-r-b-l
      strip.text=element_text(size=15),
      axis.title=element_text(size=15),
      axis.text=element_text(size=12),
      legend.text=element_text(size=12),
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)
   ) +
   facet_grid(MR~pop,scales="free_y") +
   NULL

##=================================================================##
##========================================================## map data
##=================================================================##

##========================================================## map data
usa<-getData('GADM',country="USA",level=1) ## states
wa<-usa[usa@data$NAME_1 %in% c("Washington"),] ## select WA
tol<-0.0001 ## low (0.1) to high (0.0001) resolution 
wa<-gSimplify(wa,tol=tol,topologyPreserve=F) ## resolution
##---------------------------------------------------## altitude data
alt.usa<-getData('alt',country='USA')
alt<-alt.usa[[1]] ## element 1 is continental US

##================================================## water/lake areas
counties<-c("Grays Harbor","Clallam","Jefferson","Pacific","Kitsap","Mason","Thurston","Pierce","King","Snohomish","Skagit")
areas<-area_water("WA",counties)
##---------------------------------------------## only specific areas
keep_areas<-c("Quinault Lk","Ozette Lk")
areas<-areas[areas$FULLNAME %in% keep_areas,]
ozette<-areas[areas$FULLNAME=="Ozette Lk",]
quinault<-areas[areas$FULLNAME=="Quinault Lk",]

##================================================## ocean bathymetry
west<--124.9;east<--121.9;south<-46.36;north<-48.44
bathys<-getNOAA.bathy(west,east,south,north,resolution=0.01)

##======================================================## river data

##---------------------------------------------## old river shapefile
rvrs<-readOGR("rivers_arc.shp")
rvrs<-spTransform(rvrs,CRS("+proj=longlat"))
rvrs<-rvrs[rvrs@data$WRIA_ID %in% c(1:29),] ## all west of Cascades

##----------------------------------------## WDFW Region 6 river data
R6_rvrs<-readOGR("Region6Streams.shp")
R6_rvrs<-spTransform(R6_rvrs,CRS("+proj=longlat"))
R6_rvrs<-R6_rvrs[R6_rvrs$SPECIES=="Steelhead Trout",] ## sthd rivers

##===================================================## select rivers

##-----------------------------------------------------## river names
use_rvrs<-c("Quillayute River","Calawah River","Sol Duc River","Dickey River","Bogachiel River","Hoh River","Queets River","Clearwater River","Salmon River","Sams River","Quinault River","Humptulips River","Chehalis River","Wynoochee River","Wishkah River","Satsop River","Black River","Skookumchuck River","Newaukum River")
## "Bear","Naselle","Nemah","North","Palix","Willapa"
##---------------------------------------------------## old shapefile
rivers<-rvrs
index<-NULL
for(i in 1:length(use_rvrs)){
   ind<-which(rivers$WTRCRS_NM %in% rivers$WTRCRS_NM[grepl(use_rvrs[i],rivers$WTRCRS_NM)])
   index<-c(index,ind)
}
rivers<-rivers[index,]
##----------------------------------------------## WDFW Region 6 data
sthd_rvrs<-R6_rvrs
index<-NULL
for(i in 1:length(use_rvrs)){
   ind<-which(sthd_rvrs$LLID_STRM_ %in% sthd_rvrs$LLID_STRM_[grepl(use_rvrs[i],sthd_rvrs$LLID_STRM_)])
   index<-c(index,ind)
}
use_rivers<-sthd_rvrs[index,]

##=================================================## watershed areas
hucs<-readOGR("HUC10.shp")
hucs<-spTransform(hucs,CRS("+proj=longlat"))
##-----------------------------------------------## select watersheds
huc_names<-c("Quillayute","Calawah","Sol Duc","Dickey","Bogachiel","Hoh","Queets","Clearwater","Salmon","Sams","Quinault","Humptulips","Chehalis","Wynoochee","Wishkah","Satsop","Black","Skookumchuck","Newaukum")
use_hucs<-hucs
index<-NULL
for(i in 1:length(huc_names)){
   ind<-which(hucs$Name %in% hucs$Name[grepl(huc_names[i],hucs$Name)])
   index<-c(index,ind)
}
use_hucs<-hucs[index,]
use_hucs<-use_hucs[!duplicated(use_hucs$Name),]

##========================================================## plot map
color_map<-TRUE
# pdf("map_with_rivers_and_time_series.pdf",width=6,height=6)
dev.new(width=6,height=6)
par(mar=c(4,4,3,0.5),mgp=c(1.5,0.2,0),tcl=-0.3,cex.lab=.8,cex.axis=.6)
##-------------------------------------------------## map coordinates
west<--124.9;east<--121.9;south<-46.36;north<-48.44
xlim<-c(west,east);ylim<-c(south,north)
center<-c(mean(xlim),mean(ylim))
##--------------------------------------------## elevation base map
if(color_map){ col_elev<-"gray99" }else{ col_elev<-"gray40" }
col_alt<-colorRampPalette(c("white",col_elev))(100) 
plot(alt,col=col_alt,alpha=1,legend=F,add=F,xlim=xlim,ylim=ylim,colNA=NA,ylab="",xlab="",lwd=NA,axes=F)
##-----------------------------------------## add ocean and land fill
if(color_map){
# land<-colorRampPalette(c("white","gray40"))
land<-colorRampPalette(c("white",brewer.pal(9,"YlOrBr")))
water<-colorRampPalette(c(rep("deepskyblue4",2),"deepskyblue3"))
plot.bathy(bathys,add=T,image=TRUE,land=TRUE,n=0,xlim=xlim,ylim=ylim, ylab="",xlab="",axes=F,deep=c(-2000,0),shallow=c(0,0),step=c(2000,0), lwd=c(1,1),lty=c(1,1),col=rep("deepskyblue3",2),drawlabels=c(F,F), bpal=list(c(0,max(bathys),land(100)),c(min(bathys),0,water(100))))
}
##--------------------------------------------------------## rivers
if(color_map){col_riv<-"deepskyblue3"}else{col_riv<-"gray40"}
if(color_map){col_all<-"slategray"}else{col_all<-"gray60"}
# lines(rvrs,col=col_all,lwd=0.1) ## all rivers
lines(use_rivers,col=col_riv,lwd=1) ## WDFW R6 data
lines(rivers,col=col_riv,lwd=1) ## old file (inlc. upper Chehalis)
##-------------------------------## add missing segments for Chehalis
segments(x0=-123.02,y0=46.73,x1=-123.022,y1=46.74,lwd=1,col=col_riv)
segments(x0=-123.02,y0=46.742,x1=-123.005,y1=46.746,lwd=1,col=col_riv)
segments(x0=-123.715,y0=47.26,x1=-123.715,y1=47.265,lwd=1,col=col_riv)
##---------------------------------------------------------## areas
if(color_map){col_lakes<-"deepskyblue3"}else{col_lakes<-"gray25"}
# plot(areas,col=col_lakes,border=NA,bg=NA,add=T)
# plot(ozette,col=col_lakes,border=NA,bg=NA,add=T)
plot(quinault,col=col_lakes,border=NA,bg=NA,add=T)
##-----------------------------------------------------## watersheds
if(color_map){col_huc<-alpha("deepskyblue3",0.2)}else{col_huc<-"gray10"}
plot(use_hucs,col=col_huc,border=NA,bg=NA,add=T)
# plot(use_hucs,col=col_huc,border="gray",lwd=1,bg=NA,add=T)
##----------------------------------------------------------## border
plot.bathy(bathys,add=T,image=FALSE,land=FALSE,n=0,xlim=xlim,ylim=ylim,ylab="",xlab="",axes=F,lwd=0.5,lty=1,col="black",drawlabels=FALSE)
##---------------------------------------------------## river names
rs<-c("Chehalis","Humptulips","Hoh","Queets","Quinault","Quillayute")
rlong<-c(-124.4, -124.45, -124.59, -124.55, -124.52, -124.59)
rlat<-c(46.95, 47.05, 47.75, 47.55, 47.35, 47.85)
text(rlong,rlat,rs,adj=c(0.5,0.5),cex=1.1,col="black")
##---------------------------------------------------------## scale
# map.scale(-122.5,46.45,ratio=F,relwidth=0.1,cex=0.5,lwd=0.5)
map.scale(-124.7,46.8,ratio=F,relwidth=0.12,cex=0.5,lwd=0.5)
##----------------------------------------------------------## axes
map.axes(cex.axis=0.8)
mtext("Longitude",side=1,line=2,outer=F,cex=1.2)
mtext("Latitude",side=2,line=2,outer=F,cex=1.2)
##-------------------------------------------------## add world map
x1<--190;x2<-450;y1<-10;y2<-450;par(usr=c(x1,x2,y1,y2))
map("worldHires",c("Canada","USA","Mexico"),fill=T,border=1,lwd=0.01, col="grey90",xlab="",ylab="",xlim=c(-220,-60),ylim=c(30,70),add=T)
points(center[1],center[2],pch=0,lwd=0.1,cex=0.4,col="firebrick")
rect(xleft=-180,ybottom=14,xright=-55,ytop=76,border=1,col=NA,lwd=0.1)
##------------------------------------------------------------## save
my_base<-recordPlot()  
dev.off()

##===================================================## combine plots
fig<-ggdraw() + 
   draw_plot(my_base,x=0.0,y=0.0,width=0.395,height=1) +
   draw_plot(gg_esc,x=0.353,y=0.6,width=0.588,height=0.38) +
   draw_plot(gg_age,x=0.36,y=0.0,width=0.63,height=0.6) +
   draw_plot_label(label=c("A","B","C"),x=c(0,0.37,0.37),y=c(1,1,0.6))

save_plot("../../output/map_with_rivers_plus_data.pdf",fig,ncol=2,nrow=2,base_height=4,base_width=10)

##=================================================================##
##=================================================================##
##=================================================================##