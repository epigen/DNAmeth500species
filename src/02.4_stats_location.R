source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(maps)
library(sf)

wd=file.path(analysis_dir,"02_vizStats/02.4_stats_location")
dir.create(wd)
setwd(wd)

locations=fread(file.path(meta_dir,"locations.csv"))
##stats_annot for filtering tumor and GPD samples + excluded
stats_annot=stats_annot[!grepl("_uc$",Sample_Name)&!grepl("Tumour|Cellline",Tissue)]

locations <- locations[Abbreviation %in% stats_annot$Sample_Name]
locations <- unique(locations)
locations[!grepl("USA",Location),Location:=sub(",",", ,",Location)]

locations_red=locations[,.(N_samples=length(unique(Abbreviation))),by=c("Institution","Location")]

locations_red[,City:=unlist(lapply(strsplit(Location,", "),"[[",1)),]
locations_red[,State:=unlist(lapply(strsplit(Location,", "),"[[",2)),]
locations_red[,Country:=unlist(lapply(strsplit(Location,", "),"[[",3)),]

locations_red[Institution=="fiwi",c("lat","lon"):=list(48.220704, 16.284393),]
locations_red[Institution=="University of Kentucky",c("lat","lon"):=list(38.034861, -84.504653),]
locations_red[Institution=="OGL",c("lat","lon"):=list(42.418756, -70.907325),]
locations_red[Institution=="Naschmarkt",c("lat","lon"):=list(48.198531, 16.363114),]
locations_red[Institution=="Vetmeduni",c("lat","lon"):=list(48.254559, 16.430109),]
locations_red[Institution=="MUW",c("lat","lon"):=list(48.234580, 16.349106),]
locations_red[Institution=="Biofish",c("lat","lon"):=list(48.215510, 16.334096),]
locations_red[Institution=="CCRI",c("lat","lon"):=list(48.216930, 16.344125),]
locations_red[Institution=="CeMM",c("lat","lon"):=list(48.219719, 16.349587),]
locations_red[Institution=="MPI for Evolutionary Biology",c("lat","lon"):=list(54.160409, 10.433943),]


locations_red[City=="Vienna",c("lat_c","lon_c"):=list(48.219056, 16.390189),]
locations_red[City=="Ploen",c("lat_c","lon_c"):=list(54.160409, 10.433943),]
locations_red[City=="Lexington",c("lat_c","lon_c"):=list(38.034861, -84.504653),]
locations_red[City=="Boston",c("lat_c","lon_c"):=list(42.418756, -70.907325),]

locations_red_red=locations_red[,.(N_samples=sum(N_samples),lat_c=unique(lat_c),lon_c=unique(lon_c)),by=c("City")]

#world
mapWorld <- borders("world",".", colour="gray50", fill="white",xlim =c(-90,25) ,ylim = c(45,60)) 
#Vienna ownly
shp <- st_read(file.path(meta_dir,"annotations_ext/BEZIRKSGRENZEOGD/BEZIRKSGRENZEOGDPolygon.shp"))

pdf(paste0("locations.pdf"),height=3,width=7)
ggplot(locations_red_red)+mapWorld+geom_point(aes(x=as.numeric(lon_c),y=as.numeric(lat_c),col=N_samples),alpha=0.7,size=8)+geom_text(aes(x=as.numeric(lon_c),y=as.numeric(lat_c),label=paste0(City,"\n",N_samples)))+scale_color_gradient(low="blue",high="red")+theme_bw()+coord_fixed()
ggplot(locations_red[Country=="Austria"]) + geom_sf(data = shp)+geom_point(aes(x=as.numeric(lon),y=as.numeric(lat),col=N_samples),alpha=0.7,size=8)+geom_text(aes(x=as.numeric(lon),y=as.numeric(lat),label=paste0(Institution,"\n",N_samples)))+scale_color_gradient(low="blue",high="red")+theme_bw()+coord_sf()
dev.off()

my_wt(locations_red, "locations_red.tsv")
my_wt(locations_red_red, "locations_red_red.tsv")
