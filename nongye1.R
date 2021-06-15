rm(list=ls())
gc()
library(PCICt)
library(climdex.pcic)
library(tidyverse)
library(lubridate)
library(RColorBrewer)
library(ggpubr)
#library(trend)
library(plyr)
library(maps)
library(mapdata)
library(sp)
library(maptools)
library(sf)
library(metR)
library(gstat)
library(ncdf4)
library(rgdal)
library(PCICt)
library(lubridate)
library(ggplot2)
library(sf)
library(RColorBrewer)
library(plyr)
library(gstat)
library(maptools)
library(ncdf4)
library(lubridate)
library(fields)
library(gridExtra)
library(raster)
library(rgeos)
library(sp)
library(cowplot)
#library(customLayout)
library(plyr)
library(rcolors)

#读取数据
ncdata1<-nc_open(filename = 'E:\\Documents\\R\\rdata\\data\\tmax.2020.nc')
ncdata2<-nc_open(filename = 'E:\\Documents\\R\\rdata\\data\\tmin.2020.nc')

#names(ncdata1$var) #查看变量名称 
#names(ncdata1$dim) #查看维度信息及名称
#ncdata1$dim$time$units #查看时间存储单位

lon1<-ncvar_get(nc=ncdata1,varid = 'lon')
lat1<-ncvar_get(nc=ncdata1,varid='lat')
time1<-ncvar_get(nc=ncdata1,varid = 'time')
time1<-as.PCICt(x='1900-01-01 00:00:00',cal='gregorian')+3600*time1
time1<-format.Date(x=time1,format='%Y-%m-%d')
dat1<-ncvar_get(nc=ncdata1,varid='tmax')
nc_close(ncdata1)

lon2<-ncvar_get(nc=ncdata2,varid = 'lon')
lat2<-ncvar_get(nc=ncdata2,varid='lat')
time2<-ncvar_get(nc=ncdata2,varid = 'time')
time2<-as.PCICt(x='1900-01-01 00:00:00',cal='gregorian')+3600*time2
time2<-format.Date(x=time2,format='%Y-%m-%d')
dat2<-ncvar_get(nc=ncdata2,varid='tmin')
nc_close(ncdata2)

#数据处理
lon.loc<-which(lon1>=70 & lon1<=140)
lat.loc<-which(lat1>=15 & lat1<=55)
#time.loc<-which(time=="2017-12-01"|time=="2018-01-01"|time=="2018-02-01")

tavg<-(dat1[lon.loc,lat.loc,]+dat2[lon.loc,lat.loc,])/2

ts<-array(dim = c(140,80,366)) #平滑温度
ts[,,1]<-(tavg[,,1]+tavg[,,2]+tavg[,,3])/3
ts[,,2]<-(tavg[,,1]+tavg[,,2]+tavg[,,3]+tavg[,,4])/4
ts[,,366]<-(tavg[,,364]+tavg[,,365]+tavg[,,366])/3
ts[,,365]<-(tavg[,,363]+tavg[,,364]+tavg[,,365]+tavg[,,366])/4
for (i in 3:364) {
  ts[,,i]<-(tavg[,,i-2]+tavg[,,i-1]+tavg[,,i]+tavg[,,i+1]+tavg[,,i+2])/5
}

#tt<-na.omit(ts)
tmp1=ifelse(is.na(ts),0,ts)
tmp2=ifelse(is.na(ts),10,ts)

m<-array(dim = c(140,80))
n<-array(dim=c(140,80))
for(i in 1:140){
  for(j in 1:80){
    n[i,j]=0
    for(k in 1:366){
      if(tmp1[i,j,k]>=5){
        n[i,j]=n[i,j]+1
      }
    }
  }
}

for(i in 1:140){
  for(j in 1:80){
    m[i,j]=0
    for(k in 1:366){
      if(tmp2[i,j,k]<5){
        m[i,j]=m[i,j]+1
        next
      }
    }
  }
}

#建立画图数据框

tslist<-list(lon=lon1[lon.loc],lat=lat1[lat.loc],time=time1,dat=n)
ts_mean<-apply(tslist$dat,c(1,2),mean)
lonlat<-expand.grid(lon=lon1[lon.loc],lat=lat1[lat.loc]) #生成对应格点的精度和纬度
alldata<-data.frame(lon=lonlat$lon,lat=lonlat$lat,ts_mean)

tslist1<-list(lon=lon1[lon.loc],lat=lat1[lat.loc],time=time1,dat=m)
ts_mean1<-apply(tslist1$dat,c(1,2),mean)
alldata1<-data.frame(lon=lonlat$lon,lat=lonlat$lat,ts_mean1)


#挑选在中国的点位
lonlat <- expand.grid(lon = lon, lat = lat)
grid <-SpatialPixelsDataFrame(points = lonlat,
                              data = data.frame(id = 1:nrow(lonlat)))

china <- readShapePoly('E:\\Documents\\R\\rdata\\data\\china area map\\bou1_4p.shp',
                       proj4string = CRS('+proj=longlat +ellps=WGS84'))
areas <- sapply(china@polygons, function(x) x@Polygons[[1]]@area)
china <- subset(china, areas == max(areas))
gridpoly <- as(grid, 'SpatialPolygonsDataFrame')
poly.clip <- raster::intersect(gridpoly, china)
#poly.clip@data$id??grid??Ӧ?ĸ???λ??
grid2 <- subset(grid, grid$id %in% poly.clip@data$id)


#画图
china <- readShapeLines(fn = "E:\\Documents\\R\\rdata\\data\\China boundary\\WGS84-1.shp") #本地地图 china <- fortify(china)

jpeg(filename='E:\\Documents\\R\\rdata\\save\\totaldays.png',units='mm',width=300,height=240,res=600)
p<-geom_path(data=china,mapping=aes(x=long,y=lat,group=group))
p12<-ggplot(data = alldata, mapping = aes(x=lon, y=lat))+coord_map(projection = 'mercator')+
  geom_contour_fill(mapping = aes(z = ts_mean),na.fill=TRUE) +
  scale_x_continuous(expand = c(0, 0), limits = c(70, 140)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(15, 55)) +
  scale_fill_gradientn(colors = brewer.pal(9,'YlOrBr'),breaks = seq(0,366,20)) + 
  ggtitle('Total Days - pass 5 °C')
p12 = p12 + p + guides(fill = guide_colourbar(barwidth = 1.5, barheight = 25))+
  labs(fill="days")+
  theme(axis.title = element_text(size=16),axis.text = element_text(size=14),
        legend.title = element_text(size=14),plot.title = element_text(hjust=0.5,size = 30))
p12 <- annotate_figure(p12,top = text_grob('',
                                             color = 'black', face = 'bold', size = 1),
                      fig.lab = '20181004233赵雨欣子',fig.lab.size = 15)
p12
dev.off()

jpeg(filename='E:\\Documents\\R\\rdata\\save\\startday.png',units='mm',width=300,height=240,res=600)
p<-geom_path(data=china,mapping=aes(x=long,y=lat,group=group))
p11<-ggplot(data = alldata1, mapping = aes(x=lon, y=lat))+coord_map(projection = 'mercator')+
  geom_contour_fill(mapping = aes(z = ts_mean1),na.fill=TRUE) +
  scale_x_continuous(expand = c(0, 0), limits = c(70, 140)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(15, 55)) +
  scale_fill_gradientn(colors = brewer.pal(9,'YlOrBr'),breaks = seq(0,366,20)) + 
  ggtitle('Start Day -  pass 5 °C')
p11 = p11 + p + guides(fill = guide_colourbar(barwidth = 1.5, barheight = 25))+
  labs(fill="day")+
  theme(axis.title = element_text(size=16),axis.text = element_text(size=14),
        legend.title = element_text(size=14),plot.title = element_text(hjust=0.5,size = 30))
p11 <- annotate_figure(p11,top = text_grob('',
                                           color = 'black', face = 'bold', size = 1),
                       fig.lab = '20181004233赵雨欣子',fig.lab.size = 15)
p11
dev.off()
