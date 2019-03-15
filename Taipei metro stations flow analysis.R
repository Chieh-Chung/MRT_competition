#建立最小統計區經緯度資料並畫圖
library(ggmap)
library(rgdal)
library(sp)
library(maptools)
library(dplyr)
library(ggplot2)

pathFolder1 <- "taipei_boundary" # 指定路徑到剛解壓縮後資料夾的那一層
tp <- readOGR(dsn = pathFolder1, layer = "G97_63000_U0200_2015", encoding = "UTF-8")# layer則是設定為此資料夾內的檔案名稱，此檔案名稱有多個同名不同格式的檔案
pathFolder2 <- "newtaipei_boundary"
nt<-readOGR(dsn = pathFolder2, layer = "G97_65000_U0200_2015", encoding = "UTF-8")

tp.twd <- fortify(tp,region = 'CODEBASE')
np.twd<- fortify(nt,region = 'CODEBASE')
gtp.twd<-rbind(tp.twd,np.twd)
twd97tolatlon <- function(x,y){ #轉換經緯度函數
  
  a = 6378137.0
  b = 6356752.3142451
  lon0 = 121 * pi / 180.0
  k0 = 0.9999
  dx = 250000
  dy = 0
  e = 1 - b^2 / a^2
  e2 = (1 - b^2 / a^2 ) / (b^2 / a^2)
  x = x - dx
  y = y - dy
  
  #Calculate the Meridional Arc
  M = y / k0
  
  #Calculate Footprint Latitude
  mu = M / (a * (1.0 - e/4.0 - 3*e^2 / 64.0 - 5 * e^3 / 256.0))
  e1 = (1.0 - (1.0 - e)^0.5) / (1.0 + (1.0 - e)^0.5)
  J1 = (3 * e1 / 2 - 27 * (e1^3) / 32.0)
  J2 = (21 * (e1)^2) / 16.0 - (55 * (e1)^4) / 32.0
  J3 = (151 * (e1)^3 / 96.0)
  J4 = (1097 * (e1)^4 / 512.0)
  fp = mu + J1 * sin(2*mu) + J2 * sin(4*mu) + J3 * sin(6*mu) + J4 * sin(8*mu)
  
  #Calculate Latitude and Longitude
  C1 = e2 * (cos(fp))^2
  T1 = (tan(fp))^2
  R1 = a * (1 - e) / (1 - e*(sin(fp))^2)^(3.0 / 2.0)
  N1 = a / (1 - e * (sin(fp))^2)^0.5
  
  D = x / (N1 * k0)
  
  #計算緯度
  Q1 = N1 * tan(fp) / R1
  Q2 = D^2/2.0
  Q3 = (5 + 3 * T1 + 10 * C1 - 4 * (C1)^2 - 9 * e2) * D^4 / 24.0
  Q4 = (61 + 90 * T1 + 298 * C1 + 45 * T1^2 - 3 * C1^2 - 252 * e2) * D^6 / 720.0
  lat = fp - Q1 * (Q2 - Q3 + Q4)
  
  #計算經度
  Q5 = D;
  Q6 = (1 + 2 * T1 + C1) * D^3 / 6.0
  Q7 = (5 - 2 * C1 + 28 * T1 - 3 * C1^2 + 8 * e2 + 24 * T1^2) * D^5 / 120.0
  lon = lon0 + (Q5 - Q6 + Q7) / cos(fp)
  
  lat1 = (lat * 180) / pi #緯度
  lon1 = (lon * 180) / pi #經度
  
  return(c(lat1,lon1))
  
}
latlon<-matrix(0,ncol = 2,nrow = nrow(gtp.twd))
for(i in 1:nrow(gtp.twd)){
  latlon[i,]<-twd97tolatlon(gtp.twd$long[i],gtp.twd$lat[i])
}

#將原先資料與轉換後的座標合併
gtp.latlon=cbind(lat=latlon[,1],long=latlon[,2],gtp.twd[,3:7])
#載入大台北地區底圖
map=get_googlemap(center = c(lon = 121.517532, lat = 25.056255), zoom = 12, size = c(640, 640), scale = 2)
ggmap(map,darken = c(0.3, "white")) +   #畫圖測試，map物件可在下方找到
  geom_polygon(data=gtp.latlon, aes(x = long, y = lat, group = group),color='red',
               size = 0.5,alpha = 0) +
  #scale_fill_gradient(low = "#FFFFB2", high = "#B10026",
  #limits = c(300, 2500)) +
  theme_bw(base_family = 'STHeiti') + 
  theme(axis.text = element_blank(), ## 移除經緯度
        axis.ticks = element_blank(),
        axis.title = element_blank())

#建立各捷運站中垂線經緯度資料，並將中垂線圖與googlemap底圖結合呈現
library(ggplot2)
library(ggmap)
library(deldir) 
MRT=read.csv(file("台北捷運經緯度.csv",encoding = 'big5'),header=T) #載入經緯度資料
#刪除重複的站(一個站可能同時有很多條線通過)
MRT=MRT[!duplicated(MRT$station_name_tw),]
MRT[,0]=c(1:108)
map=get_googlemap(center = c(lon = 121.517532, lat = 25.056255), zoom = 12, size = c(640, 640), scale = 2) #載入底圖
ggmap(map)
#資料建立部分
各站代號<-1:108
names(各站代號)<-unique(MRT$station_name_tw)
df=data.frame(name=MRT$station_name_tw,lat=MRT$lat,long=MRT$lon) #建立中垂線的資料
voronoi <- deldir(df$long, df$lat)
ddd<-voronoi$dirsgs[,1:6]
#dirsgs的x1x2y1y2代表邊界線的兩端點ind1和ind2代表隔開的是哪兩個點

#中垂線圖不會將邊界逼標上去，手動將邊界加上(須注意加上的線要連接好，之後會較好搜尋)
which(names(各站代號) %in% c('紅樹林','竹圍','新北投','大湖公園','東湖','南港軟體園區',
                         '南港展覽館','動物園','新店','頂埔','永寧','土城','迴龍','丹鳳',
                         '蘆洲'))#尋找邊界的站

ddd[c(which(ddd$ind1==90),which(ddd$ind2==90)),]
ddd[304,]<-c(ddd[35,3],ddd[35,4],ddd[33,3],ddd[33,4],1,0)      #手動編輯資料
ddd[305,]<-c(ddd[279,1],ddd[279,2],ddd[274,1],ddd[274,2],20,0)
ddd[306,]<-c(ddd[280,3],ddd[280,4],ddd[282,1],ddd[282,2],22,0)
ddd[307,]<-c(ddd[36,1],ddd[36,2],ddd[280,3],ddd[280,4],23,0)
ddd[308,]<-c(ddd[35,3],ddd[35,4],ddd[36,1],ddd[36,2],24,0)
ddd[309,]<-c(ddd[287,3],ddd[287,4],ddd[298,3],ddd[298,4],45,0)
ddd[310,]<-c(ddd[296,1],ddd[296,2],ddd[302,1],ddd[302,2],49,0)
ddd[311,]<-c(ddd[298,3],ddd[298,4],ddd[303,3],ddd[303,4],50,0)
ddd[312,]<-c(ddd[13,3],ddd[13,4],ddd[33,3],ddd[33,4],52,0)
ddd[313,]<-c(ddd[222,1],ddd[222,2],ddd[194,3],ddd[194,4],84,0)
ddd[314,]<-c(ddd[190,1],ddd[190,2],ddd[194,3],ddd[194,4],85,0)
ddd[315,]<-c(ddd[222,1],ddd[222,2],ddd[296,1],ddd[296,2],90,0)
ddd[316,]<-c(ddd[1,3],ddd[1,4],ddd[190,1],ddd[190,2],91,0)
ddd[317,]<-c(ddd[1,3],ddd[1,4],ddd[2,3],ddd[2,4],92,0)
ddd[318,]<-c(ddd[2,3],ddd[2,4],ddd[13,3],ddd[13,4],93,0)
ddd[319,]<-c(ddd[302,1],ddd[302,2],ddd[303,3],ddd[303,4],51,0)
ddd[279,2]<-ddd[282,2]
ddd[304,1:2]<-ddd[33,3:4]
ddd[304,3:4]<-ddd[35,3:4]


ggmap(map,darken = c(0.5, "white")) +              #畫圖
  #ggplot(data=df, aes(x=long,y=lat)) +
  #Plot the voronoi lines
  geom_segment(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    size = 1.2,
    data = ddd,
    linetype = 1,
    color= "#FFB958") +
  #Plot the points
  geom_point(
    aes(x = long, y = lat), data = df,
    fill=rgb(70,130,180,255,maxColorValue=255),
    pch=21,
    size = 3.5,
    color="#333333") 

ggplot(data=df, aes(x=long,y=lat)) +          #只劃出區域
  geom_segment(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    size = 1.5,
    data = ddd,
    linetype = 1,
    color= "#FFB958") +
  geom_point(
    aes(x = long, y = lat), data = df,
    fill=rgb(70,130,180,255,maxColorValue=255),
    pch=21,
    size = 3.5,
    color="#333333")

#將voronoi資料轉換成spatial polygon的資料型態
orderpoly<-function(x){    #將中垂線資料點排序的函數
  x<-data.matrix(x)
  cc<-as.numeric(x[1,])
  gg<-1
  for (i in 1:(nrow(x)-2)){
    x<-x[-gg[i],]
    for(j in 1:nrow(x)){
      if (sum(c(x[j,1:2]) == cc[(2*i+1):(2*i+2)])==2){cc<-c(cc,as.numeric(x[j,3:4])) 
      gg<-c(gg,j) }
      else{if(sum(c(x[j,3:4]) == cc[(2*i+1):(2*i+2)])==2){cc<-c(cc,as.numeric(x[j,1:2])) 
      gg<-c(gg,j)}}
    }
  }
  return(cc)
}
MRTpoly<-SpatialPolygons(list(Polygons(list(Polygon(matrix(orderpoly(ddd[c(which(ddd$ind1==1),
                                                                           which(ddd$ind2==1)),-5:-6]), ncol=2, byrow=TRUE))),ID=names(各站代號)[1])))
for (i in 2:108){   #開始轉換
  MRTpoly<-spRbind(MRTpoly,SpatialPolygons(list(Polygons(list(Polygon(matrix(orderpoly(ddd[c(which(ddd$ind1==i),
                                                                                             which(ddd$ind2==i)),-5:-6]), ncol=2, byrow=TRUE))),ID = names(各站代號)[i]))))
}
plot(MRTpoly)  #測試
#判斷捷運站的範圍包含哪幾個最小統計區
library(pracma)
#gtp.latlon為大台北最小統計區的邊界資訊，此步驟計算出每個統計區的中點
tp.midpoint<-gtp.latlon %>% group_by(id) %>% summarise(mid.long=mean(long),mid.lat =mean(lat))
#使用迴圈判斷各統計區的中點是否落在某個捷運涵蓋範圍內
MRTregion<-matrix(0,ncol=28767,nrow=108)
for(i in 1:108){
  s<-NULL
  s<-which(inpolygon(tp.midpoint$mid.long, tp.midpoint$mid.lat, MRTpoly@polygons[[i]]@Polygons[[1]]@coords[,1], MRTpoly@polygons[[i]]@Polygons[[1]]@coords[,2])==T)
  MRTregion[i,s]<-1
}
colnames(MRTregion)<-unique(gtp.latlon$id)


#將人口資料加上去
tpp<-read.csv("107年6月臺北市統計區人口統計_最小統計區.csv")
tpp<-tpp[-1,]
ntpp<-read.csv('107年6月新北市統計區人口統計_最小統計區.csv')
ntpp<-ntpp[-1,]

TPP<-rbind(tpp,ntpp)
TPP$P_CNT<-as.numeric(as.character(TPP$P_CNT))
g<-NULL
amount<-NULL
for( i in 1:108){
  g<-colnames(MRTregion)[which(MRTregion[i,]==1)]
  amount[i]<-sum(TPP[which(TPP$CODEBASE %in% g),5])
}
#各捷運站包含區域的總人口
MRTP<-data.frame(station_name = names(各站代號), population = amount)

#捷運流量檔
#mrtout
mrtout201601 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站出站量統計_201601.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtout201602 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站出站量統計_201602.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtout201603 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站出站量統計_201603.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtout201604 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站出站量統計_201604.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtout201605 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站出站量統計_201605.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtout201606 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站出站量統計_201606.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtout201607 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站出站量統計_201607.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtout201608 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站出站量統計_201608.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtout201609 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站出站量統計_201609.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtout201610 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站出站量統計_201610.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtout201611 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站出站量統計_201611.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtout201612 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站出站量統計_201612.csv",fileEncoding = "big5",stringsAsFactors = F ))

#mrtin
mrtin201601 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站進站量統計_201601.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtin201602 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站進站量統計_201602.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtin201603 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站進站量統計_201603.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtin201604 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站進站量統計_201604.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtin201605 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站進站量統計_201605.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtin201606 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站進站量統計_201606.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtin201607 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站進站量統計_201607.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtin201608 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站進站量統計_201608.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtin201609 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站進站量統計_201609.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtin201610 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站進站量統計_201610.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtin201611 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站進站量統計_201611.csv",fileEncoding = "big5",stringsAsFactors = F ))
mrtin201612 <- cutpaste(read.csv("/Users/leechenghsuan/Documents/政大資料競賽/臺北捷運各站進出量/臺北捷運各站進站量統計_201612.csv",fileEncoding = "big5",stringsAsFactors = F ))

for(i in 1:109){
  mrtin201601[,i] <- as.numeric(mrtin201601[,i])
  mrtin201602[,i] <- as.numeric(mrtin201602[,i])
  mrtin201603[,i] <- as.numeric(mrtin201603[,i])
  mrtin201604[,i] <- as.numeric(mrtin201604[,i])
  mrtin201605[,i] <- as.numeric(mrtin201605[,i])
  mrtin201606[,i] <- as.numeric(mrtin201606[,i])
  mrtin201607[,i] <- as.numeric(mrtin201607[,i])
  mrtin201608[,i] <- as.numeric(mrtin201608[,i])
  mrtin201609[,i] <- as.numeric(mrtin201609[,i])
  mrtin201610[,i] <- as.numeric(mrtin201610[,i])
  mrtin201611[,i] <- as.numeric(mrtin201611[,i])
  mrtin201612[,i] <- as.numeric(mrtin201612[,i])
  mrtout201601[,i] <- as.numeric(mrtout201601[,i])
  mrtout201602[,i] <- as.numeric(mrtout201602[,i])
  mrtout201603[,i] <- as.numeric(mrtout201603[,i])
  mrtout201604[,i] <- as.numeric(mrtout201604[,i])
  mrtout201605[,i] <- as.numeric(mrtout201605[,i])
  mrtout201606[,i] <- as.numeric(mrtout201606[,i])
  mrtout201607[,i] <- as.numeric(mrtout201607[,i])
  mrtout201608[,i] <- as.numeric(mrtout201608[,i])
  mrtout201609[,i] <- as.numeric(mrtout201609[,i])
  mrtout201610[,i] <- as.numeric(mrtout201610[,i])
  mrtout201611[,i] <- as.numeric(mrtout201611[,i])
  mrtout201612[,i] <- as.numeric(mrtout201612[,i])
  
}
mrtin201601 <- mrtin201601[,-1]
mrtin201602 <- mrtin201602[,-1]
mrtin201603 <- mrtin201603[,-1]
mrtin201604 <- mrtin201604[,-1]
mrtin201605 <- mrtin201605[,-1]
mrtin201606 <- mrtin201606[,-1]
mrtin201607 <- mrtin201607[,-1]
mrtin201608 <- mrtin201608[,-1]
mrtin201609 <- mrtin201609[,-1]
mrtin201610 <- mrtin201610[,-1]
mrtin201611 <- mrtin201611[,-1]
mrtin201612 <- mrtin201612[,-1]
mrtout201601 <- mrtout201601[,-1]
mrtout201602 <- mrtout201602[,-1]
mrtout201603 <- mrtout201603[,-1]
mrtout201604 <- mrtout201604[,-1]
mrtout201605 <- mrtout201605[,-1]
mrtout201606 <- mrtout201606[,-1]
mrtout201607 <- mrtout201607[,-1]
mrtout201608 <- mrtout201608[,-1]
mrtout201609 <- mrtout201609[,-1]
mrtout201610 <- mrtout201610[,-1]
mrtout201611 <- mrtout201611[,-1]
mrtout201612 <- mrtout201612[,-1]

mrtall201601 <- as.data.frame(as.matrix(mrtin201601)+as.matrix(mrtout201601))
mrtall201602 <- as.data.frame(as.matrix(mrtin201602)+as.matrix(mrtout201602))
mrtall201603 <- as.data.frame(as.matrix(mrtin201603)+as.matrix(mrtout201603))
mrtall201604 <- as.data.frame(as.matrix(mrtin201604)+as.matrix(mrtout201604))
mrtall201605 <- as.data.frame(as.matrix(mrtin201605)+as.matrix(mrtout201605))
mrtall201606 <- as.data.frame(as.matrix(mrtin201606)+as.matrix(mrtout201606))
mrtall201607 <- as.data.frame(as.matrix(mrtin201607)+as.matrix(mrtout201607))
mrtall201608 <- as.data.frame(as.matrix(mrtin201608)+as.matrix(mrtout201608))
mrtall201609 <- as.data.frame(as.matrix(mrtin201609)+as.matrix(mrtout201609))
mrtall201610 <- as.data.frame(as.matrix(mrtin201610)+as.matrix(mrtout201610))
mrtall201611 <- as.data.frame(as.matrix(mrtin201611)+as.matrix(mrtout201611))
mrtall201612 <- as.data.frame(as.matrix(mrtin201612)+as.matrix(mrtout201612))

mrtall201601_holiday <- mrtall201601[c(1,2,3,9,10,16,17,23,24,31),]
mrtall201602_holiday <- mrtall201602[c(6:14,20,21,27,28,29),]
mrtall201603_holiday <- mrtall201603[c(5,6,12,13,19,20,26,27),]
mrtall201604_holiday <- mrtall201604[c(2,3,4,5,9,10,16,17,23,24,30),]
mrtall201605_holiday <- mrtall201605[c(1,7,8,14,15,21,22,28,29),]
mrtall201606_holiday <- mrtall201606[c(5,9,10,11,12,18,19,25,26),]
mrtall201607_holiday <- mrtall201607[c(2,3,9,10,16,17,23,24,30,31),]
mrtall201608_holiday <- mrtall201608[c(6,7,13,14,20,21,27,28),]
mrtall201609_holiday <- mrtall201609[c(3,4,11,15,16,17,18,24,25),]
mrtall201610_holiday <- mrtall201610[c(1,2,8,9,10,15,16,22,23,29,30),]
mrtall201611_holiday <- mrtall201611[c(5,6,12,13,19,20,26,27),]
mrtall201612_holiday <- mrtall201612[c(3,4,10,11,17,18,24,25),]

mrtall201601_weekday <- mrtall201601[which(1-c(c(1:31)%in%c(1,2,3,9,10,16,17,23,24,31))==1),]
mrtall201602_weekday <- mrtall201602[which(1-c(c(1:29)%in%c(6:14,20,21,27,28,29))==1),]
mrtall201603_weekday <- mrtall201603[which(1-c(c(1:31)%in%c(5,6,12,13,19,20,26,27))==1),]
mrtall201604_weekday <- mrtall201604[which(1-c(c(1:30)%in%c(2,3,4,5,9,10,16,17,23,24,30))==1),]
mrtall201605_weekday <- mrtall201605[which(1-c(c(1:31)%in%c(1,7,8,14,15,21,22,28,29))==1),]
mrtall201606_weekday <- mrtall201606[which(1-c(c(1:30)%in%c(5,9,10,11,12,18,19,25,26))==1),]
mrtall201607_weekday <- mrtall201607[which(1-c(c(1:31)%in%c(2,3,9,10,16,17,23,24,30,31))==1),]
mrtall201608_weekday <- mrtall201608[which(1-c(c(1:31)%in%c(6,7,13,14,20,21,27,28))==1),]
mrtall201609_weekday <- mrtall201609[which(1-c(c(1:30)%in%c(3,4,11,15,16,17,18,24,25))==1),]
mrtall201610_weekday <- mrtall201610[which(1-c(c(1:31)%in%c(1,2,8,9,10,15,16,22,23,29,30))==1),]
mrtall201611_weekday <- mrtall201611[which(1-c(c(1:30)%in%c(5,6,12,13,19,20,26,27))==1),]
mrtall201612_weekday <- mrtall201612[which(1-c(c(1:31)%in%c(3,4,10,11,17,18,24,25,31))==1),]
mrtall2016_weekday <- rbind(mrtall201601_weekday,
                            mrtall201602_weekday,
                            mrtall201603_weekday,
                            mrtall201604_weekday,
                            mrtall201605_weekday,
                            mrtall201606_weekday,
                            mrtall201607_weekday,
                            mrtall201608_weekday,
                            mrtall201609_weekday,
                            mrtall201610_weekday,
                            mrtall201611_weekday,
                            mrtall201612_weekday)
mrtall2016_weekday_mu <- as.data.frame(apply(mrtall2016_weekday,2,mean))
rownames(mrtall2016_weekday_mu)[83] <- "台北101/世貿"
mrtall2016_weekday_mu <- as.data.frame(cbind(row.names(mrtall2016_weekday_mu),mrtall2016_weekday_mu) )
colnames(mrtall2016_weekday_mu) <- c("station_name","mean")
rownames(mrtall2016_weekday_mu) <- c(1:108)
saveRDS(mrtall2016_weekday_mu,file = "weekday1.rds")

mrtall2016_holiday <- rbind(mrtall201601_holiday,
                            mrtall201602_holiday,
                            mrtall201603_holiday,
                            mrtall201604_holiday,
                            mrtall201605_holiday,
                            mrtall201606_holiday,
                            mrtall201607_holiday,
                            mrtall201608_holiday,
                            mrtall201609_holiday,
                            mrtall201610_holiday,
                            mrtall201611_holiday,
                            mrtall201612_holiday)

mrtall2016_holiday_mu <- as.data.frame(apply(mrtall2016_holiday,2,mean))
rownames(mrtall2016_holiday_mu)[83] <- "台北101/世貿"
mrtall2016_holiday_mu <- as.data.frame(cbind(row.names(mrtall2016_holiday_mu),mrtall2016_holiday_mu) )
colnames(mrtall2016_holiday_mu) <- c("station_name","mean")
rownames(mrtall2016_holiday_mu) <- c(1:108)
dim(mrtall2016_holiday_mu)
saveRDS(mrtall2016_holiday_mu,file = "holiday1.rds")



#效率值分平日假日
holidays<-readRDS('/Users/theone/Desktop/MRT/holiday1.rds')
holidays<-data.frame(station_name = rownames(holidays),'流量.假日' = holidays$mrtall2016_holiday_mu) 
holidays<-merge(MRTP[,1:2],holidays,by='station_name',sort=F) %>% 
  mutate(efficiency=流量.假日/population)
plot(x = as.factor(1:108),y= holidays$efficiency,type = "n",ylim=c(0,2))
par(family="STKaiti") 
text(x= as.factor(1:108) ,y =holidays$efficiency,labels = holidays$station_name,cex=0.6,col='red')

weekd<-readRDS('/Users/theone/Desktop/MRT/weekday1.rds')
weekd<-data.frame(station_name = rownames(weekd),'流量.平日' = weekd$mrtall2016_weekday_mu) 
weekd<-merge(MRTP[,1:2],weekd,by='station_name',sort=F) %>% 
  mutate(efficiency=流量.平日/population)
plot(x = as.factor(1:108),y= weekd$efficiency,type = "n",ylim=c(0,2))
par(family="STKaiti") 
text(x= as.factor(1:108) ,y =weekd$efficiency,labels = weekd$station_name,cex=0.6,col='red')

#產業結構變數處理
#計算重疊區域的面積
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(maptools)


pathFolder <- "/Users/theone/Desktop/MRT/mapdata201701120616/鄉(鎮、市、區)界線檔(TWD97經緯度)1051214 2"
# layer則是設定為此資料夾內的檔案名稱，此檔案名稱有多個同名不同格式的檔案
tw3<- readOGR(dsn = pathFolder, layer = "TOWN_MOI_1051214", encoding = "UTF-8")
tw3<-tw3[tw3$COUNTYNAME %in% c('臺北市','新北市'),]
MRTpoly<-SpatialPolygonsDataFrame(MRTpoly, data.frame(MRT=names(各站代號)), match.ID=F)
projection(MRTpoly)<-projection(tw3)
pi<- intersect(tw3,MRTpoly)
plot(tw3, axes=T); plot(MRTpoly, add=T); plot(pi, add=T, col='red')
area<-data.frame(area=sapply(pi@polygons, FUN=function(x) {slot(x, 'area')}))
row.names(area) <- sapply(pi@polygons, FUN=function(x) {slot(x, 'ID')})
attArea <- spCbind(pi, area)
MRTindist<-aggregate(area~ TOWNNAME + MRT, data=attArea, FUN=sum)
MRTarea<-matrix(0,ncol=42,nrow=108) %>% data.frame()
colnames(MRTarea)<-c('MRT',as.character(tw3$TOWNNAME))
MRTarea$MRT<-names(各站代號)
for (i in 1:108){
  for (j in which(MRTindist$MRT == MRTarea$MRT[i])){
    g<-which(colnames(MRTarea) == MRTindist$TOWNNAME[j])
    MRTarea[i,g]<-MRTindist$area[j]
  }
}  
distarea <- data.frame( area = sapply(tw3@polygons, FUN=function(x){slot(x,'area')}))  
distarea$ID<-tw3$TOWNNAME

#計算所佔面積比例
MRTprop<-MRTarea
for (i in 2:42){
  MRTprop[,i]<-MRTprop[,i]/distarea$area[i-1]
}
MRTarea<-MRTarea[,-(which(apply(MRTarea[,-1],MARGIN = 2,FUN = sum)==0)+1)]
MRTprop<-MRTprop[,-(which(apply(MRTprop[,-1],MARGIN = 2,FUN = sum)==0)+1)]

#計算持有車輛數
motorcar<-read.csv(file('/Users/theone/Desktop/MRT/持有車輛data.csv',encoding = 'BIG5'),h=T)
MRTcars<-data.frame(MRT = MRTarea$MRT ,cars=0,motors=0)

for (i in 1:108){
  n<-NULL
  n1<-NULL
  n2<-NULL
  n<-colnames(MRTprop)[which(!MRTprop[i,-1]==0)+1]
  for (j in 1:length(n)){
    n1<-c(n1,motorcar[which(motorcar$Dist == n[j]),2] * MRTprop[i,which(colnames(MRTprop)==n[j])])
    n2<-c(n2,motorcar[which(motorcar$Dist == n[j]),3] * MRTprop[i,which(colnames(MRTprop)==n[j])])
  }
  MRTcars[i,2:3]<-c(sum(n1),sum(n2))
}
MRTcars


MRTdata<-cbind(MRTP,MRTcars[,-1])

#國土利用調查
國土利用<-read.csv(file('/Users/theone/Desktop/MRT/104年臺北市國土利用調查統計_最小統計區/104年臺北市國土利用調查統計_最小統計區.csv',encoding = 'BIG5'),h=T) %>%
  rbind(read.csv(file('/Users/theone/Desktop/MRT/104年新北市國土利用調查統計_最小統計區/104年新北市國土利用調查統計_最小統計區.csv',encoding = 'BIG5'),h=T))
colnames(國土利用)
國土利用<-國土利用[,c(3,4,14,17,28,29,30,33,34,35,36,39)] 
assetusing<-matrix(0,ncol=11,nrow=108)
for( i in 1:108){
  g<-NULL
  g<-colnames(MRTregion)[which(MRTregion[i,]==1)]
  for(j in 1:11){
    yy<-NULL
    yy<-sum(國土利用[which(國土利用[,1] %in% g),j+1])
    assetusing[i,j]<-yy
  }
}
colnames(assetusing)<-colnames(國土利用[,-1])
for (i in 1:108){
  assetusing[i,]<-assetusing[i,]/assetusing[i,1]
}
MRTdata<-cbind(MRTdata,data.frame(assetusing[,-1]))

#教育程度變數
education<-rbind(read.csv(file('/Users/theone/Desktop/MRT/105年12月臺北市統計區15歲以上人口教育程度統計_最小統計區/105年12月臺北市統計區15歲以上人口教育程度統計_最小統計區.csv',encoding = 'BIG5')),
                 read.csv(file('/Users/theone/Desktop/MRT/105年12月新北市統計區15歲以上人口教育程度統計_最小統計區/105年12月新北市統計區15歲以上人口教育程度統計_最小統計區.csv',encoding = 'BIG5')))
education<-education[,c(3:10)]
edu<-matrix(0,ncol=7,nrow=108)
for ( i in 1:108){
  g<-NULL
  g<-colnames(MRTregion)[which(MRTregion[i,]==1)]
  for(j in 1:7){
    yy<-NULL
    yy<-sum(education[which(education[,1] %in% g),j+1])
    edu[i,j]<-yy
  }
}
colnames(edu)<-colnames(education[,-1])
for(i in 1:108){
  edu[i,]<-edu[i,]/MRTdata$population[i]
}
MRTdata<- cbind(MRTdata,edu) 

#年齡分佈變數
age<-rbind(read.csv(file('/Users/theone/Desktop/MRT/105年6月臺北市統計區五歲年齡組性別人口統計_最小統計區/105年6月臺北市統計區五歲年齡組性別人口統計_最小統計區.csv',encoding = 'BIG5')),
           read.csv(file('/Users/theone/Desktop/MRT/105年6月新北市統計區五歲年齡組性別人口統計_最小統計區/105年6月新北市統計區五歲年齡組性別人口統計_最小統計區.csv',encoding = 'BIG5')))
age$未成年 <-apply(age[,c(4,7,10,13)],MARGIN = 1 , sum)
age$青年<-apply(age[,seq(16,40,by=3)],MARGIN = 1,sum)
age$中年<-apply(age[,seq(43,49,by=3)],MARGIN = 1,sum)
age$老年<-apply(age[,seq(52,64,by=3)],MARGIN = 1,sum)
age<-age[,c(3,68:71)]
oldman<-matrix(0,ncol = 4,nrow=108)
for(i in 1:108){
  g<-NULL
  g<-colnames(MRTregion)[which(MRTregion[i,]==1)]
  for(j in 1:4){
    yy<-NULL
    yy<-sum(age[which(age[,1] %in% g),j+1])
    oldman[i,j]<-yy
  }
}
colnames(oldman)<-colnames(age[,-1])
oldman<-data.frame(oldman)
oldman$total<-apply(oldman,MARGIN = 1,sum)
for(i in 1:108){
  oldman[i,]<-oldman[i,]/oldman[i,5]
}

MRTdata<-cbind(MRTdata,oldman[,-5])
MRTdata<-MRTdata[,-8]

###############################################
###計算google map上捷運站與最近Ubike站的距離###
###############################################

#####合併台北市與新北市ubike站資料表
#台北市ubike資料
library(jsonlite)
x=readLines("YouBikeTP.gz")
Encoding(x)="UTF-8"
x=fromJSON(x)
x=unlist(x$retVal)
#資料共有14個欄位，沒有遺漏值
ubike=data.frame(matrix(unname(x),ncol=14,byrow=T))
ubike=ubike[,-1]
#新北市ubike資料
y=readLines("ubike新北.txt")
y=fromJSON(y)
head(y)
ubike=rbind(ubike,y)
write.csv(ubike,"ubike.csv") ##檔案輸出
#整理台北捷運各站資訊
MRT$station_name_tw=as.character(MRT$station_name_tw)
Encoding(MRT$station_name_tw)="UTF-8"
MRT$address=as.character(MRT$address)
Encoding(MRT$address)="UTF-8"
location_MRT=MRT[,c(3,10,9)]
location_MRT$lat=round(location_MRT$lat,4)
#匯入台北新北ubike資訊(站名、經緯度...)
ubike=read.csv("ubike.csv")
location_ubike=ubike[,c(3,6,9,8)]
location_ubike$lat=round(location_ubike$lat,4)
location_ubike$lng=round(location_ubike$lng,4)
#避免過度取用goole API，先以經緯度計算離各捷運站直線距離最近的三個ubike站
library(geosphere)
nearest3=list()
for(i in 1:nrow(location_MRT)){
  xx=distm (c(location_MRT$lon[i],location_MRT$lat[i]),as.matrix(cbind(location_ubike$lng,location_ubike$lat)), fun = distHaversine)
  tmp=location_ubike[order(as.vector(xx),decreasing = F),][1:3,]
  nearest3[[i]]=tmp
}
###get distance from google map
library(googleway)
##set API key and try the function
APIkey="GOOGLEAPIKEY" 
test <- google_distance(origins = list(c(26.19660, -98.23591)),
                        destinations = list(c(31.62327, -94.64276),c(31.62327, -98.23276)), 
                        mode = "walking", key = APIkey, simplify = TRUE)
test$rows$elements
#整理欲輸入google map的地址名稱
area=substr(MRT$address,4,6)
location_MRT=cbind(MRT[,c(3,10,9)],area)
location_MRT$area=as.character(location_MRT$area)
location_ubike$sarea=as.character(location_ubike$sarea)
d.table=NULL
for(i in 1:121)
{d.table=c(d.table,g_distance[[i]][[1]][1]$distance$value)}
d.table=matrix(d.table,ncol=3,byrow=T)
apply(d.table,1,min)
#use address to calculate distance on gmap
#test
google_distance(origins = list("捷運木柵站"),
                destinations = list("捷運永春站"), 
                mode = "walking", key = APIkey, simplify = TRUE)
##使用地址計算google map上的距離
MRT_name=paste0("捷運",as.character(MRT[,3]))
ubike_name=paste0("youbike",as.character(ubike[,3]))
char_distance=list()
for(i in 1:nrow(MRT_location))
{ x=google_distance(origins = as.list(MRT_name)[i],
                    destinations = as.list(paste0("youbike",as.character(nearest3[[i]][,1]))), 
                    mode = "walking", key = APIkey, simplify = TRUE)
char_distance[[i]]=x$rows$elements
}
#將距離由小到大排序
distance_rank=NULL
for(i in 1:121)
{x=char_distance[[i]][[1]][1]$distance$value
distance_rank=c(distance_rank,x)
}
distance_rank=matrix(distance_rank,ncol=3,byrow=T)
##some problems happened but has been solved
distance_rank=read.csv("ubike_nearest3_gdistance.csv",header=F)
distance_rank=as.matrix(distance_rank)
distance_rank[is.na(distance_rank)]=1000000
##匯出資料(csv)
write.csv(cbind(as.character(MRT$station_name_tw),apply(distance_rank,1,min)),"捷運最近ubike站距離(公尺).csv")

######################################
########各捷運站區域內停車位推估######
######################################
#合併台北與新北公共停車場資訊
park=readOGR(dsn ='.',layer='park05',stringsAsFactors = F,encoding = "UTF-8",use_iconv=T)
library(jsonlite)
park2=readLines("新北市路外公共停車場資訊.txt")
#Encoding(park2)="UTF-8"
park2=fromJSON(park2)
head(park2)
library(rvest)
library(magrittr)
#預設座標  台北:TWD67 新北:TWD97，今欲取得經緯度
#使用"http://chihwai.blogspot.tw/2015/11/blog-post.html" 網站進行座標批次轉換
lalon.taipei=read.table("公共停車場經緯度_台北.txt",sep=",")
lalon.xinbei=read.table("公共停車場經緯度_新北.txt",sep=",")
names(lalon.taipei)=c("lon","lat")
names(lalon.xinbei)=c("lon","lat")

park.Taipei=data.frame(zone=park@data$ZONE1,name=park@data$GATENAME,car=park@data$NUM_S,motor=park@data$NUM_M,lon=lalon.taipei$lon,lat=lalon.taipei$lat)
park.Xinbei=data.frame(zone=park2$AREA,name=park2$NAME,car=park2$TOTALCAR,motor=park2$TOTALMOTOR,lon=lalon.xinbei$lon,lat=lalon.xinbei$lat)
park.Taipeicity=rbind(park.Taipei,park.Xinbei)
##匯入各捷運站劃分資訊
load("MRT.RData")
for(i in 1:108){
  Encoding(MRTpoly@polygons[[i]]@ID)="UTF-8"
}
library(pracma)
#inpolygon(要判斷的點的X軸,要判斷的點的Y軸,你選的範圍的X軸,你選的範圍的Y軸)
#判斷該點是否有在你的多邊形裡面，輸出為TRUE、FALSE判斷句
#接著把資訊指定到每個捷運站的欄位
##inpolygon(P[, 1], P[, 2], pg[, 1], pg[,2])
#extract coords
MRTpoly@polygons[[i]]@Polygons[[1]]@coords[,1]
#"motor" variable has some NA value
park.Taipeicity$motor[is.na(park.Taipeicity$motor)]=0
##計算捷運涵蓋範圍內的公共停車場數和可停放汽機車數
PARKinMRT=data.frame()
for(i in 1:108){
  x=inpolygon(park.Taipeicity$lon,park.Taipeicity$lat,MRTpoly@polygons[[i]]@Polygons[[1]]@coords[,1],MRTpoly@polygons[[i]]@Polygons[[1]]@coords[,2])
  PARKinMRT=rbind(PARKinMRT,c(sum(x),sum(as.numeric(park.Taipeicity$car[x])),sum(as.numeric(park.Taipeicity$motor[x]))))
}
names(PARKinMRT)=c("park_included","car","motor")
PARKinMRT=cbind(MRT=MRTlabel,PARKinMRT)
##輸出資料
write.csv(PARKinMRT,"各捷運站範圍內公共停車場車位數.csv")

#######變數處理和EDA_程式碼#######
#是否為轉運站
transta_var <- data.frame(station_names=c(  "動物園"       ,"木柵"         ,"萬芳社區"     ,"萬芳醫院"     ,"辛亥"         ,"麟光"        
                                            ,"六張犁"       ,"科技大樓"     ,"大安"         ,"忠孝復興"     ,"南京復興"     ,"中山國中"    
                                            ,"松山機場"     ,"大直"         ,"劍南路"       ,"西湖"         ,"港墘"         ,"文德"        
                                            ,"內湖"         ,"大湖公園"     ,"葫洲"         ,"東湖"         ,"南港軟體園區" ,"南港展覽館"  
                                            ,"象山"         ,"台北101/世貿" ,"信義安和"     ,"大安森林公園" ,"東門"         ,"中正紀念堂"  
                                            ,"中山"         ,"雙連"         ,"民權西路"     ,"圓山"         ,"劍潭"         ,"士林"        
                                            ,"芝山"         ,"明德"         ,"石牌"         ,"唭哩岸"       ,"奇岩"         ,"北投"        
                                            ,"新北投"       ,"復興崗"       ,"忠義"         ,"關渡"         ,"竹圍"         ,"紅樹林"      
                                            ,"淡水"         ,"新店"         ,"新店區公所"   ,"七張"         ,"小碧潭"       ,"大坪林"      
                                            ,"景美"         ,"萬隆"         ,"公館"         ,"台電大樓"     ,"古亭"         ,"小南門"      
                                            ,"西門"         ,"北門"         ,"松江南京"     ,"台北小巨蛋"   ,"南京三民"     ,"松山"        
                                            ,"南勢角"       ,"景安"         ,"永安市場"     ,"頂溪"         ,"忠孝新生"     ,"行天宮"      
                                            ,"中山國小"     ,"大橋頭"       ,"台北橋"       ,"菜寮"         ,"三重"         ,"先嗇宮"      
                                            ,"頭前庄"       ,"新莊"         ,"輔大"         ,"丹鳳"         ,"迴龍"         ,"三重國小"    
                                            ,"三和國中"     ,"徐匯中學"     ,"三民高中"     ,"蘆洲"         ,"頂埔"         ,"永寧"        
                                            ,"土城"         ,"海山"         ,"亞東醫院"     ,"府中"         ,"板橋"         ,"新埔"        
                                            ,"江子翠"       ,"龍山寺"       ,"善導寺"       ,"忠孝敦化"     ,"國父紀念館"   ,"市政府"      
                                            ,"永春"         ,"後山埤"       ,"昆陽"         ,"南港"     ) ,transfer_station = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1
                                                                                                                       ,0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                                                                                                                       ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

#附近有無 Youbike 站點
###############################################
###計算google map上捷運站與最近Ubike站的距離###
###############################################

#####合併台北市與新北市ubike站資料表
#台北市ubike資料
library(jsonlite)
x=readLines("YouBikeTP.gz")
Encoding(x)="UTF-8"
x=fromJSON(x)
x=unlist(x$retVal)
#資料共有14個欄位，沒有遺漏值
ubike=data.frame(matrix(unname(x),ncol=14,byrow=T))
ubike=ubike[,-1]
#新北市ubike資料
y=readLines("ubike新北.txt")
y=fromJSON(y)
head(y)
ubike=rbind(ubike,y)
write.csv(ubike,"ubike.csv") ##檔案輸出
#整理台北捷運各站資訊
MRT$station_name_tw=as.character(MRT$station_name_tw)
Encoding(MRT$station_name_tw)="UTF-8"
MRT$address=as.character(MRT$address)
Encoding(MRT$address)="UTF-8"
location_MRT=MRT[,c(3,10,9)]
location_MRT$lat=round(location_MRT$lat,4)
#匯入台北新北ubike資訊(站名、經緯度...)
ubike=read.csv("ubike.csv")
location_ubike=ubike[,c(3,6,9,8)]
location_ubike$lat=round(location_ubike$lat,4)
location_ubike$lng=round(location_ubike$lng,4)
#避免過度取用goole API，先以經緯度計算離各捷運站直線距離最近的三個ubike站
library(geosphere)
nearest3=list()
for(i in 1:nrow(location_MRT)){
  xx=distm (c(location_MRT$lon[i],location_MRT$lat[i]),as.matrix(cbind(location_ubike$lng,location_ubike$lat)), fun = distHaversine)
  tmp=location_ubike[order(as.vector(xx),decreasing = F),][1:3,]
  nearest3[[i]]=tmp
}
###get distance from google map
library(googleway)
##set API key and try the function
APIkey="GOOGLEAPIKEY" 
test <- google_distance(origins = list(c(26.19660, -98.23591)),
                        destinations = list(c(31.62327, -94.64276),c(31.62327, -98.23276)), 
                        mode = "walking", key = APIkey, simplify = TRUE)
test$rows$elements
#整理欲輸入google map的地址名稱
area=substr(MRT$address,4,6)
location_MRT=cbind(MRT[,c(3,10,9)],area)
location_MRT$area=as.character(location_MRT$area)
location_ubike$sarea=as.character(location_ubike$sarea)
d.table=NULL
for(i in 1:121)
{d.table=c(d.table,g_distance[[i]][[1]][1]$distance$value)}
d.table=matrix(d.table,ncol=3,byrow=T)
apply(d.table,1,min)
#use address to calculate distance on gmap
#test
google_distance(origins = list("捷運木柵站"),
                destinations = list("捷運永春站"), 
                mode = "walking", key = APIkey, simplify = TRUE)
##使用地址計算google map上的距離
MRT_name=paste0("捷運",as.character(MRT[,3]))
ubike_name=paste0("youbike",as.character(ubike[,3]))
char_distance=list()
for(i in 1:nrow(MRT_location))
{ x=google_distance(origins = as.list(MRT_name)[i],
                    destinations = as.list(paste0("youbike",as.character(nearest3[[i]][,1]))), 
                    mode = "walking", key = APIkey, simplify = TRUE)
char_distance[[i]]=x$rows$elements
}
#將距離由小到大排序
distance_rank=NULL
for(i in 1:121)
{x=char_distance[[i]][[1]][1]$distance$value
distance_rank=c(distance_rank,x)
}
distance_rank=matrix(distance_rank,ncol=3,byrow=T)
##some problems happened but has been solved
distance_rank=read.csv("ubike_nearest3_gdistance.csv",header=F)
distance_rank=as.matrix(distance_rank)
distance_rank[is.na(distance_rank)]=1000000
##匯出資料(csv)
write.csv(cbind(as.character(MRT$station_name_tw),apply(distance_rank,1,min)),"捷運最近ubike站距離(公尺).csv")

#公共停車場車位數
###############################################
###計算google map上捷運站與最近Ubike站的距離###
###############################################

#####合併台北市與新北市ubike站資料表
#台北市ubike資料
library(jsonlite)
x=readLines("YouBikeTP.gz")
Encoding(x)="UTF-8"
x=fromJSON(x)
x=unlist(x$retVal)
#資料共有14個欄位，沒有遺漏值
ubike=data.frame(matrix(unname(x),ncol=14,byrow=T))
ubike=ubike[,-1]
#新北市ubike資料
y=readLines("ubike新北.txt")
y=fromJSON(y)
head(y)
ubike=rbind(ubike,y)
write.csv(ubike,"ubike.csv") ##檔案輸出
#整理台北捷運各站資訊
MRT$station_name_tw=as.character(MRT$station_name_tw)
Encoding(MRT$station_name_tw)="UTF-8"
MRT$address=as.character(MRT$address)
Encoding(MRT$address)="UTF-8"
location_MRT=MRT[,c(3,10,9)]
location_MRT$lat=round(location_MRT$lat,4)
#匯入台北新北ubike資訊(站名、經緯度...)
ubike=read.csv("ubike.csv")
location_ubike=ubike[,c(3,6,9,8)]
location_ubike$lat=round(location_ubike$lat,4)
location_ubike$lng=round(location_ubike$lng,4)
#避免過度取用goole API，先以經緯度計算離各捷運站直線距離最近的三個ubike站
library(geosphere)
nearest3=list()
for(i in 1:nrow(location_MRT)){
  xx=distm (c(location_MRT$lon[i],location_MRT$lat[i]),as.matrix(cbind(location_ubike$lng,location_ubike$lat)), fun = distHaversine)
  tmp=location_ubike[order(as.vector(xx),decreasing = F),][1:3,]
  nearest3[[i]]=tmp
}
###get distance from google map
library(googleway)
##set API key and try the function
APIkey="GOOGLEAPIKEY" 
test <- google_distance(origins = list(c(26.19660, -98.23591)),
                        destinations = list(c(31.62327, -94.64276),c(31.62327, -98.23276)), 
                        mode = "walking", key = APIkey, simplify = TRUE)
test$rows$elements
#整理欲輸入google map的地址名稱
area=substr(MRT$address,4,6)
location_MRT=cbind(MRT[,c(3,10,9)],area)
location_MRT$area=as.character(location_MRT$area)
location_ubike$sarea=as.character(location_ubike$sarea)
d.table=NULL
for(i in 1:121)
{d.table=c(d.table,g_distance[[i]][[1]][1]$distance$value)}
d.table=matrix(d.table,ncol=3,byrow=T)
apply(d.table,1,min)
#use address to calculate distance on gmap
#test
google_distance(origins = list("捷運木柵站"),
                destinations = list("捷運永春站"), 
                mode = "walking", key = APIkey, simplify = TRUE)
##使用地址計算google map上的距離
MRT_name=paste0("捷運",as.character(MRT[,3]))
ubike_name=paste0("youbike",as.character(ubike[,3]))
char_distance=list()
for(i in 1:nrow(MRT_location))
{ x=google_distance(origins = as.list(MRT_name)[i],
                    destinations = as.list(paste0("youbike",as.character(nearest3[[i]][,1]))), 
                    mode = "walking", key = APIkey, simplify = TRUE)
char_distance[[i]]=x$rows$elements
}
#將距離由小到大排序
distance_rank=NULL
for(i in 1:121)
{x=char_distance[[i]][[1]][1]$distance$value
distance_rank=c(distance_rank,x)
}
distance_rank=matrix(distance_rank,ncol=3,byrow=T)
##some problems happened but has been solved
distance_rank=read.csv("ubike_nearest3_gdistance.csv",header=F)
distance_rank=as.matrix(distance_rank)
distance_rank[is.na(distance_rank)]=1000000
##匯出資料(csv)
write.csv(cbind(as.character(MRT$station_name_tw),apply(distance_rank,1,min)),"捷運最近ubike站距離(公尺).csv")


#公車路線數 
install.packages("rvest")
library(rvest)
busURL=readLines("/Users/leechenghsuan/Documents/政大資料競賽/200m內公車站牌網址.txt") #讀取各捷運站附近（200m內）站牌的網址
xx=paste0(ifelse(nchar(busURL)<20,"Y","N"),collapse = "") #區分各捷運站的站牌
ccc= unlist(strsplit(xx,"Y"))[-1]
an <- NULL
cn <- NULL
dn <- list()
for(i in 1:108){
  an <- c(an,nchar(ccc[i-1]))
  cn <- i
  bn <- sum(an)+sum(cn)
  dn[[i]] <- bn + which(unlist(strsplit(ccc[i],""))=="N")
}
a1_1   <- lapply(dn[[1  ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus"))) #抓出各捷運站的公車路線
a1_2   <- lapply(dn[[2  ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_3   <- lapply(dn[[3  ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_4   <- lapply(dn[[4  ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_5   <- lapply(dn[[5  ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_6   <- lapply(dn[[6  ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_7   <- lapply(dn[[7  ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_8   <- lapply(dn[[8  ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_9   <- lapply(dn[[9  ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_10  <- lapply(dn[[10 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_11  <- lapply(dn[[11 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_12  <- lapply(dn[[12 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_13  <- lapply(dn[[13 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_14  <- lapply(dn[[14 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_15  <- lapply(dn[[15 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_16  <- lapply(dn[[16 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_17  <- lapply(dn[[17 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_18  <- lapply(dn[[18 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_19  <- lapply(dn[[19 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_20  <- lapply(dn[[20 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_21  <- lapply(dn[[21 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_22  <- lapply(dn[[22 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_23  <- lapply(dn[[23 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_24  <- lapply(dn[[24 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_25  <- lapply(dn[[25 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_26  <- lapply(dn[[26 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_27  <- lapply(dn[[27 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_28  <- lapply(dn[[28 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_29  <- lapply(dn[[29 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_30  <- lapply(dn[[30 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_31  <- lapply(dn[[31 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_32  <- lapply(dn[[32 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_33  <- lapply(dn[[33 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_34  <- lapply(dn[[34 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_35  <- lapply(dn[[35 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_36  <- lapply(dn[[36 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_37  <- lapply(dn[[37 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_38  <- lapply(dn[[38 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_39  <- lapply(dn[[39 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_40  <- lapply(dn[[40 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_41  <- lapply(dn[[41 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_42  <- lapply(dn[[42 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_43  <- lapply(dn[[43 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_44  <- lapply(dn[[44 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_45  <- lapply(dn[[45 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_46  <- lapply(dn[[46 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_47  <- lapply(dn[[47 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_48  <- lapply(dn[[48 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_49  <- lapply(dn[[49 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_50  <- lapply(dn[[50 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_51  <- lapply(dn[[51 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_52  <- lapply(dn[[52 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_53  <- lapply(dn[[53 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_54  <- lapply(dn[[54 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_55  <- lapply(dn[[55 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_56  <- lapply(dn[[56 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_57  <- lapply(dn[[57 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_58  <- lapply(dn[[58 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_59  <- lapply(dn[[59 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_60  <- lapply(dn[[60 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_61  <- lapply(dn[[61 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_62  <- lapply(dn[[62 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_63  <- lapply(dn[[63 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_64  <- lapply(dn[[64 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_65  <- lapply(dn[[65 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_66  <- lapply(dn[[66 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_67  <- lapply(dn[[67 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_68  <- lapply(dn[[68 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_69  <- lapply(dn[[69 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_70  <- lapply(dn[[70 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_71  <- lapply(dn[[71 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_72  <- lapply(dn[[72 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_73  <- lapply(dn[[73 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_74  <- lapply(dn[[74 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_75  <- lapply(dn[[75 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_76  <- lapply(dn[[76 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_77  <- lapply(dn[[77 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_78  <- lapply(dn[[78 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_79  <- lapply(dn[[79 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_80  <- lapply(dn[[80 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_81  <- lapply(dn[[81 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_82  <- lapply(dn[[82 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_83  <- lapply(dn[[83 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_84  <- lapply(dn[[84 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_85  <- lapply(dn[[85 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_86  <- lapply(dn[[86 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_87  <- lapply(dn[[87 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_88  <- lapply(dn[[88 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_89  <- lapply(dn[[89 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_90  <- lapply(dn[[90 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_91  <- lapply(dn[[91 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_92  <- lapply(dn[[92 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_93  <- lapply(dn[[93 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_94  <- lapply(dn[[94 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_95  <- lapply(dn[[95 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_96  <- lapply(dn[[96 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_97  <- lapply(dn[[97 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_98  <- lapply(dn[[98 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_99  <- lapply(dn[[99 ]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_100 <- lapply(dn[[100]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_101 <- lapply(dn[[101]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_102 <- lapply(dn[[102]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_103 <- lapply(dn[[103]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_104 <- lapply(dn[[104]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_105 <- lapply(dn[[105]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_106 <- lapply(dn[[106]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_107 <- lapply(dn[[107]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
a1_108 <- lapply(dn[[108]], function(x) html_text(html_nodes(read_html(as.character(busURL[[x]])), ".auto-list-routelist-bus")))
t1_1   <- length(table(unlist(a1_1))) #計算各捷運站附近公車站牌的路線數
t1_2   <- length(table(unlist(a1_2)))
t1_3   <- length(table(unlist(a1_3)))
t1_4   <- length(table(unlist(a1_4)))
t1_5   <- length(table(unlist(a1_5)))
t1_6   <- length(table(unlist(a1_6)))
t1_7   <- length(table(unlist(a1_7)))
t1_8   <- length(table(unlist(a1_8)))
t1_9   <- length(table(unlist(a1_9)))
t1_10  <- length(table(unlist(a1_10)))
t1_11  <- length(table(unlist(a1_11)))
t1_12  <- length(table(unlist(a1_12)))
t1_13  <- length(table(unlist(a1_13)))
t1_14  <- length(table(unlist(a1_14)))
t1_15  <- length(table(unlist(a1_15)))
t1_16  <- length(table(unlist(a1_16 )))
t1_17  <- length(table(unlist(a1_17 )))
t1_18  <- length(table(unlist(a1_18 )))
t1_19  <- length(table(unlist(a1_19 )))
t1_20  <- length(table(unlist(a1_20 )))
t1_21  <- length(table(unlist(a1_21 )))
t1_22  <- length(table(unlist(a1_22 )))
t1_23  <- length(table(unlist(a1_23 )))
t1_24  <- length(table(unlist(a1_24 )))
t1_25  <- length(table(unlist(a1_25 )))
t1_26  <- length(table(unlist(a1_26 )))
t1_27  <- length(table(unlist(a1_27 )))
t1_28  <- length(table(unlist(a1_28 )))
t1_29  <- length(table(unlist(a1_29 )))
t1_30  <- length(table(unlist(a1_30 )))
t1_31  <- length(table(unlist(a1_31 )))
t1_32  <- length(table(unlist(a1_32 )))
t1_33  <- length(table(unlist(a1_33 )))
t1_34  <- length(table(unlist(a1_34 )))
t1_35  <- length(table(unlist(a1_35 )))
t1_36  <- length(table(unlist(a1_36 )))
t1_37  <- length(table(unlist(a1_37 )))
t1_38  <- length(table(unlist(a1_38 )))
t1_39  <- length(table(unlist(a1_39 )))
t1_40  <- length(table(unlist(a1_40 )))
t1_41  <- length(table(unlist(a1_41 )))
t1_42  <- length(table(unlist(a1_42 )))
t1_43  <- length(table(unlist(a1_43 )))
t1_44  <- length(table(unlist(a1_44 )))
t1_45  <- length(table(unlist(a1_45 )))
t1_46  <- length(table(unlist(a1_46 )))
t1_47  <- length(table(unlist(a1_47 )))
t1_48  <- length(table(unlist(a1_48 )))
t1_49  <- length(table(unlist(a1_49 )))
t1_50  <- length(table(unlist(a1_50 )))
t1_51  <- length(table(unlist(a1_51 )))
t1_52  <- length(table(unlist(a1_52 )))
t1_53  <- length(table(unlist(a1_53 )))
t1_54  <- length(table(unlist(a1_54 )))
t1_55  <- length(table(unlist(a1_55 )))
t1_56  <- length(table(unlist(a1_56 )))
t1_57  <- length(table(unlist(a1_57 )))
t1_58  <- length(table(unlist(a1_58 )))
t1_59  <- length(table(unlist(a1_59 )))
t1_60  <- length(table(unlist(a1_60 )))
t1_61  <- length(table(unlist(a1_61 )))
t1_62  <- length(table(unlist(a1_62 )))
t1_63  <- length(table(unlist(a1_63 )))
t1_64  <- length(table(unlist(a1_64 )))
t1_65  <- length(table(unlist(a1_65 )))
t1_66  <- length(table(unlist(a1_66 )))
t1_67  <- length(table(unlist(a1_67 )))
t1_68  <- length(table(unlist(a1_68 )))
t1_69  <- length(table(unlist(a1_69 )))
t1_70  <- length(table(unlist(a1_70 )))
t1_71  <- length(table(unlist(a1_71 )))
t1_72  <- length(table(unlist(a1_72 )))
t1_73  <- length(table(unlist(a1_73 )))
t1_74  <- length(table(unlist(a1_74 )))
t1_75  <- length(table(unlist(a1_75 )))
t1_76  <- length(table(unlist(a1_76 )))
t1_77  <- length(table(unlist(a1_77 )))
t1_78  <- length(table(unlist(a1_78 )))
t1_79  <- length(table(unlist(a1_79 )))
t1_80  <- length(table(unlist(a1_80 )))
t1_81  <- length(table(unlist(a1_81 )))
t1_82  <- length(table(unlist(a1_82 )))
t1_83  <- length(table(unlist(a1_83 )))
t1_84  <- length(table(unlist(a1_84 )))
t1_85  <- length(table(unlist(a1_85 )))
t1_86  <- length(table(unlist(a1_86 )))
t1_87  <- length(table(unlist(a1_87 )))
t1_88  <- length(table(unlist(a1_88 )))
t1_89  <- length(table(unlist(a1_89 )))
t1_90  <- length(table(unlist(a1_90 )))
t1_91  <- length(table(unlist(a1_91 )))
t1_92  <- length(table(unlist(a1_92 )))
t1_93  <- length(table(unlist(a1_93 )))
t1_94  <- length(table(unlist(a1_94 )))
t1_95  <- length(table(unlist(a1_95 )))
t1_96  <- length(table(unlist(a1_96 )))
t1_97  <- length(table(unlist(a1_97 )))
t1_98  <- length(table(unlist(a1_98 )))
t1_99  <- length(table(unlist(a1_99 )))
t1_100 <- length(table(unlist(a1_100)))
t1_101 <- length(table(unlist(a1_101)))
t1_102 <- length(table(unlist(a1_102)))
t1_103 <- length(table(unlist(a1_103)))
t1_104 <- length(table(unlist(a1_104)))
t1_105 <- length(table(unlist(a1_105)))
t1_106 <- length(table(unlist(a1_106)))
t1_107 <- length(table(unlist(a1_107)))
t1_108 <- length(table(unlist(a1_108)))
eachstanum <- mget(paste0("t1_",c(1:108))) #合併並儲存各捷運站的路線數
eachstanum <- do.call("rbind",eachstanum)
busnorepeat <- read.csv(file = "/Users/leechenghsuan/Documents/政大資料競賽/捷運經緯度（未重複）.csv",header = T,check.names = F,fileEncoding = "big5")

#車輛持有情形 
motorcar<-read.csv(file('/Users/theone/Desktop/MRT/持有車輛data.csv',encoding = 'BIG5'),h=T)
MRTcars<-data.frame(MRT = MRTarea$MRT ,cars=0,motors=0)

for (i in 1:108){
  n<-NULL
  n1<-NULL
  n2<-NULL
  n<-colnames(MRTprop)[which(!MRTprop[i,-1]==0)+1]
  for (j in 1:length(n)){
    n1<-c(n1,motorcar[which(motorcar$Dist == n[j]),2] * MRTprop[i,which(colnames(MRTprop)==n[j])])
    n2<-c(n2,motorcar[which(motorcar$Dist == n[j]),3] * MRTprop[i,which(colnames(MRTprop)==n[j])])
  }
  MRTcars[i,2:3]<-c(sum(n1),sum(n2))
}
MRTcars


MRTdata<-cbind(MRTP,MRTcars[,-1])

#產業結構 
#產業結構變數處理
#計算重疊區域的面積
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(maptools)


pathFolder <- "/Users/theone/Desktop/MRT/mapdata201701120616/鄉(鎮、市、區)界線檔(TWD97經緯度)1051214 2"
# layer則是設定為此資料夾內的檔案名稱，此檔案名稱有多個同名不同格式的檔案
tw3<- readOGR(dsn = pathFolder, layer = "TOWN_MOI_1051214", encoding = "UTF-8")
tw3<-tw3[tw3$COUNTYNAME %in% c('臺北市','新北市'),]
MRTpoly<-SpatialPolygonsDataFrame(MRTpoly, data.frame(MRT=names(各站代號)), match.ID=F)
projection(MRTpoly)<-projection(tw3)
pi<- intersect(tw3,MRTpoly)
plot(tw3, axes=T); plot(MRTpoly, add=T); plot(pi, add=T, col='red')
area<-data.frame(area=sapply(pi@polygons, FUN=function(x) {slot(x, 'area')}))
row.names(area) <- sapply(pi@polygons, FUN=function(x) {slot(x, 'ID')})
attArea <- spCbind(pi, area)
MRTindist<-aggregate(area~ TOWNNAME + MRT, data=attArea, FUN=sum)
MRTarea<-matrix(0,ncol=42,nrow=108) %>% data.frame()
colnames(MRTarea)<-c('MRT',as.character(tw3$TOWNNAME))
MRTarea$MRT<-names(各站代號)
for (i in 1:108){
  for (j in which(MRTindist$MRT == MRTarea$MRT[i])){
    g<-which(colnames(MRTarea) == MRTindist$TOWNNAME[j])
    MRTarea[i,g]<-MRTindist$area[j]
  }
}  
distarea <- data.frame( area = sapply(tw3@polygons, FUN=function(x){slot(x,'area')}))  
distarea$ID<-tw3$TOWNNAME

#計算所佔面積比例
MRTprop<-MRTarea
for (i in 2:42){
  MRTprop[,i]<-MRTprop[,i]/distarea$area[i-1]
}
MRTarea<-MRTarea[,-(which(apply(MRTarea[,-1],MARGIN = 2,FUN = sum)==0)+1)]
MRTprop<-MRTprop[,-(which(apply(MRTprop[,-1],MARGIN = 2,FUN = sum)==0)+1)]

#國土利用情形 
國土利用<-read.csv(file('/Users/theone/Desktop/MRT/104年臺北市國土利用調查統計_最小統計區/104年臺北市國土利用調查統計_最小統計區.csv',encoding = 'BIG5'),h=T) %>%
  rbind(read.csv(file('/Users/theone/Desktop/MRT/104年新北市國土利用調查統計_最小統計區/104年新北市國土利用調查統計_最小統計區.csv',encoding = 'BIG5'),h=T))
colnames(國土利用)
國土利用<-國土利用[,c(3,4,14,17,28,29,30,33,34,35,36,39)] 
assetusing<-matrix(0,ncol=11,nrow=108)
for( i in 1:108){
  g<-NULL
  g<-colnames(MRTregion)[which(MRTregion[i,]==1)]
  for(j in 1:11){
    yy<-NULL
    yy<-sum(國土利用[which(國土利用[,1] %in% g),j+1])
    assetusing[i,j]<-yy
  }
}
colnames(assetusing)<-colnames(國土利用[,-1])
for (i in 1:108){
  assetusing[i,]<-assetusing[i,]/assetusing[i,1]
}
MRTdata<-cbind(MRTdata,data.frame(assetusing[,-1]))

#年齡與教育程度 
#年齡分佈變數
age<-rbind(read.csv(file('/Users/theone/Desktop/MRT/105年6月臺北市統計區五歲年齡組性別人口統計_最小統計區/105年6月臺北市統計區五歲年齡組性別人口統計_最小統計區.csv',encoding = 'BIG5')),
           read.csv(file('/Users/theone/Desktop/MRT/105年6月新北市統計區五歲年齡組性別人口統計_最小統計區/105年6月新北市統計區五歲年齡組性別人口統計_最小統計區.csv',encoding = 'BIG5')))
age$未成年 <-apply(age[,c(4,7,10,13)],MARGIN = 1 , sum)
age$青年<-apply(age[,seq(16,40,by=3)],MARGIN = 1,sum)
age$中年<-apply(age[,seq(43,49,by=3)],MARGIN = 1,sum)
age$老年<-apply(age[,seq(52,64,by=3)],MARGIN = 1,sum)
age<-age[,c(3,68:71)]
oldman<-matrix(0,ncol = 4,nrow=108)
for(i in 1:108){
  g<-NULL
  g<-colnames(MRTregion)[which(MRTregion[i,]==1)]
  for(j in 1:4){
    yy<-NULL
    yy<-sum(age[which(age[,1] %in% g),j+1])
    oldman[i,j]<-yy
  }
}
colnames(oldman)<-colnames(age[,-1])
oldman<-data.frame(oldman)
oldman$total<-apply(oldman,MARGIN = 1,sum)
for(i in 1:108){
  oldman[i,]<-oldman[i,]/oldman[i,5]
}

MRTdata<-cbind(MRTdata,oldman[,-5])
MRTdata<-MRTdata[,-8]
#教育程度變數
education<-rbind(read.csv(file('/Users/theone/Desktop/MRT/105年12月臺北市統計區15歲以上人口教育程度統計_最小統計區/105年12月臺北市統計區15歲以上人口教育程度統計_最小統計區.csv',encoding = 'BIG5')),
                 read.csv(file('/Users/theone/Desktop/MRT/105年12月新北市統計區15歲以上人口教育程度統計_最小統計區/105年12月新北市統計區15歲以上人口教育程度統計_最小統計區.csv',encoding = 'BIG5')))
education<-education[,c(3:10)]
edu<-matrix(0,ncol=7,nrow=108)
for ( i in 1:108){
  g<-NULL
  g<-colnames(MRTregion)[which(MRTregion[i,]==1)]
  for(j in 1:7){
    yy<-NULL
    yy<-sum(education[which(education[,1] %in% g),j+1])
    edu[i,j]<-yy
  }
}
colnames(edu)<-colnames(education[,-1])
for(i in 1:108){
  edu[i,]<-edu[i,]/MRTdata$population[i]
}
MRTdata<- cbind(MRTdata,edu) 

#############################
######探索性資料分析#########
#############################
##觀察各個變數與效率值的關係
MRT_final=read.csv("MRT_final.csv")
View(MRT_final)
MRT_final=MRT_final[order(MRT_final$effweekday,decreasing = T),]
library(plotly)
p1=plot_ly(MRT_final)%>%
  add_trace(x=as.list(MRT_final$station_name),y=MRT_final$cars,type="bar",name="汽車持有數", marker = list(color = '#20B2AA')) %>%
  add_trace(x=as.list(MRT_final$station_name),y=MRT_final$effweekday,type="scatter",mode="lines",name="效率值",yaxis="y2",line = list(color = '#45171D')) %>%
  layout(title=F,xaxis=list(title=""),yaxis=list(side="left",title="汽車持有數",showgrid=F,zeroline=F),
         yaxis2=list(side = 'right', overlaying = "y", title = '效率值', showgrid = FALSE, zeroline = FALSE))
cor(MRT_final$effweekday,MRT_final$cars)
#汽機車持有數
p2=plot_ly(MRT_final)%>%
  add_trace(x=as.list(MRT_final$station_name),y=MRT_final$motors,type="bar",name="機車持有數", marker = list(color = '#20B2AA')) %>%
  add_trace(x=as.list(MRT_final$station_name),y=MRT_final$effweekday,type="scatter",mode="lines",name="效率值",yaxis="y2",line = list(color = '#45171D')) %>%
  layout(title=F,xaxis=list(title=""),yaxis=list(side="left",title="機車持有數",showgrid=F,zeroline=F),
         yaxis2=list(side = 'right', overlaying = "y", title = '效率值', showgrid = FALSE, zeroline = FALSE))
cor(MRT_final$effweekday,MRT_final$motors)
#停車格數量
p3=plot_ly(MRT_final)%>%
  add_trace(x=as.list(MRT_final$station_name),y=MRT_final$park_motor,type="bar",name="機車停車格數量", marker = list(color = '#20B2AA')) %>%
  add_trace(x=as.list(MRT_final$station_name),y=MRT_final$effweekday,type="scatter",mode="lines",name="效率值",yaxis="y2",line = list(color = '#45171D')) %>%
  layout(title=F,xaxis=list(title=""),yaxis=list(side="left",title="機車停車格數量",showgrid=F,zeroline=F),
         yaxis2=list(side = 'right', overlaying = "y", title = '效率值', showgrid = FALSE, zeroline = FALSE))
cor(MRT_final$effweekday,MRT_final$park_motor)
#公車路線數
p4=plot_ly(MRT_final)%>%
  add_trace(x=as.list(MRT_final$station_name),y=MRT_final$bus,type="bar",name="公車路線數", marker = list(color = '#20B2AA')) %>%
  add_trace(x=as.list(MRT_final$station_name),y=MRT_final$effweekday,type="scatter",mode="lines",name="效率值",yaxis="y2",line = list(color = '#45171D')) %>%
  layout(title=F,xaxis=list(title=""),yaxis=list(side="left",title="公車路線數",showgrid=F,zeroline=F),
         yaxis2=list(side = 'right', overlaying = "y", title = '效率值', showgrid = FALSE, zeroline = FALSE))
cor(MRT_final$effweekday,MRT_final$bus)
#製造業從業人口比例
plot_ly(MRT_final)%>%
  add_trace(x=as.list(MRT_final$station_name),y=MRT_final$industry,type="bar",name="製造業從業人口比例", marker = list(color = '#20B2AA')) %>%
  add_trace(x=as.list(MRT_final$station_name),y=MRT_final$effweekday,type="scatter",mode="lines",name="效率值",yaxis="y2",line = list(color = '#45171D')) %>%
  layout(title=F,xaxis=list(title=""),yaxis=list(side="left",title="製造業從業人口比例",showgrid=F,zeroline=F),
         yaxis2=list(side = 'right', overlaying = "y", title = '效率值', showgrid = FALSE, zeroline = FALSE))
cor(MRT_final$effweekday,MRT_final$industry)
#用地規劃
plot_ly(MRT_final)%>%
  add_trace(x=as.list(MRT_final$station_name),y=MRT_final$商業,type="bar",name="商業用地比例", marker = list(color = '#20B2AA')) %>%
  add_trace(x=as.list(MRT_final$station_name),y=MRT_final$effweekday,type="scatter",mode="lines",name="效率值",yaxis="y2",line = list(color = '#45171D')) %>%
  layout(title=F,xaxis=list(title=""),yaxis=list(side="left",title="商業用地比例",showgrid=F,zeroline=F),
         yaxis2=list(side = 'right', overlaying = "y", title = '效率值', showgrid = FALSE, zeroline = FALSE))
cor(MRT_final$effweekday,MRT_final$商業)
#人口結構
p5 <- plot_ly(MRT_final, x = as.list(MRT_final$station_name), y = MRT_final$未成年, type = 'bar', name = '未成年') %>%
  add_trace(y = MRT_final$青年, name = '青年') %>%
  add_trace(y = MRT_final$中年, name = '中年') %>%
  add_trace(y = MRT_final$老年, name = '老年') %>%
  layout(yaxis = list(title = 'percentage'), barmode = 'stack')

plot_ly(MRT_final, x = as.list(MRT_final$station_name), y = MRT_final$商業, type = 'bar', name = '商業') %>%
  add_trace(y = MRT_final$住宅, name = '住宅') %>%
  add_trace(y = MRT_final$工業, name = '工業') %>%
  add_trace(y = MRT_final$政府機關, name = '政府機關') %>%
  add_trace(y = MRT_final$學校, name = '學校') %>%
  add_trace(y = MRT_final$醫療保健, name = '醫療保健') %>%
  add_trace(y = MRT_final$社會福利設施, name = '社會福利設施') %>%
  add_trace(y = MRT_final$遊憩使用土地總數, name = '遊憩使用') %>%
  add_trace(y = MRT_final$交通使用土地總數, name = '交通使用') %>%
  layout(yaxis = list(title = 'percentage'), barmode = 'stack'
         #繪製變數間相關係數熱圖
         MRT=read.csv("MRT_final.csv")
         head(MRT)
         colnames(MRT)[2]="有無鄰近ubike站"
         colnames(MRT)[3]="是否為轉運站"
         colnames(MRT)[5]="汽車持有數"
         colnames(MRT)[6]="機車持有數"
         colnames(MRT)[27]="附近公車路線數"
         colnames(MRT)[28]="汽車停車格數"
         colnames(MRT)[29]="機車停車格數"
         colnames(MRT)[30]="交通使用土地"
         colnames(MRT)[31]="製造業從業比例"
         colnames(MRT)[33]="平日效率值"
         library(reshape2)
         library(ggplot2)
         corr <- cor(MRT[, -c(1,4,32)], method="pearson")
         ggplot(melt(corr, varnames=c("x", "y"), value.name="correlation"), 
                aes(x=x, y=y)) +
           geom_tile(aes(fill=correlation)) +
           scale_fill_gradient2(low="green", mid="yellow", high="red",
                                guide=guide_colorbar(ticks=FALSE, barheight = 5),
                                limits=c(-1,1)) + 
           theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
           labs(title="Heatmap of Correlation Matrix", 
                x=NULL, y=NULL)
         
         fitweekday <- lm(log(effweekday)~ .,data=MRTfinal5) #平日完整迴歸模型
         summary(fitweekday)
         nortest::pearson.test(fitweekday$residuals) #殘差normal 檢定
         car::ncvTest(fitweekday)#constant variance 檢定
         
         
         fitholiday <- lm(log(effholiday)~.,data=MRTfinal4)  #假日完整迴歸模型
         summary(fitholiday)
         nortest::pearson.test(fitholiday$residuals) #殘差normal 檢定
         car::ncvTest(fitholiday)#constant variance 檢定
         
         stepweekday <- stepAIC(fitweekday, direction="both") #平日逐步迴歸
         summary(stepweekday)
         nortest::pearson.test(stepweekday$residuals) #殘差normal 檢定
         car::ncvTest(stepweekday)#constant variance 檢定
         
         stepholiday <- stepAIC(fitholiday, direction="both") #假日逐步迴歸
         summary(stepholiday)
         nortest::pearson.test(stepholiday$residuals) #殘差normal 檢定
         car::ncvTest(stepholiday)#constant variance 檢定
         
         
         test<-list()
         test2<-list()
         test3<-NULL
         for(i in 1: 1000){
           test[[i]]<-sample(1:106,90,replace=T)
           test2[[i]]<-MRTfinal4[test[[i]],]
           fitholiday2 <- lm(log(effholiday)~.,data=test2[[i]])  
           stepholiday2 <- stepAIC(fitholiday2, direction="both") #平日逐步迴歸
           test3<-c(test3,names(stepholiday2$coefficients))
         }   #假日逐步回歸模型穩定度測試
         
         table(factor(test3))
         names(stepholiday$coefficients)
         
         
         
         
         test4<-list()
         test5<-list()
         test6<-NULL
         for(i in 1: 1000){
           test4[[i]]<-sample(1:106,90,replace=T)
           test5[[i]]<-MRTfinal5[test4[[i]],]
           fitweekday2 <- lm(log(effweekday)~.,data=test5[[i]])  
           stepweekday2 <- stepAIC(fitweekday2, direction="both") #平日逐步迴歸
           test6<-c(test6,names(stepweekday2$coefficients))
         }        #平日逐步回歸模型穩定度測試
         table(factor(test6))
         names(stepweekday$coefficients)
         
         
         