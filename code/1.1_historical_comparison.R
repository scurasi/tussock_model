rm(list=ls())

par(mfrow=c(1,1))
library(xlsx)
setwd("C:/Users/scurasi/Desktop/tussock_models/for_archive/data/")

curasi_2016<-read.xlsx("data_spreadsheets/curasi2016surveys.xlsx",1,startRow=2)
fetcher_1982<-read.xlsx("data_spreadsheets/fetcher1982surveys.xlsx",1,startRow=2)
curasi_2016$diam<-curasi_2016$diam_avg/2
curasi_2016$diam_avg<-curasi_2016$diam_avg/2
fetcher_1982$diam<-fetcher_1982$diam/2

past<-aggregate(diam~type,fetcher_1982,mean)
past<-cbind(past,list(err=aggregate(diam~type,fetcher_1982,function(x){sd(x)})[,2]))

present<-aggregate(diam~type,curasi_2016,mean)
present<-cbind(present,list(err=aggregate(diam~type,curasi_2016,function(x){sd(x)})[,2]))

#fig 3a
plot(c(1970,1977,2016),c(0,past$diam[1],present$diam[1]),xlim=c(1970,2020),ylim=c(0,22),font.axis=2,font.lab=2,cex.lab=1.5,xlab="year",ylab="Tussock radius (cm)",pch=15,cex=2,col="black",type="o",lwd=2,main="Eagle creek bladed")
arrows(x0=c(1977,2016),y0=c(past$diam[1]+past$err[1],present$diam[1]+present$err[1]),y1=c(past$diam[1]-past$err[1],present$diam[1]-present$err[1]),code=3,angle=90,length=0.1,lwd=2,col="black")
to_test_diam<-c(subset(curasi_2016,curasi_2016$type=="Bladed")$diam,subset(fetcher_1982,fetcher_1982$type=="Bladed")$diam)
to_test_year<-c(rep(2016,length(subset(curasi_2016,curasi_2016$type=="Bladed")$diam)),rep(1977,length(subset(fetcher_1982,fetcher_1982$type=="Bladed")$diam)))
summary(lm(to_test_diam~to_test_year))
rm(to_test_year,to_test_diam)

#fig 3b
plot(c(1977,2018),c(past$diam[2],present$diam[2]),xlim=c(1970,2020),ylim=c(2.5,22),font.axis=2,font.lab=2,cex.lab=1.5,xlab="year",ylab="Tussock radius (cm)",pch=15,cex=2,col=rgb(211/255,66/255,52/255,0.75),type="p",lwd=2,main="Cape thompson L41")
lines(c(1977,2018),c(past$diam[2],present$diam[2]),col=rgb(211/255,66/255,52/255,0.75),lwd=2,lty=2)
rect(c(1977,2018)+0.5,c(past$diam[2]+past$err[2],present$diam[2]+present$err[2]),c(1977,2018)-0.5,c(past$diam[2]-past$err[2],present$diam[2]-present$err[2]),col=rgb(211/255,66/255,52/255,0.4),border=NA)
to_test_diam<-c(subset(curasi_2016,curasi_2016$type=="L41")$diam,subset(fetcher_1982,fetcher_1982$type=="L41")$diam)
to_test_year<-c(rep(2018,length(subset(curasi_2016,curasi_2016$type=="L41")$diam)),rep(1977,length(subset(fetcher_1982,fetcher_1982$type=="L41")$diam)))
summary(lm(to_test_diam~to_test_year))
rm(to_test_year,to_test_diam)

points(c(1977,2018),c(past$diam[3],present$diam[3]),xlim=c(1970,2020),ylim=c(0,22),font.axis=2,font.lab=2,cex.lab=1.5,xlab="year",ylab="r_tussock (cm)",pch=15,cex=2,col=rgb(66/255,133/255,193/255,0.75),type="p",lwd=2,main="Cape thompson L52")
lines(c(1977,2018),c(past$diam[3],present$diam[3]),col=rgb(66/255,133/255,193/255,0.75),lwd=2,lty=2)
rect(c(1977,2018)+0.5,c(past$diam[3]+past$err[3],present$diam[3]+present$err[3]),c(1977,2018)-0.5,c(past$diam[3]-past$err[3],present$diam[3]-present$err[3]),col=rgb(66/255,133/255,193/255,0.4),border=NA)
to_test_diam<-c(subset(curasi_2016,curasi_2016$type=="L52")$diam,subset(fetcher_1982,fetcher_1982$type=="L52")$diam)
to_test_year<-c(rep(2018,length(subset(curasi_2016,curasi_2016$type=="L52")$diam)),rep(1977,length(subset(fetcher_1982,fetcher_1982$type=="L52")$diam)))
summary(lm(to_test_diam~to_test_year))
rm(to_test_year,to_test_diam)

points(c(1977,2016),c(past$diam[4],present$diam[4]),xlim=c(1970,2020),ylim=c(0,22),font.axis=2,font.lab=2,cex.lab=1.5,xlab="year",ylab="Tussock radius (cm)",pch=15,cex=2,col=rgb(93/255,178/255,104/255,0.75),type="p",lwd=2,main="Eagle creek undisturbed")
lines(c(1977,2016),c(past$diam[4],present$diam[4]),col=rgb(93/255,178/255,104/255,0.75),lwd=2,lty=2)
rect(c(1977,2016)+0.5,c(past$diam[4]+past$err[4],present$diam[4]+present$err[4]),c(1977,2016)-0.5,c(past$diam[4]-past$err[4],present$diam[4]-present$err[4]),col=rgb(93/255,178/255,104/255,0.4),border=NA)
to_test_diam<-c(subset(curasi_2016,curasi_2016$type=="Undisturbed")$diam,subset(fetcher_1982,fetcher_1982$type=="und")$diam)
to_test_year<-c(rep(2016,length(subset(curasi_2016,curasi_2016$type=="Undisturbed")$diam)),rep(1977,length(subset(fetcher_1982,fetcher_1982$type=="und")$diam)))
summary(lm(to_test_diam~to_test_year))
rm(to_test_year,to_test_diam)
