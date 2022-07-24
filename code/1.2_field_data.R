#priors for tussock model based on field data
#and plots and formatting of other data

#--------------setup---------------#
rm(list=ls())
setwd("C:/Users/scurasi/Desktop/tussock_models/for_archive/data/")
library(xlsx);library(minpack.lm)
par(mar=c(5.1,5,4.1,1.2),mfrow=c(1,1))

#----------decomposition-------------#
#mound
decomp<-read.xlsx("data_spreadsheets/all_data.xlsx",2)
hist(decomp$k*365,cex.lab=1.5,cex.axis=1.15,xlab=expression(bold(paste("Mound decomposition constant (yr"^-1,")"))) ,ylab="Probability density",freq=F,main="",cex=1.1,font.lab=2,font=2,xlim=c(0,0.35),breaks=seq(0,0.35,0.025),xaxt='n')
round(mean(decomp$k*365),2)
round(sd(decomp$k*365)/sqrt(length(decomp$k*365)),2)
max(decomp$k*365)
min(decomp$k*365)
prior_data<-list(k_mound=decomp$k*365)
rm(decomp)

#tillers
decomp<-read.xlsx("data_spreadsheets/all_data.xlsx",3)
decomp$k<--log(decomp$mass_t2/decomp$mass_t1)/decomp$deployment_days

#seperate into live and dead
decomp_d<-subset(decomp,decomp$tussok_type=="d")
decomp<-subset(decomp,decomp$tussok_type=="l")

hist(decomp$k*365,cex.lab=1.5,cex.axis=1.15,xlab=expression(bold(paste("Tiller decomposition constant (yr"^-1,")"))) ,ylab="Probability density",freq=F,main="",cex=1.1,font.lab=2,font=2,xlim=c(0,0.35),breaks=seq(0,0.45,0.025),xaxt='n')
axis(1,at=seq(0,0.45,0.05),labels=seq(0,0.45,0.05),cex=1.1,font.lab=2,font=2)
round(mean(decomp$k*365),2)
round(sd(decomp$k*365)/sqrt(length(decomp$k*365)),2)
min(decomp$k*365)
max(decomp$k*365)
rm(decomp)

#dead tillers
hist(decomp_d$k*365,cex.lab=1.5,cex.axis=1.15,xlab=expression(bold(paste("Tiller decomposition constant (yr"^-1,")"))) ,ylab="Probability density",freq=F,main="",cex=1.1,font.lab=2,font=2,breaks=seq(0,0.45,0.025),xaxt='n',add=T,col="grey")
round(mean(decomp_d$k*365),2)
round(sd(decomp_d$k*365)/sqrt(length(decomp_d$k*365)),2)
min(decomp_d$k*365)
max(decomp_d$k*365)
prior_data$k_dead<-decomp_d$k*365
rm(decomp_d)

#----------bulk_density-------------#
bulk_d<-read.xlsx("data_spreadsheets/all_data.xlsx",4)
bulk_d<-bulk_d$bd[!is.na(bulk_d$bd)]*0.000001
hist(bulk_d,cex.lab=1.5,cex.axis=1.15,xlab=expression(bold(paste("bulk density (g cm"^-3,")"))) ,ylab="Density",main="",cex=1.1,font.lab=2,font=2,ylim=c(0,10),xlim=c(0,0.3),breaks=seq(0,0.5,0.025),freq=F)
round(mean(bulk_d),2)
round(sd(bulk_d)/sqrt(length(bulk_d)),2)
prior_data$p<-bulk_d
max(bulk_d)
min(bulk_d)
rm(bulk_d)

#----------tiller radius-------------#
radius<-read.xlsx("data_spreadsheets/all_data.xlsx",5)
radius<-aggregate(radius~id,radius,mean)[,2]
hist(radius/10,cex.lab=1.5,cex.axis=1.15,xlab=expression(bold(paste("Tiller radius (cm)"))) ,ylab="Probability density",freq=F,main="",cex=1.1,font.lab=2,font=2,xlim=c(0,3)/10,ylim=c(0,15),breaks=seq(0,3,0.25)/10)
round(mean(radius),2)/10
round(sd(radius)/sqrt(length(radius)),2)/10
max(radius)/10
min(radius)/10
prior_data$r_tiller<-radius/10
rm(radius)

#----------allocation to roots-------------#
temp_pluck<-read.xlsx("data_spreadsheets/all_data.xlsx",6)
hist(temp_pluck$roots_tiller,cex.lab=1.5,cex.axis=1.15,xlab=expression(bold("root mass per tiller (g)")) ,ylab="Probability density",main="",cex=1.1,font.lab=2,font=2,ylim=c(0,15),xlim=c(0,0.2),breaks=seq(0,0.5,0.025),freq=F)
max(temp_pluck$roots_tiller)
min(temp_pluck$roots_tiller)
prior_data$a<-temp_pluck$roots_tiller
rm(temp_pluck)
saveRDS(prior_data,"r_outputs/prior_data.rds")

#---------plots based on data used to fit the model---------#
#for the pluck
pluck<-read.xlsx("data_spreadsheets/all_data.xlsx",7)[,-1]
pluck$d<-pluck$d/2
colnames(pluck)[1]<-c("r")
saveRDS(pluck,"r_outputs/pluck.rds")

#fig s2a inset
plot(pluck$M/1000,pluck$r,ylab="Tussock radius (cm)",xlab="Mass (kg)",cex.lab=1.3,cex.axis=1.1,cex=1.2,font.lab=2,ylim=c(0,17),xlim=c(0,2250/1000))
points((pi*(seq(0,22,0.01)^3)*0.1178466)/1000,seq(0,22,0.01),type="l",lwd=2,lty=4)
points(((2/3)*pi*(seq(0,22,0.01)^3)*0.1178466)/1000,seq(0,22,0.01),type="l",lty=2,lwd=2)
points(((1/3)*pi*(seq(0,22,0.01)^3)*0.1178466)/1000,seq(0,22,0.01),type="l",lty=3,lwd=2)
points((1.6*pi*(seq(0,22,0.01)^3)*0.1178466)/1000,seq(0,22,0.01),type="l",lwd=2)

#fig s2a
barplot(width=1,xlim=c(0,5),space = c(0.1,0.5,0.5,0.5),c(mean(abs(pluck$r-(pluck$M/((1/3)*pi*0.1178466))^(1/3))),mean(abs(pluck$r-(pluck$M/((2/3)*pi*0.1178466))^(1/3))),mean(abs(pluck$r-(pluck$M/(pi*0.1178466))^(1/3))),mean(abs(pluck$r-(pluck$M/(1.6*pi*0.1178466))^(1/3)))),ylim=c(0,10),ylab="MAE (cm)",cex.lab=1.3,cex.axis=1.1,cex=1.2,font.lab=2)
axis(1,at=c(0.6,2.1,3.6,5.1))
mean(abs(pluck$r-(pluck$M/(pi*0.1178466))^(1/3)))
mean(abs(pluck$r-(pluck$M/(1.6*pi*0.1178466))^(1/3)))

#fig s2b
plot(pluck$r,pluck$height_total,font.axis=2,font.lab=2,cex.lab=1.5,cex.axis=1.25,pch=16,xlim=c(0,20),ylim=c(0,40))
lm(height_total~r,pluck)
abline(lm(height_total~r,pluck))

pluck$h2<-pluck$height_total-pluck$height_live
hist(pluck$h2/pluck$r)
mean(pluck$h2/pluck$r)

#surveys
diams<-read.xlsx("data_spreadsheets/all_data.xlsx",8)

#fig s2b inset
dh<-diams[complete.cases(diams),]
dh<-dh[dh$h>0,]
plot(dh$r,dh$h,font.axis=2,font.lab=2,cex.lab=1.5,cex.axis=1.25,lwd=1,col="grey55",pch=16,cex=0.5)
summary(lm(log(h)~log(r),dh))
points(seq(0,30,0.1),exp(0.004447+log(seq(0,30,0.1))*1.030828),type="l",lwd=2)
rm(dh)

#fig s5
diams<-diams[,-3]
boxplot(diams$r~diams$coor,ylab="Tussock radius (cm)",xlab="Latitude (dd.d)",font.axis=2,font.lab=2,cex.lab=1.5,cex.axis=1.25,xaxt='n')
abline(mean(diams$r),0,lwd=2,lty=2)

#table S1
climate<-read.xlsx("data_spreadsheets/all_data.xlsx",9)
cor(climate[,c(2,4,5,6)])^2
cor.test(climate$r,climate$gs)
cor.test(climate$r,climate$t)
cor.test(climate$r,climate$pr)

#fig s5 inset
plot(1:nrow(climate),climate$t[order(as.numeric(as.character(climate$coords)))],xlim=c(0,46),ylim=c(-13,-9),ylab="Mean annual temp. ",xlab="",font.axis=2,font.lab=2,cex.lab=1.5,cex.axis=1.25,yaxt='n',xaxt='n',type="o",lwd=2,pch=16)
axis(side=2, at = c(-13,-11,-9),font.axis=2,font.lab=2,cex.lab=1.5,cex.axis=1.25)
par(new = TRUE)
plot(1:nrow(climate),climate$pr[order(as.numeric(as.character(climate$coords)))],ylim=c(107,285),axes = FALSE, bty = "n", xlab = "", ylab = "",type="o",lwd=2,pch=16,col="grey50")
par(new = F)
axis(side=4, at = c(285,196,107),font.axis=2,font.lab=2,cex.lab=1.5,cex.axis=1.25)

saveRDS(diams,"r_outputs/diams.rds")
