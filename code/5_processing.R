###########################
############setup##########
###########################
rm(list=ls())
setwd("C:/Users/scurasi/Desktop/tussock_models/for_archive/data/")
load(file="r_outputs/Step2.Rdata")
library(psych);library(FME);library(xlsx);library(fields)

#Fig S6
par(mfrow=c(1,1))
pairs.panels(param.keep, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = FALSE, # show correlation ellipses
             smooth = FALSE,
             lm = TRUE)

#table s2
cor(param.keep)

#fig 6/table 1
par(mfrow=c(3,3))
for(plot_index in 1:n.par){
  
  prior_data<-readRDS("prior_data.rds")
  
  if(names(param.keep)[plot_index] %in% names(prior_data)){
  
  prior_inter<-prior_data[which(names(param.keep)[plot_index] == names(prior_data))]  
    
  x_max<-max(c(param.keep[,plot_index],unlist(unname(prior_inter)),param.max[plot_index]))
  x_min<-min(c(param.keep[,plot_index],unlist(unname(prior_inter)),param.min[plot_index]))
  y_max<-max(density(x = param.keep[,plot_index])$y,density(x = unlist(unname(prior_inter)))$y)
  actual_name <-list(param_name=c("p","r_tiller","a","k_mound","r","d","k_dead"), 
  real_name=c("p, bulk density (g cm^-3)","r_tillers, tiller radius (cm)","a, root mass per tiller (g #^-1 yr^-1)","k_mound, mound decomposition rate (yr^-1)","r, per-capita growth rate (yr^-1)","d, per_capital death rate (yr^-1)","k_dead, dead tiller decomposition rate (yr^-1)"))
  
  param.prior = runif(1000, param.min[plot_index], param.max[plot_index])
  fraction_iqr<-((quantile(param.keep[,plot_index],0.75)-quantile(param.keep[,plot_index],0.25))/(quantile(param.prior,0.75)-quantile(param.prior,0.25)))*100
  
  plot(density(x = param.keep[,plot_index]),xlim=c(x_min,x_max),ylim=c(0,y_max),lwd=2,xlab=actual_name$real_name[which(names(param.keep)[plot_index]==actual_name$param_name)],font=2,font.lab=2,cex.lab=1.5,main=paste0("% of prior IQR =",round(fraction_iqr,2)))
  points(density(x = unlist(unname(prior_inter))),type="l",col="red",lwd=2)
  
  rm(x_max,x_min,y_max,actual_name,param.prior,fraction_iqr)
  
  }else{
    
  x_max<-param.max[plot_index]
  x_min<-param.min[plot_index]
  y_max<-max(density(x = param.keep[,plot_index])$y)
  actual_name <-list(param_name=c("p","r_tiller","a","k_mound","r","d","k_dead"), 
                     real_name=c("p, bulk density (g cm^-3)","r_tillers, tiller radius (cm)","a, root mass per tiller (g #^-1 yr^-1)","k_mound, mound decomposition rate (yr^-1)","r, per-capita growth rate (yr^-1)","d, per_capital death rate (yr^-1)","k_dead, dead tiller decomposition rate (yr^-1)"))
  param.prior = runif(1000, param.min[plot_index], param.max[plot_index])
  fraction_iqr<-((quantile(param.keep[,plot_index],0.75)-quantile(param.keep[,plot_index],0.25))/(quantile(param.prior,0.75)-quantile(param.prior,0.25)))*100
    
  plot(density(x = param.keep[,plot_index]),xlim=c(x_min,x_max),ylim=c(0,y_max),lwd=2,xlab=actual_name$real_name[which(names(param.keep)[plot_index]==actual_name$param_name)],font=2,font.lab=2,cex.lab=1.5,main=paste0("% of prior IQR =",round(fraction_iqr,2)))
  
  rm(x_max,x_min,y_max,actual_name,param.prior,fraction_iqr)
  }
}
rm(plot_index,prior_data)

############################
#####run the model once#####
############################
# library(parallel);library(doParallel);library(foreach);
# registerDoParallel(detectCores())
# par_runs<-foreach(par_loop=1:draws,.combine=rbind)%dopar%
# {
# out = data.frame(solvemodel(param.keep[par_loop,],times=seq(0,250,0.25)))
# out2 = ((seq(0,2500,1))/(1.6*pi*rep(unname(unlist(param.keep[par_loop,1])))))^(1/3)
# o1<-cbind(param.keep[par_loop,],t(out$M),t(out$N_alive),t(out$N_dead),t(out$K),t(out$r_tussock),t(out2))
# colnames(o1)<-c(colnames(param.keep[par_loop,]),paste0("M",out$time),paste0("N_alive",out$time),paste0("N_dead",out$time),paste0("K",out$time),paste0("r_tussock",out$time),paste0("a_m",seq(0,2500,1)))
# return(o1)
# rm(o1,out,out2)
# }
# stopImplicitCluster()
# saveRDS(par_runs,"Step2_run_replicate.rds")
par_runs<-readRDS("r_outputs/Step2_run.rds")

s.global3 <- par_runs[,c(1:7,8:1008)] 
s.global4 <- par_runs[,c(1:7,1009:2009)]
s.global5 <- par_runs[,c(1:7,2010:3010)]
s.global2 <- par_runs[,c(1:7,3011:4011)]
s.global <- par_runs[,c(1:7,4012:5012)]
s.allometry<- par_runs[,c(1:7,5013:7513)]
rm(par_runs)

s.global<-s.global[complete.cases(s.global),-c(1:7)]
s.global2<-s.global2[complete.cases(s.global2),-c(1:7)]
s.global3<-s.global3[complete.cases(s.global3),-c(1:7)]
s.global4<-s.global4[complete.cases(s.global4),-c(1:7)]
s.global5<-s.global5[complete.cases(s.global5),-c(1:7)]
s.allometry<-s.allometry[complete.cases(s.allometry),-c(1:7)]

s.global.summ = data.frame(x=seq(0,250,0.25),q05=unname(apply(s.global,2,function(x){quantile(x,0.05)})),q95=unname(apply(s.global,2,function(x){quantile(x,0.95)})))
s.global.summ2 = data.frame(x=seq(0,250,0.25),q05=unname(apply(s.global2,2,function(x){quantile(x,0.05)})),q95=unname(apply(s.global2,2,function(x){quantile(x,0.95)})))
s.global.summ3 = data.frame(x=seq(0,250,0.25),q05=unname(apply(s.global3,2,function(x){quantile(x,0.05)})),q95=unname(apply(s.global3,2,function(x){quantile(x,0.95)})))
s.global.summ4 = data.frame(x=seq(0,250,0.25),q05=unname(apply(s.global4,2,function(x){quantile(x,0.05)})),q95=unname(apply(s.global4,2,function(x){quantile(x,0.95)})))
s.global.summ5 = data.frame(x=seq(0,250,0.25),q05=unname(apply(s.global5,2,function(x){quantile(x,0.05)})),q95=unname(apply(s.global5,2,function(x){quantile(x,0.95)})))
s.allometry.summ = data.frame(x=seq(0,2500,1),q05=unname(apply(s.allometry,2,function(x){quantile(x,0.05)})),q95=unname(apply(s.allometry,2,function(x){quantile(x,0.95)})))

model_run<-data.frame(solvemodel(params=param.best,times = seq(0,250,0.25)))
model_allometry<-((seq(0,2500,1))/(1.6*pi*rep(unname(unlist(param.best[1])))))^(1/3)
  
#fig 7a
par(mfrow=c(2,2))
hist(data.compare1$r,xlab="Tussock radius (cm)",font.axis=2,font.lab=2,cex.lab=1.5,main="",freq=F,ylim=c(0,.25),col=rgb(0,0,1,0.25))
hist(s.global$r_tussock250,add=T,freq=F,col=rgb(1,0,0,0.25),breaks=seq(0,100,2))
points(density(data.compare1$r),type="l",lwd=2,col=rgb(0,0,1,0.75))
points(density(s.global$r_tussock250),col=rgb(1,0,0,0.75),type="l",lwd=2)
1-(length(data.compare1$r[data.compare1$r>max(s.global$r_tussock250)|data.compare1$r<min(s.global$r_tussock250)])/length(data.compare1$r))

#fig 7b
pluck<-readRDS("pluck.rds")
plot(NA,NA,font=2,font.lab=2,cex.lab=1.5,xlab="Mass (kg)",ylab="Tussock radius (cm)",ylim=c(0,max(s.allometry.summ$q95)),xlim=c(0,max(pluck$M)/1000),type="l",lwd=2)
polygon(c(rev(s.allometry.summ$x/1000),s.allometry.summ$x/1000),c(rev(s.allometry.summ$q95),s.allometry.summ$q05),col=rgb(1,0,0,0.15))
points(seq(0,2500,1)/1000,model_allometry,type="l",lwd=3,col="red")
points(pluck$M/1000,pluck$r,pch=15,col="blue")
1-sum((pluck$r-(((pluck$M)/(1.5*pi*rep(unname(unlist(param.best[1])))))^(1/3)))^2)/sum((pluck$r-mean(pluck$r))^2)
arrows(pluck$M/1000,pluck$r+1.014479469,pluck$M/1000,pluck$r-1.014479469,col="blue",code=3,angle=90,length=0.05)

#fig 7c
r_all<-round(unname(unlist(s.global)),0)
n_all<-unname(unlist(s.global4))
q1<-unname(unlist(lapply(seq(0,20,1),FUN=function(x){quantile(n_all[r_all==x],0.05)})))
q2<-unname(unlist(lapply(seq(0,20,1),FUN=function(x){quantile(n_all[r_all==x],0.95)})))
qm<-unname(unlist(lapply(seq(0,20,1),FUN=function(x){mean(n_all[r_all==x])})))
qc<-unname(unlist(lapply(round(pluck$r,0),FUN=function(x){mean(n_all[r_all==x])})))
1-sum((pluck$N_alive-qc)^2)/sum((pluck$N_alive-mean(pluck$N_alive))^2)
x<-seq(0,20,1)[!is.na(q1)]
q1<-q1[!is.na(q1)]
q2<-q2[!is.na(q2)]
qm<-qm[!is.na(qm)]
par(mfrow=c(1,1))
plot(NA,NA,font=2,font.lab=2,cex.lab=1.5,xlab="Tussock radius (cm)",ylab="# alive",xlim=c(0,max(pluck$r)),ylim=c(0,1000),type="l",lwd=2)
points(pluck$r,pluck$N_alive,pch=15,col="blue")
polygon(c(x,rev(x)),c(q1,rev(q2)),col=rgb(1,0,0,0.15))
points(x,qm,pch=15,col="red",lwd=2,type="l")
arrows(pluck$r,pluck$N_alive+(mean(sigma.obs1$N_alive,na.rm=T)*pluck$r),pluck$r,pluck$N_alive-(mean(sigma.obs1$N_alive,na.rm=T)*pluck$r),col="blue",code=3,angle=90,length=0.05)

########################
#turn off decomposition#
########################
param.best6<-param.best
param.best7<-param.best
param.best6[c(2)]<-0
param.best7[c(2)]<-param.best[c(2)]*3

model_run<-data.frame(solvemodel(params=param.best,times=seq(0,250,0.1)))
model_run6<-data.frame(solvemodel(params=param.best6,times=seq(0,250,0.1)))
model_run7<-data.frame(solvemodel(params=param.best7,times=seq(0,250,0.1)))

#fig 4
plot(model_run$time,model_run$r_tussock,xlab="Time (years)", ylab="Tussock radius (cm)",ylim=c(0,15),type="l",lwd=2,font.axis=2,font.lab=2,cex.lab=1.5,xlim=c(0,250))
points(model_run6$time,model_run6$r_tussock,type="l",lty=2,lwd=2)
points(model_run7$time,model_run7$r_tussock,type="l",lty=3,lwd=2)

#fig 5a
plot(model_run$time,model_run$M*param.best[4],type="l",ylim=c(0,20),lwd=2,font.axis=2,font.lab=2,cex.lab=1.5)
points(model_run$time,model_run$N_alive*param.best[3],type="l",lwd=2,lty=2)
polygon(c(model_run$time,rev(model_run$time)),c(model_run$M*param.best[4],rev(model_run$N_alive*param.best[3])),col=grey(0.5,0.5),border=NA)
par(new=T)
plot(model_run$time,(model_run$N_alive+model_run$N_dead)/(model_run$K+model_run$N_dead),ylim=c(0,1),lty=3,yaxt="n",typ="l",lwd=2,font.axis=2,font.lab=2,cex.lab=1.5,ylab="")
axis(4,font.axis=2,font.lab=2,cex.lab=1.5)
par(new=F)

#fig 5b
plot(model_run$time,model_run$N_alive/(model_run$N_alive+model_run$N_dead),xlab="Time (years)", ylab="Proportion of total",ylim=c(0,1),type="l",lwd=2,font.axis=2,font.lab=2,cex.lab=1.5,xlim=c(0,250))
points(model_run$time,model_run$N_dead/(model_run$N_alive+model_run$N_dead),ylim=c(0,1),type="l",lwd=2,font.axis=2,font.lab=2,cex.lab=1.5,xlim=c(0,250),lty=2)

#fig 8
out2<-s.global$r_tussock250
par(mfrow=c(3,3))
for(i in 1:7){
plot(param.keep[order(param.keep[,i]),i],predict(loess(out2[order(param.keep[,i])]~param.keep[order(param.keep[,i]),i])), col='black', lwd=2,type="l",ylab="Tussock radius (cm)",xlab=names(param.keep)[i],font.lab=2,font.axis=2,cex.lab=1.5,cex.axis=1.2,ylim=c(1,12))
plx<-predict(loess(out2[order(param.keep[,i])]~param.keep[order(param.keep[,i]),i]), se=T)
lines(param.keep[order(param.keep[,i]),i],plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(param.keep[order(param.keep[,i]),i],plx$fit + qt(0.975,plx$df)*plx$se, lty=2)
print(paste(colnames(param.keep)[i],plx$fit[length(plx$fit)]-plx$fit[1]))
rm(plx)
}

#fig 8 derived parameter
newb<-param.keep$d + (param.keep$r*(1-(s.global4[,1001]/s.global2[,1001])))
plot(newb[order(newb)],predict(loess(out2[order(newb)]~newb[order(newb)])), col='black', lwd=2,type="l",ylab="Tussock radius (cm)",xlab="new sensitivity",font.lab=2,font.axis=2,cex.lab=1.5,cex.axis=1.2,ylim=c(1,12))
plx<-predict(loess(out2[order(newb)]~newb[order(newb)]), se=T)
lines(newb[order(newb)],plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(newb[order(newb)],plx$fit + qt(0.975,plx$df)*plx$se, lty=2)
print(paste(colnames(param.keep)[i],plx$fit[length(plx$fit)]-plx$fit[1]))

#fig 9a
par(mfrow=c(1,1))
plot(pluck$root_mass/pluck$adult,pluck$r,ylim=c(0,18),xlim=c(0,0.16),font.axis=2,font.lab=2,cex.lab=1.5,cex.axis=1.25,pch=16)
o1<-pluck$root_mass/pluck$adult
o2<-pluck$r
summary(lm(o2~o1))
abline(lm(o2~o1),lwd=2,lty=2)

#fig 9b 
o1<-pluck$M/(pi*((pluck$r)^2)*(pluck$height_total-pluck$height_live))
o2<-pluck$r
plot(o1,o2,ylim=c(0,18),font.axis=2,font.lab=2,cex.lab=1.5,cex.axis=1.25,pch=16)
summary(lm(o2~o1))
abline(lm(o2~o1),lwd=2,lty=2)

#fig 9c
plot((pluck$v0+pluck$v1)/pluck$adult,pluck$r,ylim=c(0,18),font.axis=2,font.lab=2,cex.lab=1.5,cex.axis=1.25,pch=16)
o1<-(pluck$v0+pluck$v1)/pluck$adult
o2<-pluck$r
summary(lm(o2~o1))