#this code is adapted from that used in Wright and Rocha 2018
#See: https://doi.org/10.1016/j.ecolmodel.2018.05.017 and
#https://github.com/kkremers/CCaN_MS2016.git
########setup###########

rm(list=ls())
setwd("C:/Users/scurasi/Desktop/tussock_models/for_archive/data/")

#load data
diam<-readRDS("r_outputs/diams.rds")
diam_err<-diam
diam_err$r<-sd(diam$r)

pluck<-readRDS("r_outputs/pluck.rds")
diam$M<-c(pluck$r,rep(NA,length(diam$r)-length(pluck$r)))
diam_err$M<-c(rep(1.014479469,length(pluck$r)),rep(NA,length(diam$M)-length(pluck$r)))

diam$N<-c(pluck$N_alive/pluck$r,rep(NA,length(diam$r)-length(pluck$r)))
diam_err$N<-c(rep(sd(diam$N,na.rm = T),length(pluck$M)),rep(NA,length(diam$r)-length(pluck$r)))

diam<-diam[,-c(1,3)]
diam_err<-diam_err[,-c(1,3)]

#this code was written with the option for a 4th data
#stream that is unused and not fully implemented see "D" below
diam$N_dead=rep(0,length(diam$r))
diam_err$N_dead=rep(0,length(diam$r))

#########
data.compare1<-data.frame(list(r=diam$r,M=diam$M,N_alive=diam$N,N_dead=diam$N_dead))
sigma.obs1<-data.frame(list(r=diam_err$r,M=diam_err$M,N_alive=diam_err$N,N_dead=diam_err$N_dead))
#########

params<-data.frame(list(p=c(0.2563975,0.04665624),r_tiller=c(0.25,0.0992),a=c(0.18860000,0.00327027),k_mound=c(0.319287710,0.006087163),r=c(0.5,0.001),d=c(0.5,0.001),k_dead=c(0.319287710,0.006087163)))
param.max=c(unname(unlist(params[1,])))
param.min=c(unname(unlist(params[2,])))
params<-(params[1,]+params[2,])/2

#values that control the run
n.param = length(params)#number of parameters to estimate
M = 50000 #number of iterations
D = 3 #number of data types being assimilated 
n.time = nrow(data.compare1) #numer of rows in the data (the original code compares data in time, the matching here is more complex)
reject=0 #reset reject counter
t=0.5

#storage matrices
J = rep(1E100, M) #storage vector for cost function output
j = matrix(1E100, M, D) #to store error calculations for this iteration
all.draws = data.frame(matrix(1, M, n.param)) #storage for all parameter estimate iterations;
colnames(all.draws) = c(names(params))
param.est = data.frame(matrix(1, M, n.param)) #storage for accepted parameter estimate iterations;
param.est[1,]=params #change first row to current guess
all.draws[1,]=params #change first row to current guess
colnames(param.est) = c(names(params))

#load and run the model for the first time
solvemodel<-readRDS("r_outputs/model.rds")
out=data.frame(solvemodel(params=params, times=seq(0,250,0.25)))

out.compare1 = data.frame(list(r=c(rep(out$r_tussock[nrow(out)],n.time)),M=c(rep(out$M[nrow(out)],n.time)),N_alive=c(rep(out$N_alive[nrow(out)],n.time)),N_dead=c(rep(out$N_dead[nrow(out)],n.time))))[,c(1:4)]
out.compare1$N_dead =c(rep(unname(unlist(params[1])),length(na.omit(data.compare1$N_dead))),rep(NA,length(data.compare1$N_dead)-length(na.omit(data.compare1$N_dead))))
out.compare1$M<-c(((pluck$M)/(1.6*pi*rep(unname(unlist(params[1])),length(pluck$M))))^(1/3),rep(NA,length(out.compare1$M)-length(pluck$M)))
out.compare1$N_alive<-c(rep(out$N_alive[nrow(out)]/out$r_tussock[nrow(out)],length(pluck$M)),rep(NA,length(out.compare1$M)-length(pluck$M)))

error.time=matrix(0, n.time, D) #create data frame to store error calculations; want all to be "0" originally because if there is no data it will remain 0
for (d in 1:D) { #for each data type
  for (m in 1:n.time){ #for each step
    if(!is.na(data.compare1[m,d])){ #if there is data at that step for that data stream
      error.time[m,d]=((data.compare1[m,d] - out.compare1[m,d])/sigma.obs1[m,d])^2 #calculates the error at that step for that data stream
    } #end of if statement
    #if there was no data at that step, the error will remain "0" so that it will not impact the sum calculation in the next step
  } #end of step loop
  
  j[1,d] = sum(error.time[,d]) #calculate cost function for each data stream
  
} #end of data type loop

J[1] = prod(j[1,]) #calculate aggregate cost function

time<-as.numeric(as.POSIXct(Sys.time()))

set.seed(1)

#start the run
for (i in 2:M) {
  repeat{
    for(p in 1:n.param){ #for each parameter
      param.est[i,p] = param.est[i-1,p] + rnorm(1, 0, t*(param.max[p]-param.min[p]))
      all.draws[i,p] = param.est[i,p]
      parms = as.numeric(param.est[i,]) #parameters for model run
      names(parms) = names(params) #fix names
      }
      if(all(param.est[i,]>=param.min) & all(param.est[i,]<=param.max)){
      break
      }}#end of parameter loop
    
    out=data.frame(solvemodel(params=parms, times=seq(0,250,0.25)))#run the model
    
    if(any(is.na(out))){
      reject = reject+1 #reject parameter set
      param.est[i,] = param.est[i-1,] #set current parameter set to previous parameter set
      J[i] = J[i-1] #set current J to previous J (the minimum J so far)
      j[i,] = j[i-1,]
      #break
    }else if((out$r_tussock[nrow(out)] - out$r_tussock[nrow(out)-50])>=0.01){ 
      reject = reject+1 #reject parameter set
      param.est[i,] = param.est[i-1,] #set current parameter set to previous parameter set
      J[i] = J[i-1] #set current J to previous J (the minimum J so far)
      j[i,] = j[i-1,]
      #break
    }else if((out$r_tussock[nrow(out)] - out$r_tussock[nrow(out)-50])<0.01 & !any(is.na(out))){#end of if loop
  
  out.compare1 = data.frame(list(r=c(rep(out$r_tussock[nrow(out)],n.time)),M=c(rep(out$M[nrow(out)],n.time)),N_alive=c(rep(out$N_alive[nrow(out)],n.time)),N_dead=c(rep(out$N_dead[nrow(out)],n.time))))[,c(1:4)]
  out.compare1$N_dead =c(rep(unname(unlist(parms[1])),length(na.omit(data.compare1$N_dead))),rep(NA,length(data.compare1$N_dead)-length(na.omit(data.compare1$N_dead))))
  out.compare1$M<-c(((pluck$M)/(1.6*pi*rep(unname(unlist(parms[1])),length(pluck$M))))^(1/3),rep(NA,length(out.compare1$M)-length(pluck$M)))
  out.compare1$N_alive<-c(rep(out$N_alive[nrow(out)]/out$r_tussock[nrow(out)],length(pluck$M)),rep(NA,length(out.compare1$M)-length(pluck$M)))
  
  #determine if parameter set is accepted or rejected
  error.time=((data.compare1 - out.compare1)/sigma.obs1)^2 #calculates the error at that timestep for that data stream
  
  j[i,1:D] = colSums(error.time[,1:D],na.rm = T)
  
  J[i] = prod(j[i,]) #calculate aggregate cost function
  
  diff = J[i-1]/J[i] #calculate the acceptance ratio
  
  if(is.na(diff)){
    reject = reject+1 #reject parameter set
    param.est[i,] = param.est[i-1,] #set current parameter set to previous parameter set
    J[i] = J[i-1] #set current J to previous J
    j[i,] = j[i-1,]
  } else { 
    
    if(diff<1){ #if difference is < 1 (or if the current J is greater than the previous J)
      
      u=runif(1, 0, 1) #draw random number between 0 and 1
      
      if(u>=diff){   
        reject = reject+1 #reject parameter set
        param.est[i,] = param.est[i-1,] #set current parameter set to previous parameter set
        J[i] = J[i-1] #set current J to previous J (the minimum J so far)
        j[i,] = j[i-1,]
      } #end of if loop
    } #end of if loop
    
  } #end of else loop
  }
  
  acceptance = 1 - (reject / i) #calculate proportion of accepted iterations
  
  if(acceptance>0.30){
    t = 1.01*t
  }
  
  if(acceptance<0.25){
    t = 0.99*t
  }
  
  if(t>0.5){
    t=0.5
  }
  
  if(i == round(M*0.001,0) | i== round(M*0.01,0) |  i %% round(M*0.1,0) == 0){
    print(paste0("itterations = ",i," / rejected = ", reject," / t = ",t," / prop accepted = ",acceptance, " / elapsed (hrs) = ",round(((as.numeric(as.POSIXct(Sys.time()))-time)/60)/60,2)," / remaining (hrs) = ",round((((M-i)*((as.numeric(as.POSIXct(Sys.time()))-time)/i))/60)/60,2)))
  }
  
} #end of exploration

steps=seq(1:i) #create a vector that represents the number of steps or iterations run
J1=data.frame(steps, J[1:i]) #create a dataframe that has "steps" as the first column and "J" as the second column
step.best = J1[which.min(J1[,2]),1] #determine which step has the minimum value of J and store as "step.best"
param.best = as.numeric(param.est[step.best,]) #store that parameter set as param.best
names(param.best) = c(names(params)) #change the names to match params
j.best = j[step.best,] #pull out the minimum j

par(mfrow=c(3,3))
for(plots in 1:n.param){
  plot(all.draws[1:i,plots],main=paste0("param=",names(param.est)[plots]," / est=",round(param.est[step.best,plots],3)) ) 
  lines(param.est[1:i,plots], col="red")
}

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.75, paste0("itterations = ", i),cex = 1.6, col = "black")
text(x = 0.5, y = 0.25, paste0("acceptance rate = ", acceptance),cex = 1.6, col = "black")
rm(plots)

model_run<-data.frame(solvemodel(params=param.best,times=seq(0,250,0.25)))

hist(diam$r)
abline(v=model_run$r_tussock[nrow(model_run)],col="red")

#plot(model_run$time,model_run$r_tussock)

save.image(file="r_outputs/Step1_replicate.Rdata")