#load packages
rm(list=ls())
require(deSolve)

setwd("C:/Users/scurasi/Desktop/tussock_models/for_archive/data/")
load(file="r_outputs/Step1_replicate.Rdata")

draws=10000

save_seed<-.Random.seed

#######STEP 2: ESTIMATE PARAMETER UNCERTAINTY
out= data.frame(solvemodel(param.best,times=seq(0,250,0.25))) #creates table of model output

out.compare1 = data.frame(list(r=c(rep(out$r_tussock[nrow(out)],n.time)),M=c(rep(out$M[nrow(out)],n.time)),N_alive=c(rep(out$N_alive[nrow(out)],n.time)),N_dead=c(rep(out$N_dead[nrow(out)],n.time))))[,c(1:4)]
out.compare1$N_dead =c(rep(unname(unlist(param.best[1])),length(na.omit(data.compare1$N_dead))),rep(NA,length(data.compare1$N_dead)-length(na.omit(data.compare1$N_dead))))
out.compare1$M<-c(((pluck$M)/(1.6*pi*rep(unname(unlist(param.best[1])),length(pluck$M))))^(1/3),rep(NA,length(out.compare1$M)-length(pluck$M)))
out.compare1$N_alive<-c(rep(out$N_alive[nrow(out)]/out$r_tussock[nrow(out)],length(pluck$M)),rep(NA,length(out.compare1$M)-length(pluck$M)))

var.jbest = rep(0, D)
error.jbest=matrix(NA, n.time, D) #create data frame to store error calculations; want all to be "0" originally because if there is no data it will remain 0
for (d in 1:D) { #for each data type
  for (m in 1:n.time){ #for each timestep
    if(!is.na(data.compare1[m,d])){ #if there is data at that timestep for that data stream
      error.jbest[m,d]=((data.compare1[m,d] - out.compare1[m,d])/sigma.obs1[m,d])^2 #calculates the error at that timestep for that data stream
    } #end of if statement
  } #end of time step loop
  
  var.jbest[d] = var(error.jbest[!is.na(data.compare1[,d]),d]) #calculate variance of the errors (excludes NAs)
  
} #end of data type loop

#storage matrices for Monte Carlo reps
j = rep(0, D)
param.keep = data.frame(matrix(1, draws, n.param)) #storage for parameter estimate iterations; 
colnames(param.keep) = c(names(param.best))
param.keep[1,]=param.best

#also need to know degrees of freedom for chi square test
n.par = n.param #number of parameters predicted by each data stream
df = rep(0, D)
for (d in 1:D) { #for each data type
  df[d] = sum(!is.na(data.compare1[,d])) - n.par
} #end of data loop

#set initial values
param.est = param.best #set initial values for parameters
reject=0 #reset reject counter
num.accepted = 0 #counter for number of accepted parameters - when this gets to 1000, loop will stop
num.reps = 0 #counter for number of repititions - calculates acceptance rate

time<-as.numeric(as.POSIXct(Sys.time()))

#start loop
repeat { #repeat until desired number of parameter sets are accepted
  
  num.reps=num.reps+1 #add to number of reps counter
  
  repeat{
    repeat{
      for(p in 1:n.param){ #for each parameter
        param.est[p] = param.best[p] + rnorm(1, 0, t*(param.max[p]-param.min[p]))
        parms = as.numeric(param.est) #parameters for model run
        names(parms) = names(params) #fix names
      }
      if(all(param.est>=param.min) & all(param.est<=param.max)){
        break
      }}#end of parameter loop
    
    out=data.frame(solvemodel(params=parms, times=seq(0,250,0.25)))#run the model
    
    if((out$r_tussock[nrow(out)] - out$r_tussock[nrow(out)-50])<0.01 & !any(is.na(out))){
      break
    } #end of if loop
  } #end of repeat
  
  out.compare1 = data.frame(list(r=c(rep(out$r_tussock[nrow(out)],n.time)),M=c(rep(out$M[nrow(out)],n.time)),N_alive=c(rep(out$N_alive[nrow(out)],n.time)),N_dead=c(rep(out$N_dead[nrow(out)],n.time))))[,c(1:4)]
  out.compare1$N_dead =c(rep(unname(unlist(parms[1])),length(na.omit(data.compare1$N_dead))),rep(NA,length(data.compare1$N_dead)-length(na.omit(data.compare1$N_dead))))
  out.compare1$M<-c(((pluck$M)/(1.6*pi*rep(unname(unlist(parms[1])),length(pluck$M))))^(1/3),rep(NA,length(out.compare1$M)-length(pluck$M)))
  out.compare1$N_alive<-c(rep(out$N_alive[nrow(out)]/out$r_tussock[nrow(out)],length(pluck$M)),rep(NA,length(out.compare1$M)-length(pluck$M)))
  
  #determine if parameter set is accepted or rejected
  error = matrix(NA, n.time, D)
  var.error=rep(0,D)  
  
  error=((data.compare1 - out.compare1)/sigma.obs1)^2 #calculates the error at that timestep for that data stream
  
  for (d in 1:D) { #for each data type
    
    var.error[d] = var(error[!is.na(data.compare1[,d]),d]) #calculate variance of the errors (excludes NAs)
    
    for (m in 1:n.time){ #for each timestep
      error[m,d] = (error[m,d]*sqrt(var.jbest[d]))/sqrt(var.error[d]) #variance normalization
    } #end of time step loop  
    
    j[d] = sum(error[!is.na(data.compare1[,d]),d]) #calculate cost function for each data stream after variance normalizaiton
  } #end of data type loop
  
  #chi-square test
  accept = rep (0, D) #vector to keep track of if each j has been accepted or rejected; 1=accept, 0=reject
  for (d in 1:D) { #for each data type  
    
    if(j[d]-j.best[d] <= qchisq(0.9, df[d])) { #conduct chi square test
      accept[d] = 1} #if accepted, change value in accept vector to 1
  } #end of data type loop
  
  d.accept = sum(accept) #calculate the number of j's accepted
  
  if(d.accept==D) { #if all j's are accepted
    num.accepted = num.accepted+1 #add to number of parameter sets accepted
    param.keep[num.accepted,]=param.est #store the parameter set in the storage dataframe
  } #end of if loop
  if(d.accept<D) { #if any j's rejected
    reject = reject+1 #reject parameter set
  } #end of if loop

  acceptance = 1 - (reject / num.reps) #calculate proportion of accepted iterations
  
  #print information during the run
  if(num.accepted == round(draws*0.001,0) | num.accepted== round(draws*0.01,0) |  num.accepted %% round(draws*0.1,0) == 0){
    print(paste0("itterations = ",num.reps," / rejected = ", reject," / prop accepted = ",acceptance, " / elapsed (hrs) = ",round(((as.numeric(as.POSIXct(Sys.time()))-time)/60)/60,2)," / remaining (hrs) = ",round((((draws-num.accepted)*((as.numeric(as.POSIXct(Sys.time()))-time)/num.accepted))/60)/60,2)))
  }
  
  if (num.accepted==draws) { #if you have accepted the number of parameter sets you want (i.e., 1000)
    break  #break repeat loop
  } 
  
} #end of repeat

save.image(file="r_outputs/Step2_replicate.Rdata")
