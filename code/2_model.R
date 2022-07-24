#code for tussock model difeqs

#--------------setup---------------#
rm(list=ls())
setwd("C:/Users/scurasi/Desktop/tussock_models/for_archive/data/")
if(!"deSolve" %in% installed.packages()[,"Package"]) install.packages("deSolve") #requires deSolve be installed

#--------------model function ---------------#
solvemodel<- function(params, state = c(M=1,N_alive=1,N_dead=0), times = seq(0,250,0.25)) {
  
  model<-function(t,state,params)
  { 
    with(as.list(c(state, params)),{
      
      #model does not include forcings, or constants
      
      #calculate K and r
      K=(((pi*(((M)/(1.6*pi*p))^(2/3)))/((r_tiller^2)*(sqrt(12)))) - N_dead)
      r_tussock=((M)/(1.6*pi*p))^(1/3)
        
      #differential equations
      dM = a*N_alive - M*k_mound
      dN_alive= (r)*N_alive*(1-(N_alive/K))
      dN_dead = d*N_alive - N_dead*k_dead
      
      #output
      list(c(dM,dN_alive,dN_dead),c(K=K,r_tussock=r_tussock))
    })
  }
  return(deSolve::ode(y=state,times=times,func=model,parms = params, method="rk4")) #solve with rk4
}

#--------------test---------------#
#use mean of bounds
params<-data.frame(list(p=c(0.2563975,0.04665624),r_tiller=c(0.25,0.0992),a=c(0.18860000,0.00327027),k_mound=c(0.319287710,0.006087163),r=c(0.5,0.001),d=c(0.5,0.001),k_dead=c(0.319287710,0.006087163)))
params<-(params[1,]+params[2,])/2

#test/plot
par(mfrow=c(3,3))
out = data.frame(solvemodel(params=params)) #creates table of model output
plot(out$time,out$N_alive,ylab="Alive #",xlab="time yrs.",type="l",lwd=2)
plot(out$time,out$N_dead,ylab="Dead #",xlab="time yrs.",type="l",lwd=2)
plot(out$time,out$M,ylab="Mass g",xlab="time yrs.",type="l",lwd=2)
plot(out$time,out$r,ylab="radius cm",xlab="time yrs.",type="l",lwd=2)
plot(out$time,out$K,ylab="K",xlab="time yrs.",type="l",lwd=2)

#--------------save the function---------------#
saveRDS(solvemodel,"r_outputs/model.rds")
