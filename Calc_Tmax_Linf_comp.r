library(freeR)
library(ggplot2)
library(rfishbase)
library(truncnorm)
library(fishmethods)
library(FSAsim)
library(FSA)
library(EnvStats)
library(reshape2)
library(viridis)
library(dplyr)
### Functions ###
VBGF<-function(Linf, k, t0, ages){ 
  Lts<-Linf * (1 - exp(-k * (ages - t0)))
  Lts
}

rand.VBGF<-function(Linf, k, t0, ages,CV){ 
  Lts<-Linf * (1 - exp(-k * (ages - t0)))
  rnorm(length(Lts),Lts,CV*Lts) 
}

vbgf.inflect<-function(Linf,k,t0)
{
  (1/k)*(exp(3*(1-(t0/Linf)^(1/3))))
}

VBGF.age<-function(Linf,k,t0,lt){ 
  t0 - (log(1 - (lt / Linf)) / k) 
} 

trunc.rand.VBGF<-function(Linf, k, t0, ages,CV,maxlt){ 
  Lts<-Linf * (1 - exp(-k * (ages - t0)))
  rlnormTrunc(length(Lts), Lts, CV*Lts, min = 0,max = maxlt)
  #rnorm(length(Lts),Lts,CV*Lts) 
}

Comp.Lmax.Linf<-function(n=100,Nsim=500,Linf.in,K.in,t0.in,CV.in,Amax.in,lowage,maxage,maxlt)
{
  Linf_Lmax_ratios<-data.frame(Lmax_TLinf=NA,Lmax_estLinf=NA,estLinf_TLinf=NA,Linf.true=NA,Linf_est=NA,Lmax=NA, Amin=NA, Amax=NA,Kest=NA,K_jit=NA,toest=NA,t0_jit=NA)
  for(i in 1:Nsim)
  {
    #Create data
    rand.ages<-rlnormTrunc(n, meanlog = log(Amax.in), sdlog = 1, min = lowage,max = maxage)
    rand.lts<-rand.VBGF(Linf.in,K.in,t0.in,rand.ages,CV.in)
    age.lt.dat<-data.frame(Ages=round(rand.ages,0),Lengths=round(rand.lts,0))
    #Filter max length
    age.lt.dat<-age.lt.dat[age.lt.dat$Lengths<maxlt,]
    
    
    
    #  Age_Lts<-data.frame(Ages=rand.ages,Lts=rand.lts)
    #  plot(rand.ages,rand.lts)
    #Fit VBGF curve
    k_t0_jits<-c(jitter(K.in,amount=0.05),jitter(t0.in,amount=0.05))
    growth.out<-try(nls(age.lt.dat$Lengths~Linf * (1 - exp(-k * (rand.ages - t0))),
                        data=age.lt.dat,start=list(Linf=max(rand.lts),k=k_t0_jits[1],t0=0)))
    #growth.out<-try(growth(size=rand.lts,age=round(rand.ages,0),Sinf=max(rand.lts),K=k_t0_jits[1],t0=0,
    #                       graph=FALSE,control=list(maxiter=10000,tol=1e-5))) #k_t0_jits[2]
    #vbmod <- Lts ~ Linf * (1 - exp(-K * (Ages - t0)))
    #vbgf.parms<-list(Linf=20,K=0.1,t0=0)
    #vbgf.fit<-nls(vbmod,data=Age_Lts,start=vbgf.parms)
    if(class(growth.out)=="try-error")
    {Linf_Lmax_ratios[i,]=NA}
    
    if(class(growth.out)=="nls")
    {
      growth.out.parms<-growth.out$m$getPars()
      print(growth.out.parms)
      Linf_Lmax_ratios[i,1]<-max(age.lt.dat$Lengths)/Linf.in
      Linf_Lmax_ratios[i,2]<-max(age.lt.dat$Lengths)/growth.out.parms[1]
      Linf_Lmax_ratios[i,3]<-growth.out.parms[1]/Linf.in
      Linf_Lmax_ratios[i,4]<-Linf.in
      Linf_Lmax_ratios[i,5]<-growth.out.parms[1]
      Linf_Lmax_ratios[i,6]<-max(age.lt.dat$Lengths)
      Linf_Lmax_ratios[i,7]<-min(age.lt.dat$Ages)
      Linf_Lmax_ratios[i,8]<-max(age.lt.dat$Ages)
      Linf_Lmax_ratios[i,9]<-growth.out.parms[2]
      Linf_Lmax_ratios[i,10]<-k_t0_jits[1]
      Linf_Lmax_ratios[i,11]<-growth.out.parms[3]
      Linf_Lmax_ratios[i,12]<-k_t0_jits[2]
      
      #Linf_Lmax_ratios[i,4]<-(growth.out.parms[1]/max(rand.lts))/(growth.out.parms[1]/Linf.in)    
    }
  }
  return(Linf_Lmax_ratios)
}

Linf_Lmax_plots<-function(Linf_Lmax_ratios)
{
  ratios<-Linf_Lmax_ratios[,1:3]
  gg_Linf_Lmax_ratios<-melt(ratios)
  colnames(gg_Linf_Lmax_ratios)<-c("Label","Ratio")
  ggplot(gg_Linf_Lmax_ratios,aes(Ratio,color=Label))+ 
    geom_density(lwd=2)+
    geom_vline(xintercept=c(1,median(Linf_Lmax_ratios$estLinf_TLinf,na.rm=TRUE)), color=c("red","black"))+
    xlim(0,2)
}

Linf_plots<-function(Linf_Lmax_ratios)
{
  ratios<-Linf_Lmax_ratios[,3]
  #gg_Linf_Lmax_ratios<-melt(ratios)
  #colnames(gg_Linf_Lmax_ratios)<-c("Label","Ratio")
  ggplot(as.data.frame(Linf_Lmax_ratios[,3]),aes(Linf_Lmax_ratios[,3]))+ 
    geom_density(lwd=2)+
    geom_vline(xintercept=1, color=c("red"))+
    xlim(0,2)+
    xlab("Ratio of Est. Linf to True Linf")
}

#############

#Upload VBGF and tmax values from FishBase
#These are only the species confirmed to have growth estiamtes 
fb.spp.parms<-read.csv("C:/Users/Jason.Cope/Documents/Current Action/Publications/tmax vs Linf/Fishbase_spp_names.csv")
#This command would pull all fish in FishBase
#fb.spp.parms<-all_fish()

#Use species names from FishBase that have growth parameters and extract Linf, K, and tmax values using FishLife. 
#Add t0 = 0 values too.
VBGF.fb.out<-data.frame(Species=NA,Linf=NA,K=NA,tmax=NA,t0=0)
for(i in 1:length(fb.spp.parms$Species))
{
    #try(fishlife(fb.spp.parms$Species[i]),silent=TRUE)
    VBGF.fb.out[i,]<-try(c(fb.spp.parms$Species[i],fishlife(fb.spp.parms$Species[i])[c(1,2,4)],0))
}

#Calculate the Lt value at tmax
VBGF.fb.out$Linf_tmax<-as.numeric(VBGF.fb.out$Linf)*(1-exp(-as.numeric(VBGF.fb.out$K)*(as.numeric(VBGF.fb.out$tmax)-as.numeric(VBGF.fb.out$t0))))
#Calculate the ratio of Linf to L_tmax
VBGF.fb.out$Linf_tmax_Linf<-VBGF.fb.out$Linf_tmax/as.numeric(VBGF.fb.out$Linf)
#Remove the NAs
VBGF.fb.out_noNA<-VBGF.fb.out[ !is.na(VBGF.fb.out$Linf_tmax_Linf),]
#Convert the rest of the inputs to numeric
VBGF.fb.out_noNA$Linf<-as.numeric(VBGF.fb.out_noNA$Linf)
VBGF.fb.out_noNA$K<-as.numeric(VBGF.fb.out_noNA$K)
VBGF.fb.out_noNA$tmax<-as.numeric(VBGF.fb.out_noNA$tmax)
#Save as R object
save(VBGF.fb.out_noNA,file="C:/Users/Jason.Cope/Documents/Current Action/Publications/tmax vs Linf/VBGF_fb_out_noNA.rds")
load("C:/Users/Jason.Cope/Documents/Github/Linf_Lmax_Amax/VBGF_fb_out_noNA.rds")

#Extract Lmax values from FishBase
Lmax.spp.fb<-popchar(fb.spp.parms$Species)
Lmax.spp.fb$Lmax<-as.numeric(Lmax.spp.fb$Lmax)
Lmax.spp.fb<-Lmax.spp.fb[!is.na(Lmax.spp.fb$Lmax),]
Lmax.spp.max<-aggregate(Lmax~Species,Lmax.spp.fb,FUN=max)
Lmax.spp.mean<-aggregate(Lmax~Species,Lmax.spp.fb,FUN=mean)
names(Lmax.spp.mean)[2]<-"Lmax.mean"
save(Lmax.spp.max,file="C:/Users/Jason.Cope/Documents/Current Action/Publications/tmax vs Linf/Lmax_spp_max.rds")
save(Lmax.spp.mean,file="C:/Users/Jason.Cope/Documents/Current Action/Publications/tmax vs Linf/Lmax_spp_mean.rds")

#Merge the vbgf and lmax objects
VBGF_tmax_Lmax<-merge(VBGF.fb.out_noNA,Lmax.spp.mean,by="Species")
VBGF_tmax_Lmax<-merge(VBGF_tmax_Lmax,Lmax.spp.max,by="Species")
VBGF_tmax_Lmax$Lmax_Linf<-VBGF_tmax_Lmax$Lmax/VBGF_tmax_Lmax$Linf
VBGF_tmax_Lmax$Linflect<-VBGF_tmax_Lmax$Linf*(8/27)
save(VBGF_tmax_Lmax,file="C:/Users/Jason.Cope/Documents/Current Action/Publications/tmax vs Linf/VBGF_tmax_Lmax.rds")
load("C:/Users/Jason.Cope/Documents/Current Action/Publications/tmax vs Linf/VBGF_tmax_Lmax.rds")

#Plots to compare different parameter to Linf:L_Amax
#ggplot(VBGF_tmax_Lmax,aes(Linf__Linf_tmax,Linf))+
#    geom_point()+
#    geom_vline(xintercept=1, color="red")

#ggplot(VBGF_tmax_Lmax,aes(Linf_tmax_Linf,K))+
#    geom_point()

#ggplot(VBGF_tmax_Lmax,aes(Linf_tmax_Linf,tmax))+
#    geom_point()

#Compare to Lmax
#ggplot(VBGF_tmax_Lmax,aes(Linf_Lmax,Linf))+
#    geom_point()+
#    geom_vline(xintercept=1, color="red")

#ggplot(VBGF_tmax_Lmax,aes(Linf_Lmax,Linflect))+
#  geom_point()+
#  geom_vline(xintercept=1, color="red")


ggplot(VBGF_tmax_Lmax,aes(Lmax_Linf,Linf_tmax_Linf))+
    geom_point()+
    geom_hline(yintercept=c(0.9,0.95,1),color=c("red","red","blue"))+
    geom_vline(xintercept=c(0.9,1),color=c("red","blue"))+
    xlim(0,2)+
    ylim(0,2)+
    xlab("Lmax:Linf")+
    ylab("Lt@Tmax:Linf")

ggplot(VBGF_tmax_Lmax,aes(Lmax_Linf))+
  geom_density(lwd=2)+
  xlim(0,2.5)+
  geom_vline(xintercept=c(median(VBGF_tmax_Lmax$Lmax_Linf),0.9), color=c("blue","red"),linetype = c("longdash","solid"))+
  xlab("Lmax:Linf")

ggplot(VBGF_tmax_Lmax,aes(Linf_tmax_Linf))+
  geom_density(lwd=2)+
  xlim(0.6,1.1)+
  geom_vline(xintercept=c(median(VBGF_tmax_Lmax$Linf_tmax_Linf),0.9,0.95,1), color=c("blue","red","red","red"),linetype = c("longdash","solid","solid","solid"),lwd=c(0.5,0.5,0.5,1.25))+
  xlab("Lt@Tmax:Linf")


#Calculate ratio categories
#High Linf vs Amax; Linf>Lmax WORST LINF
#High Linf vs Amax; Linf<Lmax BAD Amax
#Linf vs Amax; Linf>Lmax BAD LINF
#Linf vs Amax; Linf<Lmax Good

Linf_Amax_ratio_hi<-0.99
Lmax_Linf_ratio_low<-0.9

#ratio.cat<-c(dim(VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf_tmax_Linf<Linf_Amax_ratio & VBGF_tmax_Lmax$Lmax_Linf<Lmax_Linf_ratio_low])[1],
#dim(VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf__Linf_tmax>Linf_Amax_ratio & VBGF_tmax_Lmax$Linf_Lmax<Linf_Lmax_ratio,])[1],
#dim(VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf__Linf_tmax<Linf_Amax_ratio & VBGF_tmax_Lmax$Linf_Lmax>Linf_Lmax_ratio,])[1],
#dim(VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf__Linf_tmax<Linf_Amax_ratio & VBGF_tmax_Lmax$Linf_Lmax<Linf_Lmax_ratio,])[1])

#Linf to Amax
dim(VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf_tmax_Linf<Linf_Amax_ratio_hi,])[1]/dim(VBGF_tmax_Lmax)[1]
#Linf to Lmax
dim(VBGF_tmax_Lmax[VBGF_tmax_Lmax$Lmax_Linf<Lmax_Linf_ratio_low,])[1]/dim(VBGF_tmax_Lmax)[1]

#Linf to Amax and Lmax
dim(VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf_tmax_Linf<Linf_Amax_ratio_hi&VBGF_tmax_Lmax$Lmax_Linf<Lmax_Linf_ratio_low,])[1]/dim(VBGF_tmax_Lmax)[1]


ratio_Linfs<-ratio.cat/sum(ratio.cat)
names(ratio_Linfs)<-c("Linf_ratio>LAmax_Lmax","Linf>LAmax","Linf>Lmax","Linf<LAmax_Lmax")

ratio_Linfs_combo<-c(sum(ratio.cat[c(1,3)]),ratio.cat[2],ratio.cat[4])/sum(ratio.cat)
names(ratio_Linfs_combo)<-c("Linf>LAmax_Lmax_Lmax","Linf>LAmax","Linf<LAmax_Lmax")

Lmax_less<-VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf__Linf_tmax>Linf_Amax_ratio & VBGF_tmax_Lmax$Linf_Lmax>Linf_Lmax_ratio,]

#############################################################

Linf.in<-50
K.in<-0.1
t0.in<--1
vbgf.inflect(K.in,t0.in,10)

Linf.in<-Lmax_less$Linf[1]
K.in<-Lmax_less$K[1]
t0.in<-0

CV.in<-0.1
rand.ages<-round(runif(1000,1,round(Lmax_less$tmax)),0)
rand.ages<-round(rlnorm(1000,log(round(Lmax_less$tmax),0.05)),0)


rand.ages<-round(runif(1000,1,50),0)
Age_sample_mat<-data.frame(Ages=c(1:50),
            Samples=c(5,10,50,100,100,100,100,200,200,200,
                    200,200,150,100,100,100,100,100,100,100,
                    100,100,50,50,50,50,50,50,50,50,
                    20,20,10,10,10,10,10,5,5,5,
                    5,5,5,5,5,5,5,4,2,1))


rand.ages.5<-round(rtruncnorm(0.15*1000,1,5,4,1),0)
rand.ages.15<-round(rtruncnorm(0.8*1000,6,15,12,3),0)
rand.ages.25<-round(rtruncnorm(0.0001*1000,16,25,18,4),0)
rand.ages.35<-round(rtruncnorm(0.0001*1000,26,35,30,6),0)
rand.ages.50<-round(rtruncnorm(0.0001*1000,36,50,36,8),0)
rand.ages<-c(rand.ages.5,rand.ages.15,rand.ages.25,rand.ages.35,rand.ages.50)

#rand.ages<-round(rtruncnorm(1000,2,10,10,10),0)

#############################
### Plot VBGF and samples ###
#############################
Linf.in<-60
M<-0.1
K.in<-M/0.5
t0.in<--1
CV.in<-0.1

lowage<-0.1
maxage<-5.4/M
maxlt<-Linf.in*0.9

ages.full<-rlnormTrunc(1000, meanlog = log(Amax.in), sdlog = 1, min = lowage,max = maxage)
rand.lts.full<-rand.VBGF(Linf.in,K.in,t0.in,ages.full,CV.in)

age.lt.dat<-data.frame(Ages=round(ages.full,0),Lengths=round(rand.lts.full,0),Sample="Full")

age.lt.dat.Alim<-age.lt.dat[age.lt.dat$Ages<maxage*0.5,]
age.lt.dat.Alim$Sample<-paste0("Age ",round(maxage*0.5,0))
age.lt.dat.Llim<-age.lt.dat[age.lt.dat$Lengths<maxlt,]
age.lt.dat.Llim$Sample<-paste0("Lt ",round(maxlt,0))
age.lt.dat.AlimLlim<-age.lt.dat[age.lt.dat$Ages<maxage*0.5&age.lt.dat$Lengths<maxlt,]
age.lt.dat.AlimLlim$Sample<-paste0("Age ",maxage*0.5,", Lt ",maxlt)
age.lt.dat.Alow<-age.lt.dat[age.lt.dat$Ages<maxage*0.2,]
age.lt.dat.Alow$Sample<-"20%Amax"

age.lt.dat.all<-rbind(age.lt.dat,age.lt.dat.Alim,age.lt.dat.Llim,age.lt.dat.AlimLlim)
age.lt.dat.all$Sample<-factor(age.lt.dat.all$Sample,levels=c("Full",paste0("Age ",round(maxage*0.5,0)),paste0("Lt ",round(maxlt,0)),paste0("Age ",maxage*0.5,", Lt ",maxlt)))

vbgf.parms.all<-nls(age.lt.dat$Lengths~Linf * (1 - exp(-k * (age.lt.dat$Ages - t0))),data=age.lt.dat,start=list(Linf=max(age.lt.dat$Lengths),k=0.1,t0=0))
vbgf.parms.Alim<-nls(age.lt.dat.Alim$Lengths~Linf * (1 - exp(-k * (age.lt.dat.Alim$Ages - t0))),data=age.lt.dat.Alim,start=list(Linf=max(age.lt.dat.Alim$Lengths),k=0.1,t0=0))
vbgf.parms.Llim<-nls(age.lt.dat.Llim$Lengths~Linf * (1 - exp(-k * (age.lt.dat.Llim$Ages - t0))),data=age.lt.dat.Llim,start=list(Linf=max(age.lt.dat.Llim$Lengths),k=0.1,t0=0))
vbgf.parms.AlimLlim<-nls(age.lt.dat.AlimLlim$Lengths~Linf * (1 - exp(-k * (age.lt.dat.AlimLlim$Ages - t0))),data=age.lt.dat.AlimLlim,start=list(Linf=max(age.lt.dat.AlimLlim$Lengths),k=0.1,t0=0))
vbgf.parms.Alow<-nls(age.lt.dat.Alow$Lengths~Linf * (1 - exp(-k * (age.lt.dat.Alow$Ages - t0))),data=age.lt.dat.Alow,start=list(Linf=max(age.lt.dat.Alow$Lengths),k=0.1,t0=0))

#Fitting all samples
ggplot(age.lt.dat,aes(Ages,Lengths,colour = Sample))+
  geom_point()+
  scale_colour_manual(values = c("red"))+
  theme_bw()+
  geom_function(fun = function(x) Linf.in * (1 - exp(-K.in * (x - t0.in))), colour = "black",lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.all$m$getPars()[1] * (1 - exp(-vbgf.parms.all$m$getPars()[2] * (x - vbgf.parms.all$m$getPars()[3]))), colour = c("red"),lwd=1.25)+
  xlim(0,50)  

#Example with only 30 as maximum age
Linf.in<-60
M<-0.05
K.in<-M/1.5
t0.in<--1
CV.in<-0.1

lowage<-0.1
maxage<-30
maxlt<-100

ages.full.30<-rlnormTrunc(1000, meanlog = log(Amax.in), sdlog = 1, min = lowage,max = maxage)
rand.lts.full.30<-rand.VBGF(Linf.in,K.in,t0.in,ages.full,CV.in)
age.lt.dat.30<-data.frame(Ages=round(ages.full.30,0),Lengths=round(rand.lts.full.30,0),Sample="Full")

ggplot(age.lt.dat.30,aes(Ages,Lengths))+
  geom_point()+
  scale_colour_manual(values = c("red"))+
  theme_bw()+
  geom_function(fun = function(x) Linf.in * (1 - exp(-K.in * (x - t0.in))), colour = "black",lwd=1.25)+
  #  geom_function(fun = function(x) vbgf.parms.all$m$getPars()[1] * (1 - exp(-vbgf.parms.all$m$getPars()[2] * (x - vbgf.parms.all$m$getPars()[3]))), colour = c("red"),lwd=1.25)+
  xlim(0,50)  


#Fitting A50%
ggplot(subset(age.lt.dat.all,Sample %in% c("Full",paste0("Age ",round(maxage*0.5,0)))),aes(Ages,Lengths,colour = Sample))+
  geom_point()+
  scale_colour_manual(values = c("red","pink"))+
  theme_bw()+
  geom_function(fun = function(x) Linf.in * (1 - exp(-K.in * (x - t0.in))), colour = c("black"),lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.all$m$getPars()[1] * (1 - exp(-vbgf.parms.all$m$getPars()[2] * (x - vbgf.parms.all$m$getPars()[3]))), colour = c("red"),lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.Alim$m$getPars()[1] * (1 - exp(-vbgf.parms.Alim$m$getPars()[2] * (x - vbgf.parms.Alim$m$getPars()[3]))), colour = c("pink"),lwd=1.25)
#Fitting Linf90%
ggplot(subset(age.lt.dat.all,Sample %in% c("Full",paste0("Lt ",round(maxlt,0)))),aes(Ages,Lengths,colour = Sample))+
  geom_point()+
  scale_colour_manual(values = c("red","darkgreen"))+
  theme_bw()+
  geom_function(fun = function(x) Linf.in * (1 - exp(-K.in * (x - t0.in))), colour = c("black"),lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.all$m$getPars()[1] * (1 - exp(-vbgf.parms.all$m$getPars()[2] * (x - vbgf.parms.all$m$getPars()[3]))), colour = c("red"),lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.Llim$m$getPars()[1] * (1 - exp(-vbgf.parms.Llim$m$getPars()[2] * (x - vbgf.parms.Llim$m$getPars()[3]))), colour = c("darkgreen"),lwd=1.25)
#Fitting A50% and Linf90%
ggplot(subset(age.lt.dat.all,Sample %in% c("Full",paste0("Age ",maxage*0.5,", Lt ",maxlt))),aes(Ages,Lengths,colour = Sample))+
  geom_point()+
  scale_colour_manual(values = c("red","blue"))+
  theme_bw()+
  geom_function(fun = function(x) Linf.in * (1 - exp(-K.in * (x - t0.in))), colour = c("black"),lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.all$m$getPars()[1] * (1 - exp(-vbgf.parms.all$m$getPars()[2] * (x - vbgf.parms.all$m$getPars()[3]))), colour = c("red"),lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.AlimLlim$m$getPars()[1] * (1 - exp(-vbgf.parms.AlimLlim$m$getPars()[2] * (x - vbgf.parms.AlimLlim$m$getPars()[3]))), colour = c("blue"),lwd=1.25)
#Fit all scenarios
ggplot(age.lt.dat.all,aes(Ages,Lengths,colour = Sample))+
  geom_point()+
  scale_colour_manual(values = c("red","pink" , "darkgreen","blue"))+
  theme_bw()+
  geom_function(fun = function(x) Linf.in * (1 - exp(-K.in * (x - t0.in))), colour = "black",lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.all$m$getPars()[1] * (1 - exp(-vbgf.parms.all$m$getPars()[2] * (x - vbgf.parms.all$m$getPars()[3]))), colour = "red",lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.Alim$m$getPars()[1] * (1 - exp(-vbgf.parms.Alim$m$getPars()[2] * (x - vbgf.parms.Alim$m$getPars()[3]))), colour = "pink",lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.Llim$m$getPars()[1] * (1 - exp(-vbgf.parms.Llim$m$getPars()[2] * (x - vbgf.parms.Llim$m$getPars()[3]))), colour ="darkgreen",lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.AlimLlim$m$getPars()[1] * (1 - exp(-vbgf.parms.AlimLlim$m$getPars()[2] * (x - vbgf.parms.AlimLlim$m$getPars()[3]))), colour = "blue",lwd=1.25)
#Low ages and lengths  
age.lt.dat.comp.Alow<-rbind(age.lt.dat,age.lt.dat.Alow)
age.lt.dat.comp.Alow$Sample<-factor(age.lt.dat.comp.Alow$Sample,levels=c("Full","20%Amax"))
ggplot(age.lt.dat.comp.Alow,aes(Ages,Lengths,colour = Sample))+
  geom_point()+
  scale_colour_manual(values = c("red","orange"))+
  theme_bw()+
  geom_function(fun = function(x) Linf.in * (1 - exp(-K.in * (x - t0.in))), colour = c("black"),lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.all$m$getPars()[1] * (1 - exp(-vbgf.parms.all$m$getPars()[2] * (x - vbgf.parms.all$m$getPars()[3]))), colour = c("red"),lwd=1.25)+
  geom_function(fun = function(x) vbgf.parms.Alow$m$getPars()[1] * (1 - exp(-vbgf.parms.Alow$m$getPars()[2] * (x - vbgf.parms.Alow$m$getPars()[3]))), colour = c("orange"),lwd=1.25)



###################
### Simulations ###
###################

Linf.in<-60
M<-0.1
K.in.05<-M/0.5
K.in.1.5<-M/1.5
K.in.2<-M/2
t0.in<--1
CV.in<-0.1
CV.in20<-0.1

lowage<-0.1
maxage<-5.4/M
maxage.50<-5.4/M/2
maxage.25<-5.4/M/5
maxlt<-100
maxlt.90<-Linf.in*0.9

#M/k=0.5
Linf_Lmax_ratios_Mk05.N100<-Comp.Lmax.Linf(n=100,Nsim=500,Linf.in=Linf.in,K.in=K.in.05,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage,maxlt)
Linf_Lmax_ratios_Mk05.N200<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.05,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage,maxlt)
Linf_Lmax_ratios_Mk05.N500<-Comp.Lmax.Linf(n=500,Nsim=500,Linf.in=Linf.in,K.in=K.in.05,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage,maxlt)
Linf_Lmax_ratios_Mk05.N200_90Linf<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.05,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage,maxlt.90)
Linf_Lmax_ratios_Mk05.N200_50Amax<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.05,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage.50,maxlt)
Linf_Lmax_ratios_Mk05.N200_90Linf_50Amax<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.05,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage.50,maxlt.90)
Linf_Lmax_ratios_Mk05.N200_25Amax<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.05,t0.in=t0.in,CV.in=CV.in20,Amax.in=Amax.in,lowage=lowage,maxage=maxage.25,maxlt)
ratios100<-data.frame(Ratio=Linf_Lmax_ratios_Mk05.N100[,3],Label="Mk_05_N100")
ratios200<-data.frame(Ratio=Linf_Lmax_ratios_Mk05.N200[,3],Label="Mk_05_N200")
ratios500<-data.frame(Ratio=Linf_Lmax_ratios_Mk05.N500[,3],Label="Mk_05_N500")
gg_Linf_ratios<-rbind(ratios100,ratios200,ratios500)
ggplot(gg_Linf_ratios,aes(Ratio,color=Label))+ 
  geom_density(lwd=2)+
  geom_vline(xintercept=c(1), color="red")+
  xlim(0,2)

#M/k=1.5
Linf_Lmax_ratios_Mk15<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.1.5,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage,maxlt)
Linf_Lmax_ratios_Mk15_90Linf<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.1.5,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage,maxlt.90)
Linf_Lmax_ratios_Mk15_50Amax<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.1.5,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage.50,maxlt)
Linf_Lmax_ratios_Mk15_50Amax_90Linf<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.1.5,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage.50,maxlt.90)
Linf_Lmax_ratios_Mk15_25Amax<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.1.5,t0.in=t0.in,CV.in=CV.in20,Amax.in=Amax.in,lowage=lowage,maxage=maxage.25,maxlt)

#M/k=2.5
Linf_Lmax_ratios_Mk20<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.2,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage,maxlt)
Linf_Lmax_ratios_Mk20_90Linf<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.2,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage,maxlt.90)
Linf_Lmax_ratios_Mk20_50Amax<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.2,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage.50,maxlt)
Linf_Lmax_ratios_Mk20_50Amax_90Linf<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.2,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage.50,maxlt.90)
Linf_Lmax_ratios_Mk20_25Amax<-Comp.Lmax.Linf(n=200,Nsim=500,Linf.in=Linf.in,K.in=K.in.2,t0.in=t0.in,CV.in=CV.in20,Amax.in=Amax.in,lowage=lowage,maxage=maxage.25,maxlt)


#Full sampling
ratiosmk05<-data.frame(Ratio=Linf_Lmax_ratios_Mk05.N200[,3],Label="Mk_0.5")
ratiosmk15<-data.frame(Ratio=Linf_Lmax_ratios_Mk15[,3],Label="Mk_1.5")
ratiosmk20<-data.frame(Ratio=Linf_Lmax_ratios_Mk20[,3],Label="Mk_2")
gg_Linf_ratios<-rbind(ratiosmk05,ratiosmk15,ratiosmk20)
ggplot(gg_Linf_ratios,aes(Ratio,color=Label))+ 
  geom_density(lwd=2)+
  geom_vline(xintercept=c(1), color="red")+
  xlim(0,2)

#90% Linf
ratiosmk05_90Linf<-data.frame(Ratio=Linf_Lmax_ratios_Mk05.N200_90Linf[,3],Label="Mk_0.5")
ratiosmk15_90Linf<-data.frame(Ratio=Linf_Lmax_ratios_Mk15_90Linf[,3],Label="Mk_1.5")
ratiosmk20_90Linf<-data.frame(Ratio=Linf_Lmax_ratios_Mk20_90Linf[,3],Label="Mk_2")
gg_Linf_ratios_90Linf<-rbind(ratiosmk05_90Linf,ratiosmk15_90Linf,ratiosmk20_90Linf)
ggplot(gg_Linf_ratios_90Linf,aes(Ratio,color=Label))+ 
  geom_density(lwd=2)+
  geom_vline(xintercept=c(1), color="red")+
  xlim(0,2)

#50% Amax
ratiosmk05_50Amax<-data.frame(Ratio=Linf_Lmax_ratios_Mk05.N200_50Amax[,3],Label="Mk_0.5")
ratiosmk15_50Amax<-data.frame(Ratio=Linf_Lmax_ratios_Mk15_50Amax[,3],Label="Mk_1.5")
ratiosmk20_50Amax<-data.frame(Ratio=Linf_Lmax_ratios_Mk20_50Amax[,3],Label="Mk_2")
gg_Linf_ratios_50Amax<-rbind(ratiosmk05_50Amax,ratiosmk15_50Amax,ratiosmk20_50Amax)
ggplot(gg_Linf_ratios_50Amax,aes(Ratio,color=Label))+ 
  geom_density(lwd=2)+
  geom_vline(xintercept=c(1), color="red")+
  xlim(0,2)

#90% Linf, 50% Amax
ratiosmk05_50Amax_90Linf<-data.frame(Ratio=Linf_Lmax_ratios_Mk05.N200_90Linf_50Amax[,3],Label="Mk_0.5")
ratiosmk15_50Amax_90Linf<-data.frame(Ratio=Linf_Lmax_ratios_Mk15_50Amax_90Linf[,3],Label="Mk_1.5")
ratiosmk20_50Amax_90Linf<-data.frame(Ratio=Linf_Lmax_ratios_Mk20_50Amax_90Linf[,3],Label="Mk_2")
gg_Linf_ratios_50Amax_90Linf<-rbind(ratiosmk05_50Amax_90Linf,ratiosmk15_50Amax_90Linf,ratiosmk20_50Amax_90Linf)
ggplot(gg_Linf_ratios_50Amax_90Linf,aes(Ratio,color=Label))+ 
  geom_density(lwd=2)+
  geom_vline(xintercept=c(1), color="red")+
  xlim(0,2)

#20% Amax
ratiosmk05_25Amax<-data.frame(Ratio=Linf_Lmax_ratios_Mk05.N200_25Amax[,3],Label="Mk_0.5")
ratiosmk15_25Amax<-data.frame(Ratio=Linf_Lmax_ratios_Mk15_25Amax[,3],Label="Mk_1.5")
ratiosmk20_25Amax<-data.frame(Ratio=Linf_Lmax_ratios_Mk20_25Amax[,3],Label="Mk_2")
gg_Linf_ratios_25Amax<-rbind(ratiosmk05_25Amax,ratiosmk15_25Amax,ratiosmk20_25Amax)
ggplot(gg_Linf_ratios_25Amax,aes(Ratio,color=Label))+ 
  geom_density(lwd=2)+
  geom_vline(xintercept=c(1), color="red")+
  xlim(0,2)

#3 indicator plots
Linf_Lmax_plots(Linf_Lmax_ratios_Mk05.N200) #3-plot N200, m/k=0.5
Linf_Lmax_plots(Linf_Lmax_ratios_Mk15) #3-plot N200, m/k=1.5
Linf_Lmax_plots(Linf_Lmax_ratios_Mk15_90Linf) #3-plot N200, m/k=1.5 90Linf
Linf_Lmax_plots(Linf_Lmax_ratios_Mk15_50Amax) #3-plot N200, m/k=1.5 50Linf
Linf_Lmax_plots(Linf_Lmax_ratios_Mk15_50Amax_90Linf) #3-plot N200, m/k=1.5 50Linf
Linf_Lmax_plots(Linf_Lmax_ratios_Mk20_50Amax_90Linf) #3-plot N200, m/k=1.5 50Linf
Linf_Lmax_plots(Linf_Lmax_ratios_Mk20_25Amax) #3-plot N200, m/k=1.5 20%Amax




Linf_plots(Linf_Lmax_ratios)
Linf_Lmax_plots(Linf_Lmax_ratios)
summary(Linf_Lmax_ratios[])  










Linf_Lmax_grid<-data.frame(TLinf_Lmax =NA,estLinf_TLinf=NA,K=NA)
K_vec<-seq(0.05,1,0.05)
for(i in 1:length(K_vec))
{
  Linf_Lmax_ratios<-Comp.Lmax.Linf(n=500,Nsim=500,Linf.in=Linf.in,K.in=K_vec[i],t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage)
  Linf_Lmax_grid[i,1]<-summary(Linf_Lmax_ratios[,1])[3] 
  Linf_Lmax_grid[i,2]<-summary(Linf_Lmax_ratios[,3])[3] 
  Linf_Lmax_grid[i,3]<-K_vec[i]
}

Linf_Lmax_plots(Linf_Lmax_ratios)
summary(Linf_Lmax_ratios[])  




abs(Linf_Lmax_ratios$estLinf_TLinf-1)>0.05









