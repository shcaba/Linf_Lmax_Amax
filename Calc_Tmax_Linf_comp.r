library(freeR)
library(ggplot2)
library(rfishbase)
library(truncnorm)
library(fishmethods)
library(FSAsim)
library(FSA)
library(EnvStats)
library(reshape2)

### Functions ###
VBGF<-function(Linf, k, t0, ages){ 
  Lts<-Linf * (1 - exp(-k * (ages - t0)))
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
#############

#Upload VBGF and tmax values from FishBase
fb.spp.parms<-read.csv("C:/Users/Jason.Cope/Documents/Current Action/Publications/tmax vs Linf/Fishbase_spp_names.csv")

#Use species names from FishBase that have growth parameters and extract Linf, K, and tmax values using FishLIfe. 
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
VBGF.fb.out$Linf__Linf_tmax<-as.numeric(VBGF.fb.out$Linf)/VBGF.fb.out$Linf_tmax
#Remove the NAs
VBGF.fb.out_noNA<-VBGF.fb.out[ !is.na(VBGF.fb.out$Linf__Linf_tmax),]
#Convert the rest of the inputs to numeric
VBGF.fb.out_noNA$Linf<-as.numeric(VBGF.fb.out_noNA$Linf)
VBGF.fb.out_noNA$K<-as.numeric(VBGF.fb.out_noNA$K)
VBGF.fb.out_noNA$tmax<-as.numeric(VBGF.fb.out_noNA$tmax)
#Save as R object
save(VBGF.fb.out_noNA,file="C:/Users/Jason.Cope/Documents/Current Action/Publications/tmax vs Linf/VBGF_fb_out_noNA.rds")
load("C:/Users/Jason.Cope/Documents/Current Action/Publications/tmax vs Linf/VBGF_fb_out_noNA.rds")

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
VBGF_tmax_Lmax$Linf_Lmax<-VBGF_tmax_Lmax$Linf/VBGF_tmax_Lmax$Lmax
VBGF_tmax_Lmax$Linflect<-VBGF_tmax_Lmax$Linf*(8/27)
save(VBGF_tmax_Lmax,file="C:/Users/Jason.Cope/Documents/Current Action/Publications/tmax vs Linf/VBGF_tmax_Lmax.rds")
load("C:/Users/Jason.Cope/Documents/Current Action/Publications/tmax vs Linf/VBGF_tmax_Lmax.rds")

#Plots to compare different parameter to Linf:L_Amax
ggplot(VBGF_tmax_Lmax,aes(Linf__Linf_tmax,Linf))+
    geom_point()+
    geom_vline(xintercept=1, color="red")

ggplot(VBGF_tmax_Lmax,aes(Linf__Linf_tmax,K))+
    geom_point()

ggplot(VBGF_tmax_Lmax,aes(Linf__Linf_tmax,tmax))+
    geom_point()

#Compare to Lmax
ggplot(VBGF_tmax_Lmax,aes(Linf_Lmax,Linf))+
    geom_point()+
    geom_vline(xintercept=1, color="red")

ggplot(VBGF_tmax_Lmax,aes(Linf_Lmax,Linflect))+
  geom_point()+
  geom_vline(xintercept=1, color="red")


ggplot(VBGF_tmax_Lmax,aes(Linf_Lmax,Linf__Linf_tmax))+
    geom_point()+
    geom_hline(yintercept=1.1,color="red")+
    geom_vline(xintercept=1, color="red")+
    xlim(0,2)


#Calculate ratio categories
#High Linf vs Amax; Linf>Lmax WORST LINF
#High Linf vs Amax; Linf<Lmax BAD Amax
#Linf vs Amax; Linf>Lmax BAD LINF
#Linf vs Amax; Linf<Lmax Good

Linf_Amax_ratio<-1.1
Linf_Lmax_ratio<-1
ratio.cat<-c(dim(VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf__Linf_tmax>Linf_Amax_ratio & VBGF_tmax_Lmax$Linf_Lmax>Linf_Lmax_ratio,])[1],
dim(VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf__Linf_tmax>Linf_Amax_ratio & VBGF_tmax_Lmax$Linf_Lmax<Linf_Lmax_ratio,])[1],
dim(VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf__Linf_tmax<Linf_Amax_ratio & VBGF_tmax_Lmax$Linf_Lmax>Linf_Lmax_ratio,])[1],
dim(VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf__Linf_tmax<Linf_Amax_ratio & VBGF_tmax_Lmax$Linf_Lmax<Linf_Lmax_ratio,])[1])

ratio_Linfs<-ratio.cat/sum(ratio.cat)
names(ratio_Linfs)<-c("Linf>LAmax_Lmax","Linf>LAmax","Linf>Lmax","Linf<LAmax_Lmax")

ratio_Linfs_combo<-c(sum(ratio.cat[c(1,3)]),ratio.cat[2],ratio.cat[4])/sum(ratio.cat)
names(ratio_Linfs_combo)<-c("Linf>LAmax_Lmax_Lmax","Linf>LAmax","Linf<LAmax_Lmax")

Lmax_less<-VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf__Linf_tmax>Linf_Amax_ratio & VBGF_tmax_Lmax$Linf_Lmax>Linf_Lmax_ratio,]



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



#rand.VBGF(Linf.in,K.in,t0.in,Amax.in,0.1)
#summary(rand.VBGF(Linf.in,K.in,t0.in,rep(Amax.in,100),0.1))

#Linf.in<-Lmax_less$Linf[1]
#K.in<-Lmax_less$K[1]
#t0.in<-0
#CV.in<-0.1
#Amax.in<-Lmax_less$tmax[1]


Comp.Lmax.Linf<-function(n=100,Nsim=500,Linf.in,K.in,t0.in,CV.in,Amax.in,lowage,maxage)
{
  Linf_Lmax_ratios<-data.frame(TLinf_Lmax=NA,estLinf_Lmax=NA,estLinf_TLinf=NA,Linf.true=NA,Linf_est=NA,Lmax=NA, Amin=NA, Amax=NA,Kest=NA,K_jit=NA,toest=NA,t0_jit=NA)
  for(i in 1:Nsim)
  {
    rand.ages<-rlnormTrunc(n, meanlog = log(Amax.in), sdlog = 1, min = lowage,max = maxage)
    #rand.ages.trunc<-rand.ages[rand.ages>lowage]
    rand.lts<-rand.VBGF(Linf.in,K.in,t0.in,rand.ages,CV.in)
    #  Age_Lts<-data.frame(Ages=rand.ages,Lts=rand.lts)
    #  plot(rand.ages,rand.lts)
    #Fit VBGF curve
    k_t0_jits<-c(jitter(K.in,amount=0.05),jitter(t0.in,amount=0.05))
    growth.out<-try(growth(size=rand.lts,age=round(rand.ages,0),Sinf=max(rand.lts),K=k_t0_jits[1],t0=0,graph=FALSE,control=list(maxiter=10000,minFactor=100,tol=1e-5))$vout) #k_t0_jits[2]
    #vbmod <- Lts ~ Linf * (1 - exp(-K * (Ages - t0)))
    #vbgf.parms<-list(Linf=20,K=0.1,t0=0)
    #vbgf.fit<-nls(vbmod,data=Age_Lts,start=vbgf.parms)
    if(any(growth.out=="Fit failed"))
    {Linf_Lmax_ratios[i,]=NA}
    
    if(any(growth.out!="Fit failed"))
    {
      growth.out.parms<-growth.out$m$getPars()
      Linf_Lmax_ratios[i,1]<-Linf.in/max(rand.lts)
      Linf_Lmax_ratios[i,2]<-growth.out.parms[1]/max(rand.lts)
      Linf_Lmax_ratios[i,3]<-growth.out.parms[1]/Linf.in
      Linf_Lmax_ratios[i,4]<-Linf.in
      Linf_Lmax_ratios[i,5]<-growth.out.parms[1]
      Linf_Lmax_ratios[i,6]<-max(rand.lts)
      Linf_Lmax_ratios[i,7]<-min(rand.ages.trunc)
      Linf_Lmax_ratios[i,8]<-max(rand.ages.trunc)
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
      geom_vline(xintercept=1, color="red")+
      xlim(0,2)
}


###########
Linf.in<-75
K.in<-0.05
t0.in<--3
CV.in<-0.1

Amax.in<-10
lowage<-2
maxage<-10
##########
Linf_Lmax_ratios<-Comp.Lmax.Linf(n=500,Nsim=500,Linf.in=Linf.in,K.in=K.in,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage)
Linf_Lmax_plots(Linf_Lmax_ratios)
summary(Linf_Lmax_ratios[])  


##########
Linf.in<-50
K.in<-0.1
t0.in<--2
CV.in<-0.1

Amax.in<-10
lowage<-1
maxage<-10
###########

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

###########
Linf.in<-20
K.in<-0.3
t0.in<--1
CV.in<-0.1

Amax.in<-5
lowage<-1
maxage<-3
###########
Linf_Lmax_ratios<-Comp.Lmax.Linf(n=500,Nsim=500,Linf.in=Linf.in,K.in=K.in,t0.in=t0.in,CV.in=CV.in,Amax.in=Amax.in,lowage=lowage,maxage=maxage)
Linf_Lmax_plots(Linf_Lmax_ratios)
summary(Linf_Lmax_ratios[])  




#plot(Linf_Lmax_ratios[,2],Linf_Lmax_ratios[,3])



abs(Linf_Lmax_ratios$estLinf_TLinf-1)>0.05










#Reduce sampling
ageless<-25
ages.less<-rand.ages[rand.ages<ageless]
Lts.less<-Lts.rand[rand.ages<ageless]
growth(size=Lts.less,age=ages.less,Sinf=40,K=0.15,t0=0)$vout

#Max length < Linf
ltsless<-40
Lts.lessLinf<-Lts.rand[rand.lts<ltsless]
ages.lessLinf<-rand.ages[rand.lts<ltsless]
growth(size=Lts.lessLinf,age=ages.lessLinf,Sinf=50,K=0.1,t0=0)$vout

#Max length < Linf
Lts.lessLinfage<-Lts.rand[rand.lts<ltsless&rand.ages<ageless]
ages.lessLinfage<-rand.ages[rand.lts<ltsless&rand.ages<ageless]
growth(size=Lts.lessLinfage,age=ages.lessLinfage,Sinf=50,K=0.1,t0=0)$vout


#Compare distribution of parameters when Linf ratios > 1
xx<-VBGF_tmax_Lmax[VBGF_tmax_Lmax$Linf__Linf_tmax>Linf_Amax_ratio &VBGF_tmax_Lmax$Linf_Lmax>Linf_Lmax_ratio,]
