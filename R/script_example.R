# empty array
param=array(0,c(3,2,7,6))

param[1,,,1] <- rbind(c(1,3,3,1,1.5,-0.1,2), c(1,3,3,1,1.5,-0.1,2))

# beta1 for both sites and all vital rates (site 2 responds stonger to environment)
param[1,,,2] <- rbind(rep(0.08,7),c(rep(0.08,5),-0.08,0.1))

# beta2 for both sites and all vital rates (no differences between sites in repsonse to intraspec. density)
param[1,,,3] <- rbind(c(rep(-0.07,5),0.12,-0.1), c(rep(-0.07,5),0.12,-0.1))

# beta3 for both sites and all vital rates (in site 1 response to species 2 is stonger than in site 2)
param[1,,,4] <- rbind(c(rep(-0.1,3),0.05,rep(-0.1,3)),c(rep(-0.1,3),0.02,rep(-0.1,3)))

# beta4 for both sites and all vital rates (in site 1 response to species 3 is weaker than in site 2)
param[1,,,5] <- rbind(c(rep(-0.02,3),0.02,rep(-0.02,3)), c(rep(-0.03,3),0.05,rep(-0.03,3)))

# beta5 for both sites and all vital rates (equal for all)
param[1,,,6] <- rbind(rep(0.01,7),rep(0.01,7))

param[2,,,1] <- rbind(c(0.8,3.1,3.5,0.8,2,0.1,0.5), c(0.8,3.1,3.5,0.8,2,0.3,0.5))

# beta1 for both sites and all vital rates (site 2 responds stonger to environment)
param[2,,,2] <- rbind(rep(0.07,7),c(rep(0.07,5),-0.07,0.07))

# beta2 for both sites and all vital rates (no differences between sites in repsonse to intraspec. density)
param[2,,,3] <- rbind(c(rep(-0.2,5),0.17,-0.2),c(rep(-0.2,5),0.17,-0.2))

# beta3 for both sites and all vital rates (in site 1 response to species 1 is stonger than in site 2)
param[2,,,4] <- rbind(c(rep(0.04,3),-0.05,rep(0.04,3)),c(rep(0.02,3),-0.02,rep(0.02,3)))

# beta4 for both sites and all vital rates (in site 1 response to species 3 is weaker than in site 2)
param[2,,,5] <- rbind(c(-0.1,rep(-0.02,),0.01,rep(-0.02,3)),c(-0.1,rep(-0.02,2),0.05,rep(-0.02,3)))

# beta5 for both sites and all vital rates (equal for all)
param[2,,,6] <- rbind(rep(0.01,7),rep(0.01,7))

param[3,,,1] <- rbind(c(1.0,3.8,4.5,0.8,2,0.1,0.8), c(1.0,3.8,4.5,0.8,2,0.1,0.8))

# beta1 for both sites and all vital rates 
param[3,,,2] <- rbind(rep(0.05,7),c(rep(0.05,5),-0.05,0.05))

# beta2 for both sites and all vital rates (no differences between sites in repsonse to intraspec. density)
param[3,,,3] <- rbind(c(rep(-0.18,5),0.2,-0.18),c(rep(-0.18,5),0.2,-0.18))

# beta3 for both sites and all vital rates 
param[3,,,4] <- rbind(c(rep(0.12,3),-0.15,rep(0.12,3)),c(rep(0.12,3),-0.15,rep(0.12,3)))

# beta4 for both sites and all vital rates (in site 1 response to species 2 is weaker than in site 2)
param[3,,,5] <- rbind(c(rep(-0.01,3),0.01,rep(-0.01,3)),c(rep(-0.01,3),0.01,rep(-0.01,3)))

# beta5 for both sites and all vital rates (equal for all)
param[3,,,6] <- rbind(rep(0.01,7),rep(0.01,7))

# vital rate function
vr.mpm <- function(species,site,vr,parameters,env,D1,D2,D3){
  
  # for rates
  if(vr!=7){
    
    mean=exp(sum(parameters[species,site,vr,]*c(1,env,D1,D2,D3,env*D3)))
    output=mean/(1+mean)
    
  }else{ # for reproductive output
    
    mean=exp(sum(parameters[species,site,vr,]*c(1,env,D1,D2,D3,env*D3)))
    output=mean
  }
  
  return(output)
}

# define number of stages in MPM:

n=3

#define number of sites:
n.site=2
# Basic structure needed (for each species: s1, s2, s3) 

# vec permutation approach
s1_B=s2_B=s3_B=matrix(0,n*n.site,n*n.site) # placeholder block-diagonal demography matrix
s1_M=s2_M=s3_M=matrix(0,n.site*n,n.site*n) # placeholder block-diagonal dispersal matrix

diag(s1_M)=diag(s2_M)=diag(s3_M)=1

P=matrix(0,n.site*n,n.site*n) # vec‐permutation matrix

for(i in 1:n){
  for(j in 1:n.site){
    E=matrix(0,n,n.site)
    E[i,j]=1
    P=P+(E%x%t(E))
  }
} 

s1_P=s2_P=s3_P=P

# Initial abundances
years=100 # years (or time steps) of simulations

# simulate environment from normal distribution

set.seed(123)
env=rnorm(years, mean=0, sd=1)

## Initial S1

densS1.1=c(10,10,15) # site 1
densS1.2=c(10,10,13) # site 2

densMet_S1=c(densS1.1,densS1.2) # density metapopulation

dens.biS1.1=sum(densS1.1)
dens.biS1.2=sum(densS1.2)

## Initial S2

densS2.1=c(15,5,3) # site 1
densS2.2=c(15,5,3) # site 2

densMet_S2=c(densS2.1,densS2.2) # density metapopulation

dens.biS2.1=sum(densS2.1)
dens.biS2.2=sum(densS2.2)

## Initial S3

densS3.1=c(5,5,2) # site 1
densS3.2=c(3,4,2) # site 2

densMet_S3=c(densS3.1,densS3.2) # density metapopulation

dens.biS3.1=sum(densS3.1)
dens.biS3.2=sum(densS3.2)

sim.data=array(0,c(3,n,n.site,years)) # stage- specific abundances for species and sites

for(i in 1:years){
  
  # Update abundance
  
  sim.data[1,,,i]= cbind(densMet_S1[1:n],densMet_S1[(n+1):length(densMet_S1)])
  sim.data[2,,,i]= cbind(densMet_S2[1:n],densMet_S2[(n+1):length(densMet_S2)])
  sim.data[3,,,i]= cbind(densMet_S3[1:n],densMet_S3[(n+1):length(densMet_S3)])
  
  # LOCAL DEMOGRAPHY
  # Local demography species 1
  
  # Site 1
  matS1.1=matrix(c(0,vr.mpm(1,1,1,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1),0,
                   0,vr.mpm(1,1,2,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1)*(1-vr.mpm(1,1,4,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1)),
                   vr.mpm(1,1,2,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1)*vr.mpm(1,1,4,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1),
                   vr.mpm(1,1,3,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1)*vr.mpm(1,1,7,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1),vr.mpm(1,1,3,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1)*(1-vr.mpm(1,1,5,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1)),
                   vr.mpm(1,1,3,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1)*vr.mpm(1,1,5,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1)),n,n) # I fill all the matrices by column (the default)
  
  # Site 2
  matS1.2=matrix(c(0,vr.mpm(1,2,1,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2),0,
                   0,vr.mpm(1,2,2,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2)*(1-vr.mpm(1,2,4,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2)),
                   vr.mpm(1,2,2,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2)*vr.mpm(1,2,4,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2),
                   vr.mpm(1,2,3,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2)*vr.mpm(1,2,7,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2),vr.mpm(1,2,3,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2)*(1-vr.mpm(1,2,5,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2)),
                   vr.mpm(1,2,3,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2)*vr.mpm(1,2,5,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2)),n,n)
  
  ############################################  
  # Local demography species 2
  
  # Site 1
  matS2.1=matrix(c(0,vr.mpm(2,1,1,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1),0,
                   0,vr.mpm(2,1,2,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1)*(1-vr.mpm(2,1,4,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1)),
                   vr.mpm(2,1,2,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1)*vr.mpm(2,1,4,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1),
                   vr.mpm(2,1,3,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1)*vr.mpm(2,1,7,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1),vr.mpm(2,1,3,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1)*(1-vr.mpm(2,1,5,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1)),
                   vr.mpm(2,1,3,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1)*vr.mpm(2,1,5,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1)),n,n) # I fill all the matrices by column (the default)
  
  # Site 2
  matS2.2=matrix(c(0,vr.mpm(2,2,1,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2),0,
                   0,vr.mpm(2,2,2,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2)*(1-vr.mpm(2,2,4,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2)),
                   vr.mpm(2,2,2,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2)*vr.mpm(2,2,4,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2),
                   vr.mpm(2,2,3,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2)*vr.mpm(2,2,7,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2),vr.mpm(2,2,3,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2)*(1-vr.mpm(2,2,5,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2)),
                   vr.mpm(2,2,3,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2)*vr.mpm(2,2,5,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2)),n,n)
  
  ############################################  
  # Local demography species 3
  
  # Site 1
  matS3.1=matrix(c(0,vr.mpm(3,1,1,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1),0,
                   0,vr.mpm(3,1,2,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1)*(1-vr.mpm(3,1,4,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1)),
                   vr.mpm(3,1,2,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1)*vr.mpm(3,1,4,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1),
                   vr.mpm(3,1,3,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1)*vr.mpm(3,1,7,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1),vr.mpm(3,1,3,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1)*(1-vr.mpm(3,1,5,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1)),
                   vr.mpm(3,1,3,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1)*vr.mpm(3,1,5,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1)),n,n) # I fill all the matrices by column (the default)
  
  # Site 2
  matS3.2=matrix(c(0,vr.mpm(3,2,1,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2),0,
                   0,vr.mpm(3,2,2,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2)*(1-vr.mpm(3,2,4,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2)),
                   vr.mpm(3,2,2,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2)*vr.mpm(3,2,4,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2),
                   vr.mpm(3,2,3,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2)*vr.mpm(3,2,7,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2),vr.mpm(3,2,3,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2)*(1-vr.mpm(3,2,5,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2)),
                   vr.mpm(3,2,3,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2)*vr.mpm(3,2,5,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2)),n,n)
  
  # METAPOPULATION
  
  s1_B[1:n,1:n]=matS1.1
  s1_B[(n+1):(2*n),(n+1):(2*n)]=matS1.2
  
  s2_B[1:n,1:n]=matS2.1
  s2_B[(n+1):(2*n),(n+1):(2*n)]=matS2.2
  
  s3_B[1:n,1:n]=matS3.1
  s3_B[(n+1):(2*n),(n+1):(2*n)]=matS3.2
  
  s1_M[(n.site+1),(n.site+1)]=1-vr.mpm(1,1,6,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1) # staying in site 1
  s1_M[(n.site+2),(n.site+2)]=1-vr.mpm(1,2,6,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2) # staying in site 2
  s1_M[(n.site*n),(n.site+1)] = vr.mpm(1,1,6,param,env[i],dens.biS1.1,dens.biS2.1,dens.biS3.1) # moving to reproductive in site 2
  s1_M[(n.site*n-1),(n.site+2)] = vr.mpm(1,2,6,param,env[i],dens.biS1.2,dens.biS2.2,dens.biS3.2) # moving to reproductive in site 1
  
  s2_M[(n.site+1),(n.site+1)]=1-vr.mpm(2,1,6,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1) # staying in site 1
  s2_M[(n.site+2),(n.site+2)]=1-vr.mpm(2,2,6,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2) # staying in site 2
  s2_M[(n.site*n),(n.site+1)] = vr.mpm(2,1,6,param,env[i],dens.biS2.1,dens.biS1.1,dens.biS3.1) # moving to reproductive in site 2
  s2_M[(n.site*n-1),(n.site+2)] = vr.mpm(2,2,6,param,env[i],dens.biS2.2,dens.biS1.2,dens.biS3.2) # moving to reproductive in site 1
  
  s3_M[(n.site+1),(n.site+1)]=1-vr.mpm(3,1,6,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1) # staying in site 1
  s3_M[(n.site+2),(n.site+2)]=1-vr.mpm(3,2,6,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2) # staying in site 2
  s3_M[(n.site*n),(n.site+1)] = vr.mpm(3,1,6,param,env[i],dens.biS3.1,dens.biS1.1,dens.biS2.1) # moving to reproductive insite 2
  s3_M[(n.site*n-1),(n.site+2)] = vr.mpm(3,2,6,param,env[i],dens.biS3.2,dens.biS1.2,dens.biS2.2) # moving to reproductive in site 1
  
  # UPDATE DENSITIES
  
  densMet_S1=t(s1_P)%*%s1_M%*%s1_P%*%s1_B%*%densMet_S1
  
  densMet_S2=t(s2_P)%*%s2_M%*%s2_P%*%s2_B%*%densMet_S2
  
  densMet_S3=t(s3_P)%*%s1_M%*%s3_P%*%s3_B%*%densMet_S3
  
  dens.biS1.1=sum(densMet_S1[1:n])
  dens.biS1.2=sum(densMet_S1[(n+1):length(densMet_S1)])
  
  dens.biS2.1=sum(densMet_S2[1:n])
  dens.biS2.2=sum(densMet_S2[(n+1):length(densMet_S2)])
  
  dens.biS3.1=sum(densMet_S3[1:n])
  dens.biS3.2=sum(densMet_S3[(n+1):length(densMet_S3)])
  
}

library(plyr)
library(ggplot2)

df=adply(sim.data,c(1,2,3,4))

colnames(df)=c("species","stage","site","year","density")

levels(df$stage)=c("J","N","R")
levels(df$species)=c("S1","S2","S3")
df$year=as.numeric(df$year)


ggplot(df,aes(year,density,col=species))+
  geom_line()+
  facet_grid(stage~site,scales = "free")+
  scale_color_manual(name="",values=c("darkgreen","orange","darkred"))+
  xlab("Simulation year")+ylab("Total density")+theme_bw(base_size=20)+
  theme(panel.grid = element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))



