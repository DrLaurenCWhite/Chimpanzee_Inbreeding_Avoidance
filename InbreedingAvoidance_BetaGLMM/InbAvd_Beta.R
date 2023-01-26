rm(list = ls())
library(data.table)
library(rethinking)

#import and format data
D=fread("./InbAvd_Input.csv")

#Remove data point(s) in which sire only occurs once (not usable with unidirectional indexing)
D[Male %in% D[,.N, by=Male][N<=1]$Male] #There are none here

#For unidirectional indexing, each individual needs to occur at least once in each column:
D$idA=as.character()
D$idB=as.character()
set.seed(42)
for (i in 1:nrow(D)){
  ida=sample(c(D[i,]$Male, D[i,]$Female), 1)
  idb=ifelse(ida==D[i,]$Male, D[i,]$Female, D[i,]$Male)
  D[i,]$idA=ida
  D[i,]$idB=idb
}
#Check that worked. The below should be empty data tables. (Run the for loop again if not)
D[!idA %in% idB]
D[!idB %in% idA]

#Sample sizes
D[, .N] #data points
D[, .N, by=Actual] #By parent status (actual(true) or potential(false))
D[Actual==TRUE, .N, by=F.ImmiOrNatal] #By female natality status
D[Actual==FALSE, .N, by=F.ImmiOrNatal] #By female natality status


dat=list(
  rel=D$GeneticR+0.00001,
  sire=ifelse(D$Actual==TRUE, 1, 0),
  Nat_I=ifelse(D$F.ImmiOrNatal=="Immigrant",1,0),
  Nat_U=ifelse(D$F.ImmiOrNatal=="Unknown",1,0),
  ida_index=as.integer(as.factor(D$idA)),
  idb_index=as.integer(as.factor(D$idB)),
  off_id=as.integer(as.factor(D$Offspring))
)

#If you do not want to fit the model, load it from this object file:
load("InbAvd_ModelObject.rdata")

#Fit model
InbAvd=map2stan(
  alist(
    rel ~ dbeta2(p, theta),
    logit(p) <- alpha + bs*sire + #fixed effects and intercept
      bnI*Nat_I + bsXnI*sire*Nat_I + 
      bnU*Nat_U + bsXnU*sire*Nat_U + 
      
      alpha_id[ida_index] + alpha_id[idb_index] + alpha_off[off_id] + #varying intercepts
      
      bs_id[ida_index]*sire + 
      bnI_id[ida_index]*Nat_I + bsXnI_id[ida_index]*Nat_I*sire + #varying slopes
      bnU_id[ida_index]*Nat_U + bsXnU_id[ida_index]*Nat_U*sire +
      
      bs_id[idb_index]*sire + 
      bnI_id[idb_index]*Nat_I + bsXnI_id[idb_index]*Nat_I*sire + #varying slopes
      bnU_id[idb_index]*Nat_U + bsXnU_id[idb_index]*Nat_U*sire +
      
      bs_off[off_id]*sire + 
      bnI_off[off_id]*Nat_I + bsXnI_off[off_id]*Nat_I*sire + #varying slopes
      bnU_off[off_id]*Nat_U + bsXnU_off[off_id]*Nat_U*sire,
    
    
    c(alpha_id, bs_id, bnI_id, bsXnI_id, bnU_id, bsXnU_id)[ida_index] ~ dmvnormNC(sigma_id, Rho_id),
    Rho_id ~ dlkjcorr(3),
    sigma_id ~ dcauchy(0,2),
    
    c(alpha_off, bs_off, bnI_off, bsXnI_off, bnU_off, bsXnU_off)[off_id] ~ dmvnormNC(sigma_off, Rho_off),
    Rho_off ~ dlkjcorr(3),
    sigma_off ~ dcauchy(0,2),
    
    c(alpha,bs,bnI,bsXnI, bnU, bsXnU) ~ normal(0,1.5),
    theta ~ dcauchy ( 2, 1 )
  ), data=dat, cores=2, iter=2000, chains=2
)


#Save model object
#save(InbAvd, file="InbAvd_ModelObject.rdata")

#Fixed effects coefficients, number of effective samples and R-hat
precis(InbAvd) #R-hat and n_eff look fine

#Fixed effects+varying effects coefficients, number of effective samples and R-hat
precis(InbAvd, depth=3) #almost all of R-hat and n_eff are fine, some varying effects are fit poorly, but should be OK

#Traceplots
plot(InbAvd) #almost all are fine, some varying effects are fit poorly, but should be OK


#Extract posterior medians and 89% credbility intervals for the different combinations of categories

#Dummy code set up for natality:
#Nat_I==0 & Nat_U==0 <-intercept=Natal females
#Nat_I==1 & Nat_U==0 <-intercept+Nat_I = Immigrant females
#Nat_I==0 & Nat_U==1 <-intercept+Nat_U = Unknown females
pred=list(
  sire=c(0,1,0,1,0,1),
  Nat_I=c(0,0,1,1,0,0),
  Nat_U=c(0,0,0,0,1,1),
  ida_index=c(1,1,1,1,1,1),
  idb_index=c(1,1,1,1,1,1),
  off_id=c(1,1,1,1,1,1)
)

a_id <- matrix(0,1000,length(unique(dat$ida_index)))
of_id <- matrix(0,1000,length(unique(dat$off_id)))
P=link(InbAvd,data=pred,n=1000,
       replace=list(alpha_id=a_id,bs_id=a_id, bnI_id=a_id, bsXnI_id=a_id, bnU_id=a_id, bsXnU_id=a_id, 
                    alpha_off=of_id, bs_off=of_id, bnI_off=of_id, bsXnI_off=of_id, bnU_off=of_id, bsXnU_off=of_id))		

round(median(P[,1]), 3) #natal potential
round(HPDI(P[,1]), 3) #natal potential
round(median(P[,2]), 3) #natal actual
round(HPDI(P[,2]), 3) #natal actual

round(median(P[,3]), 3) #immigrant potential
round(HPDI(P[,3]), 3) #immigrant potential
round(median(P[,4]), 3) #immigrant actual
round(HPDI(P[,4]), 3) #immigrant actual

round(median(P[,5]), 3) #unknown potential
round(HPDI(P[,5]), 3) #unknown potential
round(median(P[,6]), 3) #unknown actual
round(HPDI(P[,6]), 3) #unknown actual

