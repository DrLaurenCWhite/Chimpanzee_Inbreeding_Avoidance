rm(list = ls())
library(ggplot2)
library(data.table)
library(rethinking)
library(ggpubr)

#Import and format data
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


#Plot the raw data
D[is.na(PedigreeCategory), PedigreeCategory:="Unknown"]
D$PedigreeCategory <- factor(D$PedigreeCategory, levels = c("PO", "HS", "HPN", "Unknown"))
D$F.ImmiOrNatal <- factor(D$F.ImmiOrNatal , levels = c("Natal", "Immigrant", "Unknown"))

ggplot(D, aes(x=Actual, y=GeneticR+0.00001)) + geom_jitter(aes(col=PedigreeCategory, alpha=PedigreeCategory), size=2, width=0.3) +
  facet_grid(.~F.ImmiOrNatal) + 
  ylab("Genetic Relatedness") + scale_x_discrete(labels=c("TRUE"="Actual\n", "FALSE"="Potential\n")) +
  scale_color_manual(name = "Pedigree\nCategory", 
                     values=c("#CC79A7", "#E69F00",  "#009E73", "#999999"),
                     labels = c("Parent-Offspring", "Half Siblings",  "Half Avuncular", "Unknown")) +
  scale_alpha_manual(values = c(1, 1,1, 0.4)) + guides(alpha=FALSE) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.title = element_text(size=15),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_blank(),
        strip.text = element_text(size=15, colour="black"))



#Load model
load("InbAvd_ModelObject.rdata")
dat=list(
  rel=D$GeneticR+0.00001,
  sire=ifelse(D$Actual==TRUE, 1, 0),
  Nat_I=ifelse(D$F.ImmiOrNatal=="Immigrant",1,0),
  Nat_U=ifelse(D$F.ImmiOrNatal=="Unknown",1,0),
  ida_index=as.integer(as.factor(D$idA)),
  idb_index=as.integer(as.factor(D$idB)),
  off_id=as.integer(as.factor(D$Offspring))
)


#Extract model predictions

#zeros for varying effects to just plot main effect
a_id <- matrix(0,1000,length(unique(dat$ida_index)))
of_id <- matrix(0,1000,length(unique(dat$off_id)))

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

P=link(InbAvd,data=pred,n=1000,
       replace=list(alpha_id=a_id,bs_id=a_id, bnI_id=a_id, bsXnI_id=a_id, bnU_id=a_id, bsXnU_id=a_id, 
                    alpha_off=of_id, bs_off=of_id, bnI_off=of_id, bsXnI_off=of_id, bnU_off=of_id, bsXnU_off=of_id))		


#Plot poster distributions with parent-type seperately -> This is the version presented in the paper
par(mfrow=c(2,1), mar=c(2,2,0.5,0.5), oma=c(2,0,0,0))

dens(P[,1], xlim=c(0.005,0.03) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE, 
     ylim=c(-11,400),mar=c(0,0,0,0), frame.plot=FALSE)
axis(1, cex=1,labels=NA, at=seq(0.005,0.03,by=0.005))
axis(1, cex.axis=1.1,at=seq(0.005,0.03,by=0.005),labels=seq(0.005,0.03,by=0.005),line=0.4,col=NA)
mtext(side=2,line=0,cex=1.1,text="Density - Potential Parents")
shade( density(P[,1]) , lim= as.vector(HPDI(P[,1], prob=0.999999)) , col = col.alpha("orange", 0.5)) #Natal potential
shade( density(P[,3]) , lim= as.vector(HPDI(P[,3], prob=0.999999)) , col = col.alpha("cyan3", 0.5)) #Immigrant potential
shade( density(P[,5]) , lim= as.vector(HPDI(P[,5], prob=0.999999)) , col = col.alpha("red", 0.5)) #Unknown potential
abline(v=median(P[,1]) , lty=1,lwd=1.5)
abline(v=median(P[,3]) , lty=2,lwd=1.5)
abline(v=median(P[,5]) , lty=3,lwd=1.5)

dens(P[,1], xlim=c(0.005,0.03) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE, 
     ylim=c(-11,400),mar=c(0,0,0,0), frame.plot=FALSE)
axis(1, cex=1,labels=NA, at=seq(0.005,0.03,by=0.005))
axis(1, cex.axis=1.1,at=seq(0.005,0.03,by=0.005),labels=seq(0.005,0.03,by=0.005),line=0.4,col=NA)
mtext(side=2,line=0,cex=1.1,text="Density - Actual Parents")
shade( density(P[,2]) , lim= as.vector(HPDI(P[,2], prob=0.999999)) , col = col.alpha("orange", 0.5)) #Natal actual
shade( density(P[,4]) , lim= as.vector(HPDI(P[,4], prob=0.999999)) , col = col.alpha("cyan3", 0.5)) #immigrant actual
shade( density(P[,6]) , lim= as.vector(HPDI(P[,6], prob=0.999999)) , col = col.alpha("red", 0.5)) #unknown actual
abline(v=median(P[,2]) , lty=1,lwd=1.5)
abline(v=median(P[,4]) , lty=2,lwd=1.5)
abline(v=median(P[,6]) , lty=3,lwd=1.5)
legend(x=0.022,y=400, legend = c("Natal", "Immigrant", "Unknown"), 
       col=c(col.alpha("orange", 0.5), col.alpha("cyan3", 0.5), col.alpha("red", 0.5)) , 
       pch=c(15,15,15),
       pt.cex=c(5,5,5) , bty="n", y.intersp=2,x.intersp=0.4, lty=c(0,0,0) ,cex=1.1, lwd=c(1,1))
legend(x=0.022,y=400, legend = c("","",""),
       col=1 , pch=c(NA,NA, NA),
       pt.cex=c(5,5,5), bty="n",y.intersp=2,x.intersp=0.4, lty=c(1,2,3) ,cex=1.1, lwd=c(1,1))
mtext(side=1,line=0,cex=1.1,text="Relatedness", outer=TRUE)

#Plot poster distributions with natality-type seperately
par(mfrow=c(3,1), mar=c(2,2,0.5,0.5), oma=c(2,0,0,0))

dens(P[,1], xlim=c(0.005,0.03) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE, 
     ylim=c(-11,400),mar=c(0,0,0,0), frame.plot=FALSE)
axis(1, cex=1,labels=NA, at=seq(0.005,0.03,by=0.005))
axis(1, cex.axis=1.3, at=seq(0.005,0.03,by=0.005),labels=seq(0.005,0.03,by=0.005),line=0.4,col=NA)
mtext(side=2,line=0,cex=1.1,text="Density - Natal Females")
shade( density(P[,1]) , lim= as.vector(HPDI(P[,1], prob=0.999999)) , col = col.alpha("forestgreen", 0.5)) #Natal potential
shade( density(P[,2]) , lim= as.vector(HPDI(P[,2], prob=0.999999)) , col = col.alpha("purple", 0.5)) #Natal actual
abline(v=median(P[,1]) , lty=1,lwd=1.5)
abline(v=median(P[,2]) , lty=2,lwd=1.5)

dens(P[,1], xlim=c(0.005,0.03) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE, 
     ylim=c(-11,400),mar=c(0,0,0,0), frame.plot=FALSE)
axis(1, cex=1,labels=NA, at=seq(0.005,0.03,by=0.005))
axis(1, cex.axis=1.3, at=seq(0.005,0.03,by=0.005),labels=seq(0.005,0.03,by=0.005),line=0.4,col=NA)
mtext(side=2,line=0,cex=1.1,text="Density - Immigrant Females")
shade( density(P[,3]) , lim= as.vector(HPDI(P[,3], prob=0.999999)) , col = col.alpha("forestgreen", 0.5)) #Immigrant potential
shade( density(P[,4]) , lim= as.vector(HPDI(P[,4], prob=0.999999)) , col = col.alpha("purple", 0.5)) #immigrant actual
abline(v=median(P[,3]) , lty=1,lwd=1.5)
abline(v=median(P[,4]) , lty=2,lwd=1.5)

dens(P[,1], xlim=c(0.005,0.03) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE, 
     ylim=c(-11,400),mar=c(0,0,0,0), frame.plot=FALSE)
axis(1, cex=1,labels=NA, at=seq(0.005,0.03,by=0.005))
axis(1, cex.axis=1.3, at=seq(0.005,0.03,by=0.005),labels=seq(0.005,0.03,by=0.005),line=0.4,col=NA)
mtext(side=2,line=0,cex=1.1,text="Density - Unknown Females")
shade( density(P[,5]) , lim= as.vector(HPDI(P[,5], prob=0.999999)) , col = col.alpha("forestgreen", 0.5)) #unknown potential
shade( density(P[,6]) , lim= as.vector(HPDI(P[,6], prob=0.999999)) , col = col.alpha("purple", 0.5)) #unknown actual
abline(v=median(P[,5]) , lty=1,lwd=1.5)
abline(v=median(P[,6]) , lty=2,lwd=1.5)

legend(x=0.023,y=350, legend = c("Actual", "Potential"), 
       col=c(col.alpha("purple", 0.5), col.alpha("forestgreen", 0.5)), 
       pch=c(15,15), pt.cex=c(5,5) , 
       bty="n", y.intersp=1,x.intersp=0.4, lty=c(0,0,0) ,cex=2, lwd=c(1,1))
legend(x=0.023,y=350, legend = c("",""),
       col=1 , pch=c(NA,NA), pt.cex=c(5,5), 
       bty="n",y.intersp=1,x.intersp=0.4, lty=c(2, 1) ,cex=2, lwd=c(1,1))
mtext(side=1,line=0.9,cex=1.1,text="Relatedness", outer=TRUE)


#Plot poster distributions all on one plot
par(mfrow=c(1,2), mar=c(4,2,0.5,0), oma=c(0,0,0,0))
layout(mat=matrix(c(1,1,2), ncol=3, nrow=1))
dens(P[,1], xlim=c(0.005,0.03) , xlab=NA , ylab=NA , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE, 
     ylim=c(-11,400),mar=c(0,0,0,0), frame.plot=FALSE)

axis(1, cex=1,labels=NA, at=seq(0.005,0.03,by=0.005))
axis(1, cex.axis=1.5,at=seq(0.005,0.03,by=0.005),labels=seq(0.005,0.03,by=0.005),line=0.4,col=NA)
mtext(side=1,line=2.5,cex=1.1,text="Relatedness")
mtext(side=2,line=0,cex=1.1,text="Density")

shade( density(P[,1]) , lim= as.vector(HPDI(P[,1], prob=0.999999)) , col = col.alpha("orange", 0.5)) #Natal potential
shade( density(P[,2]) , lim= as.vector(HPDI(P[,2], prob=0.999999)) , col = col.alpha("cyan3", 0.5)) #Natal actual
shade( density(P[,3]) , lim= as.vector(HPDI(P[,3], prob=0.999999)) , col = col.alpha("red", 0.5)) #Immigrant potential
shade( density(P[,4]) , lim= as.vector(HPDI(P[,4], prob=0.999999)) , col = col.alpha("slateblue", 0.5)) #immigrant actual
shade( density(P[,5]) , lim= as.vector(HPDI(P[,5], prob=0.999999)) , col = col.alpha("forestgreen", 0.5)) #unknown potential
shade( density(P[,6]) , lim= as.vector(HPDI(P[,6], prob=0.999999)) , col = col.alpha("deeppink", 0.5)) #unknown actual

abline(v=median(P[,1]) , lty=1,lwd=1.5)
abline(v=median(P[,2]) , lty=2,lwd=1.5)
abline(v=median(P[,3]) , lty=3,lwd=1.5)
abline(v=median(P[,4]) , lty=4,lwd=1.5)
abline(v=median(P[,5]) , lty=5,lwd=1.5)
abline(v=median(P[,6]) , lty=6,lwd=1.5)

par(mar=c(0, 0 ,0, 0))
plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), 
     xaxt='n' ,yaxt='n', frame.plot=FALSE)
legend(x=0,y=10, legend = c("Natal - Potential","Natal - Actual","Immigrant - Potential","Immigrant - Actual", "Unknown - Potential", "Unknown - Actual"), 
       col=c(col.alpha("orange", 0.5) ,col.alpha("cyan3", 0.5) ,col.alpha("red", 0.5),col.alpha("slateblue", 0.5),col.alpha("forestgreen", 0.5),col.alpha("deeppink", 0.5)) , 
       pch=c(15,15,15,15,15,15),
       pt.cex=c(5,5,5,5,5,5) , bty="n", y.intersp=1.5,x.intersp=0.2, lty=c(0,0,0,0) ,cex=1.6,lwd=c(1,1))
legend(x=0,y=10, legend = c("","","","", "",""),
       col=1 , pch=c(NA,NA),
       pt.cex=c(2,2,2,2,2,2) , bty="n",y.intersp=1.5,x.intersp=0.2, lty=c(1,2,3,4,5,6) ,cex=1.6,lwd=c(1,1))










