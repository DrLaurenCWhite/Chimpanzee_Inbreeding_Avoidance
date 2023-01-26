rm(list=ls())
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)

pairs=fread("RelatednessComparison_AllSNPS_SubsetMAF.csv")
pairsLD=fread("RelatednessComparison_LDPrunedSNPS_SubsetMAF.csv")
pairsAF=fread("RelatednessComparison_AllSNPS_AllIndsMAF.csv")
pairsLDAF=fread("RelatednessComparison_LDPrunedSNPS_AllIndsMAF.csv")

#Compare SNP set and different MAF calclations
fit1 <- lm(GeneticR ~ PedigreeR, data = pairs)
summary(fit1) #p<<0.001 AdjustedR^2=0.95.8

fit2 <- lm(GeneticR ~ PedigreeR, data = pairsLD)
summary(fit2) #pvalue=<<0.001 adjusted R^2=0.95.6

fit3 <- lm(GeneticR ~ PedigreeR, data = pairsAF)
summary(fit3) #pvalue=<<0.001 adjusted R^2=0.95.6

fit4 <- lm(GeneticR ~ PedigreeR, data = pairsLDAF)
summary(fit4) #pvalue=<<0.001 adjusted R^2=0.95.8


#Compare depth of coverage cutoffs
#Sample sizes
nrow(pairs)
nrow(pairs[AvDepth.a>=8 & AvDepth.b>=8])
nrow(pairs[AvDepth.a>=12 & AvDepth.b>=12])
nrow(pairs[AvDepth.a>=16 & AvDepth.b>=16])

fitq1 <- lm(GeneticR ~ PedigreeR, data = pairs) #Same as fit1
summary(fitq1) #p<<0.001 AdjustedR^2=0.95.8

fitq2 <- lm(GeneticR ~ PedigreeR, data = pairs[AvDepth.a>=8 & AvDepth.b>=8])
summary(fitq2) #pvalue=<<0.001 adjusted R^2=0.95.6

fitq3 <- lm(GeneticR ~ PedigreeR, data = pairs[AvDepth.a>=12 & AvDepth.b>=12])
summary(fitq3) #pvalue=<<0.001 adjusted R^2=0.95.1

fitq4 <- lm(GeneticR ~ PedigreeR, data = pairs[AvDepth.a>=16 & AvDepth.b>=16])
summary(fitq4) #pvalue=<<0.001 adjusted R^2=0.95.2


#PLOTS

#Plot comparisons of SNP set and different MAF calclations
#red line is 1:1, blue line is linear regression/correlation. adjusted R^2=0.9577
p1=ggplot(pairs, aes(x=PedigreeR, y=GeneticR)) + geom_jitter(alpha=0.5) + 
  geom_abline(slope=1, intercept=0, col="red", lwd=1) +
  annotate("text", label="Entire SNP set\nMAF calculated from subset", hjust=0, x=0.02, y=0.58, size=5) +
  stat_smooth(method="lm") + scale_y_continuous(limits=c(0.0, 0.6)) +
  xlab("Pedigree R") + ylab("Genetic R") + labs(col="Dyad\nCategory") +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

p2=ggplot(pairsLD, aes(x=PedigreeR, y=GeneticR)) + geom_jitter(alpha=0.5) + 
  annotate("text", label="LD pruned SNP set\nMAF calculated from subset", hjust=0, x=0.02, y=0.58, size=5) +
  geom_abline(slope=1, intercept=0, col="red", lwd=1) +
  stat_smooth(method="lm") + scale_y_continuous(limits=c(0.0, 0.6)) +
  xlab("Pedigree R") + ylab("Genetic R") + labs(col="Dyad\nCategory") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

p3=ggplot(pairsAF, aes(x=PedigreeR, y=GeneticR)) + geom_jitter(alpha=0.5) + 
  annotate("text", label="Entire SNP set\nMAF calculated from entire sample", hjust=0, x=0.02, y=0.58, size=5) +
  geom_abline(slope=1, intercept=0, col="red", lwd=1) +
  stat_smooth(method="lm") + scale_y_continuous(limits=c(0.0, 0.6)) +
  xlab("Pedigree R") + ylab("Genetic R") + labs(col="Dyad\nCategory") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

p4=ggplot(pairsLDAF, aes(x=PedigreeR, y=GeneticR)) + geom_jitter(alpha=0.5) + 
  annotate("text", label="LD pruned SNP set\nMAF calculated from entire sample", hjust=0, x=0.02, y=0.58, size=5) +
  geom_abline(slope=1, intercept=0, col="red", lwd=1) +
  stat_smooth(method="lm") + scale_y_continuous(limits=c(0.0, 0.6)) +
  xlab("Pedigree R") + ylab("Genetic R") + labs(col="Dyad\nCategory") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

ggarrange(p1, p2, p3, p4)



##Plot Average Depth check
q1=ggplot(pairs, aes(x=PedigreeR, y=GeneticR)) + geom_jitter(alpha=0.5) + 
  geom_abline(slope=1, intercept=0, col="red", lwd=1) +
  annotate("text", label="Average Depth of Coverage >4x", hjust=0, x=0.02, y=0.58, size=5) +
  stat_smooth(method="lm") + scale_y_continuous(limits=c(0.0, 0.6)) +
  xlab("Pedigree R") + ylab("Genetic R") + labs(col="Dyad\nCategory") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))
q2=ggplot(pairs[AvDepth.a>=8 & AvDepth.b>=8], aes(x=PedigreeR, y=GeneticR)) + geom_jitter(alpha=0.5) + 
  geom_abline(slope=1, intercept=0, col="red", lwd=1) +
  annotate("text", label="Average Depth of Coverage >8x", hjust=0, x=0.02, y=0.58, size=5) +
  stat_smooth(method="lm") + scale_y_continuous(limits=c(0.0, 0.6)) +
  xlab("Pedigree R") + ylab("Genetic R") + labs(col="Dyad\nCategory") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))
q3=ggplot(pairs[AvDepth.a>=12 & AvDepth.b>=12], aes(x=PedigreeR, y=GeneticR)) + geom_jitter(alpha=0.5) + 
  geom_abline(slope=1, intercept=0, col="red", lwd=1) +
  annotate("text", label="Average Depth of Coverage >12x", hjust=0, x=0.02, y=0.58, size=5) +
  stat_smooth(method="lm") + scale_y_continuous(limits=c(0.0, 0.6)) +
  xlab("Pedigree R") + ylab("Genetic R") + labs(col="Dyad\nCategory") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))
q4=ggplot(pairs[AvDepth.a>=16 & AvDepth.b>=16], aes(x=PedigreeR, y=GeneticR)) + geom_jitter(alpha=0.5) + 
  geom_abline(slope=1, intercept=0, col="red", lwd=1) +
  annotate("text", label="Average Depth of Coverage >16x", hjust=0, x=0.02, y=0.58, size=5) +
  stat_smooth(method="lm") + scale_y_continuous(limits=c(0.0, 0.6)) +
  xlab("Pedigree R") + ylab("Genetic R") + labs(col="Dyad\nCategory") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))


ggarrange(q1, q2, q3, q4)


