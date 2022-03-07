
###########################################################################
###########################################################################
###                                                                     ###
###         GXE AND STABILITY ANALYSIS USING METAN PACKAGE IN R         ###
###                                                                     ###
###########################################################################
###########################################################################



# Speaker: Oluwafemi Oyedele - MSc Candidate at University of Ibadan
# Email: oluwafemioyedele908@gmail.com

#-------------------------------------------------------------
rm(list = ls())
library(readxl)


data <- read.csv("Lecture.csv",header = T)
head(data)
str(data)
data$Hybrid <-  as.factor(data$Hybrid)
data$ME     <-  as.factor(data$ME)
data$Rep    <-  as.factor(data$Rep)
data$Env    <-  as.factor(data$Env)

str(data)


#Checking for NAs
colSums(is.na(data))

# Removing NA (trait dependent)
#install.packages("tidyverse")
library(tidyverse)


data.clean <- data %>% drop_na("Yield")

##### Exploratory analysis #####

summary(data.clean)
hist(data.clean$Yield, main = "Yield", xlab = "t/ha", col="lightblue")


##### Statistical analysis #####
#install.packages("lme4")
#install.packages("lmerTest")
library(lme4)
library(lmerTest)

# Hybrid as fixed
model <- lmer(Yield ~ Hybrid + (1|Env) + (1|Hybrid:Env) + (1|Rep:Env), data=data.clean)

anova(model,type = "III") #Wald test for fixed effects

#####################
##### Mean test ##### (Note: When hybrids are fix)
#####################

# Extracting BLUES
#summary(emmeans::emmeans(model, "Hybrid"))
#install.packages("emmeans")
library(emmeans)
BLUES <- emmeans(model, "Hybrid", lmer.df = "satterthwaite")

BLUES
# Mean separation
#install.packages("multcompView")
library(multcompView)

marginal = emmeans(model, "Hybrid")
multcomp::cld(marginal,   # compact letter displays (not recommended) # use pwpm() or pwpp()
              reversed = T,
              alpha=0.05,
              Letters=letters,      ### Use lower-case letters for .group
              adjust="none")        ### adjust='none'(Fisher LSD); adjust='Tukey'(Tukey HSD); adjust='bon'(Bonferroni correction)

# Variances for heritability (Quality of the experiment)
model1 <- lmer(Yield ~ (1|Hybrid) + (1|Env) + (1|Hybrid:Env) + (1|Env:Rep), data=data.clean)

ranova(model1) ## Anova-like table of random-effect terms using likelihood ratio tests
summary(model1)

# Normality test of the residuals (h0 = normal, p-value > 0.05 "fail to reject h0", indication of errors are normal)
# REMEMBER: if P-VALUE IS LOW... NULL MUST GO!
shapiro.test(resid(model1))
library(fBasics)
qqnormPlot(resid(model1))

#--------- Parameters for Environment ---------#
# Extracting Variance Components and mean

#Repeatability & CVs

varcomp <- VarCorr(model1)
blup <- ranef(model1, condVar = TRUE) # "SCA"
VG   <- attr(varcomp$Hybrid, "stddev")^2
VGE  <- attr(varcomp$"Hybrid:Env", "stddev")^2
VE   <- attr(varcomp, "sc")^2
mu   <- mean(coef(model1)$Hybrid[,,2])
nR   <- length(unique(data.clean$Rep))
nEnv <- length(unique(data.clean$Env))# 2 Environment!

#Repeatability & CVs
h2  <- VG/(VG+(VGE/nEnv)+(VE/(nR*nEnv)))
CVe <- sqrt(VE)/mu
CVg <- sqrt(VG)/mu
CVr <- CVg/CVe
print(matrix(c(mu,VG,VGE,VE,CVe,CVg,CVr,h2)))

############
# Exercise #
############

#Do the enalysis for each environment.
#Change the trait (PHT, DA)
#Compare the outcomes -> across loc x combined
#Does the heritability magnitude increase or decrease for the combined analysis?

#--------------- Filter for each env or ME ---------------
#install.packages("dplyr")
library(dplyr)
data.clean <- data.clean %>%
  filter(Env %in% c("E1"))
model_env <- lmer(Yield ~ (1|Hybrid) + (1|Rep), data=data.clean)
#---------------------------------------------

##################################
########## GGE Analysis ##########
##################################
#install.packages("gge")
#install.packages("agridat")
library(agridat)
library(gge)
#center = TRUE, environment-centered (G+GE)
#scale = If TRUE, scale values for each environment
#hull = If TRUE, show a which-won-where polygon.
#origin	= If "auto", the plotting window is centered on genotypes, otherwise the origin is at the middle of the window.


GGE <- gge(data.clean, Yield~Hybrid*Env, center = T, scale=T)
plot(GGE)
biplot(GGE, main="GxE Interaction",scale=T,
       cex.gen = 0.8,
       cex.env = 1.3,
       col.gen = "grey50",
       col.env = "orange4",
       zoom.gen = 1,
       comps=(1:2),flip=c(0,1), origin=0, hull=T)

GGB <- gge(data.clean, Yield~Hybrid:Env,
           env.group=ME, center=T,scale=T,ggb = T)
plot(GGB)
biplot(GGB, main ="Genotype-by-Block-of-Environment", #"GGB biplot
       cex.gen = 0.8,
       cex.env = 1.3,
       col.gen = "grey50",
       col.env = c("darkred","purple"),
       comps=(1:2),flip=c(2,1), origin=0, #origin If "auto", the plotting window is centered on genotypes, otherwise the origin is at the middle of the window.
       zoom.gen = .98)

#--------------- Filter for each ME ----------
#install.packages("dplyr")
library(dplyr)
data.clean2 <- data %>%
  filter(Env %in% c("E1","E3"))
#---------------------------------------------

GGB <- gge(data.clean2, Yield~Hybrid:Env,
           env.group=ME, center=T,scale=T,ggb = T)
plot(GGB)
biplot(GGB, main ="Genotype-by-Block-of-Environment", #"Mean vs stability biplot" for ME
       cex.gen = 0.8,
       cex.env = 1.3,
       col.gen = "grey50",
       col.env = c("darkred","purple"),
       comps=(1:2),flip=c(2,1), origin=0, #origin If "auto", the plotting window is centered on genotypes, otherwise the origin is at the middle of the window.
       zoom.gen = .98)

#--------------- Different data set less genotypes ------------
# Example 1.  Data is a data.frame in 'matrix' format
B <- matrix(c(50, 67, 90, 98, 120,
              55, 71, 93, 102, 129,
              65, 76, 95, 105, 134,
              50, 80, 102, 130, 138,
              60, 82, 97, 135, 151,
              65, 89, 106, 137, 153,
              75, 95, 117, 133, 155), ncol=5, byrow=TRUE)
rownames(B) <- c("G1","G2","G3","G4","G5","G6","G7")
colnames(B) <- c("E1","E2","E3","E4","E5")
m1 = gge(B)
plot(m1)
biplot(m1, main="Example biplot")

#Multi-environment: 18 wheat genotypes across 25 locations
data(crossa.wheat)
dat2 <- crossa.wheat
m2 <- gge(dat2, yield~gen*loc, env.group=locgroup, scale=FALSE)
plot(m2)
biplot(m2, lab.env=TRUE, main="crossa.wheat")
biplot3d(m2) #requires 3+ env.

################################
##### Stability Analysis ####### Should be done within each ME (Yan & tinker, 2006)
################################
#install.packages("agricolae")
library(agricolae)
TAMMI.mod<-with(data.clean,AMMI(ENV = Env, GEN =  Hybrid, REP = Rep, Y = Yield, MSE = 0,console=F,PC=T))

TAMMI.mod$ANOVA
TAMMI.mod$analysis #IPC & F-test

plot.AMMI(TAMMI.mod,first=1,second=2,third=3,
          type=1,number=F, # type=1 produce graphs biplot. type=2 produce graphs triplot
          gcol="lightgrey",
          ecol="darkred",
          angle=25,lwd=1.5,length=0.1,
          xlab=NULL,ylab=NULL,
          xlim=NULL,ylim=NULL)# biplot PC2 vs PC1 (AMMI2)

plot(TAMMI.mod,0,1,gcol="lightgrey",
     ecol="darkred",number=T) ## plot PC1 vs Yield (AMMI1)

library(ammistability)
# FA.AMMI = Equivalent to Wrick's ecovalance when all PCA are considered (default)
FA.AMMI(TAMMI.mod) # |ssi.method default = farshadfar|alternive ssi.method = "rao"

out <- ammistability(TAMMI.mod, AMGE = TRUE, ASI = TRUE, ASV = TRUE, ASTAB = TRUE,
                     AVAMGE = TRUE, DA = TRUE, DZ = TRUE, EV = TRUE,
                     FA = TRUE, MASI = TRUE, MASV = TRUE, SIPC = TRUE,
                     ZA = TRUE)
out$`SP Correlogram`

# Additional packages for GGE and Stability analysis
#install.packages(metan)
library("metan")
#https://github.com/TiagoOlivoto/metan
#VIGNETTE: https://tiagoolivoto.github.io/metan/articles/vignettes_gge.html

#---------- Model specification ----------# (Yan et al. 2007; Yan and Kang 2003)
###### Centering methods available
#0 or "none" for no centering;
#1 or "global" for global centered (E+G+GE);
#2 or "environment" (default), for environment-centered (G+GE);
#3 or "double" for double centred (GE). A biplot cannot be produced with models produced without centering.

#Scaling methods available
#0 or "none" (default) for no scaling;
#1 or "sd" where each value is divided by the standard deviation of its corresponding environment (column). This will put all testers roughly the same range of values.

###### Singular Value Partitioning methods available (svp)
#1 or "genotype" The singular value is entirely partitioned into the genotype eigenvectors, also called row metric preserving;
#2 or "environment" (default) the singular value is entirely partitioned into the environment eigenvectors, also called column metric preserving;
#3 or "symmetrical" The singular value is symmetrically partitioned into the genotype and the environment eigenvectors This SVP is most often used in AMMI analysis and other biplot analysis, but it is not ideal for visualizing either the relationship among entries or that among the testers.

###### Visualizing the Biplot (Yan and Kang (2003))
#type = 1 A basic biplot.
#type = 2 Mean performance vs. stability.
#type = 3 Which-won-where.
#type = 4 Discriminativeness vs. representativeness.
#type = 5 Examine an environment.
#type = 6 Ranking environments.
#type = 7 Examine a genotype.
#type = 8 Ranking genotypes
#type = 9 Compare two genotypes.
#type = 10 Relationship among environments.

model <- gge(.data = dat2, #Wheat example
             env = loc,
             gen = gen,
             resp = yield,
             centering = "environment",
             svp = "genotype",
             axis.expand = 1.5)
plot(model)

meanXstab <- plot(model, type = 2, #should be done within ME
                  cex = 0.4,
                  col.gen = "Blue",
                  col.env = "darkred",
                  size.text.gen = 2.2,
                  size.text.env = 4,
                  plot_theme = theme_metan_minimal())
meanXstab

ExamineEnv <- plot(model, type = 6, #should be done within ME
                   cex = 0.4,
                   col.gen = "Blue",
                   col.env = "darkred",
                   size.text.gen = 2.2,
                   size.text.env = 4,
                   plot_theme = theme_metan_minimal())
ExamineEnv

############
# Exercise #
############
# Use metan package to analyze the "data.clean" dataset

#There are addicional ways to assess GxE and Stability
#assuming random effects -> Factor analytics models
#Bayesian models -> Bayesian Finley-Wilkinson Regression

# Added agricolae stability page

https://myaseen208.github.io/agricolae/articles/StabilityAnalysis.html