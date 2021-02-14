# Environment Preparation and Data Import
## Import of relevant Libraries.

library(readr)
library(tidyverse)
library(readxl)
library(meta)
library(metafor)
library(esc)
library(dplyr)
library(meta)
library(metafor)
library(metaviz)
library(onehot)
library(mltools)
library(caret)
library(plyr)
library(ggplot2)
library(devtools)
library(dmetar)
library(car)
library(rmarkdown)
library(writexl)
library(rmdformats)
library(extrafont)
library(forestplot)
library(rmeta)
library(jtools)
library(PerformanceAnalytics)
library(psych)
library(GGally)
font_install("fontcm")

if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("MathiasHarrer/dmetar")
remotes::install_github("crsh/papaja")

options(scipen=999)


# Generate functions to be used throughout the analysis.
## Function to automatically calculate effect sizes.

effectsizecalculator <-  function (i) {
  esc_mean_sd(grp1m = data$minttc[effectsizes$id_ctrl[i] == data$id],
              grp1sd = data$minttc_sd[effectsizes$id_ctrl[i] == data$id],
              grp1n =  data$sample_size[effectsizes$id_ctrl[i] == data$id],
              grp2m = data$minttc[effectsizes$id_exp[i] == data$id],
              grp2sd = data$minttc_sd[effectsizes$id_exp[i] == data$id],
              grp2n = data$sample_size[effectsizes$id_exp[i] == data$id],
              es.type = 'g',
              study = studies$authors_meta[effectsizes$study_id[i] == studies$id])
}

## Function to locate authors to studies
authors <- function (a) {
  studies$authors[effectsizes$study_id[a] == studies$id]
}

## Function to calculate variance distribution
mlm.variance.distribution = function(x){
  
  m = x
  
  # Check class
  if (!(class(m)[1] %in% c("rma.mv", "rma"))){
    stop("x must be of class 'rma.mv'.")
  }
  
  # Check for three level model
  if (m$sigma2s != 2){
    stop("The model you provided does not seem to be a three-level model. This function can only be used for three-level models.")
  }
  
  # Get variance diagonal and calculate total variance
  n = m$k.eff
  vector.inv.var = 1/(diag(m$V))
  sum.inv.var = sum(vector.inv.var)
  sum.sq.inv.var = (sum.inv.var)^2
  vector.inv.var.sq = 1/(diag(m$V)^2)
  sum.inv.var.sq = sum(vector.inv.var.sq)
  num = (n-1)*sum.inv.var
  den = sum.sq.inv.var - sum.inv.var.sq
  est.samp.var = num/den
  
  # Calculate variance proportions
  level1=((est.samp.var)/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
  level2=((m$sigma2[2])/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
  level3=((m$sigma2[1])/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
  
  # Prepare df for return
  Level=c("Level 1", "Level 2", "Level 3")
  Variance=c(level1, level2, level3)
  df.res=data.frame(Variance)
  colnames(df.res) = c("% of total variance")
  rownames(df.res) = Level
  I2 = c("---", round(Variance[2:3], 2))
  df.res = as.data.frame(cbind(df.res, I2))
  
  totalI2 = Variance[2] + Variance[3]
  
  
  # Generate plot
  df1 = data.frame("Level" = c("Sampling Error", "Total Heterogeneity"),
                   "Variance" = c(df.res[1,1], df.res[2,1]+df.res[3,1]),
                   "Type" = rep(1,2))
  
  df2 = data.frame("Level" = rownames(df.res),
                   "Variance" = df.res[,1],
                   "Type" = rep(2,3))
  
  df = as.data.frame(rbind(df1, df2))
  
  
  g = ggplot(df, aes(fill=Level, y=Variance, x=as.factor(Type))) +
    coord_cartesian(ylim = c(0,1), clip = "off") +
    geom_bar(stat="identity", position="fill", width = 1, color="black") +
    scale_y_continuous(labels = scales::percent)+
    theme(axis.title.x=element_blank(),
          axis.text.y = element_text(color="black"),
          axis.line.y = element_blank(),
          axis.title.y=element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_line(lineend = "round"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.background = element_rect(linetype="solid",
                                           colour ="black"),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"cm"),
          axis.ticks.length=unit(.25, "cm"),
          plot.margin = unit(c(1,3,1,1), "lines")) +
    scale_fill_manual(values = c("darkseagreen3", "deepskyblue3", "darkseagreen2",
                                 "deepskyblue1", "deepskyblue2")) +
    
    # Add Annotation
    
    # Total Variance
    annotate("text", x = 1.5, y = 1.05,
             label = paste("Total Variance:",
                           round(m$sigma2[1]+m$sigma2[2]+est.samp.var, 3))) +
    
    # Sampling Error
    annotate("text", x = 1, y = (df[1,2]/2+df[2,2])/100,
             label = paste("Sampling Error Variance: \n", round(est.samp.var, 3)), size = 3) +
    
    # Total I2
    annotate("text", x = 1, y = ((df[2,2])/100)/2-0.02,
             label = bquote("Total"~italic(I)^2*":"~.(round(df[2,2],2))*"%"), size = 3) +
    annotate("text", x = 1, y = ((df[2,2])/100)/2+0.05,
             label = paste("Variance not attributable \n to sampling error: \n", round(m$sigma2[1]+m$sigma2[2],3)), size = 3) +
    
    # Level 1
    annotate("text", x = 2, y = (df[1,2]/2+df[2,2])/100, label = paste("Level 1: \n",
                                                                       round(df$Variance[3],2), "%", sep=""), size = 3) +
    
    # Level 2
    annotate("text", x = 2, y = (df[5,2]+(df[4,2]/2))/100,
             label = bquote(italic(I)[Level2]^2*":"~.(round(df[4,2],2))*"%"), size = 3) +
    
    # Level 3
    annotate("text", x = 2, y = (df[5,2]/2)/100,
             label = bquote(italic(I)[Level3]^2*":"~.(round(df[5,2],2))*"%"), size = 3)
  
  print(df.res)
  cat("Total I2: ", round(totalI2, 2), "% \n", sep="")
  suppressWarnings(print(g))
  invisible(df.res)
}

## Automatically spot outliers
spot.outliers.random<-function(data){
  data<-data
  Author<-data$studlab
  lowerci<-data$lower
  upperci<-data$upper 
  m.outliers<-data.frame(Author,lowerci,upperci) 
  te.lower<-data$lower.random 
  te.upper<-data$upper.random 
  dplyr::filter(m.outliers,upperci < te.lower) 
  dplyr::filter(m.outliers,lowerci > te.upper)
}

## Plot hat values
hat.plot <-function(fit) {
  p <- length(coefficients(fit))
  n <- length(fitted(fit))
  plot(hatvalues(fit), main="Index Plot of Hat Values")
  abline(h=c(2,3)*p/n, col="red", lty=2)
  identify(1:n, hatvalues(fit), names(hatvalues(fit)))
}


# Import files.
studies <- read_excel('~/Documents/1TU Dresden/4. Semester/Master-Thesis/Automatisiertes Fahren - Metaanalyse/00 Arbeitsmaterialien/Daten_final/studies.xlsx')
data <- read_excel('~/Documents/1TU Dresden/4. Semester/Master-Thesis/Automatisiertes Fahren - Metaanalyse/00 Arbeitsmaterialien/Daten_final/dataset.xlsx')
data_num <- read_excel('~/Documents/1TU Dresden/4. Semester/Master-Thesis/Automatisiertes Fahren - Metaanalyse/00 Arbeitsmaterialien/Daten_final/dataset_num.xlsx')
effectsizes <- read_excel('~/Documents/1TU Dresden/4. Semester/Master-Thesis/Automatisiertes Fahren - Metaanalyse/00 Arbeitsmaterialien/Daten_final/effsize_calc.xlsx')

head(studies)
head(data)
head(data_num)
head(effectsizes)
```



# Data Preprocessing

##Change datatypes.
data$minttc <- sapply(data$minttc, as.numeric)
data$minttc_sd <- sapply(data$minttc_sd, as.numeric)
data$effectsize_id <- sapply(data$effectsize_id, as.numeric)

##Convert files to dataframe.
data <-  as.data.frame(data)
studies <- as.data.frame(studies)
effectsizes <- as.data.frame(effectsizes)

##Convert moderators to factor.
data$sim <- factor(data$sim)
data$lad <-factor(data$lad)
data$ndt_v <- factor(data$ndt_v)
data$ndt_a <- factor(data$ndt_a)
data$ndt_m <- factor(data$ndt_m)
data$ndt_c <- factor(data$ndt_c)
data$hand <- factor(data$hand)
data$ndt_p <- factor(data$ndt_p)
data$tor_p <- factor(data$tor_p)
data$dre <- factor(data$dre)
data$iru <-factor(data$iru)
data$urg <- factor(data$urg)


# Effectsize Calculation
## Calculate effect sizes for the studies.
for (i in 1:35) {
  effectsizes[i, 5:15] <- effectsizecalculator(i)
}


## Change effectsizes columns to numeric.
effectsizes$effect_size <- sapply(effectsizes$effect_size, as.numeric)
effectsizes$standard_error <- sapply(effectsizes$standard_error, as.numeric)
effectsizes$variance <- sapply(effectsizes$variance,as.numeric)
effectsizes$lower_CI <- sapply(effectsizes$lower_CI, as.numeric)
effectsizes$upper_CI <-  sapply(effectsizes$upper_CI, as.numeric)
effectsizes$weight <- sapply(effectsizes$weight, as.numeric)


## Merge effectsizs with studies for details.
effectsizes <- merge(effectsizes, studies, by.x = "study_id", by.y = "id")

## Merge effectsizes with data for moderators.
effectsizes <- merge(effectsizes, data, by.x= "id_exp", by.y = "id")

## Remove unnecessary columns from dataframe.
drop <- c("effectsize_id.y","authors_meta","authors","id_exp","id_ctrl","info","es","studynr","exp_cond","k/e","subject_group","")
effectsizes = effectsizes[,!(names(effectsizes) %in% drop)]

## Rename column(s).
effectsizes <- rename(effectsizes, c("effectsize_id.x" = "effectsize_id"))

## Show output of new effectsizes dataframe.
str(effectsizes)

## Print a copy of the effectsizes dataframe.
write_xlsx(effectsizes,"/Users/timcerta/Documents/Programming/R/Projekte/automated driving_metaanalysis/effectsizes.xlsx")



# Descriptive Statistics

## Descriptive statstics summary 
summary(data)

## Histograms
### Histogram of minTTC
ggplot(data=data, aes(minttc)) + 
  geom_histogram(breaks=seq(0, 15, by=0.4),
                 col="black", 
                 fill="black", 
                 alpha = .2,
                 binwidth = 5) + 
  labs(title="Histogram of minTTC") +
  labs(x="minTTC (s)", y="frequency") + 
  theme_apa()

### Histogram of minTTC standard deviation
ggplot(data=data, aes(minttc_sd)) + 
  geom_histogram(breaks=seq(0, 10, by=0.3),
                 col="black", 
                 fill="black", 
                 alpha = .2) + 
  labs(title="Histogram of minTTC standard deviance") +
  labs(x="standard deviance of minTTC (s)", y="frequency") + 
  theme_apa()

### Histogram of effect sizes
ggplot(data=effectsizes, aes(x=effect_size)) + 
  geom_histogram(breaks=seq(-4, 4, by=0.2),
                 col="black", 
                 fill="black", 
                 alpha = .2) + 
  labs(title="Histogram of calculated effect sizes") +
  labs(x="effect size", y="frequency") + 
  theme_apa()

### Histogram of mTOT
ggplot(data=data, aes(x=mtot)) + 
  geom_histogram(breaks=seq(0, 5, by=0.2),
                 col="black", 
                 fill="black", 
                 alpha = .2) + 
  labs(title="Histogram of meanTOT") +
  labs(x="meanTOT (s)", y="frequency") + 
  theme_apa()

### Histogram of mTOT standard deviation
ggplot(data=data, aes(x=mtot_sd)) + 
  geom_histogram(breaks=seq(0, 5, by=0.2),
                 col="black", 
                 fill="black", 
                 alpha = .2) + 
  labs(title="Histogram of meanTOT standard deviance") +
  labs(x="standard deviance of meanTOT (s)", y="frequency") + 
  theme_apa()

### Histogram of age
ggplot(data=data, aes(age)) + 
  geom_histogram(breaks=seq(0, 75, by=1),
                 col="black", 
                 fill="black", 
                 alpha = .2) + 
  geom_density(alpha=.2, fill="#FF6666") +
  labs(title="Histogram of age") +
  labs(x="age", y="frequency") +
  theme_apa()

### Histogram of age standard deviation
ggplot(data=data, aes(age_sd)) + 
  geom_histogram(breaks=seq(0, 25, by=1),
                 col="black", 
                 fill="black", 
                 alpha = .2) + 
  geom_density(alpha=.2, fill="#FF6666") +
  labs(title="Histogram for age_sd") +
  labs(x="age_sd", y="frequency") +
  theme_apa()

### Histogram of time budget to collision
ggplot(data=data, aes(tbtc)) + 
  geom_histogram(breaks=seq(0, 25, by=1),
                 col="black", 
                 fill="black", 
                 alpha = .2) + 
  geom_density(alpha=.2, fill="#FF6666") +
  labs(title="Histogram for tbtc") +
  labs(x="time budget to collision", y="frequency") +
  theme_apa()


## Correlation matrices 
### Pearson correlation
vcols <- data_num %>% select(minttc, ndt_v, ndt_a, ndt_m, ndt_c, hand, ndt_p, dre, iru, urg)
var.cor <- cor(vcols, method = "pearson")
print(var.cor)

### Spearman correlation
var.cor <- cor(vcols, method = "spearman")
print(var.cor)

## Scatterplots
### Scatterplot of minTTC and minTTC standard deviation
ggplot(data=data, aes(minttc, minttc_sd)) +
  geom_point(alpha = 0.5,) + 
  xlab("minTTC (s)") + 
  ylab("minTTC (s) standard deviance") +
  theme_apa()

### Scatterplot of and mTOT and mTOT standard deviation
ggplot(data=data, aes(mtot, mtot_sd)) +
  geom_point(alpha = 0.5,) + 
  xlab("mTOT (s)") + 
  ylab("mTOT (s) standard deviance") +
  theme_apa()

### Scatterplot of meanTOT and minTTC
ggplot(data=data, aes(minttc, mtot)) +
  geom_point(alpha = 0.5,) + 
  xlab("minTTC (s)") + 
  ylab("meanTOT (s)") +
  theme_apa()

### Scatterplot of meanTOT and minTTC with a fitted least squares line
ggplot(data,aes(y=minttc,x=mtot))+geom_point()+geom_smooth(method="lm") + 
  theme_apa()

### Scatterplot of TBTB and minTTC
ggplot(data=data, aes(tbtc, minttc)) +
  geom_point(alpha = 0.5) +
  theme_apa()

### Scatterplot of TBTC and minTTC with a fitted least squares line
ggplot(data,aes(y=minttc,x=tbtc))+geom_point()+geom_smooth(method="lm") +
  theme_apa()

## Scatterplot Matrix
###Scattterplot Matrix for continuous variables
mcols <- data %>% select(minttc, minttc_sd, mtot, mtot_sd, tbtc)

ggpairs(mcols, axisLabels="show",
        upper = list(continuous = wrap("cor", method="spearman", size = 3)), 
        lower = list(continuous = "points"))+ 
  theme_apa()



## Group Comparisons for the dependent variable minTTC
### Group comparison minTTC - presence of NDT
count_ndt_p <- data %>% count("ndt_p")
minttc_compar_ndt_p <- data %>% group_by(ndt_p) %>% summarise_at(vars(minttc, minttc_sd), list(mean))
freq_ndt_p <- merge(minttc_compar_ndt_p,count_ndt_p)
print(freq_ndt_p)

### Group comparison minTTC - visual task
count_ndt_v <- data %>% count("ndt_v")
minttc_compar_ndt_v <- data %>% group_by(ndt_v) %>% summarise_at(vars(minttc, minttc_sd), list(mean))
freq_ndt_v <- merge(minttc_compar_ndt_v,count_ndt_v)
print(freq_ndt_v)

### Group comparison minTTC - acoustic task
count_ndt_a <- data %>% count("ndt_a")
minttc_compar_ndt_a <- data %>% group_by(ndt_a) %>% summarise_at(vars(minttc, minttc_sd), list(mean))
freq_ndt_a <- merge(minttc_compar_ndt_a,count_ndt_a)
print(freq_ndt_a)

### Group comparison minTTC - manual task
count_ndt_m <- data %>% count("ndt_m")
minttc_compar_ndt_m <- data %>% group_by(ndt_m) %>% summarise_at(vars(minttc, minttc_sd), list(mean))
freq_ndt_m <- merge(minttc_compar_ndt_m,count_ndt_m)
print(freq_ndt_m)

### Group comparison minTTC - cognitive task
count_ndt_c <- data %>% count("ndt_c")
minttc_compar_ndt_c <- data %>% group_by(ndt_c) %>% summarise_at(vars(minttc, minttc_sd), list(mean))
freq_ndt_c <- merge(minttc_compar_ndt_c,count_ndt_c)
print(freq_ndt_c)

### Group comparison minTTC - handheld task
count_hand <- data %>% count("hand")
minttc_compar_hand <- data %>% group_by(hand) %>% summarise_at(vars(minttc, minttc_sd), list(mean))
freq_hand <- merge(minttc_compar_hand,count_hand)
print(freq_hand)

### Group comparison minTTC - urgency(urg)
count_urg <- data %>% count("urg")
minttc_compar_urg <- data %>% group_by(urg) %>% summarise_at(vars(minttc, minttc_sd), list(mean))
freq_urg <- merge(minttc_compar_urg,count_urg)
print(freq_urg)

### Group comparison minTTC - driver response (dre)
count_dre <- data %>% count("dre")
minttc_compar_dre <- data %>% group_by(dre) %>% summarise_at(vars(minttc, minttc_sd), list(mean))
freq_dre <- merge(minttc_compar_dre,count_dre)
print(freq_dre)

### Group comparison minTTC - interaction with other road users (iru)
count_iru <- data %>% count("iru")
minttc_compar_iru <- data %>% group_by(iru) %>% summarise_at(vars(minttc, minttc_sd), list(mean))
freq_iru <- merge(minttc_compar_iru,count_iru)
print(freq_iru)


## Group comparisons mTOT
### Group Comparisons for mTOT - ndt_p
count_ndt_p <- data %>% count("ndt_p")
mtot_compar_ndt_p <- data %>% group_by(ndt_p) %>% summarise_at(vars(mtot, mtot_sd), list(mean), na.rm=TRUE)
freq_ndt_p <- merge(mtot_compar_ndt_p,count_ndt_p)
print(freq_ndt_p)

### Group comparison mTOT - visual task
count_ndt_v <- data %>% count("ndt_v")
mtot_compar_ndt_v <- data %>% group_by(ndt_v) %>% summarise_at(vars(mtot, mtot_sd), list(mean), na.rm=TRUE)
freq_ndt_v <- merge(mtot_compar_ndt_v,count_ndt_v)
print(freq_ndt_v)

### Group comparison mTOT - acoustic task
count_ndt_a <- data %>% count("ndt_a")
mtot_compar_ndt_a <- data %>% group_by(ndt_a) %>% summarise_at(vars(mtot, mtot_sd), list(mean), na.rm=TRUE)
freq_ndt_a <- merge(mtot_compar_ndt_a,count_ndt_a)
print(freq_ndt_a)

### Group comparison mTOT - manual task
count_ndt_m <- data %>% count("ndt_m")
mtot_compar_ndt_m <- data %>% group_by(ndt_m) %>% summarise_at(vars(mtot, mtot_sd), list(mean), na.rm=TRUE)
freq_ndt_m <- merge(mtot_compar_ndt_m,count_ndt_m)
print(freq_ndt_m)

### Group comparison mTOT - cognitive task
count_ndt_c <- data %>% count("ndt_c")
mtot_compar_ndt_c <- data %>% group_by(ndt_c) %>% summarise_at(vars(mtot, mtot_sd), list(mean), na.rm=TRUE)
freq_ndt_c <- merge(mtot_compar_ndt_c,count_ndt_c)
print(freq_ndt_c)

### Group comparison mTOT - handheld task
count_hand <- data %>% count("hand")
mtot_compar_hand <- data %>% group_by(hand) %>% summarise_at(vars(mtot, mtot_sd), list(mean), na.rm=TRUE)
freq_hand <- merge(mtot_compar_hand,count_hand)
print(freq_hand)

### Group comparison mTOT - urgency(urg)
count_urg <- data %>% count("urg")
mtot_compar_urg <- data %>% group_by(urg) %>% summarise_at(vars(mtot, mtot_sd), list(mean), na.rm=TRUE)
freq_urg <- merge(mtot_compar_urg,count_urg)
print(freq_urg)

### Group comparison mTOT - driver response (dre)
count_dre <- data %>% count("dre")
mtot_compar_dre <- data %>% group_by(dre) %>% summarise_at(vars(mtot, mtot_sd), list(mean), na.rm=TRUE)
freq_dre <- merge(mtot_compar_dre,count_dre)
print(freq_dre)

### Group comparison mTOT - interaction with other road users (iru)
count_iru <- data %>% count("iru")
mtot_compar_iru <- data %>% group_by(iru) %>% summarise_at(vars(mtot, mtot_sd), list(mean), na.rm=TRUE)
freq_iru <- merge(mtot_compar_iru,count_iru)
print(freq_iru)


## Mean differences
### Mean difference in minTTC (s) between "NDT" and "no NDT" condition 
means <- data %>% group_by(`ndt_p`) %>% summarise_at(vars(minttc, minttc_sd), list(mean))
means$upper <- means$minttc + (0.5 * means$minttc_sd)
means$lower <- means$minttc - (0.5 * means$minttc_sd)

compar_minttc <- data %>% group_by(`ndt_p`) %>% summarise_at(vars(minttc, minttc_sd), list(mean))
compar_minttc$upper <- compar_minttc$minttc + (0.5 * means$minttc_sd)
compar_minttc$lower <- compar_minttc$minttc - (0.5 * means$minttc_sd)
ggplot(compar_minttc, aes(x = `ndt_p`, y=`minttc`)) +
  geom_col(width=0.6, fill=c("grey", "white"), col="black") +
  geom_errorbar(data=compar_minttc, mapping=aes(x=`ndt_p`, ymin=upper, ymax=lower, width=0.1)) + 
  geom_text(aes(label = format(minttc, digits = 3), y = 0.7)) +
  scale_x_discrete(labels=c("NDT", "No NDT")) +
  xlab("Engaged vs not engaged in a NDT while driving highly automated") + 
  ylab("minTTC (s)") +
  annotate("text", x = 1.5, y=5, label = "mean difference = 1.09") +
  theme_apa()

### Mean difference in mTOT (s) between "NDT" and "no NDT" condition 
compar_mtot <- data %>% group_by(`ndt_p`) %>% summarise_at(vars(mtot, mtot_sd), list(mean), na.rm=TRUE)
compar_mtot$upper <- compar_mtot$mtot + (0.5 * compar_mtot$mtot_sd)
compar_mtot$lower <- compar_mtot$mtot - (0.5 * compar_mtot$mtot_sd)

ggplot(compar_mtot, aes(x = `ndt_p`, y=`mtot`)) +
  geom_col(width=0.6, fill=c("grey", "white"), col="black") +
  geom_errorbar(data=compar_mtot, mapping=aes(x=`ndt_p`, ymin=upper, ymax=lower, width=0.1)) + 
  geom_text(aes(label = format(mtot, digits = 3), y = 0.7)) +
  scale_x_discrete(labels=c("NDT", "No NDT")) +
  xlab("Engaged vs not engaged in a NDT while driving highly automated") + 
  ylab("mTOT (s)") +
  annotate("text", x = 1.5, y=2.8, label = "mean difference = 0.3") +
  theme_apa()

### Summary of the datasett
summary(data)


# Meta-analytc modeling


# Two-Level meta-analytic model
# Build a two-level model whose results then can be compared to the three-level model.
twolevel <- rma(effect_size, variance, data=effectsizes)
print(twolevel, digits=3)
confint(twolevel, digits=3)

# Three-Level Model
## Intercept-Only Model
##Estimate the overall effect by fitting an intercept-only model.
overall <- rma.mv(effect_size, variance, random = list(~ 1 | effectsize_id, ~ 1 | study_id), tdist=
                    TRUE, data=effectsizes)

## Request a print of the results stored in the object # ‘‘overall’’ in three digits.
summary(overall, digits=3)

## The estimated effect is the same as for the two-level model, but the confidence intervals differ.


## Variance Distribution
## Determine how the total variance is distributed over the three levels of the meta-analytic model.
## Print the results in percentages on screen.
n <- length(effectsizes$variance)
list.inverse.variances <- 1 / (effectsizes$variance) 
sum.inverse.variances <- sum(list.inverse.variances)
squared.sum.inverse.variances <- (sum.inverse.variances) ^ 2 
list.inverse.variances.square <- 1 / (effectsizes$variance^2) 
sum.inverse.variances.square <- sum(list.inverse.variances.square) 
numerator <- (n - 1) * sum.inverse.variances 
denominator <- squared.sum.inverse.variances - sum.inverse.variances.square
estimated.sampling.variance <- numerator / denominator
I2_1 <- (estimated.sampling.variance) / (overall$sigma2[1] + overall$sigma2[2] + estimated.sampling.variance)
I2_2 <- (overall$sigma2[1]) / (overall$sigma2[1] + overall$sigma2[2] + estimated.sampling.variance)
I2_3 <- (overall$sigma2[2]) / (overall$sigma2[1] + overall$sigma2[2] + estimated.sampling.variance)

amountvariancelevel1 <- I2_1 * 100
amountvariancelevel2 <- I2_2 * 100
amountvariancelevel3 <- I2_3 * 100

amountvariancelevel1
amountvariancelevel2
amountvariancelevel3

## Prepare df for return
Level=c("Level 1", "Level 2", "Level 3")
Variance=c(amountvariancelevel1, amountvariancelevel2, amountvariancelevel3)
df.res=data.frame(Variance)
colnames(df.res) = c("% of total variance")
rownames(df.res) = Level
I2 = c("---", round(Variance[2:3], 2))
df.res = as.data.frame(cbind(df.res, I2))

totalI2 = Variance[2] + Variance[3]


## Generate plot
df1 = data.frame("Level" = c("Sampling Error", "Total Heterogeneity"),
                 "Variance" = c(df.res[1,1], df.res[2,1]+df.res[3,1]),
                 "Type" = rep(1,2))

df2 = data.frame("Level" = rownames(df.res),
                 "Variance" = df.res[,1],
                 "Type" = rep(2,3))

df = as.data.frame(rbind(df1, df2))


g = ggplot(df, aes(fill=Level, y=Variance, x=as.factor(Type))) +
  coord_cartesian(ylim = c(0,1), clip = "off") +
  geom_bar(stat="identity", position="fill", width = 1, color="black") +
  scale_y_continuous(labels = scales::percent)+
  theme(axis.title.x=element_blank(),
        axis.text.y = element_text(color="black"),
        axis.line.y = element_blank(),
        axis.title.y=element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(lineend = "round"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.background = element_rect(linetype="solid",
                                         colour ="black"),
        legend.title = element_blank(),
        legend.key.size = unit(0.75,"cm"),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = unit(c(1,3,1,1), "lines")) +
  scale_fill_manual(values = c("darkseagreen3", "deepskyblue3", "darkseagreen2",
                               "deepskyblue1", "deepskyblue2")) +
  
  # Add Annotation
  
  # Total Variance
  annotate("text", x = 1.5, y = 1.05,
           label = paste("Total Variance:",
                         round(effectsizes$sigma2[1]+effectsizes$sigma2[2]+estimated.sampling.variance, 3))) +
  
  # Sampling Error
  annotate("text", x = 1, y = (df[1,2]/2+df[2,2])/100,
           label = paste("Sampling Error Variance: \n", round(estimated.sampling.variance, 3)), size = 3) +
  
  # Total I2
  annotate("text", x = 1, y = ((df[2,2])/100)/2-0.02,
           label = bquote("Total"~italic(I)^2*":"~.(round(df[2,2],2))*"%"), size = 3) +
  annotate("text", x = 1, y = ((df[2,2])/100)/2+0.05,
           label = paste("Variance not attributable \n to sampling error: \n", round(effectsizes$sigma2[1]+effectsizes$sigma2[2],3)), size = 3) +
  
  # Level 1
  annotate("text", x = 2, y = (df[1,2]/2+df[2,2])/100, label = paste("Level 1: \n",
                                                                     round(df$Variance[3],2), "%", sep=""), size = 3) +
  
  # Level 2
  annotate("text", x = 2, y = (df[5,2]+(df[4,2]/2))/100,
           label = bquote(italic(I)[Level2]^2*":"~.(round(df[4,2],2))*"%"), size = 3) +
  
  # Level 3
  annotate("text", x = 2, y = (df[5,2]/2)/100,
           label = bquote(italic(I)[Level3]^2*":"~.(round(df[5,2],2))*"%"), size = 3)

print(df.res)
cat("Total I2: ", round(totalI2, 2), "% \n", sep="")
suppressWarnings(print(g))
invisible(df.res)

## Calculation of the Intraclass Correlation (ICC)
## ICC describes how strongly units in the same group resemble each other.
round(overall$sigma2[1] / sum(overall$sigma2), 3)

## ICC = 0 indicates, that the values within clusters (=studies) are not similar. Because of the underlying research question, this was to be expected.

## Profile Likelihood Plot
σ_1^2 (variance at the effectsize level) and σ_2^2 (variance at the study level) were fixed at different values. For each value of σ_1^2 and σ_2^2, the likelihood over the remaining model parameters, such as the fixed effects, was estimated. This means it was estimated how likely the values of these parameters are given the observed data. Less negative values of the logarithm indicate a higher likelihood than more negative values. 


par(mfrow=c(2,1))
profile1 <-profile(overall, sigma2=1)
profile2 <-profile(overall, sigma2=2)


## The likelihood are highest for the values of σ_1^2 and σ_2^2 that had been estimated in the original model. Since profiling for variance components is done for non-negative value of the variance, the line is highest at this point in the second plot. We can be “fairly confident” that our meta-analytic models could identify the variance components (Viechtbauer, 2017).

## Confidence Intervals for Variance
confint(overall)

## The confindence intervals are very wide. For 𝜎22 hence, while the most likely value is 0, other values above 0 cannot be rejected either.

## Residuals --- QQ-Plot
qqnorm(residuals(overall,type="pearson"),main="QQ plot: residuals")
qqline(residuals(overall,type="pearson"),col="red")

## Two-Level Model without within-study variance
Build a two-level model without within-study variance.
```{r}
modelnovar2 <- rma.mv(effect_size, variance, random = list(~ 1 | effectsize_id, ~ 1 | study_id),
                      sigma2=c(0,NA), tdist=TRUE, data=effectsizes)

## Request a print of the results stored in the object # ‘‘movelnovar2’’ in three digits.
summary(modelnovar2, digits=3)

## Perform a likelihood-ratio-test to determine the significance of the within-study (level 2) variance.
anova(overall,modelnovar2)

## The model does not perform better than the three-level model.

## Two-level model without between-study variance
Build a two-level model without between-study variance.

modelnovar3 <- rma.mv(effect_size, variance, random = list(~ 1 | effectsize_id, ~ 1 | study_id),
                      sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 

## Request a print of the results stored in the object # ‘‘modelnovar3’’ in three digits.
summary(modelnovar3, digits=3)

## Perform a likelihood-ratio-test to determine the significance of the between-study (level 3) variance.
anova(overall,modelnovar3)

## This model does perform very similar to the three-level model. A third model for modeling differences between studies is not needed.


# Moderator Analysis

## NDT Modality
## Test of different modalities independently
### acoustic NDT
ndt_a <- rma.mv(effect_size, variance, mods = ~ ndt_a, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(ndt_a, digits=3)

### manual NDT
ndt_m <- rma.mv(effect_size, variance, mods = ~ ndt_m, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(ndt_m, digits=3)

### cognitive NDT
ndt_c <- rma.mv(effect_size, variance, mods = ~ ndt_c, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(ndt_c, digits=3)

### visual NDT
ndt_cv<- rma.mv(effect_size, variance, mods = ~ ndt_v, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(ndt_cv digits=3)

### handheld NDT
ndt_hand <- rma.mv(effect_size, variance, mods = ~ hand, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(ndt_hand, digits=3)

### All task modalities in one model.
ndt_modality <- rma.mv(effect_size, variance, mods = ~ ndt_a + ndt_m + ndt_c + ndt_v + hand, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(ndt_modality, digits=3)

## Urgency

### One-hot encoding of urgency.
effectsizes <- effectsizes %>% 
  mutate(urg_low = ifelse(urg=='low', 1, 0),
         urg_medium = ifelse(urg=='medium', 1, 0),
         urg_high = ifelse(urg=='high', 1, 0))

### Determine the potential moderating effect of urgency as categorical moderator.
###Low urgency is chosen as the reference category.
urgency <- rma.mv(effect_size, variance, mods = ~ urg_medium + urg_high, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(urgency, digits=3)

### Determine the potential moderating effect of urgency as continuous moderator.
### High urgency is chosen as intercept category.
urgency_cont <- rma.mv(effect_size, variance, mods = ~ urg, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(urgency_cont, digits=3)

### Determine the potential moderating effect of time budget to collision.
tbtc <- rma.mv(effect_size, variance, mods = ~ tbtc, random = list(~ 1 | effectsize_id, ~1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(tbtc, digits=3)

## Exploratory Moderator Analyses
### Simulator Fidelity
sim_fidelity <- rma.mv(effect_size, variance, mods = ~ sim, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(sim_fidelity, digits=3)

### Driver response complexity
dre <- rma.mv(effect_size, variance, mods = ~ dre, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(dre, digits=3)

### Interaction with other road users
iru <- rma.mv(effect_size, variance, mods = ~ iru, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(iru, digits=3)

### Level of Automated Driving
lad <- rma.mv(effect_size, variance, mods = ~ lad, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes) 
summary(lad, digits=3)

# Forest Plot
Forest Plot for individual effect sizes.

## Sort effectsizes
effectsizes_sorted <- effectsizes[with(effectsizes, order(year, study, effect_size)),]
effectsize_num <- c(1:35)

## Forest Plot sorted
forest_sorted <- viz_forest(x = effectsizes_sorted[, c("effect_size", "standard_error")],
                            variant = "classic",
                            type = "standard",
                            study_table = effectsizes_sorted$study,
                            method = "REML",
                            annotate_CI = T,
                            table_headers = c("Study"),
                            xlab = "Hedges' g",
                            summary_table="Overall effect",
                            col = "darkblue",
                            summary_col = "firebrick",
                            table_layout = matrix(c(1, 2, 2, 3), nrow = 1))
print(forest_sorted)


# Publication Bias
## Visualization: Funnel Plot
## Simple funnel plot.
funnel(overall, xlab = "Hedges' g", main="Funnel Plot") 

## Contour-enhanced funnel plot.
viz_funnel(overall, method="REML", contours_col="Greys", xlab="Hedges' g")

## Contour-enhanced sunset funnel plot estimating the power of studies.
viz_sunset(overall, true_effect=.3, power_contours="continuous")

## Trim-and-fill method
## This is a non-parametric approach, where pseudo-studies are being imputed until the funnel plot symmetry is restored. 
trimfill <- trimfill(twolevel)
summary(trimfill)
funnel(trimfill, xlab = "Hedges' g", main="Funnel plot after imputing pseudo-studies") 

## Testing the results of the the trim-and-fill method for the three-level model against the two-level model
viz_funnel(overall,
           method="REML", 
           contours_col="Greys", 
           xlab="Hedges' g",
           trim_and_fill = TRUE,
           trim_and_fill_side = "left",
           egger = FALSE)

## No differences in the results can be detected.


## Egger Regression Symmetry Test 
## Test for asymmetry of the funnel plot. It tests for the Y intercept = 0 from a linear regression of normalized effect estimate.
egger <- rma.mv(effect_size, variance, mods =  standard_error, random = list(~ 1 | study_id, ~ 1 | effectsize_id), tdist= TRUE, data=effectsizes, method="REML")
summary(egger , digits=3)

## The p-value is > 0.05, which indicates non-existence of publication bias.

## Rank-Correlation Test 
## Correlates the standardized treatment effect with the variance of the treatment effect using Kendall's tau as the measure of association. Since it is not available for three-level meta-analysis, it was calculated with the two-level model.
ranktest(twolevel)

## The p-value is > 0.05, which indicates non-existence of publication bias.


# Outlier & influence analysis
## Leave-one-out analysis
##Shows, how much Hedge's g change after removing the studies individually.
leaveoneout <- viz_forest(x = effectsizes_sorted[, c("effect_size", "standard_error", "weight")],
                          variant = "classic",
                          type = "sensitivity",
                          study_table = effectsizes_sorted$study,
                          method = "REML",
                          annotate_CI = T,
                          table_headers = c("Study"),
                          xlab = "Hedges' g",
                          summary_table="Overall effect",
                          col = "darkblue",
                          summary_col = "firebrick",
                          table_layout = matrix(c(1, 2, 2, 3), nrow = 1))
print(leaveoneout)

## Cook's Distances 
## Result for each effect size or study shows how much the regression model would change when the value would be removed.

## Cook's Distances for each observed Outcome.
cook_ind <- cooks.distance(overall)
plot(cook_ind, type="o", pch=19, xlab="Effect sizes", ylab="Cook's Distance")



## Cook's Distances for each Cluster.
cook_cluster <- cooks.distance(overall, cluster=effectsizes$study_id)
plot(cook_cluster, type="o", pch=19, xlab="Studies", ylab="Cook's Distance", xaxt="n")
axis(side=1, at=seq_along(cook_cluster), labels=as.numeric(names(cook_cluster)))

## Hat values
## The hat values are the fitted values, the predictions made by the model for each observation.
hat <- hatvalues(overall)
print(hat)

## Plot values
plot(hat, type='o',
main="Hat Values",
col="black",
pch=20,
ylim = c(0, 0.030))


# Meta-analysis without outliers on effectsize level


# Two-Level Model without Outliers on effectsize level
## Removal of Outliers
effectsizes_filtered <- effectsizes %>% filter(between(effect_size, -1, 1.5))

## Build a two-level model whose robustness then can be compared to the three-level model.
twolevel_wo <- rma(effect_size, variance, data=effectsizes_filtered)
print(twolevel_wo, digits=3)
confint(twolevel_wo, digits=3)

# Three-Level-Model without Outliers
## Intercept-Only Model
## Estimate the overall effect by fitting an intercept-only model.
overall_wo <- rma.mv(effect_size, variance, random = list(~ 1 | effectsize_id, ~ 1 | study_id), tdist=TRUE, data=effectsizes_filtered)

## Request a print of the results stored in the object ‘‘overall’’ in three digits.
summary(overall_wo, digits=3)

## Variance Distribution
## Determine how the total variance is distributed over the three levels of the meta-analytic model.
## Print the results in percentages on screen. 

n <- length(effectsizes_filtered$variance)
list.inverse.variances <- 1 / (effectsizes_filtered$variance) 
sum.inverse.variances <- sum(list.inverse.variances)
squared.sum.inverse.variances <- (sum.inverse.variances) ^ 2 
list.inverse.variances.square <- 1 / (effectsizes_filtered$variance^2) 
sum.inverse.variances.square <- sum(list.inverse.variances.square) 
numerator <- (n - 1) * sum.inverse.variances 
denominator <- squared.sum.inverse.variances - sum.inverse.variances.square
estimated.sampling.variance <- numerator / denominator
I2_1 <- (estimated.sampling.variance) / (overall_wo$sigma2[1] + overall_wo$sigma2[2] + estimated.sampling.variance)
I2_2 <- (overall_wo$sigma2[1]) / (overall_wo$sigma2[1] + overall_wo$sigma2[2] + estimated.sampling.variance)
I2_3 <- (overall_wo$sigma2[2]) / (overall_wo$sigma2[1] + overall_wo$sigma2[2] + estimated.sampling.variance)

amountvariancelevel1 <- I2_1 * 100
amountvariancelevel2 <- I2_2 * 100
amountvariancelevel3 <- I2_3 * 100

amountvariancelevel1
amountvariancelevel2
amountvariancelevel3

# Prepare df for return
Level=c("Level 1", "Level 2", "Level 3")
Variance=c(amountvariancelevel1, amountvariancelevel2, amountvariancelevel3)
df.res=data.frame(Variance)
colnames(df.res) = c("% of total variance")
rownames(df.res) = Level
I2 = c("---", round(Variance[2:3], 2))
df.res = as.data.frame(cbind(df.res, I2))

totalI2 = Variance[2] + Variance[3]


# Generate plot
df1 = data.frame("Level" = c("Sampling Error", "Total Heterogeneity"),
                 "Variance" = c(df.res[1,1], df.res[2,1]+df.res[3,1]),
                 "Type" = rep(1,2))

df2 = data.frame("Level" = rownames(df.res),
                 "Variance" = df.res[,1],
                 "Type" = rep(2,3))

df = as.data.frame(rbind(df1, df2))


g = ggplot(df, aes(fill=Level, y=Variance, x=as.factor(Type))) +
  coord_cartesian(ylim = c(0,1), clip = "off") +
  geom_bar(stat="identity", position="fill", width = 1, color="black") +
  scale_y_continuous(labels = scales::percent)+
  theme(axis.title.x=element_blank(),
        axis.text.y = element_text(color="black"),
        axis.line.y = element_blank(),
        axis.title.y=element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(lineend = "round"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.background = element_rect(linetype="solid",
                                         colour ="black"),
        legend.title = element_blank(),
        legend.key.size = unit(0.75,"cm"),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = unit(c(1,3,1,1), "lines")) +
  scale_fill_manual(values = c("darkseagreen3", "deepskyblue3", "darkseagreen2",
                               "deepskyblue1", "deepskyblue2")) +
  
  # Add Annotation
  
  # Total Variance
  annotate("text", x = 1.5, y = 1.05,
           label = paste("Total Variance:",
                         round(effectsizes_filtered$sigma2[1]+effectsizes_filtered$sigma2[2]+estimated.sampling.variance, 3))) +
  
  # Sampling Error
  annotate("text", x = 1, y = (df[1,2]/2+df[2,2])/100,
           label = paste("Sampling Error Variance: \n", round(estimated.sampling.variance, 3)), size = 3) +
  
  # Total I2
  annotate("text", x = 1, y = ((df[2,2])/100)/2-0.02,
           label = bquote("Total"~italic(I)^2*":"~.(round(df[2,2],2))*"%"), size = 3) +
  annotate("text", x = 1, y = ((df[2,2])/100)/2+0.05,
           label = paste("Variance not attributable \n to sampling error: \n", round(effectsizes_filtered$sigma2[1]+effectsizes_filtered$sigma2[2],3)), size = 3) +
  
  # Level 1
  annotate("text", x = 2, y = (df[1,2]/2+df[2,2])/100, label = paste("Level 1: \n",
                                                                     round(df$Variance[3],2), "%", sep=""), size = 3) +
  
  # Level 2
  annotate("text", x = 2, y = (df[5,2]+(df[4,2]/2))/100,
           label = bquote(italic(I)[Level2]^2*":"~.(round(df[4,2],2))*"%"), size = 3) +
  
  # Level 3
  annotate("text", x = 2, y = (df[5,2]/2)/100,
           label = bquote(italic(I)[Level3]^2*":"~.(round(df[5,2],2))*"%"), size = 3)

print(df.res)
cat("Total I2: ", round(totalI2, 2), "% \n", sep="")
suppressWarnings(print(g))
invisible(df.res)



## Profile Likelihood Plots
par(mfrow=c(2,1))
profile1wo <- profile(overall_wo, sigma2=1)
profile2wo <- profile(overall_wo, sigma2=2)

## Confidence intervals
confint(overall_wo)

## Calculation of Intraclass Coefficient
round(overall_wo$sigma2[2] / sum(overall_wo$sigma2), 3)

## Total sum of heterogeneity
round(sum(overall_wo$sigma2), 3)


## Residuals - QQ-plot
qqnorm(residuals(overall_wo,type="pearson"),main="QQ plot: residuals")
qqline(residuals(overall_wo,type="pearson"),col="red")

## Two-Level Model without within-study variance
## Build a two-level model without within-study variance.
modelnovar2_wo <- rma.mv(effect_size, variance, random = list(~ 1 | effectsize_id, ~ 1 | study_id),
                      sigma2=c(0,NA), tdist=TRUE, data=effectsizes_filtered)

## Request a print of the results stored in the object # ‘‘movelnovar2_wo’’ in three digits.
summary(modelnovar2_wo, digits=3)

## Perform a likelihood-ratio-test to determine the significance of the within-study (level2) variance.
anova(overall_wo,modelnovar2_wo)

## Two-Level Model without between-study variance
Build a two-level model without between-study variance

modelnovar3_wo <- rma.mv(effect_size, variance, random = list(~ 1 | effectsize_id, ~ 1 | study_id),
                      sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered) 

## Request a print of the results stored in the object # ‘‘modelnovar3’’ in three digits.
summary(modelnovar3_wo, digits=3)

## Perform a likelihood-ratio-test to determine the significance of the between-study (level3) variance.
anova(overall_wo,modelnovar3_wo)


# Moderator Analysis without Outliers

## NDT Modality
## Determine the potential moderating effect of NDT modality.

##Test of different modalities independently.
### Acoustic NDT
ndt_a_wo <- rma.mv(effect_size, variance, mods = ~ ndt_a, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered) 
summary(ndt_a_wo, digits=3)

### Manual NDT
ndt_m_wo <- rma.mv(effect_size, variance, mods = ~ ndt_m, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered) 
summary(ndt_m_wo, digits=3)

### Cognitive NDT
ndt_c_wo <- rma.mv(effect_size, variance, mods = ~ ndt_c, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered) 
summary(ndt_c_wo, digits=3)

### Handheld NDT
ndt_hand_wo <- rma.mv(effect_size, variance, mods = ~ hand, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered) 
summary(ndt_hand_wo, digits=3)

### All modalities in one model.
ndt_modality_wo <- rma.mv(effect_size, variance, mods = ~ ndt_a + ndt_m + ndt_c + hand, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered) 
summary(ndt_modality_wo, digits=3)


## Urgency
### One-hot encoding of hand-coded urgency
effectsizes_filtered <- effectsizes_filtered %>% 
  mutate(urg_low = ifelse(urg=='low', 1, 0),
         urg_medium = ifelse(urg=='medium', 1, 0),
         urg_high = ifelse(urg=='high', 1, 0))

### Determine the potential moderating effect of hand-coded urgency as categorical moderator.
### No urgency is chosen as the reference category
urgency_wo <- rma.mv(effect_size, variance, mods = ~ urg_medium + urg_high, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered) 
summary(urgency_wo, digits=3)

### Determine the potential moderating effect of urgency as continuous moderator.
### High urgency is chosen as inttercept category
urgency_cont_wo <- rma.mv(effect_size, variance, mods = ~ urg, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered) 
summary(urgency_cont_wo, digits=3)

### Determine the potential moderating effect of time budget to collision.
tbtc_wo <- rma.mv(effect_size, variance, mods = ~ tbtc, random = list(~ 1 | effectsize_id, ~1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered) 
summary(tbtc_wo, digits=3)

## Further Moderator Analyses
### Simulator fidelity
sim_fidelity_wo <- rma.mv(effect_size, variance, mods = ~ sim, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered) 
summary(sim_fidelity_wo, digits=3)

### Interaction with other road users
iru_wo <- rma.mv(effect_size, standard_error, mods = ~ iru, random = list(~1 | effectsize_id, ~ 1 | study_id), tdist=TRUE, data=effectsizes_filtered) 
summary(iru_wo, digits=3)

### Level of automated driving
lad_wo <- rma.mv(effect_size, standard_error, mods = ~ lad, random = list(~1 | effectsize_id, ~ 1 | study_id), tdist=TRUE, data=effectsizes_filtered) 
summary(lad_wo, digits=3)


# Publication Bias without Outliers

## Visualization: Funnel Plot
funnel(overall_wo, main = "Standard Error", xlab = "Hedges' g")

## Trim-and-fill method
trimfill <- trimfill(twolevel_wo)
summary(trimfill)
funnel(trimfill,main = "Standard Error", xlab = "Hedges' g") 

## Egger Regression Symmetry Test 
##(eggers.test not working, bc 'overall' is no 'meta'-type. Regtest is carried out with the two-level model)

## eggers.test(x = overall)
regtest(twolevel_wo,  model="rma")

## The p-value is > 0.05, publication bias can be rejected.

## Another method for calculating Egger's regression test - same result though.
egger <- rma.mv(effect_size, variance, mods =  standard_error, random = list(~ 1 | study_id, ~ 1 | effectsize_id), tdist= TRUE, data=effectsizes_filtered, method="REML")
summary(egger , digits=3)

## Rank-Correlation Test 
ranktest(twolevel_wo)


# Two-Level Model without Lin et al., 2020
## Removal of Outliers
effectsizes_filtered2 <- effectsizes %>% filter(study_id != 8)

# Two-Level Model without Lin et al., 2020
# Build a two-level model whose robustness then can be compared to the three-level model.

twolevel_w8 <- rma(effect_size, variance, data=effectsizes_filtered2)
print(twolevel_w8, digits=3)
confint(twolevel_w8, digits=3)

# Three-Level-Model without Lin et al., 2020
## Intercept-Only Model
## Estimate the overall effect by fitting an intercept-only model.
overall_w8 <- rma.mv(effect_size, variance, random = list(~ 1 | effectsize_id, ~ 1 | study_id), tdist=TRUE, data=effectsizes_filtered2)

##mRequest a print of the results stored in the object ‘‘overall’’ in three digits.
summary(overall_w8, digits=3)

## Variance Distribution
## Determine how the total variance is distributed over the three levels of the meta-analytic model.
## Print the results in percentages on screen. 
n <- length(effectsizes_filtered2$variance)
list.inverse.variances <- 1 / (effectsizes_filtered2$variance) 
sum.inverse.variances <- sum(list.inverse.variances)
squared.sum.inverse.variances <- (sum.inverse.variances) ^ 2 
list.inverse.variances.square <- 1 / (effectsizes_filtered2$variance^2) 
sum.inverse.variances.square <- sum(list.inverse.variances.square) 
numerator <- (n - 1) * sum.inverse.variances 
denominator <- squared.sum.inverse.variances - sum.inverse.variances.square
estimated.sampling.variance <- numerator / denominator
I2_1 <- (estimated.sampling.variance) / (overall_w8$sigma2[1] + overall_w8$sigma2[2] + estimated.sampling.variance)
I2_2 <- (overall_w8$sigma2[1]) / (overall_w8$sigma2[1] + overall_w8$sigma2[2] + estimated.sampling.variance)
I2_3 <- (overall_w8$sigma2[2]) / (overall_w8$sigma2[1] + overall_w8$sigma2[2] + estimated.sampling.variance)

amountvariancelevel1 <- I2_1 * 100
amountvariancelevel2 <- I2_2 * 100
amountvariancelevel3 <- I2_3 * 100

amountvariancelevel1
amountvariancelevel2
amountvariancelevel3

# Prepare df for return
Level=c("Level 1", "Level 2", "Level 3")
Variance=c(amountvariancelevel1, amountvariancelevel2, amountvariancelevel3)
df.res=data.frame(Variance)
colnames(df.res) = c("% of total variance")
rownames(df.res) = Level
I2 = c("---", round(Variance[2:3], 2))
df.res = as.data.frame(cbind(df.res, I2))

totalI2 = Variance[2] + Variance[3]


# Generate plot
df1 = data.frame("Level" = c("Sampling Error", "Total Heterogeneity"),
                 "Variance" = c(df.res[1,1], df.res[2,1]+df.res[3,1]),
                 "Type" = rep(1,2))

df2 = data.frame("Level" = rownames(df.res),
                 "Variance" = df.res[,1],
                 "Type" = rep(2,3))

df = as.data.frame(rbind(df1, df2))


g = ggplot(df, aes(fill=Level, y=Variance, x=as.factor(Type))) +
  coord_cartesian(ylim = c(0,1), clip = "off") +
  geom_bar(stat="identity", position="fill", width = 1, color="black") +
  scale_y_continuous(labels = scales::percent)+
  theme(axis.title.x=element_blank(),
        axis.text.y = element_text(color="black"),
        axis.line.y = element_blank(),
        axis.title.y=element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(lineend = "round"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.background = element_rect(linetype="solid",
                                         colour ="black"),
        legend.title = element_blank(),
        legend.key.size = unit(0.75,"cm"),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = unit(c(1,3,1,1), "lines")) +
  scale_fill_manual(values = c("darkseagreen3", "deepskyblue3", "darkseagreen2",
                               "deepskyblue1", "deepskyblue2")) +
  
  # Add Annotation
  
  # Total Variance
  annotate("text", x = 1.5, y = 1.05,
           label = paste("Total Variance:",
                         round(effectsizes_filtered2$sigma2[1]+effectsizes_filtered2$sigma2[2]+estimated.sampling.variance, 3))) +
  
  # Sampling Error
  annotate("text", x = 1, y = (df[1,2]/2+df[2,2])/100,
           label = paste("Sampling Error Variance: \n", round(estimated.sampling.variance, 3)), size = 3) +
  
  # Total I2
  annotate("text", x = 1, y = ((df[2,2])/100)/2-0.02,
           label = bquote("Total"~italic(I)^2*":"~.(round(df[2,2],2))*"%"), size = 3) +
  annotate("text", x = 1, y = ((df[2,2])/100)/2+0.05,
           label = paste("Variance not attributable \n to sampling error: \n", round(effectsizes_filtered2$sigma2[1]+effectsizes_filtered2$sigma2[2],3)), size = 3) +
  
  # Level 1
  annotate("text", x = 2, y = (df[1,2]/2+df[2,2])/100, label = paste("Level 1: \n",
                                                                     round(df$Variance[3],2), "%", sep=""), size = 3) +
  
  # Level 2
  annotate("text", x = 2, y = (df[5,2]+(df[4,2]/2))/100,
           label = bquote(italic(I)[Level2]^2*":"~.(round(df[4,2],2))*"%"), size = 3) +
  
  # Level 3
  annotate("text", x = 2, y = (df[5,2]/2)/100,
           label = bquote(italic(I)[Level3]^2*":"~.(round(df[5,2],2))*"%"), size = 3)

print(df.res)
cat("Total I2: ", round(totalI2, 2), "% \n", sep="")
suppressWarnings(print(g))
invisible(df.res)


### Profile Likelihood Plots
par(mfrow=c(2,1))
profile1wo <- profile(overall_w8, sigma2=1)
profile2wo <- profile(overall_w8, sigma2=2)

### Confidence intervals
confint(overall_w8)

### Calculation of Intraclass Coefficient
round(overall_w8$sigma2[2] / sum(overall_w8$sigma2), 3)

### Total sum of heterogeneity
round(sum(overall_w8$sigma2), 3)

### Residuals # QQ-plot
qqnorm(residuals(overall_w8,type="pearson"),main="QQ plot: residuals")
qqline(residuals(overall_w8,type="pearson"),col="red")

## Two-Level Model without within-study variance
## Build a two-level model without within-study variance.
modelnovar2_w8 <- rma.mv(effect_size, variance, random = list(~ 1 | effectsize_id, ~ 1 | study_id),
                         sigma2=c(0,NA), tdist=TRUE, data=effectsizes_filtered2)

## Request a print of the results stored in the object # ‘‘movelnovar2’’ in three digits.
summary(modelnovar2_w8, digits=3)

## Perform a likelihood-ratio-test to determine the significance of the within-study (level2) variance.
anova(overall_w8,modelnovar2_w8)

## Two-Level Model without between-study variance
## Build a two-level model without between-study variance

modelnovar3_w8 <- rma.mv(effect_size, variance, random = list(~ 1 | effectsize_id, ~ 1 | study_id),
                         sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered2) 

## Request a print of the results stored in the object # ‘‘modelnovar3’’ in three digits.
summary(modelnovar3_w8, digits=3)

## Perform a likelihood-ratio-test to determine the significance of the between-study (level3) variance.
anova(overall_w8,modelnovar3_w8)


# Moderator Analysis without Outliers
## NDT Modality
## Determine the potential moderating effect of ndt modality.

## Test of different modalities independently
ndt_a_w8 <- rma.mv(effect_size, variance, mods = ~ ndt_a, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered2) 
summary(ndt_a_wo, digits=3)

### manual
ndt_m_w8 <- rma.mv(effect_size, variance, mods = ~ ndt_m, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered2) 
summary(ndt_m_w8, digits=3)

### cognitive
ndt_c_w8 <- rma.mv(effect_size, variance, mods = ~ ndt_c, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered2) 
summary(ndt_c_w8, digits=3)

### handheld
ndt_hand_w8 <- rma.mv(effect_size, variance, mods = ~ hand, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered2) 
summary(ndt_hand_w8, digits=3)

### All modalities in one model
ndt_modality_w8 <- rma.mv(effect_size, variance, mods = ~ ndt_a + ndt_m + ndt_c + hand, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered2) 
summary(ndt_modality_w8, digits=3)

## Urgency
###One-hot encoding of urgency
effectsizes_filtered2 <- effectsizes_filtered2 %>% 
  mutate(urg_low = ifelse(urg=='low', 1, 0),
         urg_medium = ifelse(urg=='medium', 1, 0),
         urg_high = ifelse(urg=='high', 1, 0))

### Determine the potential moderating effect of urgency as categorical moderator.
### No urgency is chosen as the reference category
urgency_w8 <- rma.mv(effect_size, variance, mods = ~ urg_medium + urg_high, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered2) 
summary(urgency_w8, digits=3)

### Determine the potential moderating effect of urgency as continuous moderator.
### High urgency is chosen as inttercept category
urgency_cont_w8 <- rma.mv(effect_size, variance, mods = ~ urg, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered2) 
summary(urgency_cont_w8, digits=3)

### Determine the potential moderating effect of time budget to collision.
tbtc_w8 <- rma.mv(effect_size, variance, mods = ~ tbtc, random = list(~ 1 | effectsize_id, ~1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered2) 
summary(tbtc_w8, digits=3)

## Further Moderator Analyses
### Simulator fidelity
sim_fidelity_w8 <- rma.mv(effect_size, variance, mods = ~ sim, random = list(~1 | effectsize_id, ~ 1 | study_id), sigma2=c(NA,0), tdist=TRUE, data=effectsizes_filtered2) 
summary(sim_fidelity_w8, digits=3)

### Interaction with other road users
iru_w8 <- rma.mv(effect_size, standard_error, mods = ~ iru, random = list(~1 | effectsize_id, ~ 1 | study_id), tdist=TRUE, data=effectsizes_filtered2) 
summary(iru_w8, digits=3)

### Level of automated driving
lad_w8 <- rma.mv(effect_size, standard_error, mods = ~ lad, random = list(~1 | effectsize_id, ~ 1 | study_id), tdist=TRUE, data=effectsizes_filtered2) 
summary(lad_w8, digits=3)



# Publication Bias without Lin et al., 2020
## Visualization: Funnel Plot
funnel(overall_w8,main = "Standard Error", xlab = "Hedges' g") 

## Trim-and-fill method
trimfill <- trimfill(twolevel_w8)
summary(trimfill)
funnel(trimfill,main = "Standard Error", xlab = "Hedges' g") 

## Egger Regression Symmetry Test 
## (eggers.test not working, bc 'overall' is no 'meta'-type. Regtest is carried out with the two-level model)

## eggers.test(x = overall)
regtest(twolevel_w8,  model="rma")

## The p-value is > 0.05, publication bias can be rejected.


##Another method for calculating Egger's regression test - same result though.
egger <- rma.mv(effect_size, variance, mods =  standard_error, random = list(~ 1 | study_id, ~ 1 | effectsize_id), tdist= TRUE, data=effectsizes_filtered2, method="REML")
summary(egger , digits=3)

## Rank-Correlation Test 
ranktest(twolevel_w8)
  
