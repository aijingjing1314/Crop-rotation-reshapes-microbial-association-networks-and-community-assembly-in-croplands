
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Carbon_degradation <- read.csv("Carbon_degradation.csv", fileEncoding = "latin1")
# Check data
head(Carbon_degradation)

# 1. The number of Obversation
total_number <- nrow(Carbon_degradation)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 275  

# 2. The number of Study
unique_studyid_number <- length(unique(Carbon_degradation$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 



#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Carbon_degradation_filteredLongitude_Sub <- subset(Carbon_degradation, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Carbon_degradation_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Carbon_degradation_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredLongitude_Sub %>%
  group_by(Longitude_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Sub Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 LongitudeDa0           117             41
# 2 LongitudeXy0           158              9

overall_model_Carbon_degradation_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Carbon_degradation_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -120.8848   241.7697   247.7697   258.5981   247.8589   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0031  0.0556     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 3126.2922, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 5.7883, p-val = 0.0553
# 
# Model Results:
# 
#                            estimate      se    zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0    0.0194  0.0095  2.0420  0.0412   0.0008  0.0380  * 
# Longitude_SubLongitudeXy0    0.0242  0.0190  1.2722  0.2033  -0.0131  0.0615        

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredLongitude_Sub)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# Longitude_SubLongitudeDa0 Longitude_SubLongitudeXy0 
#                       "a"                       "a" 



### 8.2 MAPmean2_Sub
Carbon_degradation_filteredMAPmean2_Sub <- subset(Carbon_degradation, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Carbon_degradation_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Carbon_degradation_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredMAPmean2_Sub %>%
  group_by(MAPmean2_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MAPmean2_Sub  Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 MAP600Dao1200           72             19
# 2 MAPDy1200               99             11
# 3 MAPXy600               104             21

overall_model_Carbon_degradation_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Carbon_degradation_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredMAPmean2_Sub)
# ultivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -118.0371   236.0741   244.0741   258.4973   244.2239   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0030  0.0552     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 3007.7400, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 15.3014, p-val = 0.0016
# 
# Model Results:
# 
#                            estimate      se    zval    pval    ci.lb   ci.ub     
# MAPmean2_SubMAP600Dao1200    0.0241  0.0108  2.2338  0.0255   0.0029  0.0452   * 
# MAPmean2_SubMAPDy1200        0.0540  0.0188  2.8669  0.0041   0.0171  0.0909  ** 
# MAPmean2_SubMAPXy600         0.0013  0.0105  0.1240  0.9013  -0.0192  0.0218     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredMAPmean2_Sub)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# MAPmean2_SubMAP600Dao1200     MAPmean2_SubMAPDy1200      MAPmean2_SubMAPXy600 
#                       "a"                       "a"                       "b" 




### 8.3 MATmean_Sub
Carbon_degradation_filteredMATmean_Sub <- subset(Carbon_degradation, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Carbon_degradation_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Carbon_degradation_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredMATmean_Sub %>%
  group_by(MATmean_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Sub Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 MAT8Dao15             46             11
# 2 MATDy15              108             17
# 3 MATXy8               121             23

overall_model_Carbon_degradation_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Carbon_degradation_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredMATmean_Sub)
# ltivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -119.5278   239.0556   247.0556   261.4788   247.2054   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0032  0.0565     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 3021.6184, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 12.5646, p-val = 0.0057
# 
# Model Results:
# 
#                       estimate      se    zval    pval    ci.lb   ci.ub     
# MATmean_SubMAT8Dao15    0.0264  0.0125  2.1085  0.0350   0.0019  0.0510   * 
# MATmean_SubMATDy15      0.0396  0.0153  2.5903  0.0096   0.0096  0.0696  ** 
# MATmean_SubMATXy8       0.0046  0.0109  0.4257  0.6703  -0.0167  0.0260       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredMATmean_Sub)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# MATmean_SubMAT8Dao15   MATmean_SubMATDy15    MATmean_SubMATXy8 
#                 "a"                 "ab"                  "b" 




### 8.4 LegumeNonlegume
Carbon_degradation_filteredLegumeNonlegume <- subset(Carbon_degradation, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Carbon_degradation_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Carbon_degradation_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredLegumeNonlegume %>%
  group_by(LegumeNonlegume) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   LegumeNonlegume          Observations Unique_StudyID
#   <fct>                           <int>          <int>
# 1 Legume to Non-legume               36             14
# 2 Non-legume to Legume              151             20
# 3 Non-legume to Non-legume           88             25

overall_model_Carbon_degradation_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Carbon_degradation_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -32.0272   64.0543   72.0543   86.4776   72.2042   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0052  0.0723     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 2976.1411, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 193.2438, p-val < .0001
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb    ci.ub      
# LegumeNonlegumeLegume to Non-legume       -0.0016  0.0121  -0.1285  0.8977  -0.0252   0.0221      
# LegumeNonlegumeNon-legume to Legume       -0.0273  0.0114  -2.4008  0.0164  -0.0495  -0.0050    * 
# LegumeNonlegumeNon-legume to Non-legume    0.0678  0.0116   5.8677  <.0001   0.0452   0.0905  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredLegumeNonlegume)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
    # LegumeNonlegumeLegume to Non-legume     LegumeNonlegumeNon-legume to Legume LegumeNonlegumeNon-legume to Non-legume 
    #                                 "a"                                     "b"                                     "c" 




### 8.5 AMnonAM
Carbon_degradation_filteredAMnonAM <- subset(Carbon_degradation, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Carbon_degradation_filteredAMnonAM$AMnonAM <- droplevels(factor(Carbon_degradation_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             209             37
# 2 AM to nonAM           30              8
# 3 nonAM to AM           36              8

overall_model_Carbon_degradation_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Carbon_degradation_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredAMnonAM)
# ltivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -122.3757   244.7514   252.7514   267.1746   252.9012   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0029  0.0540     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 2538.6054, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 9.8585, p-val = 0.0198
# 
# Model Results:
# 
#                     estimate      se    zval    pval    ci.lb   ci.ub     
# AMnonAMAM to AM       0.0237  0.0087  2.7274  0.0064   0.0067  0.0406  ** 
# AMnonAMAM to nonAM    0.0080  0.0103  0.7781  0.4365  -0.0122  0.0283     
# AMnonAMnonAM to AM    0.0146  0.0127  1.1468  0.2514  -0.0103  0.0395     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredAMnonAM)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
   # AMnonAMAM to AM AMnonAMAM to nonAM AMnonAMnonAM to AM 
   #             "a"                "a"                "a" 


### 8.6 C3C4
Carbon_degradation_filteredC3C4 <- subset(Carbon_degradation, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Carbon_degradation_filteredC3C4$C3C4 <- droplevels(factor(Carbon_degradation_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredC3C4 %>%
  group_by(C3C4) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   C3C4     Observations Unique_StudyID
#   <fct>           <int>          <int>
# 1 C3 to C3           73             18
# 2 C3 to C4           69             26
# 3 C4 to C3          132             14

overall_model_Carbon_degradation_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Carbon_degradation_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredC3C4)
# ltivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -117.4782   234.9564   242.9564   257.3649   243.1068   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0026  0.0512     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 3039.7272, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 22.8820, p-val < .0001
# 
# Model Results:
# 
#               estimate      se    zval    pval    ci.lb   ci.ub      
# C3C4C3 to C3    0.0085  0.0094  0.9055  0.3652  -0.0099  0.0268      
# C3C4C3 to C4    0.0328  0.0085  3.8681  0.0001   0.0162  0.0494  *** 
# C3C4C4 to C3    0.0158  0.0089  1.7720  0.0764  -0.0017  0.0332    . 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredC3C4)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredC3C4)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# C3C4C3 to C3 C3C4C3 to C4 C3C4C4 to C3 
#          "a"          "b"          "a" 


### 8.7 Annual_Pere
Carbon_degradation_filteredAnnual_Pere <- subset(Carbon_degradation, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
Carbon_degradation_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Carbon_degradation_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredAnnual_Pere %>%
  group_by(Annual_Pere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Annual_Pere         Observations Unique_StudyID
#   <fct>                      <int>          <int>
# 1 Annual to Annual             241             40
# 2 Annual to Perennial           13              8
# 3 Perennial to Annual           21              7

overall_model_Carbon_degradation_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Carbon_degradation_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredAnnual_Pere)
# ltivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -124.2439   248.4879   256.4879   270.9111   256.6377   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0030  0.0549     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 3036.7216, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 6.3744, p-val = 0.0947
# 
# Model Results:
# 
#                                 estimate      se    zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual       0.0185  0.0088  2.1005  0.0357   0.0012  0.0357  * 
# Annual_PereAnnual to Perennial    0.0251  0.0130  1.9303  0.0536  -0.0004  0.0505  . 
# Annual_PerePerennial to Annual    0.0291  0.0148  1.9641  0.0495   0.0001  0.0582  * 

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredAnnual_Pere)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
   # Annual_PereAnnual to Annual Annual_PereAnnual to Perennial Annual_PerePerennial to Annual 
   #                         "a"                           "a"                            "a" 



### 8.8 Plant_stage
Carbon_degradation_filteredPlant_stage <- subset(Carbon_degradation, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
Carbon_degradation_filteredPlant_stage$Plant_stage <- droplevels(factor(Carbon_degradation_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredPlant_stage %>%
  group_by(Plant_stage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Plant_stage        Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Harvest                      98             23
# 2 Reproductive stage           61             11
# 3 Vegetative stage             81              8

overall_model_Carbon_degradation_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Carbon_degradation_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredPlant_stage)
# ultivariate Meta-Analysis Model (k = 240; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -155.5632   311.1263   319.1263   332.9986   319.2987   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0039  0.0628     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 237) = 2769.8030, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 55.2046, p-val < .0001
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb    ci.ub     
# Plant_stageHarvest               0.0310  0.0113   2.7472  0.0060   0.0089   0.0531  ** 
# Plant_stageReproductive stage   -0.0387  0.0129  -2.9945  0.0027  -0.0640  -0.0134  ** 
# Plant_stageVegetative stage      0.0303  0.0117   2.6044  0.0092   0.0075   0.0532  ** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredPlant_stage)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# Plant_stageHarvest Plant_stageReproductive stage   Plant_stageVegetative stage 
# "a"                           "b"                           "a" 


### 8.9 Rotation_cycles2
Carbon_degradation_filteredRotation_cycles2 <- subset(Carbon_degradation, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
Carbon_degradation_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Carbon_degradation_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredRotation_cycles2 %>%
  group_by(Rotation_cycles2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Rotation_cycles2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                         80             19
# 2 D1-3                       87             12
# 3 D10                        11              7
# 4 D3-5                       35              9
# 5 D5-10                      62             15

overall_model_Carbon_degradation_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Carbon_degradation_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredRotation_cycles2)
# ultivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -102.8705   205.7410   217.7410   239.3315   218.0604   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0028  0.0529     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 2719.3274, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 57.7078, p-val < .0001
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1       0.0266  0.0109   2.4419  0.0146   0.0052  0.0479  * 
# Rotation_cycles2D1-3     0.0254  0.0111   2.2884  0.0221   0.0036  0.0471  * 
# Rotation_cycles2D10      0.0097  0.0175   0.5521  0.5809  -0.0247  0.0441    
# Rotation_cycles2D3-5    -0.0103  0.0133  -0.7721  0.4401  -0.0363  0.0158    
# Rotation_cycles2D5-10    0.0263  0.0131   2.0081  0.0446   0.0006  0.0520  * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredRotation_cycles2)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
#    Rotation_cycles2D1  Rotation_cycles2D1-3   Rotation_cycles2D10  Rotation_cycles2D3-5 Rotation_cycles2D5-10 
#             "a"                   "a"                  "ab"                   "b"                   "a"                    "a" 


### 8.10 Duration2
Carbon_degradation_filteredDuration2 <- subset(Carbon_degradation, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Carbon_degradation_filteredDuration2$Duration2 <- droplevels(factor(Carbon_degradation_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredDuration2 %>%
  group_by(Duration2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Duration2 Observations Unique_StudyID
#   <fct>            <int>          <int>
# 1 D1-5               147             26
# 2 D11-20              58             11
# 3 D20-40              37              8
# 4 D6-10               33             11

overall_model_Carbon_degradation_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Carbon_degradation_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredDuration2)
# ltivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -123.3555   246.7110   256.7110   274.7216   256.9374   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0031  0.0561     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 2790.2412, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 9.3734, p-val = 0.0524
# 
# Model Results:
# 
#                  estimate      se    zval    pval    ci.lb   ci.ub    
# Duration2D1-5      0.0220  0.0118  1.8629  0.0625  -0.0011  0.0452  . 
# Duration2D11-20    0.0143  0.0136  1.0520  0.2928  -0.0124  0.0410    
# Duration2D20-40    0.0058  0.0145  0.4000  0.6891  -0.0227  0.0343    
# Duration2D6-10     0.0318  0.0145  2.1890  0.0286   0.0033  0.0603  * 
 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredDuration2)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredDuration2)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
  # Duration2D1-5 Duration2D11-20 Duration2D20-40  Duration2D6-10 
            # "a"             "a"             "a"             "a" 


### 8.11 SpeciesRichness2
Carbon_degradation_filteredSpeciesRichness2 <- subset(Carbon_degradation, SpeciesRichness2 %in% c("R2", "R3"))
#
Carbon_degradation_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Carbon_degradation_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredSpeciesRichness2 %>%
  group_by(SpeciesRichness2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   SpeciesRichness2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 R2                        215             44
# 2 R3                         51              8

overall_model_Carbon_degradation_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Carbon_degradation_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredSpeciesRichness2)
# variate Meta-Analysis Model (k = 266; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -76.3108  152.6216  158.6216  169.3494  158.7139   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0032  0.0567     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 264) = 3063.4215, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 120.4009, p-val < .0001
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb    ci.ub      
# SpeciesRichness2R2    0.0242  0.0089   2.7226  0.0065   0.0068   0.0416   ** 
# SpeciesRichness2R3   -0.0388  0.0103  -3.7742  0.0002  -0.0590  -0.0187  *** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredSpeciesRichness2)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# SpeciesRichness2R2 SpeciesRichness2R3 
#                "a"                "b"          



### 8.12 Bulk_Rhizosphere
Carbon_degradation_filteredBulk_Rhizosphere <- subset(Carbon_degradation, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Carbon_degradation_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Carbon_degradation_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
 # Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non_Rhizosphere           210             36
# 2 Rhizosphere                65             24

overall_model_Carbon_degradation_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Carbon_degradation_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredBulk_Rhizosphere)
# ultivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -110.3367   220.6735   226.6735   237.5019   226.7627   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0029  0.0540     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 2828.1820, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 30.5537, p-val < .0001
# 
# Model Results:
# 
#                                  estimate      se    zval    pval    ci.lb   ci.ub      
# Bulk_RhizosphereNon_Rhizosphere    0.0124  0.0084  1.4696  0.1417  -0.0041  0.0289      
# Bulk_RhizosphereRhizosphere        0.0345  0.0088  3.9401  <.0001   0.0173  0.0517  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredBulk_Rhizosphere)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# Bulk_RhizosphereNon-Rhizosphere     Bulk_RhizosphereRhizosphere 
#                             "a"                             "b" 


### 8.13 Soil_texture
Carbon_degradation_filteredSoil_texture <- subset(Carbon_degradation, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Carbon_degradation_filteredSoil_texture$Soil_texture <- droplevels(factor(Carbon_degradation_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredSoil_texture %>%
  group_by(Soil_texture) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Soil_texture Observations Unique_StudyID
#   <fct>               <int>          <int>
# 1 Coarse                 30              7
# 2 Fine                   20              8
# 3 Medium                175             17

overall_model_Carbon_degradation_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Carbon_degradation_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredSoil_texture)
# ltivariate Meta-Analysis Model (k = 225; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -187.0406   374.0812   382.0812   395.6919   382.2655   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0031  0.0558     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 222) = 2927.9468, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.0456, p-val = 0.2566
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse   -0.0109  0.0218  -0.5014  0.6161  -0.0537  0.0318    
# Soil_textureFine      0.0441  0.0244   1.8088  0.0705  -0.0037  0.0919  . 
# Soil_textureMedium    0.0102  0.0141   0.7228  0.4698  -0.0175  0.0379    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredSoil_texture)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# Soil_textureCoarse   Soil_textureFine Soil_textureMedium 
#                "a"                "a"                "a" 


### 8.14 Tillage
Carbon_degradation_filteredTillage <- subset(Carbon_degradation, Tillage %in% c("Tillage", "No_tillage"))
#
Carbon_degradation_filteredTillage$Tillage <- droplevels(factor(Carbon_degradation_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredTillage %>%
  group_by(Tillage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Tillage    Observations Unique_StudyID
#   <fct>             <int>          <int>
# 1 No_tillage          120              7
# 2 Tillage              31              9

overall_model_Carbon_degradation_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Carbon_degradation_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredTillage)
# te Meta-Analysis Model (k = 151; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  134.7432  -269.4864  -263.4864  -254.4746  -263.3209   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0024  0.0488     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 149) = 568.8023, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 3.9291, p-val = 0.1402
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0098  0.0153  -0.6404  0.5219  -0.0397  0.0202    
# TillageTillage       0.0042  0.0151   0.2752  0.7832  -0.0255  0.0338    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredTillage)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredTillage)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# TillageNo_tillage    TillageTillage 
#               "a"               "b" 


### 8.15 Straw_retention
Carbon_degradation_filteredStraw_retention <- subset(Carbon_degradation, Straw_retention %in% c("Retention", "No_retention"))
#
Carbon_degradation_filteredStraw_retention$Straw_retention <- droplevels(factor(Carbon_degradation_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredStraw_retention %>%
  group_by(Straw_retention) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Straw_retention Observations Unique_StudyID
#   <fct>                  <int>          <int>
# 1 No_retention              12              5
# 2 Retention                 43              8

overall_model_Carbon_degradation_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Carbon_degradation_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 55; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -73.3136  146.6272  152.6272  158.5380  153.1170   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0052  0.0719     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 759.7178, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.8805, p-val = 0.2369
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention    0.0131  0.0240   0.5470  0.5844  -0.0339  0.0602    
# Straw_retentionRetention      -0.0102  0.0229  -0.4450  0.6563  -0.0551  0.0347     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredStraw_retention)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# Straw_retentionNo_retention    Straw_retentionRetention 
#                         "a"                         "a" 

### 8.16 Primer
Carbon_degradation_filteredPrimer <- subset(Carbon_degradation, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
Carbon_degradation_filteredPrimer$Primer <- droplevels(factor(Carbon_degradation_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Carbon_degradation_filteredPrimer %>%
  group_by(Primer) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Primer Observations Unique_StudyID
#   <fct>         <int>          <int>
# 1 V3-V4           126             28
# 2 V4              100              9
# 3 V4-V5            39              9

overall_model_Carbon_degradation_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Carbon_degradation_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Carbon_degradation_filteredPrimer)
# ultivariate Meta-Analysis Model (k = 265; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -111.0811   222.1622   230.1622   244.4356   230.3179   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0036  0.0597     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 262) = 2863.4845, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.4860, p-val = 0.2135
# 
# Model Results:
# 
#              estimate      se    zval    pval    ci.lb   ci.ub    
# PrimerV3-V4    0.0164  0.0119  1.3738  0.1695  -0.0070  0.0397    
# PrimerV4       0.0241  0.0224  1.0736  0.2830  -0.0199  0.0681    
# PrimerV4-V5    0.0263  0.0218  1.2025  0.2292  -0.0165  0.0690
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_degradation_filteredPrimer)
vcov_rotation <- vcov(overall_model_Carbon_degradation_filteredPrimer)
# define
pairwise_comparison <- function(coefs, vcovs, group1, group2) {
  diff <- coefs[group1] - coefs[group2]  # 
  se_diff <- sqrt(vcovs[group1, group1] + vcovs[group2, group2] - 2 * vcovs[group1, group2]) 
  z <- diff / se_diff  # Z 
  p <- 2 * (1 - pnorm(abs(z)))  # 
  return(p)
}
# Compare
group_names <- names(coef_rotation)
p_matrix <- matrix(NA, nrow = length(group_names), ncol = length(group_names),
                   dimnames = list(group_names, group_names))
for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i < j) {
      p_matrix[i, j] <- pairwise_comparison(coef_rotation, vcov_rotation, group_names[i], group_names[j])
    }
  }
}
# Convert to letter
p_matrix[lower.tri(p_matrix)] <- t(p_matrix)[lower.tri(p_matrix)] 
significance_letters <- multcompLetters(p_matrix)$Letters
# Output
print(significance_letters)
# PrimerV3-V4    PrimerV4 PrimerV4-V5 
#         "a"         "a"         "a"









#### 9. Linear Mixed Effect Model
# 
Carbon_degradation$Wr <- 1 / Carbon_degradation$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Carbon_degradation)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Carbon_degradation)

# Define
models <- list(
  Model1 = Model1,
  Model2 = Model2,
  Model3 = Model3,
  Model4 = Model4,
  Model5 = Model5,
  Model6 = Model6,
  Model7 = Model7,
  Model8 = Model8,
  Model9 = Model9,
  Model10 = Model10
)
# 创建一个数据框来存储每个模型的AIC、BIC和logLik
model_comparison <- data.frame(
  Model = names(models),
  AIC = sapply(models, AIC),
  BIC = sapply(models, BIC),
  logLik = sapply(models, logLik)
)
# 查看结果
print(model_comparison)
#           Model       AIC       BIC   logLik
# Model1   Model1 -482.7387 -461.0381 247.3694
# Model2   Model2 -484.4591 -462.7585 248.2296###########################
# Model3   Model3 -482.8943 -461.1936 247.4471
# Model4   Model4 -483.8760 -462.1754 247.9380
# Model5   Model5 -482.7597 -461.0591 247.3799
# Model6   Model6 -483.9401 -462.2395 247.9701
# Model7   Model7 -483.3079 -461.6073 247.6540
# Model8   Model8 -483.9987 -462.2980 247.9993
# Model9   Model9 -455.0328 -418.8651 237.5164
# Model10 Model10 -459.7119 -423.5442 239.8559

##### Model 2 is the best model
summary(Model2)
# Number of obs: 275, groups:  StudyID, 50
anova(Model2)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles))  0.087   0.087     1 122.74  0.0151 0.902471   
# scale(log(SpeciesRichness)) 61.981  61.981     1 261.02 10.7714 0.001171 **
# scale(log(Duration))         0.090   0.090     1 116.53  0.0156 0.900860  


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  4.677   4.677     1 21.758  0.5863 0.45209  
# scale(log(SpeciesRichness))  2.207   2.207     1 68.122  0.2766 0.60063  
# scale(log(Duration))         1.212   1.212     1 30.360  0.1519 0.69943  
# scale(pHCK)                 38.996  38.996     1 25.251  4.8885 0.03631 *

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  3.4044  3.4044     1 28.983  0.4349 0.5148
# scale(log(SpeciesRichness))  0.8912  0.8912     1 65.669  0.1138 0.7369
# scale(log(Duration))         0.6805  0.6805     1 39.777  0.0869 0.7697
# scale(SOCCK)                10.5474 10.5474     1 12.991  1.3473 0.2666

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 5.2157  5.2157     1 13.167  2.0791 0.17269  
# scale(log(SpeciesRichness)) 1.5856  1.5856     1 39.001  0.6321 0.43140  
# scale(log(Duration))        7.5258  7.5258     1 19.581  3.0000 0.09899 .
# scale(TNCK)                 1.2280  1.2280     1 11.901  0.4895 0.49759  


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.92157 1.92157     1 22.799  1.2747 0.2706
# scale(log(SpeciesRichness)) 0.35090 0.35090     1 28.761  0.2328 0.6331
# scale(log(Duration))        1.97265 1.97265     1 32.344  1.3086 0.2610
# scale(NO3CK)                0.00803 0.00803     1 20.780  0.0053 0.9425

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 2.04041 2.04041     1 24.068  1.3976 0.2487
# scale(log(SpeciesRichness)) 0.35849 0.35849     1 29.050  0.2455 0.6240
# scale(log(Duration))        1.93847 1.93847     1 32.276  1.3278 0.2577
# scale(NH4CK)                0.64773 0.64773     1 12.398  0.4437 0.5176

# #### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(APCK) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# > anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  4.8023  4.8023     1 12.833  0.5243 0.4820
# scale(log(SpeciesRichness))  0.0398  0.0398     1 50.749  0.0043 0.9477
# scale(log(Duration))         8.2539  8.2539     1 16.147  0.9011 0.3565
# scale(APCK)                 13.6797 13.6797     1 36.774  1.4935 0.2295

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  7.1307  7.1307     1  8.0859  0.4902 0.5035
# scale(log(SpeciesRichness))  0.0719  0.0719     1 28.3255  0.0049 0.9444
# scale(log(Duration))        24.2754 24.2754     1  9.6131  1.6689 0.2266
# scale(AKCK)                 12.3140 12.3140     1 19.5334  0.8466 0.3687

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  8.2364  8.2364     1 12.8151  0.4633 0.5082
# scale(log(SpeciesRichness))  3.7789  3.7789     1 24.4029  0.2126 0.6489
# scale(log(Duration))         4.7023  4.7023     1 11.8456  0.2645 0.6165
# scale(ANCK)                 21.9869 21.9869     1  6.5773  1.2367 0.3051

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
summary(ModelLatitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles))  0.104   0.104     1 126.198  0.0180 0.893481   
# scale(log(SpeciesRichness)) 62.182  62.182     1 260.323 10.7981 0.001156 **
# scale(log(Duration))         0.104   0.104     1 120.491  0.0181 0.893136   
# scale(Latitude)              0.260   0.260     1  87.927  0.0451 0.832337  

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
summary(ModelLongitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles))  0.023   0.023     1 122.164  0.0040 0.949960   
# scale(log(SpeciesRichness)) 63.346  63.346     1 260.304 11.0600 0.001009 **
# scale(log(Duration))         0.329   0.329     1 116.802  0.0574 0.811063   
# scale(Longitude)             6.389   6.389     1  27.043  1.1155 0.300231 

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
summary(ModelMAPmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)    
# scale(log(Rotation_cycles))  0.203   0.203     1 134.386  0.0366 0.848621    
# scale(log(SpeciesRichness)) 69.268  69.268     1 264.812 12.4526 0.000492 ***
# scale(log(Duration))         0.168   0.168     1 122.361  0.0301 0.862502    
# scale(MAPmean)              32.097  32.097     1  64.298  5.7703 0.019198 * 

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Carbon_degradation)
summary(ModelMATmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)    
# scale(log(Rotation_cycles))  0.289   0.289     1 129.69  0.0503 0.822923    
# scale(log(SpeciesRichness)) 64.234  64.234     1 259.68 11.1897 0.000944 ***
# scale(log(Duration))         0.311   0.311     1 120.90  0.0542 0.816333    
# scale(MATmean)               2.707   2.707     1  43.04  0.4716 0.495914  



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Carbon_degradation$SpeciesRichness)) ## n = 275
p1 <- ggplot(Carbon_degradation, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnCarbon_degradation275")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Carbon_degradation$Duration)) ## n = 275
p2 <- ggplot(Carbon_degradation, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnCarbon_degradation275")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Carbon_degradation$Rotation_cycles)) ## n = 275
p3 <- ggplot(Carbon_degradation, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnCarbon_degradation275")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Carbon_degradation$Latitude)) ## n = 275
p5 <- ggplot(Carbon_degradation, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnCarbon_degradation275")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Carbon_degradation$Longitude)) ## n = 275
p6 <- ggplot(Carbon_degradation, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnCarbon_degradation275")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Carbon_degradation$MAPmean)) ## n = 275
p7 <- ggplot(Carbon_degradation, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnCarbon_degradation275")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Carbon_degradation$MATmean)) ## n = 275
p8 <- ggplot(Carbon_degradation, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnCarbon_degradation275")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Carbon_degradation$pHCK)) ## n = 86
p9 <- ggplot(Carbon_degradation, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnCarbon_degradation86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Carbon_degradation$SOCCK)) ## n = 83
p10 <- ggplot(Carbon_degradation, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnCarbon_degradation83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Carbon_degradation$TNCK)) ## n = 52
p11 <- ggplot(Carbon_degradation, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnCarbon_degradation52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Carbon_degradation$NO3CK)) ## n = 46
p12 <- ggplot(Carbon_degradation, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnCarbon_degradation46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Carbon_degradation$NH4CK)) ## n = 45
p13<- ggplot(Carbon_degradation, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnCarbon_degradation45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Carbon_degradation$APCK)) ## n = 64
p14 <- ggplot(Carbon_degradation, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnCarbon_degradation64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Carbon_degradation$AKCK)) ## n = 39
p15 <- ggplot(Carbon_degradation, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnCarbon_degradation39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Carbon_degradation$ANCK)) ## n = 31
p16 <- ggplot(Carbon_degradation, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnCarbon_degradation31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Carbon_degradation$RRpH)) ## n = 86
p17 <- ggplot(Carbon_degradation, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="RRpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRpH" , y="RR86")
p17
pdf("RRpH.pdf",width=8,height=8)
p17
dev.off() 

## RRSOC
sum(!is.na(Carbon_degradation$RRSOC)) ## n = 83
p18 <- ggplot(Carbon_degradation, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="RRSOC")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRSOC" , y="RR83")
p18
pdf("RRSOC.pdf",width=8,height=8)
p18
dev.off() 

## RRTN
sum(!is.na(Carbon_degradation$RRTN)) ## n = 52
p19 <- ggplot(Carbon_degradation, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="RRTN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRTN" , y="RR52")
p19
pdf("RRTN.pdf",width=8,height=8)
p19
dev.off() 

## RRNO3
sum(!is.na(Carbon_degradation$RRNO3)) ## n = 44
p20 <- ggplot(Carbon_degradation, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="RRNO3")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNO3" , y="RR44")
p20
pdf("RRNO3.pdf",width=8,height=8)
p20
dev.off() 

## RRNH4
sum(!is.na(Carbon_degradation$RRNH4)) ## n = 45
p21 <- ggplot(Carbon_degradation, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="RRNH4")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNH4" , y="RR45")
p21
pdf("RRNH4.pdf",width=8,height=8)
p21
dev.off() 

## RRAP
sum(!is.na(Carbon_degradation$RRAP)) ## n = 64
p22 <- ggplot(Carbon_degradation, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="RRAP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAP" , y="RR64")
p22
pdf("RRAP.pdf",width=8,height=8)
p22
dev.off() 

## RRAK
sum(!is.na(Carbon_degradation$RRAK)) ## n = 39
p23 <- ggplot(Carbon_degradation, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="RRAK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAK" , y="RR39")
p23
pdf("RRAK.pdf",width=8,height=8)
p23
dev.off() 

## RRAN
sum(!is.na(Carbon_degradation$RRAN)) ## n = 31
p24 <- ggplot(Carbon_degradation, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_degradation", x="RRAN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAN" , y="RR31")
p24
pdf("RRAN.pdf",width=8,height=8)
p24
dev.off() 

## RRYield
sum(!is.na(Carbon_degradation$RRYield)) ## n = 49
p25 <- ggplot(Carbon_degradation, aes(x=RR, y=RRYield)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm, color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") +
  theme_bw() +
  theme(text = element_text(family = "serif", size=20)) +
  geom_vline(aes(xintercept=0), colour="black", linewidth=0.5, linetype="dashed") +
  labs(x="RR", y="RRYield49") +
  theme(panel.grid=element_blank()) + 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = x ~ y,  
    parse = TRUE,
    color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85
  ) +
  stat_cor(method = "spearman", size = 5) +
  scale_x_continuous(limits=c(-0.34, 0.2), expand=c(0, 0)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) 
p25
pdf("RRYield_swapped_limited.pdf", width=8, height=8)
p25
dev.off()

library(patchwork)

combined_plot <- (
  p17 | p18 | p19
) / (
  p20 | p21 | p22
) / (
  p23 | p24 | p25
)

combined_plot

ggsave(
  "Combined_RR.pdf",
  combined_plot,
  width = 10,
  height = 10
)


