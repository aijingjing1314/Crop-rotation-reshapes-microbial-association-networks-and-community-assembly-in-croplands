
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Denitrification <- read.csv("Denitrification.csv", fileEncoding = "latin1")
# Check data
head(Denitrification)

# 1. The number of Obversation
total_number <- nrow(Denitrification)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 275  

# 2. The number of Study
unique_studyid_number <- length(unique(Denitrification$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 





#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Denitrification_filteredLongitude_Sub <- subset(Denitrification, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Denitrification_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Denitrification_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredLongitude_Sub %>%
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

overall_model_Denitrification_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Denitrification_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   52.1828  -104.3656   -98.3656   -87.5372   -98.2764   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0040  0.0633     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 1469.2150, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.8128, p-val = 0.2450
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0   -0.0173  0.0110  -1.5661  0.1173  -0.0389  0.0043    
# Longitude_SubLongitudeXy0    0.0130  0.0217   0.6002  0.5484  -0.0295  0.0555       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Denitrification_filteredLongitude_Sub)
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
Denitrification_filteredMAPmean2_Sub <- subset(Denitrification, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Denitrification_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Denitrification_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredMAPmean2_Sub %>%
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

overall_model_Denitrification_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Denitrification_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredMAPmean2_Sub)
# ultivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   64.5919  -129.1839  -121.1839  -106.7607  -121.0341   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0063  0.0792     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1453.5337, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 32.2455, p-val < .0001
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb    ci.ub      
# MAPmean2_SubMAP600Dao1200   -0.0599  0.0153  -3.9046  <.0001  -0.0900  -0.0298  *** 
# MAPmean2_SubMAPDy1200        0.0174  0.0302   0.5776  0.5635  -0.0417   0.0765      
# MAPmean2_SubMAPXy600         0.0209  0.0148   1.4124  0.1578  -0.0081   0.0499   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Denitrification_filteredMAPmean2_Sub)
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
#           "a"                       "b"                       "b" 




### 8.3 MATmean_Sub
Denitrification_filteredMATmean_Sub <- subset(Denitrification, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Denitrification_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Denitrification_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredMATmean_Sub %>%
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

overall_model_Denitrification_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Denitrification_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredMATmean_Sub)
# ultivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   66.2661  -132.5323  -124.5323  -110.1091  -124.3825   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0061  0.0783     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1436.4199, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 35.9014, p-val < .0001
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb    ci.ub      
# MATmean_SubMAT8Dao15   -0.0790  0.0177  -4.4685  <.0001  -0.1136  -0.0443  *** 
# MATmean_SubMATDy15      0.0028  0.0223   0.1252  0.9003  -0.0409   0.0465      
# MATmean_SubMATXy8       0.0109  0.0149   0.7322  0.4640  -0.0183   0.0400   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Denitrification_filteredMATmean_Sub)
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
#      "a"                  "b"                  "b" 




### 8.4 LegumeNonlegume
Denitrification_filteredLegumeNonlegume <- subset(Denitrification, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Denitrification_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Denitrification_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredLegumeNonlegume %>%
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

overall_model_Denitrification_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Denitrification_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   52.1668  -104.3336   -96.3336   -81.9104   -96.1838   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0039  0.0622     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1434.6709, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 10.3148, p-val = 0.0161
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume       -0.0192  0.0112  -1.7175  0.0859  -0.0410  0.0027  . 
# LegumeNonlegumeNon-legume to Legume       -0.0038  0.0104  -0.3627  0.7168  -0.0242  0.0166    
# LegumeNonlegumeNon-legume to Non-legume   -0.0135  0.0108  -1.2493  0.2115  -0.0348  0.0077    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Denitrification_filteredLegumeNonlegume)
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
    #                                 "a"                                     "b"                                    "ab" 



### 8.5 AMnonAM
Denitrification_filteredAMnonAM <- subset(Denitrification, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Denitrification_filteredAMnonAM$AMnonAM <- droplevels(factor(Denitrification_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredAMnonAM %>%
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

overall_model_Denitrification_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Denitrification_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredAMnonAM)
# ultivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   61.3121  -122.6242  -114.6242  -100.2010  -114.4744   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0044  0.0661     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1447.6675, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 25.7916, p-val < .0001
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM      -0.0205  0.0110  -1.8655  0.0621  -0.0421  0.0010  . 
# AMnonAMAM to nonAM    0.0332  0.0138   2.4020  0.0163   0.0061  0.0603  * 
# AMnonAMnonAM to AM   -0.0070  0.0233  -0.3009  0.7635  -0.0526  0.0386 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Denitrification_filteredAMnonAM)
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
   #     "a"                "b"               "ab" 


### 8.6 C3C4
Denitrification_filteredC3C4 <- subset(Denitrification, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Denitrification_filteredC3C4$C3C4 <- droplevels(factor(Denitrification_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredC3C4 %>%
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

overall_model_Denitrification_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Denitrification_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredC3C4)
# ltivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   65.1375  -130.2751  -122.2751  -107.8666  -122.1247   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0029  0.0541     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1425.4274, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 29.8019, p-val < .0001
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb    ci.ub     
# C3C4C3 to C3    0.0214  0.0108   1.9737  0.0484   0.0001   0.0426   * 
# C3C4C3 to C4   -0.0280  0.0094  -2.9701  0.0030  -0.0465  -0.0095  ** 
# C3C4C4 to C3   -0.0105  0.0098  -1.0708  0.2843  -0.0298   0.0087     



# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredC3C4)
vcov_rotation <- vcov(overall_model_Denitrification_filteredC3C4)
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
#          "a"          "b"          "c" 


### 8.7 Annual_Pere
Denitrification_filteredAnnual_Pere <- subset(Denitrification, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
Denitrification_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Denitrification_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredAnnual_Pere %>%
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

overall_model_Denitrification_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Denitrification_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  48.7855  -97.5709  -89.5709  -75.1477  -89.4211   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0042  0.0644     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1453.3338, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.4152, p-val = 0.7020
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual      -0.0098  0.0104  -0.9416  0.3464  -0.0303  0.0106    
# Annual_PereAnnual to Perennial   -0.0124  0.0158  -0.7812  0.4347  -0.0433  0.0186    
# Annual_PerePerennial to Annual   -0.0211  0.0250  -0.8437  0.3988  -0.0700  0.0279   

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Denitrification_filteredAnnual_Pere)
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
Denitrification_filteredPlant_stage <- subset(Denitrification, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
Denitrification_filteredPlant_stage$Plant_stage <- droplevels(factor(Denitrification_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredPlant_stage %>%
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

overall_model_Denitrification_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Denitrification_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredPlant_stage)
# Multivariate Meta-Analysis Model (k = 240; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  16.2349  -32.4698  -24.4698  -10.5975  -24.2974   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0033  0.0576     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 237) = 1228.6087, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 17.5876, p-val = 0.0005
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb    ci.ub      
# Plant_stageHarvest               0.0056  0.0111   0.5009  0.6165  -0.0162   0.0273      
# Plant_stageReproductive stage   -0.0481  0.0141  -3.4042  0.0007  -0.0758  -0.0204  *** 
# Plant_stageVegetative stage      0.0007  0.0115   0.0586  0.9533  -0.0218   0.0232      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Denitrification_filteredPlant_stage)
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
#      "a"                           "b"                           "a" 


### 8.9 Rotation_cycles2
Denitrification_filteredRotation_cycles2 <- subset(Denitrification, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
Denitrification_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Denitrification_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredRotation_cycles2 %>%
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

overall_model_Denitrification_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Denitrification_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredRotation_cycles2)
# ultivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   52.4063  -104.8126   -92.8126   -71.2220   -92.4932   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0042  0.0647     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 1416.0006, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 17.6110, p-val = 0.0035
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb    ci.ub    
# Rotation_cycles2D1      -0.0213  0.0134  -1.5959  0.1105  -0.0475   0.0049    
# Rotation_cycles2D1-3    -0.0187  0.0135  -1.3908  0.1643  -0.0451   0.0077    
# Rotation_cycles2D10     -0.0426  0.0217  -1.9617  0.0498  -0.0851  -0.0000  * 
# Rotation_cycles2D3-5     0.0043  0.0165   0.2607  0.7944  -0.0281   0.0367    
# Rotation_cycles2D5-10    0.0217  0.0163   1.3306  0.1833  -0.0103   0.0537    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Denitrification_filteredRotation_cycles2)
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
#           "ab"                 "abc"                   "a"                   "b"                   "c"


### 8.10 Duration2
Denitrification_filteredDuration2 <- subset(Denitrification, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Denitrification_filteredDuration2$Duration2 <- droplevels(factor(Denitrification_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredDuration2 %>%
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

overall_model_Denitrification_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Denitrification_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredDuration2)
# ultivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   57.7110  -115.4220  -105.4220   -87.4114  -105.1956   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0040  0.0632     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1439.8911, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 21.4694, p-val = 0.0003
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb    ci.ub     
# Duration2D1-5     -0.0279  0.0136  -2.0479  0.0406  -0.0545  -0.0012   * 
# Duration2D11-20   -0.0071  0.0156  -0.4516  0.6515  -0.0377   0.0236     
# Duration2D20-40    0.0495  0.0177   2.7979  0.0051   0.0148   0.0842  ** 
# Duration2D6-10    -0.0072  0.0167  -0.4300  0.6672  -0.0399   0.0255     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredDuration2)
vcov_rotation <- vcov(overall_model_Denitrification_filteredDuration2)
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
            # "a"             "a"             "b"             "a"


### 8.11 SpeciesRichness2
Denitrification_filteredSpeciesRichness2 <- subset(Denitrification, SpeciesRichness2 %in% c("R2", "R3"))
#
Denitrification_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Denitrification_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredSpeciesRichness2 %>%
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

overall_model_Denitrification_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Denitrification_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredSpeciesRichness2)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  44.5260  -89.0520  -83.0520  -72.3242  -82.9597   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0040  0.0634     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 264) = 1433.6529, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.0814, p-val = 0.3532
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2   -0.0145  0.0101  -1.4334  0.1517  -0.0344  0.0053    
# SpeciesRichness2R3   -0.0123  0.0124  -0.9881  0.3231  -0.0366  0.0121  


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Denitrification_filteredSpeciesRichness2)
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
#                "a"                "a"          



### 8.12 Bulk_Rhizosphere
Denitrification_filteredBulk_Rhizosphere <- subset(Denitrification, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Denitrification_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Denitrification_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredBulk_Rhizosphere %>%
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

overall_model_Denitrification_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Denitrification_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   50.1064  -100.2129   -94.2129   -83.3845   -94.1237   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0040  0.0635     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 1469.3684, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.4167, p-val = 0.4925
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere   -0.0120  0.0101  -1.1814  0.2374  -0.0318  0.0079    
# Bulk_RhizosphereRhizosphere       -0.0091  0.0110  -0.8268  0.4083  -0.0306  0.0125   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Denitrification_filteredBulk_Rhizosphere)
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
#                             "a"                             "a" 


### 8.13 Soil_texture
Denitrification_filteredSoil_texture <- subset(Denitrification, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Denitrification_filteredSoil_texture$Soil_texture <- droplevels(factor(Denitrification_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredSoil_texture %>%
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

overall_model_Denitrification_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Denitrification_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredSoil_texture)
# ltivariate Meta-Analysis Model (k = 225; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#   2.0335   -4.0671    3.9329   17.5436    4.1173   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0033  0.0571     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 222) = 1114.3657, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.0742, p-val = 0.5571
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse    0.0246  0.0224   1.1020  0.2705  -0.0192  0.0684    
# Soil_textureFine     -0.0042  0.0264  -0.1588  0.8739  -0.0559  0.0475    
# Soil_textureMedium   -0.0138  0.0151  -0.9136  0.3609  -0.0434  0.0158  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Denitrification_filteredSoil_texture)
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
Denitrification_filteredTillage <- subset(Denitrification, Tillage %in% c("Tillage", "No_tillage"))
#
Denitrification_filteredTillage$Tillage <- droplevels(factor(Denitrification_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredTillage %>%
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

overall_model_Denitrification_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Denitrification_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredTillage)
# Multivariate Meta-Analysis Model (k = 151; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  13.4956  -26.9912  -20.9912  -11.9794  -20.8257   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0055  0.0742     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 149) = 660.1280, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 11.4321, p-val = 0.0033
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0061  0.0217  -0.2799  0.7795  -0.0487  0.0365    
# TillageTillage      -0.0258  0.0216  -1.1928  0.2329  -0.0681  0.0166   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredTillage)
vcov_rotation <- vcov(overall_model_Denitrification_filteredTillage)
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
Denitrification_filteredStraw_retention <- subset(Denitrification, Straw_retention %in% c("Retention", "No_retention"))
#
Denitrification_filteredStraw_retention$Straw_retention <- droplevels(factor(Denitrification_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredStraw_retention %>%
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

overall_model_Denitrification_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Denitrification_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 55; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -11.2452   22.4904   28.4904   34.4013   28.9802   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0020  0.0443     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 379.9511, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.5065, p-val = 0.7763
# 
# Model Results:
# 
#                              estimate      se    zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention    0.0059  0.0168  0.3488  0.7273  -0.0271  0.0389    
# Straw_retentionRetention       0.0104  0.0150  0.6976  0.4855  -0.0189  0.0398    
#    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Denitrification_filteredStraw_retention)
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
Denitrification_filteredPrimer <- subset(Denitrification, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
Denitrification_filteredPrimer$Primer <- droplevels(factor(Denitrification_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Denitrification_filteredPrimer %>%
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

overall_model_Denitrification_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Denitrification_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Denitrification_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 265; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  47.0081  -94.0162  -86.0162  -71.7428  -85.8605   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0040  0.0633     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 262) = 1296.9872, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 0.5117, p-val = 0.9163
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4    0.0004  0.0128   0.0315  0.9748  -0.0247  0.0255    
# PrimerV4      -0.0098  0.0260  -0.3769  0.7063  -0.0607  0.0411    
# PrimerV4-V5   -0.0143  0.0236  -0.6072  0.5437  -0.0606  0.0320    
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Denitrification_filteredPrimer)
vcov_rotation <- vcov(overall_model_Denitrification_filteredPrimer)
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
Denitrification$Wr <- 1 / Denitrification$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Denitrification)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Denitrification)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Denitrification)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Denitrification)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Denitrification)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Denitrification)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Denitrification)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Denitrification)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Denitrification)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Denitrification)

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
# print(model_comparison)
#           Model       AIC       BIC   logLik
# Model1   Model1 -485.5199 -463.8193 248.7599##############################
# Model2   Model2 -481.3547 -459.6541 246.6773
# Model3   Model3 -482.4460 -460.7454 247.2230
# Model4   Model4 -485.4951 -463.7945 248.7476
# Model5   Model5 -482.2532 -460.5526 247.1266
# Model6   Model6 -482.1728 -460.4722 247.0864
# Model7   Model7 -481.4559 -459.7552 246.7279
# Model8   Model8 -482.3630 -460.6623 247.1815
# Model9   Model9 -455.6379 -419.4702 237.8189
# Model10 Model10 -451.9411 -415.7734 235.9706

##### Model 1 is the best model
summary(Model1)
# Number of obs: 275, groups:  StudyID, 50
anova(Model1) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 16.6754 16.6754     1  25.935  3.8374 0.06095 .
# scale(SpeciesRichness)  0.3062  0.3062     1 219.477  0.0705 0.79090  
# scale(Duration)        26.9544 26.9544     1  43.195  6.2029 0.01668 *


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Denitrification)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 3.1648  3.1648     1 15.516  0.6331 0.4382
# scale(SpeciesRichness) 0.1392  0.1392     1 57.230  0.0278 0.8681
# scale(Duration)        2.6110  2.6110     1 28.769  0.5223 0.4757
# scale(pHCK)            0.0630  0.0630     1 30.102  0.0126 0.9114

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Denitrification)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 6.4201  6.4201     1 11.533  1.1507 0.3053
# scale(SpeciesRichness) 1.5222  1.5222     1 77.527  0.2728 0.6029
# scale(Duration)        7.1358  7.1358     1 32.376  1.2789 0.2664
# scale(SOCCK)           0.2834  0.2834     1  9.552  0.0508 0.8264

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Denitrification)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 3.7138  3.7138     1 31.797  1.5909 0.2164
# scale(SpeciesRichness) 0.0153  0.0153     1 29.618  0.0065 0.9361
# scale(Duration)        5.3358  5.3358     1 45.472  2.2857 0.1375
# scale(TNCK)            1.9211  1.9211     1 10.625  0.8229 0.3844


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Denitrification)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# scale(Rotation_cycles) 13.7369 13.7369     1 40.362 13.0871 0.0008184 ***
# scale(SpeciesRichness)  0.0317  0.0317     1 26.658  0.0302 0.8634548    
# scale(Duration)        10.7691 10.7691     1 37.646 10.2597 0.0027661 ** 
# scale(NO3CK)            0.2932  0.2932     1 35.931  0.2794 0.6003701    

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Denitrification)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# scale(Rotation_cycles) 15.1964 15.1964     1 37.349 15.3703 0.0003643 ***
# scale(SpeciesRichness)  0.0465  0.0465     1 27.915  0.0471 0.8298273    
# scale(Duration)        11.1645 11.1645     1 34.601 11.2922 0.0019083 ** 
# scale(NH4CK)            1.1547  1.1547     1 14.701  1.1679 0.2972513

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = Denitrification)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# > anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 3.9825  3.9825     1 24.662  0.6536 0.4266
# scale(SpeciesRichness) 0.1733  0.1733     1 43.209  0.0284 0.8669
# scale(Duration)        3.4817  3.4817     1 32.510  0.5714 0.4552
# scale(APCK)            1.1916  1.1916     1 49.184  0.1956 0.6603

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Denitrification)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 44.115  44.115     1  5.500  5.3022 0.06480 .
# scale(SpeciesRichness)  0.430   0.430     1 33.944  0.0517 0.82156  
# scale(Duration)        50.125  50.125     1  9.486  6.0244 0.03521 *
# scale(AKCK)             0.204   0.204     1 20.481  0.0245 0.87721 

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Denitrification)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.77160 0.77160     1  4.9213  0.0639 0.8107
# scale(SpeciesRichness) 2.37763 2.37763     1 19.6545  0.1968 0.6621
# scale(Duration)        0.81031 0.81031     1  5.0097  0.0671 0.8059
# scale(ANCK)            1.40950 1.40950     1  2.6044  0.1167 0.7583

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Denitrification)
summary(ModelLatitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 16.7750 16.7750     1  26.517  3.8664 0.05981 .
# scale(SpeciesRichness)  0.2452  0.2452     1 228.402  0.0565 0.81231  
# scale(Duration)        27.2372 27.2372     1  44.239  6.2778 0.01598 *
# scale(Latitude)         0.2570  0.2570     1  52.554  0.0592 0.80867  
# ---

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Denitrification)
summary(ModelLongitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 16.4105 16.4105     1  27.316  3.7972 0.06168 .
# scale(SpeciesRichness)  0.2974  0.2974     1 218.897  0.0688 0.79331  
# scale(Duration)        23.0134 23.0134     1  49.532  5.3250 0.02524 *
# scale(Longitude)        0.0214  0.0214     1  13.920  0.0050 0.94485  

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Denitrification)
summary(ModelMAPmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 13.1474 13.1474     1  22.392  3.0156 0.09621 .
# scale(SpeciesRichness)  0.1567  0.1567     1 201.168  0.0359 0.84982  
# scale(Duration)        24.3171 24.3171     1  35.844  5.5776 0.02375 *
# scale(MAPmean)          3.7803  3.7803     1  22.661  0.8671 0.36158 

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Denitrification)
summary(ModelMATmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 14.4204 14.4204     1  24.926  3.3260 0.08021 .
# scale(SpeciesRichness)  0.2667  0.2667     1 217.812  0.0615 0.80435  
# scale(Duration)        25.1874 25.1874     1  40.978  5.8093 0.02051 *
# scale(MATmean)          0.4607  0.4607     1  19.925  0.1063 0.74784  



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Denitrification$SpeciesRichness)) ## n = 275
p1 <- ggplot(Denitrification, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnDenitrification275")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Denitrification$Duration)) ## n = 275
p2 <- ggplot(Denitrification, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnDenitrification275")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Denitrification$Rotation_cycles)) ## n = 275
p3 <- ggplot(Denitrification, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnDenitrification275")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Denitrification$Latitude)) ## n = 275
p5 <- ggplot(Denitrification, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnDenitrification275")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Denitrification$Longitude)) ## n = 275
p6 <- ggplot(Denitrification, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnDenitrification275")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Denitrification$MAPmean)) ## n = 275
p7 <- ggplot(Denitrification, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnDenitrification275")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Denitrification$MATmean)) ## n = 275
p8 <- ggplot(Denitrification, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnDenitrification275")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Denitrification$pHCK)) ## n = 86
p9 <- ggplot(Denitrification, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnDenitrification86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Denitrification$SOCCK)) ## n = 83
p10 <- ggplot(Denitrification, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnDenitrification83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Denitrification$TNCK)) ## n = 52
p11 <- ggplot(Denitrification, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnDenitrification52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Denitrification$NO3CK)) ## n = 46
p12 <- ggplot(Denitrification, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnDenitrification46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Denitrification$NH4CK)) ## n = 45
p13<- ggplot(Denitrification, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnDenitrification45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Denitrification$APCK)) ## n = 64
p14 <- ggplot(Denitrification, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnDenitrification64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Denitrification$AKCK)) ## n = 39
p15 <- ggplot(Denitrification, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnDenitrification39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Denitrification$ANCK)) ## n = 31
p16 <- ggplot(Denitrification, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnDenitrification31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Denitrification$RRpH)) ## n = 86
p17 <- ggplot(Denitrification, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="RRpH")+
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
sum(!is.na(Denitrification$RRSOC)) ## n = 83
p18 <- ggplot(Denitrification, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="RRSOC")+
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
sum(!is.na(Denitrification$RRTN)) ## n = 52
p19 <- ggplot(Denitrification, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="RRTN")+
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
sum(!is.na(Denitrification$RRNO3)) ## n = 44
p20 <- ggplot(Denitrification, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="RRNO3")+
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
sum(!is.na(Denitrification$RRNH4)) ## n = 45
p21 <- ggplot(Denitrification, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="RRNH4")+
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
sum(!is.na(Denitrification$RRAP)) ## n = 64
p22 <- ggplot(Denitrification, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="RRAP")+
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
sum(!is.na(Denitrification$RRAK)) ## n = 39
p23 <- ggplot(Denitrification, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="RRAK")+
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
sum(!is.na(Denitrification$RRAN)) ## n = 31
p24 <- ggplot(Denitrification, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDenitrification", x="RRAN")+
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
sum(!is.na(Denitrification$RRYield)) ## n = 49
p25 <- ggplot(Denitrification, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-0.3, 0.35), expand=c(0, 0)) + 
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



