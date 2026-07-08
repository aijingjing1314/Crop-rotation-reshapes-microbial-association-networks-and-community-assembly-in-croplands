
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
P_starvation_regulation <- read.csv("P_starvation_regulation.csv", fileEncoding = "latin1")
# Check data
head(P_starvation_regulation)

# 1. The number of Obversation
total_number <- nrow(P_starvation_regulation)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 275  

# 2. The number of Study
unique_studyid_number <- length(unique(P_starvation_regulation$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 





#### 8. Subgroup analysis
### 8.1 Longitude_Sub
P_starvation_regulation_filteredLongitude_Sub <- subset(P_starvation_regulation, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
P_starvation_regulation_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(P_starvation_regulation_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredLongitude_Sub %>%
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

overall_model_P_starvation_regulation_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  14.8798  -29.7596  -23.7596  -12.9312  -23.6704   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0372     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 2154.4220, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.5163, p-val = 0.7725
# 
# Model Results:
# 
#                            estimate      se    zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0    0.0001  0.0064  0.0199  0.9841  -0.0125  0.0128    
# Longitude_SubLongitudeXy0    0.0092  0.0128  0.7183  0.4726  -0.0158  0.0342      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredLongitude_Sub)
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
P_starvation_regulation_filteredMAPmean2_Sub <- subset(P_starvation_regulation, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
P_starvation_regulation_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(P_starvation_regulation_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredMAPmean2_Sub %>%
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

overall_model_P_starvation_regulation_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  14.2101  -28.4202  -20.4202   -5.9970  -20.2704   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0378     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 2068.2558, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.8702, p-val = 0.4121
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean2_SubMAP600Dao1200    0.0059  0.0082   0.7173  0.4732  -0.0102  0.0220    
# MAPmean2_SubMAPDy1200        0.0153  0.0142   1.0795  0.2804  -0.0125  0.0431    
# MAPmean2_SubMAPXy600        -0.0062  0.0077  -0.8106  0.4176  -0.0213  0.0088    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredMAPmean2_Sub)
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
#                       "a"                       "a"                       "a" 




### 8.3 MATmean_Sub
P_starvation_regulation_filteredMATmean_Sub <- subset(P_starvation_regulation, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
P_starvation_regulation_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(P_starvation_regulation_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredMATmean_Sub %>%
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

overall_model_P_starvation_regulation_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  13.8378  -27.6756  -19.6756   -5.2523  -19.5257   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0380     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 2136.8611, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.5551, p-val = 0.4654
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubMAT8Dao15    0.0070  0.0100   0.7048  0.4810  -0.0125  0.0266    
# MATmean_SubMATDy15      0.0109  0.0109   1.0043  0.3152  -0.0104  0.0322    
# MATmean_SubMATXy8      -0.0053  0.0076  -0.6985  0.4849  -0.0201  0.0096 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredMATmean_Sub)
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
#                  "a"                 "a"                  "a" 




### 8.4 LegumeNonlegume
P_starvation_regulation_filteredLegumeNonlegume <- subset(P_starvation_regulation, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
P_starvation_regulation_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(P_starvation_regulation_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredLegumeNonlegume %>%
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

overall_model_P_starvation_regulation_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  256.3599  -512.7197  -504.7197  -490.2965  -504.5699   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0019  0.0440     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 2076.2008, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 494.8848, p-val < .0001
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb    ci.ub      
# LegumeNonlegumeLegume to Non-legume       -0.0277  0.0071  -3.8821  0.0001  -0.0417  -0.0137  *** 
# LegumeNonlegumeNon-legume to Legume        0.0102  0.0070   1.4544  0.1458  -0.0036   0.0240      
# LegumeNonlegumeNon-legume to Non-legume    0.0116  0.0073   1.5833  0.1133  -0.0028   0.0260   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredLegumeNonlegume)
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
    #                                 "a"                                     "b"                                     "b" 




### 8.5 AMnonAM
P_starvation_regulation_filteredAMnonAM <- subset(P_starvation_regulation, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
P_starvation_regulation_filteredAMnonAM$AMnonAM <- droplevels(factor(P_starvation_regulation_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredAMnonAM %>%
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

overall_model_P_starvation_regulation_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  11.8689  -23.7378  -15.7378   -1.3145  -15.5879   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0371     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 2232.4865, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 0.6055, p-val = 0.8952
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM       0.0023  0.0062   0.3708  0.7108  -0.0099  0.0145    
# AMnonAMAM to nonAM   -0.0018  0.0082  -0.2239  0.8228  -0.0179  0.0142    
# AMnonAMnonAM to AM    0.0040  0.0109   0.3646  0.7154  -0.0174  0.0253 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredAMnonAM)
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
P_starvation_regulation_filteredC3C4 <- subset(P_starvation_regulation, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
P_starvation_regulation_filteredC3C4$C3C4 <- droplevels(factor(P_starvation_regulation_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredC3C4 %>%
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

overall_model_P_starvation_regulation_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  283.6262  -567.2523  -559.2523  -544.8438  -559.1019   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0013  0.0367     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1978.5553, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 540.6539, p-val < .0001
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb   ci.ub      
# C3C4C3 to C3    0.0001  0.0071   0.0141  0.9888  -0.0138  0.0140      
# C3C4C3 to C4   -0.0076  0.0062  -1.2425  0.2140  -0.0197  0.0044      
# C3C4C4 to C3    0.0329  0.0062   5.2774  <.0001   0.0207  0.0452  *** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredC3C4)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredC3C4)
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
#       "a"          "a"          "b" 


### 8.7 Annual_Pere
P_starvation_regulation_filteredAnnual_Pere <- subset(P_starvation_regulation, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
P_starvation_regulation_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(P_starvation_regulation_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredAnnual_Pere %>%
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

overall_model_P_starvation_regulation_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  12.4861  -24.9722  -16.9722   -2.5490  -16.8224   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0013  0.0366     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 2271.3008, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.7250, p-val = 0.4360
# 
# Model Results:
# 
#                                 estimate      se    zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual       0.0000  0.0060  0.0022  0.9983  -0.0117  0.0118    
# Annual_PereAnnual to Perennial    0.0138  0.0094  1.4577  0.1449  -0.0047  0.0323    
# Annual_PerePerennial to Annual    0.0047  0.0119  0.3912  0.6957  -0.0187  0.0281   

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredAnnual_Pere)
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
   #              "a"                            "a"                            "a" 



### 8.8 Plant_stage
P_starvation_regulation_filteredPlant_stage <- subset(P_starvation_regulation, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
P_starvation_regulation_filteredPlant_stage$Plant_stage <- droplevels(factor(P_starvation_regulation_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredPlant_stage %>%
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

overall_model_P_starvation_regulation_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredPlant_stage)
# Multivariate Meta-Analysis Model (k = 240; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -36.7710   73.5419   81.5419   95.4142   81.7143   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0017  0.0411     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 237) = 2048.3289, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.0040, p-val = 0.8003
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb   ci.ub    
# Plant_stageHarvest               0.0016  0.0075   0.2127  0.8316  -0.0131  0.0163    
# Plant_stageReproductive stage   -0.0052  0.0088  -0.5849  0.5586  -0.0224  0.0121    
# Plant_stageVegetative stage      0.0020  0.0077   0.2582  0.7962  -0.0132  0.0171    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredPlant_stage)
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
# "a"                           "a"                           "a" 


### 8.9 Rotation_cycles2
P_starvation_regulation_filteredRotation_cycles2 <- subset(P_starvation_regulation, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
P_starvation_regulation_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(P_starvation_regulation_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredRotation_cycles2 %>%
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

overall_model_P_starvation_regulation_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#   8.7931  -17.5862   -5.5862   16.0043   -5.2668   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0370     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 2156.7113, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 5.3131, p-val = 0.3789
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1      -0.0000  0.0076  -0.0025  0.9980  -0.0150  0.0150    
# Rotation_cycles2D1-3    -0.0056  0.0077  -0.7184  0.4725  -0.0207  0.0096    
# Rotation_cycles2D10      0.0046  0.0140   0.3269  0.7438  -0.0229  0.0321    
# Rotation_cycles2D3-5     0.0097  0.0096   1.0085  0.3132  -0.0091  0.0285    
# Rotation_cycles2D5-10    0.0076  0.0093   0.8231  0.4105  -0.0106  0.0259    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredRotation_cycles2)
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
#         "a"                   "b"                  "ab"                  "ab"                  "ab"                    "a" 


### 8.10 Duration2
P_starvation_regulation_filteredDuration2 <- subset(P_starvation_regulation, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
P_starvation_regulation_filteredDuration2$Duration2 <- droplevels(factor(P_starvation_regulation_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredDuration2 %>%
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

overall_model_P_starvation_regulation_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  11.4497  -22.8994  -12.8994    5.1111  -12.6730   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0373     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 2152.9588, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 2.2836, p-val = 0.6838
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5     -0.0050  0.0080  -0.6229  0.5333  -0.0207  0.0107    
# Duration2D11-20    0.0077  0.0101   0.7662  0.4436  -0.0120  0.0274    
# Duration2D20-40    0.0050  0.0112   0.4458  0.6558  -0.0169  0.0269    
# Duration2D6-10     0.0146  0.0108   1.3476  0.1778  -0.0066  0.0358  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredDuration2)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredDuration2)
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
P_starvation_regulation_filteredSpeciesRichness2 <- subset(P_starvation_regulation, SpeciesRichness2 %in% c("R2", "R3"))
#
P_starvation_regulation_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(P_starvation_regulation_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredSpeciesRichness2 %>%
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

overall_model_P_starvation_regulation_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredSpeciesRichness2)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  -3.2647    6.5294   12.5294   23.2573   12.6218   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0368     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 264) = 2180.4192, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.0245, p-val = 0.9878
# 
# Model Results:
# 
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2    0.0003  0.0059  0.0451  0.9640  -0.0113  0.0118    
# SpeciesRichness2R3    0.0010  0.0073  0.1342  0.8933  -0.0133  0.0153  


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredSpeciesRichness2)
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
P_starvation_regulation_filteredBulk_Rhizosphere <- subset(P_starvation_regulation, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
P_starvation_regulation_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(P_starvation_regulation_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredBulk_Rhizosphere %>%
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

overall_model_P_starvation_regulation_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  15.9180  -31.8360  -25.8360  -15.0076  -25.7468   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0368     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 2212.5831, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 5.5030, p-val = 0.0638
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere    0.0050  0.0058   0.8513  0.3946  -0.0065  0.0164    
# Bulk_RhizosphereRhizosphere       -0.0040  0.0063  -0.6409  0.5216  -0.0163  0.0082     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredBulk_Rhizosphere)
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
P_starvation_regulation_filteredSoil_texture <- subset(P_starvation_regulation, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
P_starvation_regulation_filteredSoil_texture$Soil_texture <- droplevels(factor(P_starvation_regulation_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredSoil_texture %>%
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

overall_model_P_starvation_regulation_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 225; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -64.4847  128.9694  136.9694  150.5801  137.1538   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0370     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 222) = 1610.0196, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.4652, p-val = 0.4816
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse    0.0160  0.0145   1.1031  0.2700  -0.0124  0.0444    
# Soil_textureFine     -0.0176  0.0158  -1.1121  0.2661  -0.0486  0.0134    
# Soil_textureMedium    0.0010  0.0094   0.1076  0.9143  -0.0175  0.0195    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredSoil_texture)
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
P_starvation_regulation_filteredTillage <- subset(P_starvation_regulation, Tillage %in% c("Tillage", "No_tillage"))
#
P_starvation_regulation_filteredTillage$Tillage <- droplevels(factor(P_starvation_regulation_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredTillage %>%
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

overall_model_P_starvation_regulation_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredTillage)
# Multivariate Meta-Analysis Model (k = 151; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  250.1705  -500.3409  -494.3409  -485.3291  -494.1754   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0025  0.0496     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 149) = 326.4178, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.5274, p-val = 0.4659
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0157  0.0145  -1.0825  0.2790  -0.0443  0.0128    
# TillageTillage      -0.0175  0.0145  -1.2080  0.2270  -0.0458  0.0109   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredTillage)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredTillage)
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
#               "a"               "a" 


### 8.15 Straw_retention
P_starvation_regulation_filteredStraw_retention <- subset(P_starvation_regulation, Straw_retention %in% c("Retention", "No_retention"))
#
P_starvation_regulation_filteredStraw_retention$Straw_retention <- droplevels(factor(P_starvation_regulation_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredStraw_retention %>%
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

overall_model_P_starvation_regulation_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 55; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -273.0209   546.0418   552.0418   557.9527   552.5316   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0024  0.0485     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 1439.5121, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.6436, p-val = 0.7248
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention   -0.0004  0.0162  -0.0239  0.9809  -0.0321  0.0314    
# Straw_retentionRetention      -0.0073  0.0154  -0.4747  0.6350  -0.0376  0.0229    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredStraw_retention)
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
P_starvation_regulation_filteredPrimer <- subset(P_starvation_regulation, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
P_starvation_regulation_filteredPrimer$Primer <- droplevels(factor(P_starvation_regulation_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- P_starvation_regulation_filteredPrimer %>%
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

overall_model_P_starvation_regulation_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = P_starvation_regulation_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_P_starvation_regulation_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 265; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  -9.1940   18.3880   26.3880   40.6614   26.5437   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0016  0.0402     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 262) = 2014.4516, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 0.2982, p-val = 0.9604
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4    0.0041  0.0082   0.5014  0.6161  -0.0119  0.0201    
# PrimerV4       0.0026  0.0152   0.1729  0.8628  -0.0271  0.0323    
# PrimerV4-V5   -0.0019  0.0149  -0.1298  0.8967  -0.0311  0.0273   
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_starvation_regulation_filteredPrimer)
vcov_rotation <- vcov(overall_model_P_starvation_regulation_filteredPrimer)
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
P_starvation_regulation$Wr <- 1 / P_starvation_regulation$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = P_starvation_regulation)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = P_starvation_regulation)

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
# Model1   Model1 -668.3871 -646.6865 340.1935
# Model2   Model2 -668.6908 -646.9902 340.3454
# Model3   Model3 -668.6911 -646.9905 340.3455
# Model4   Model4 -668.3529 -646.6522 340.1764
# Model5   Model5 -668.0354 -646.3347 340.0177
# Model6   Model6 -668.0074 -646.3068 340.0037
# Model7   Model7 -668.7243 -647.0237 340.3622###################################
# Model8   Model8 -668.6558 -646.9552 340.3279
# Model9   Model9 -629.8735 -593.7058 324.9368
# Model10 Model10 -629.4934 -593.3257 324.7467

##### Model 7 is the best model
summary(Model7)
# Number of obs: 275, groups:  StudyID, 50
anova(Model7)
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq  Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.022491 0.022491     1  98.997  0.0039 0.9504
# scale(SpeciesRichness)      0.174757 0.174757     1 137.314  0.0302 0.8622
# scale(log(Duration))        0.013392 0.013392     1  94.242  0.0023 0.9617


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(pHCK) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 5.8983  5.8983     1 38.889  0.4315 0.5151
# scale(SpeciesRichness)      0.0004  0.0004     1 79.815  0.0000 0.9956
# scale(log(Duration))        5.1476  5.1476     1 44.526  0.3766 0.5426
# scale(pHCK)                 0.0013  0.0013     1 49.953  0.0001 0.9921

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.6179  0.6179     1 25.066  0.0417 0.8398
# scale(SpeciesRichness)      4.9652  4.9652     1 44.112  0.3354 0.5654
# scale(log(Duration))        0.0197  0.0197     1 35.261  0.0013 0.9711
# scale(SOCCK)                0.9659  0.9659     1  3.429  0.0653 0.8129

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(TNCK) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.92335 0.92335     1 23.058  0.5257 0.4757
# scale(SpeciesRichness)      1.09164 1.09164     1 29.032  0.6215 0.4369
# scale(log(Duration))        0.76737 0.76737     1 24.814  0.4369 0.5147
# scale(TNCK)                 0.00691 0.00691     1 10.223  0.0039 0.9512


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.02693 0.02693     1  9.8758  0.0155 0.9035
# scale(SpeciesRichness)      1.14789 1.14789     1 30.2333  0.6595 0.4231
# scale(log(Duration))        0.37603 0.37603     1 14.0119  0.2160 0.6492
# scale(NO3CK)                0.29121 0.29121     1 17.7294  0.1673 0.6874

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.25663 0.25663     1 10.1920  0.1502 0.7063
# scale(SpeciesRichness)      1.19776 1.19776     1 30.0695  0.7011 0.4090
# scale(log(Duration))        0.86755 0.86755     1 13.2910  0.5078 0.4884
# scale(NH4CK)                0.43727 0.43727     1  8.8801  0.2560 0.6252

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(APCK) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  2.145   2.145     1 21.637  0.1247 0.7274
# scale(SpeciesRichness)       1.588   1.588     1 58.173  0.0924 0.7623
# scale(log(Duration))         9.948   9.948     1 31.268  0.5786 0.4526
# scale(APCK)                 40.609  40.609     1 24.742  2.3619 0.1370

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(AKCK) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 16.620  16.620     1    34  0.6144 0.4386
# scale(SpeciesRichness)       0.562   0.562     1    34  0.0208 0.8863
# scale(log(Duration))        16.984  16.984     1    34  0.6278 0.4337
# scale(AKCK)                 58.710  58.710     1    34  2.1702 0.1499

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(ANCK) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  7.063   7.063     1    26  0.2236 0.6403
# scale(SpeciesRichness)       2.795   2.795     1    26  0.0885 0.7685
# scale(log(Duration))         6.894   6.894     1    26  0.2182 0.6443
# scale(ANCK)                 59.719  59.719     1    26  1.8905 0.1809

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(Latitude) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
summary(ModelLatitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq  Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.027304 0.027304     1  96.119  0.0047 0.9454
# scale(SpeciesRichness)      0.168951 0.168951     1 139.827  0.0292 0.8646
# scale(log(Duration))        0.015580 0.015580     1  91.346  0.0027 0.9588
# scale(Latitude)             0.027727 0.027727     1  71.487  0.0048 0.9450

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(Longitude) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
summary(ModelLongitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.03552 0.03552     1  99.971  0.0062 0.9376
# scale(SpeciesRichness)      0.15037 0.15037     1 137.828  0.0260 0.8720
# scale(log(Duration))        0.11276 0.11276     1  91.864  0.0195 0.8892
# scale(Longitude)            0.55977 0.55977     1  22.496  0.0970 0.7584

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
summary(ModelMAPmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.02508 0.02508     1  97.171  0.0043 0.9477
# scale(SpeciesRichness)      0.19921 0.19921     1 133.910  0.0344 0.8531
# scale(log(Duration))        0.01022 0.01022     1  92.996  0.0018 0.9666
# scale(MAPmean)              0.32045 0.32045     1  39.572  0.0553 0.8152

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(MATmean) + (1 | StudyID), weights = Wr, data = P_starvation_regulation)
summary(ModelMATmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.07857 0.07857     1  97.284  0.0136 0.9075
# scale(SpeciesRichness)      0.24834 0.24834     1 134.454  0.0429 0.8362
# scale(log(Duration))        0.06545 0.06545     1  93.233  0.0113 0.9155
# scale(MATmean)              0.71133 0.71133     1  33.590  0.1229 0.7281



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(P_starvation_regulation$SpeciesRichness)) ## n = 275
p1 <- ggplot(P_starvation_regulation, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnP_starvation_regulation275")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(P_starvation_regulation$Duration)) ## n = 275
p2 <- ggplot(P_starvation_regulation, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnP_starvation_regulation275")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(P_starvation_regulation$Rotation_cycles)) ## n = 275
p3 <- ggplot(P_starvation_regulation, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnP_starvation_regulation275")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(P_starvation_regulation$Latitude)) ## n = 275
p5 <- ggplot(P_starvation_regulation, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnP_starvation_regulation275")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(P_starvation_regulation$Longitude)) ## n = 275
p6 <- ggplot(P_starvation_regulation, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnP_starvation_regulation275")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(P_starvation_regulation$MAPmean)) ## n = 275
p7 <- ggplot(P_starvation_regulation, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnP_starvation_regulation275")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(P_starvation_regulation$MATmean)) ## n = 275
p8 <- ggplot(P_starvation_regulation, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnP_starvation_regulation275")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(P_starvation_regulation$pHCK)) ## n = 86
p9 <- ggplot(P_starvation_regulation, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnP_starvation_regulation86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(P_starvation_regulation$SOCCK)) ## n = 83
p10 <- ggplot(P_starvation_regulation, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnP_starvation_regulation83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(P_starvation_regulation$TNCK)) ## n = 52
p11 <- ggplot(P_starvation_regulation, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnP_starvation_regulation52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(P_starvation_regulation$NO3CK)) ## n = 46
p12 <- ggplot(P_starvation_regulation, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnP_starvation_regulation46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(P_starvation_regulation$NH4CK)) ## n = 45
p13<- ggplot(P_starvation_regulation, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnP_starvation_regulation45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(P_starvation_regulation$APCK)) ## n = 64
p14 <- ggplot(P_starvation_regulation, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnP_starvation_regulation64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(P_starvation_regulation$AKCK)) ## n = 39
p15 <- ggplot(P_starvation_regulation, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnP_starvation_regulation39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(P_starvation_regulation$ANCK)) ## n = 31
p16 <- ggplot(P_starvation_regulation, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnP_starvation_regulation31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(P_starvation_regulation$RRpH)) ## n = 86
p17 <- ggplot(P_starvation_regulation, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="RRpH")+
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
sum(!is.na(P_starvation_regulation$RRSOC)) ## n = 83
p18 <- ggplot(P_starvation_regulation, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="RRSOC")+
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
sum(!is.na(P_starvation_regulation$RRTN)) ## n = 52
p19 <- ggplot(P_starvation_regulation, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="RRTN")+
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
sum(!is.na(P_starvation_regulation$RRNO3)) ## n = 44
p20 <- ggplot(P_starvation_regulation, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="RRNO3")+
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
sum(!is.na(P_starvation_regulation$RRNH4)) ## n = 45
p21 <- ggplot(P_starvation_regulation, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="RRNH4")+
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
sum(!is.na(P_starvation_regulation$RRAP)) ## n = 64
p22 <- ggplot(P_starvation_regulation, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="RRAP")+
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
sum(!is.na(P_starvation_regulation$RRAK)) ## n = 39
p23 <- ggplot(P_starvation_regulation, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="RRAK")+
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
sum(!is.na(P_starvation_regulation$RRAN)) ## n = 31
p24 <- ggplot(P_starvation_regulation, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_starvation_regulation", x="RRAN")+
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
sum(!is.na(P_starvation_regulation$RRYield)) ## n = 49
p25 <- ggplot(P_starvation_regulation, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-0.34, 0.22), expand=c(0, 0)) + 
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


