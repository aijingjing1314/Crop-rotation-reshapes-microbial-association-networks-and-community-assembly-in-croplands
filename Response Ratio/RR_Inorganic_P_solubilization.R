
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Inorganic_P_solubilization <- read.csv("Inorganic_P_solubilization.csv", fileEncoding = "latin1")
# Check data
head(Inorganic_P_solubilization)

# 1. The number of Obversation
total_number <- nrow(Inorganic_P_solubilization)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 275  

# 2. The number of Study
unique_studyid_number <- length(unique(Inorganic_P_solubilization$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 






#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Inorganic_P_solubilization_filteredLongitude_Sub <- subset(Inorganic_P_solubilization, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Inorganic_P_solubilization_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Inorganic_P_solubilization_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredLongitude_Sub %>%
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

overall_model_Inorganic_P_solubilization_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  191.0703  -382.1406  -376.1406  -365.3122  -376.0514   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0011  0.0326     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 1828.4100, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.6230, p-val = 0.4442
# 
# Model Results:
# 
#                            estimate      se    zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0    0.0074  0.0058  1.2729  0.2031  -0.0040  0.0187    
# Longitude_SubLongitudeXy0    0.0006  0.0113  0.0533  0.9575  -0.0216  0.0228       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredLongitude_Sub)
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
Inorganic_P_solubilization_filteredMAPmean2_Sub <- subset(Inorganic_P_solubilization, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Inorganic_P_solubilization_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Inorganic_P_solubilization_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredMAPmean2_Sub %>%
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

overall_model_Inorganic_P_solubilization_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  191.9625  -383.9251  -375.9251  -361.5019  -375.7753   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0012  0.0342     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1806.9160, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 7.4736, p-val = 0.0582
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean2_SubMAP600Dao1200   -0.0055  0.0072  -0.7609  0.4467  -0.0197  0.0087    
# MAPmean2_SubMAPDy1200        0.0128  0.0131   0.9773  0.3284  -0.0129  0.0384    
# MAPmean2_SubMAPXy600         0.0131  0.0068   1.9210  0.0547  -0.0003  0.0265  .  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredMAPmean2_Sub)
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
#              "a"                      "ab"                       "b" 




### 8.3 MATmean_Sub
Inorganic_P_solubilization_filteredMATmean_Sub <- subset(Inorganic_P_solubilization, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Inorganic_P_solubilization_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Inorganic_P_solubilization_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredMATmean_Sub %>%
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

overall_model_Inorganic_P_solubilization_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  193.6940  -387.3880  -379.3880  -364.9648  -379.2382   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0324     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1763.6419, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 11.4378, p-val = 0.0096
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubMAT8Dao15   -0.0146  0.0084  -1.7428  0.0814  -0.0309  0.0018  . 
# MATmean_SubMATDy15      0.0136  0.0094   1.4398  0.1499  -0.0049  0.0321    
# MATmean_SubMATXy8       0.0099  0.0065   1.5135  0.1302  -0.0029  0.0227    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredMATmean_Sub)
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
#             "a"                  "b"                  "b" 




### 8.4 LegumeNonlegume
Inorganic_P_solubilization_filteredLegumeNonlegume <- subset(Inorganic_P_solubilization, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Inorganic_P_solubilization_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Inorganic_P_solubilization_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredLegumeNonlegume %>%
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

overall_model_Inorganic_P_solubilization_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredLegumeNonlegume)
# ultivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  192.2941  -384.5882  -376.5882  -362.1650  -376.4384   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0323     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1793.1429, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 13.3013, p-val = 0.0040
# 
# Model Results:
# 
#                                          estimate      se    zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume        0.0140  0.0059  2.3912  0.0168   0.0025  0.0255  * 
# LegumeNonlegumeNon-legume to Legume        0.0024  0.0053  0.4538  0.6500  -0.0080  0.0128    
# LegumeNonlegumeNon-legume to Non-legume    0.0047  0.0054  0.8599  0.3898  -0.0060  0.0153    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredLegumeNonlegume)
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
#     LegumeNonlegumeLegume to Non-legume     LegumeNonlegumeNon-legume to Legume LegumeNonlegumeNon-legume to Non-legume 
#                                     "a"                                     "b"                                     "b" 



### 8.5 AMnonAM
Inorganic_P_solubilization_filteredAMnonAM <- subset(Inorganic_P_solubilization, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Inorganic_P_solubilization_filteredAMnonAM$AMnonAM <- droplevels(factor(Inorganic_P_solubilization_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredAMnonAM %>%
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

overall_model_Inorganic_P_solubilization_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  189.6947  -379.3894  -371.3894  -356.9662  -371.2396   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0011  0.0329     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1705.7172, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.5032, p-val = 0.2120
# 
# Model Results:
# 
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM       0.0040  0.0057  0.7069  0.4796  -0.0071  0.0151    
# AMnonAMAM to nonAM    0.0158  0.0078  2.0396  0.0414   0.0006  0.0310  * 
# AMnonAMnonAM to AM    0.0070  0.0124  0.5671  0.5707  -0.0173  0.0313 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredAMnonAM)
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
Inorganic_P_solubilization_filteredC3C4 <- subset(Inorganic_P_solubilization, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Inorganic_P_solubilization_filteredC3C4$C3C4 <- droplevels(factor(Inorganic_P_solubilization_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredC3C4 %>%
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

overall_model_Inorganic_P_solubilization_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  193.2170  -386.4340  -378.4340  -364.0255  -378.2836   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0310     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1808.7979, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 10.2435, p-val = 0.0166
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb   ci.ub    
# C3C4C3 to C3    0.0113  0.0059   1.8929  0.0584  -0.0004  0.0229  . 
# C3C4C3 to C4    0.0081  0.0053   1.5315  0.1257  -0.0023  0.0185    
# C3C4C4 to C3   -0.0011  0.0057  -0.1897  0.8495  -0.0123  0.0101  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredC3C4)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredC3C4)
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
#         "a"          "a"          "b" 


### 8.7 Annual_Pere
Inorganic_P_solubilization_filteredAnnual_Pere <- subset(Inorganic_P_solubilization, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
Inorganic_P_solubilization_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Inorganic_P_solubilization_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredAnnual_Pere %>%
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

overall_model_Inorganic_P_solubilization_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  206.6217  -413.2434  -405.2434  -390.8202  -405.0936   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0012  0.0352     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1733.0821, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 40.2365, p-val < .0001
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb    ci.ub     
# Annual_PereAnnual to Annual       0.0146  0.0057   2.5516  0.0107   0.0034   0.0258   * 
# Annual_PereAnnual to Perennial   -0.0217  0.0076  -2.8512  0.0044  -0.0365  -0.0068  ** 
# Annual_PerePerennial to Annual   -0.0352  0.0124  -2.8303  0.0047  -0.0595  -0.0108  ** 

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredAnnual_Pere)
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
   #              "a"                            "b"                            "b" 



### 8.8 Plant_stage
Inorganic_P_solubilization_filteredPlant_stage <- subset(Inorganic_P_solubilization, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
Inorganic_P_solubilization_filteredPlant_stage$Plant_stage <- droplevels(factor(Inorganic_P_solubilization_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredPlant_stage %>%
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

overall_model_Inorganic_P_solubilization_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredPlant_stage)
# Multivariate Meta-Analysis Model (k = 240; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  131.1780  -262.3561  -254.3561  -240.4838  -254.1836   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0008  0.0289     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 237) = 1529.1451, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 9.0302, p-val = 0.0289
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb   ci.ub    
# Plant_stageHarvest               0.0023  0.0058   0.3926  0.6946  -0.0091  0.0136    
# Plant_stageReproductive stage   -0.0007  0.0080  -0.0906  0.9278  -0.0164  0.0149    
# Plant_stageVegetative stage      0.0090  0.0060   1.5104  0.1310  -0.0027  0.0208   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredPlant_stage)
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
#  "a"                          "ab"                           "b" 

### 8.9 Rotation_cycles2
Inorganic_P_solubilization_filteredRotation_cycles2 <- subset(Inorganic_P_solubilization, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
Inorganic_P_solubilization_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Inorganic_P_solubilization_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredRotation_cycles2 %>%
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

overall_model_Inorganic_P_solubilization_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  185.7538  -371.5075  -359.5075  -337.9170  -359.1881   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0011  0.0335     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 1805.1146, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 8.4514, p-val = 0.1330
# 
# Model Results:
# 
#                        estimate      se    zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1       0.0040  0.0070  0.5752  0.5651  -0.0097  0.0177    
# Rotation_cycles2D1-3     0.0079  0.0070  1.1293  0.2588  -0.0058  0.0217    
# Rotation_cycles2D10      0.0114  0.0138  0.8289  0.4071  -0.0156  0.0384    
# Rotation_cycles2D3-5     0.0004  0.0090  0.0430  0.9657  -0.0173  0.0181    
# Rotation_cycles2D5-10    0.0087  0.0088  0.9928  0.3208  -0.0085  0.0259    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredRotation_cycles2)
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
#             "a"                   "a"                   "a"                   "a"                   "a"                  "a" 


### 8.10 Duration2
Inorganic_P_solubilization_filteredDuration2 <- subset(Inorganic_P_solubilization, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Inorganic_P_solubilization_filteredDuration2$Duration2 <- droplevels(factor(Inorganic_P_solubilization_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredDuration2 %>%
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

overall_model_Inorganic_P_solubilization_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  195.8056  -391.6113  -381.6113  -363.6007  -381.3849   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0009  0.0305     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1657.0266, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 20.4269, p-val = 0.0004
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5      0.0104  0.0067   1.5508  0.1209  -0.0027  0.0234    
# Duration2D11-20   -0.0153  0.0087  -1.7536  0.0795  -0.0324  0.0018  . 
# Duration2D20-40    0.0197  0.0098   2.0052  0.0449   0.0004  0.0389  * 
# Duration2D6-10     0.0066  0.0099   0.6638  0.5068  -0.0128  0.0260    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredDuration2)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredDuration2)
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
#            "a"             "b"             "a"            "ab" 


### 8.11 SpeciesRichness2
Inorganic_P_solubilization_filteredSpeciesRichness2 <- subset(Inorganic_P_solubilization, SpeciesRichness2 %in% c("R2", "R3"))
#
Inorganic_P_solubilization_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Inorganic_P_solubilization_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredSpeciesRichness2 %>%
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

overall_model_Inorganic_P_solubilization_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredSpeciesRichness2)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  190.9843  -381.9686  -375.9686  -365.2407  -375.8763   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0312     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 264) = 1775.8583, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 23.6911, p-val < .0001
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb    ci.ub    
# SpeciesRichness2R2    0.0069  0.0051   1.3678  0.1714  -0.0030   0.0169    
# SpeciesRichness2R3   -0.0130  0.0062  -2.0843  0.0371  -0.0252  -0.0008  * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredSpeciesRichness2)
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
Inorganic_P_solubilization_filteredBulk_Rhizosphere <- subset(Inorganic_P_solubilization, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Inorganic_P_solubilization_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Inorganic_P_solubilization_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredBulk_Rhizosphere %>%
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

overall_model_Inorganic_P_solubilization_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  205.4503  -410.9005  -404.9005  -394.0721  -404.8113   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0011  0.0334     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 1818.2552, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 33.1506, p-val < .0001
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere    0.0135  0.0054   2.4889  0.0128   0.0029  0.0241  * 
# Bulk_RhizosphereRhizosphere       -0.0076  0.0058  -1.3224  0.1860  -0.0190  0.0037    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredBulk_Rhizosphere)
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
Inorganic_P_solubilization_filteredSoil_texture <- subset(Inorganic_P_solubilization, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Inorganic_P_solubilization_filteredSoil_texture$Soil_texture <- droplevels(factor(Inorganic_P_solubilization_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredSoil_texture %>%
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

overall_model_Inorganic_P_solubilization_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredSoil_texture)
# ultivariate Meta-Analysis Model (k = 225; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  113.5406  -227.0811  -219.0811  -205.4704  -218.8968   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0006  0.0250     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 222) = 1452.8714, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 8.1005, p-val = 0.0440
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse    0.0146  0.0099   1.4826  0.1382  -0.0047  0.0340    
# Soil_textureFine      0.0271  0.0127   2.1306  0.0331   0.0022  0.0521  * 
# Soil_textureMedium   -0.0077  0.0066  -1.1675  0.2430  -0.0206  0.0052      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredSoil_texture)
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
#               "ab"                "a"                "b"  


### 8.14 Tillage
Inorganic_P_solubilization_filteredTillage <- subset(Inorganic_P_solubilization, Tillage %in% c("Tillage", "No_tillage"))
#
Inorganic_P_solubilization_filteredTillage$Tillage <- droplevels(factor(Inorganic_P_solubilization_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredTillage %>%
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

overall_model_Inorganic_P_solubilization_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredTillage)
# Multivariate Meta-Analysis Model (k = 151; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   70.6853  -141.3705  -135.3705  -126.3587  -135.2050   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0007  0.0258     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 149) = 908.8312, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 24.9859, p-val < .0001
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0155  0.0083  -1.8550  0.0636  -0.0318  0.0009  . 
# TillageTillage       0.0041  0.0082   0.4972  0.6191  -0.0120  0.0201    

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredTillage)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredTillage)
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
Inorganic_P_solubilization_filteredStraw_retention <- subset(Inorganic_P_solubilization, Straw_retention %in% c("Retention", "No_retention"))
#
Inorganic_P_solubilization_filteredStraw_retention$Straw_retention <- droplevels(factor(Inorganic_P_solubilization_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredStraw_retention %>%
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

overall_model_Inorganic_P_solubilization_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 55; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  39.0209  -78.0418  -72.0418  -66.1310  -71.5520   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0015  0.0392     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 661.1192, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.2956, p-val = 0.8626
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention   -0.0004  0.0132  -0.0289  0.9769  -0.0263  0.0255    
# Straw_retentionRetention       0.0030  0.0127   0.2344  0.8146  -0.0219  0.0279    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredStraw_retention)
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
Inorganic_P_solubilization_filteredPrimer <- subset(Inorganic_P_solubilization, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
Inorganic_P_solubilization_filteredPrimer$Primer <- droplevels(factor(Inorganic_P_solubilization_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Inorganic_P_solubilization_filteredPrimer %>%
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

overall_model_Inorganic_P_solubilization_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Inorganic_P_solubilization_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Inorganic_P_solubilization_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 265; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  219.5519  -439.1038  -431.1038  -416.8304  -430.9482   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0012  0.0344     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 262) = 1690.8856, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.8972, p-val = 0.4078
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4    0.0109  0.0072   1.5233  0.1277  -0.0031  0.0250    
# PrimerV4       0.0079  0.0129   0.6136  0.5395  -0.0173  0.0331    
# PrimerV4-V5   -0.0059  0.0132  -0.4476  0.6544  -0.0318  0.0200  
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Inorganic_P_solubilization_filteredPrimer)
vcov_rotation <- vcov(overall_model_Inorganic_P_solubilization_filteredPrimer)
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
Inorganic_P_solubilization$Wr <- 1 / Inorganic_P_solubilization$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)

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
# Model1   Model1 -734.1836 -712.4829 373.0918
# Model2   Model2 -733.3690 -711.6684 372.6845
# Model3   Model3 -733.3665 -711.6659 372.6833
# Model4   Model4 -734.5144 -712.8138 373.2572###############################
# Model5   Model5 -732.3013 -710.6006 372.1506
# Model6   Model6 -732.5774 -710.8768 372.2887
# Model7   Model7 -733.0805 -711.3798 372.5402
# Model8   Model8 -733.6837 -711.9831 372.8418
# Model9   Model9 -694.7780 -658.6103 357.3890
# Model10 Model10 -701.0246 -664.8569 360.5123

##### Model 4 is the best model
summary(Model4)
# Number of obs: 275, groups:  StudyID, 50
anova(Model4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(Rotation_cycles)       4.7238  4.7238     1  58.901  1.1817 0.281440   
# scale(log(SpeciesRichness)) 28.1718 28.1718     1 266.909  7.0473 0.008415 **
# scale(Duration)              8.1293  8.1293     1 100.247  2.0336 0.156964  


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles)       0.4552  0.4552     1 25.909  0.0877 0.76954  
# scale(log(SpeciesRichness)) 19.1381 19.1381     1 70.041  3.6852 0.05897 .
# scale(Duration)              0.0242  0.0242     1 45.398  0.0047 0.94586  
# scale(pHCK)                 12.6590 12.6590     1 25.472  2.4376 0.13080 

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# scale(Rotation_cycles)       1.928   1.928     1 24.058  0.5271 0.4748182    
# scale(log(SpeciesRichness)) 46.056  46.056     1 74.135 12.5945 0.0006756 ***
# scale(Duration)              0.051   0.051     1 56.484  0.0140 0.9063124    
# scale(SOCCK)                 1.664   1.664     1 14.173  0.4549 0.5108560  

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)      0.62642 0.62642     1 16.189  0.1428 0.7104
# scale(log(SpeciesRichness)) 0.00000 0.00000     1 44.018  0.0000 0.9999
# scale(Duration)             0.06364 0.06364     1 35.638  0.0145 0.9048
# scale(TNCK)                 0.02606 0.02606     1 25.901  0.0059 0.9391


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles)      4.6565  4.6565     1 10.281  3.3441 0.09657 .
# scale(log(SpeciesRichness)) 0.9753  0.9753     1 38.615  0.7004 0.40779  
# scale(Duration)             3.3360  3.3360     1 23.086  2.3958 0.13526  
# scale(NO3CK)                0.1340  0.1340     1 10.533  0.0962 0.76246

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
# scale(Rotation_cycles)       5.0321  5.0321     1    40  3.2196 0.08032 .
# scale(log(SpeciesRichness))  0.4980  0.4980     1    40  0.3187 0.57557  
# scale(Duration)              0.8848  0.8848     1    40  0.5661 0.45621  
# scale(NH4CK)                11.1500 11.1500     1    40  7.1338 0.01088 *

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# > anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)      0.40431 0.40431     1 31.897  0.1112 0.7409
# scale(log(SpeciesRichness)) 1.97815 1.97815     1 44.682  0.5443 0.4645
# scale(Duration)             0.25344 0.25344     1 39.988  0.0697 0.7931
# scale(APCK)                 0.00070 0.00070     1 55.714  0.0002 0.9890

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)      2.4881  2.4881     1 14.126  0.5218 0.4819
# scale(log(SpeciesRichness)) 1.2377  1.2377     1 28.866  0.2596 0.6143
# scale(Duration)             8.3223  8.3223     1 17.780  1.7454 0.2032
# scale(AKCK)                 5.2686  5.2686     1 21.802  1.1050 0.3047

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)       1.6160  1.6160     1  4.0600  0.1532 0.7152
# scale(log(SpeciesRichness)) 18.8900 18.8900     1 11.3103  1.7907 0.2071
# scale(Duration)              0.9416  0.9416     1  3.9346  0.0893 0.7802
# scale(ANCK)                 28.7426 28.7426     1  2.2016  2.7247 0.2291

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
summary(ModelLatitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles)       6.2140  6.2140     1  54.825  1.5494 0.21851  
# scale(log(SpeciesRichness)) 25.7789 25.7789     1 266.411  6.4279 0.01181 *
# scale(Duration)              8.0822  8.0822     1  94.003  2.0153 0.15903  
# scale(Latitude)             10.2893 10.2893     1  89.974  2.5656 0.11272  

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
summary(ModelLongitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles)       4.4867  4.4867     1  58.420  1.1260 0.29299  
# scale(log(SpeciesRichness)) 24.9692 24.9692     1 265.042  6.2664 0.01291 *
# scale(Duration)             11.9638 11.9638     1 112.902  3.0025 0.08587 .
# scale(Longitude)             5.0519  5.0519     1  31.664  1.2679 0.26863  

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
summary(ModelMAPmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(Rotation_cycles)       4.4129  4.4129     1  55.710  1.1029 0.298171   
# scale(log(SpeciesRichness)) 27.6672 27.6672     1 263.903  6.9146 0.009051 **
# scale(Duration)              7.9348  7.9348     1  88.038  1.9831 0.162590   
# scale(MAPmean)               0.0029  0.0029     1  59.791  0.0007 0.978650  

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Inorganic_P_solubilization)
summary(ModelMATmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(Rotation_cycles)       7.5748  7.5748     1  50.430  1.8866 0.175661   
# scale(log(SpeciesRichness)) 29.3324 29.3324     1 262.899  7.3056 0.007322 **
# scale(Duration)             10.5744 10.5744     1  84.501  2.6337 0.108344   
# scale(MATmean)               7.4773  7.4773     1  37.851  1.8623 0.180416   


############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Inorganic_P_solubilization$SpeciesRichness)) ## n = 275
p1 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnInorganic_P_solubilization275")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Inorganic_P_solubilization$Duration)) ## n = 275
p2 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnInorganic_P_solubilization275")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Inorganic_P_solubilization$Rotation_cycles)) ## n = 275
p3 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnInorganic_P_solubilization275")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Inorganic_P_solubilization$Latitude)) ## n = 275
p5 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnInorganic_P_solubilization275")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Inorganic_P_solubilization$Longitude)) ## n = 275
p6 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnInorganic_P_solubilization275")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Inorganic_P_solubilization$MAPmean)) ## n = 275
p7 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnInorganic_P_solubilization275")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Inorganic_P_solubilization$MATmean)) ## n = 275
p8 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnInorganic_P_solubilization275")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Inorganic_P_solubilization$pHCK)) ## n = 86
p9 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnInorganic_P_solubilization86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Inorganic_P_solubilization$SOCCK)) ## n = 83
p10 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnInorganic_P_solubilization83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Inorganic_P_solubilization$TNCK)) ## n = 52
p11 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnInorganic_P_solubilization52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Inorganic_P_solubilization$NO3CK)) ## n = 46
p12 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnInorganic_P_solubilization46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Inorganic_P_solubilization$NH4CK)) ## n = 45
p13<- ggplot(Inorganic_P_solubilization, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnInorganic_P_solubilization45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Inorganic_P_solubilization$APCK)) ## n = 64
p14 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnInorganic_P_solubilization64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Inorganic_P_solubilization$AKCK)) ## n = 39
p15 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnInorganic_P_solubilization39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Inorganic_P_solubilization$ANCK)) ## n = 31
p16 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnInorganic_P_solubilization31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Inorganic_P_solubilization$RRpH)) ## n = 86
p17 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="RRpH")+
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
sum(!is.na(Inorganic_P_solubilization$RRSOC)) ## n = 83
p18 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="RRSOC")+
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
sum(!is.na(Inorganic_P_solubilization$RRTN)) ## n = 52
p19 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="RRTN")+
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
sum(!is.na(Inorganic_P_solubilization$RRNO3)) ## n = 44
p20 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="RRNO3")+
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
sum(!is.na(Inorganic_P_solubilization$RRNH4)) ## n = 45
p21 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="RRNH4")+
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
sum(!is.na(Inorganic_P_solubilization$RRAP)) ## n = 64
p22 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="RRAP")+
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
sum(!is.na(Inorganic_P_solubilization$RRAK)) ## n = 39
p23 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="RRAK")+
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
sum(!is.na(Inorganic_P_solubilization$RRAN)) ## n = 31
p24 <- ggplot(Inorganic_P_solubilization, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnInorganic_P_solubilization", x="RRAN")+
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
sum(!is.na(Inorganic_P_solubilization$RRYield)) ## n = 49
p25 <- ggplot(Inorganic_P_solubilization, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-0.19, 0.22), expand=c(0, 0)) + 
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




