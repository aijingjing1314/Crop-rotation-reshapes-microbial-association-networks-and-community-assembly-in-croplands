
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Nitrogen_fixation <- read.csv("Nitrogen_fixation.csv", fileEncoding = "latin1")
# Check data
head(Nitrogen_fixation)

# 1. The number of Obversation
total_number <- nrow(Nitrogen_fixation)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 275  

# 2. The number of Study
unique_studyid_number <- length(unique(Nitrogen_fixation$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 





#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Nitrogen_fixation_filteredLongitude_Sub <- subset(Nitrogen_fixation, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Nitrogen_fixation_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Nitrogen_fixation_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredLongitude_Sub %>%
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

overall_model_Nitrogen_fixation_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  332.0427  -664.0854  -658.0854  -647.2570  -657.9962   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0034  0.0584     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 878.4211, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.2892, p-val = 0.5249
# 
# Model Results:
# 
#                            estimate      se    zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0    0.0098  0.0100  0.9833  0.3255  -0.0098  0.0295    
# Longitude_SubLongitudeXy0    0.0113  0.0198  0.5678  0.5702  -0.0276  0.0501       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredLongitude_Sub)
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
Nitrogen_fixation_filteredMAPmean2_Sub <- subset(Nitrogen_fixation, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Nitrogen_fixation_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Nitrogen_fixation_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredMAPmean2_Sub %>%
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

overall_model_Nitrogen_fixation_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  333.3895  -666.7791  -658.7791  -644.3559  -658.6293   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0035  0.0589     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 835.2157, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 7.7973, p-val = 0.0504
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean2_SubMAP600Dao1200    0.0195  0.0111   1.7558  0.0791  -0.0023  0.0414  . 
# MAPmean2_SubMAPDy1200        0.0218  0.0207   1.0550  0.2914  -0.0187  0.0623    
# MAPmean2_SubMAPXy600        -0.0029  0.0108  -0.2673  0.7892  -0.0241  0.0183    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredMAPmean2_Sub)
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
#             "a"                      "ab"                       "b" 




### 8.3 MATmean_Sub
Nitrogen_fixation_filteredMATmean_Sub <- subset(Nitrogen_fixation, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Nitrogen_fixation_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Nitrogen_fixation_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredMATmean_Sub %>%
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

overall_model_Nitrogen_fixation_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  333.0092  -666.0184  -658.0184  -643.5952  -657.8686   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0035  0.0592     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 842.5104, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 7.5061, p-val = 0.0574
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubMAT8Dao15    0.0211  0.0126   1.6767  0.0936  -0.0036  0.0459  . 
# MATmean_SubMATDy15      0.0203  0.0164   1.2361  0.2164  -0.0119  0.0524    
# MATmean_SubMATXy8      -0.0011  0.0112  -0.0984  0.9216  -0.0231  0.0209    
   
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredMATmean_Sub)
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
#     "a"                 "ab"                  "b"




### 8.4 LegumeNonlegume
Nitrogen_fixation_filteredLegumeNonlegume <- subset(Nitrogen_fixation, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Nitrogen_fixation_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Nitrogen_fixation_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredLegumeNonlegume %>%
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

overall_model_Nitrogen_fixation_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  332.9457  -665.8914  -657.8914  -643.4682  -657.7416   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0034  0.0582     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 887.3331, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 10.6066, p-val = 0.0141
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume       -0.0050  0.0102  -0.4880  0.6256  -0.0250  0.0150    
# LegumeNonlegumeNon-legume to Legume        0.0084  0.0096   0.8733  0.3825  -0.0104  0.0271    
# LegumeNonlegumeNon-legume to Non-legume    0.0207  0.0100   2.0612  0.0393   0.0010  0.0404  * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredLegumeNonlegume)
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
    #             "a"                                     "b"                                     "b" 



### 8.5 AMnonAM
Nitrogen_fixation_filteredAMnonAM <- subset(Nitrogen_fixation, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Nitrogen_fixation_filteredAMnonAM$AMnonAM <- droplevels(factor(Nitrogen_fixation_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredAMnonAM %>%
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

overall_model_Nitrogen_fixation_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  329.6967  -659.3934  -651.3934  -636.9702  -651.2436   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0034  0.0584     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 897.9190, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.1854, p-val = 0.5348
# 
# Model Results:
# 
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM       0.0125  0.0095  1.3241  0.1855  -0.0060  0.0311    
# AMnonAMAM to nonAM    0.0051  0.0114  0.4425  0.6581  -0.0173  0.0274    
# AMnonAMnonAM to AM    0.0025  0.0184  0.1345  0.8930  -0.0337  0.0386 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredAMnonAM)
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
Nitrogen_fixation_filteredC3C4 <- subset(Nitrogen_fixation, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Nitrogen_fixation_filteredC3C4$C3C4 <- droplevels(factor(Nitrogen_fixation_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredC3C4 %>%
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

overall_model_Nitrogen_fixation_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  335.0163  -670.0326  -662.0326  -647.6241  -661.8822   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0033  0.0572     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 839.5005, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 14.3642, p-val = 0.0024
# 
# Model Results:
# 
#               estimate      se    zval    pval    ci.lb   ci.ub     
# C3C4C3 to C3    0.0046  0.0102  0.4466  0.6552  -0.0154  0.0245     
# C3C4C3 to C4    0.0107  0.0094  1.1323  0.2575  -0.0078  0.0291     
# C3C4C4 to C3    0.0262  0.0097  2.6961  0.0070   0.0071  0.0452  ** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredC3C4)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredC3C4)
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
Nitrogen_fixation_filteredAnnual_Pere <- subset(Nitrogen_fixation, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
Nitrogen_fixation_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Nitrogen_fixation_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredAnnual_Pere %>%
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

overall_model_Nitrogen_fixation_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  333.7010  -667.4019  -659.4019  -644.9787  -659.2521   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0032  0.0565     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 888.8118, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 10.9662, p-val = 0.0119
# 
# Model Results:
# 
#                                 estimate      se    zval    pval    ci.lb   ci.ub     
# Annual_PereAnnual to Annual       0.0036  0.0092  0.3938  0.6937  -0.0145  0.0218     
# Annual_PereAnnual to Perennial    0.0465  0.0149  3.1303  0.0017   0.0174  0.0756  ** 
# Annual_PerePerennial to Annual    0.0222  0.0172  1.2919  0.1964  -0.0115  0.0559     

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredAnnual_Pere)
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
   #              "a"                            "b"                           "ab"



### 8.8 Plant_stage
Nitrogen_fixation_filteredPlant_stage <- subset(Nitrogen_fixation, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
Nitrogen_fixation_filteredPlant_stage$Plant_stage <- droplevels(factor(Nitrogen_fixation_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredPlant_stage %>%
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

overall_model_Nitrogen_fixation_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredPlant_stage)
# Multivariate Meta-Analysis Model (k = 240; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  293.3770  -586.7540  -578.7540  -564.8818  -578.5816   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0033  0.0571     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 237) = 667.2632, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 12.6708, p-val = 0.0054
# 
# Model Results:
# 
#                                estimate      se    zval    pval    ci.lb   ci.ub    
# Plant_stageHarvest               0.0097  0.0105  0.9261  0.3544  -0.0109  0.0303    
# Plant_stageReproductive stage    0.0098  0.0115  0.8570  0.3915  -0.0127  0.0324    
# Plant_stageVegetative stage      0.0247  0.0110  2.2452  0.0248   0.0031  0.0463  * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredPlant_stage)
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
#    "a"                          "ab"                           "b" 


### 8.9 Rotation_cycles2
Nitrogen_fixation_filteredRotation_cycles2 <- subset(Nitrogen_fixation, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
Nitrogen_fixation_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Nitrogen_fixation_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredRotation_cycles2 %>%
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

overall_model_Nitrogen_fixation_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  329.3895  -658.7790  -646.7790  -625.1884  -646.4596   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0033  0.0577     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 902.1689, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 11.2831, p-val = 0.0460
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1       0.0122  0.0117   1.0421  0.2974  -0.0107  0.0351    
# Rotation_cycles2D1-3     0.0016  0.0118   0.1345  0.8930  -0.0216  0.0248    
# Rotation_cycles2D10     -0.0049  0.0202  -0.2411  0.8095  -0.0444  0.0347    
# Rotation_cycles2D3-5     0.0229  0.0146   1.5671  0.1171  -0.0057  0.0516    
# Rotation_cycles2D5-10    0.0152  0.0144   1.0549  0.2915  -0.0131  0.0435    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredRotation_cycles2)
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
#            "a"                   "b"                  "ab"                  "ab"                  "ab"


### 8.10 Duration2
Nitrogen_fixation_filteredDuration2 <- subset(Nitrogen_fixation, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Nitrogen_fixation_filteredDuration2$Duration2 <- droplevels(factor(Nitrogen_fixation_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredDuration2 %>%
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

overall_model_Nitrogen_fixation_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  333.9820  -667.9640  -657.9640  -639.9534  -657.7376   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0034  0.0579     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 881.7846, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 12.9339, p-val = 0.0116
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub     
# Duration2D1-5      0.0031  0.0124   0.2541  0.7995  -0.0211  0.0274     
# Duration2D11-20    0.0073  0.0150   0.4884  0.6253  -0.0221  0.0367     
# Duration2D20-40   -0.0121  0.0155  -0.7767  0.4373  -0.0425  0.0184     
# Duration2D6-10     0.0473  0.0171   2.7614  0.0058   0.0137  0.0808  **       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredDuration2)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredDuration2)
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
#        "ab"             "a"             "b"             "c" 


### 8.11 SpeciesRichness2
Nitrogen_fixation_filteredSpeciesRichness2 <- subset(Nitrogen_fixation, SpeciesRichness2 %in% c("R2", "R3"))
#
Nitrogen_fixation_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Nitrogen_fixation_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredSpeciesRichness2 %>%
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

overall_model_Nitrogen_fixation_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredSpeciesRichness2)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  316.8080  -633.6159  -627.6159  -616.8881  -627.5236   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0035  0.0595     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 264) = 889.4284, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.4590, p-val = 0.2924
# 
# Model Results:
# 
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2    0.0068  0.0094  0.7223  0.4701  -0.0116  0.0252    
# SpeciesRichness2R3    0.0168  0.0115  1.4571  0.1451  -0.0058  0.0393    



# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredSpeciesRichness2)
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
Nitrogen_fixation_filteredBulk_Rhizosphere <- subset(Nitrogen_fixation, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Nitrogen_fixation_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Nitrogen_fixation_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredBulk_Rhizosphere %>%
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

overall_model_Nitrogen_fixation_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  339.4803  -678.9606  -672.9606  -662.1322  -672.8714   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0037  0.0606     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 900.5520, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 19.1224, p-val < .0001
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere    0.0184  0.0094   1.9499  0.0512  -0.0001  0.0369  . 
# Bulk_RhizosphereRhizosphere       -0.0074  0.0101  -0.7299  0.4655  -0.0272  0.0125    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredBulk_Rhizosphere)
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
Nitrogen_fixation_filteredSoil_texture <- subset(Nitrogen_fixation, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Nitrogen_fixation_filteredSoil_texture$Soil_texture <- droplevels(factor(Nitrogen_fixation_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredSoil_texture %>%
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

overall_model_Nitrogen_fixation_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 225; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  262.8875  -525.7751  -517.7751  -504.1644  -517.5907   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0033  0.0575     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 222) = 664.4879, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.3454, p-val = 0.2265
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse   -0.0112  0.0227  -0.4923  0.6225  -0.0557  0.0333    
# Soil_textureFine      0.0461  0.0246   1.8709  0.0614  -0.0022  0.0943  . 
# Soil_textureMedium    0.0116  0.0149   0.7764  0.4375  -0.0177  0.0409 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredSoil_texture)
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
Nitrogen_fixation_filteredTillage <- subset(Nitrogen_fixation, Tillage %in% c("Tillage", "No_tillage"))
#
Nitrogen_fixation_filteredTillage$Tillage <- droplevels(factor(Nitrogen_fixation_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredTillage %>%
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

overall_model_Nitrogen_fixation_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredTillage)
# Multivariate Meta-Analysis Model (k = 151; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  154.2537  -308.5074  -302.5074  -293.4956  -302.3419   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0080  0.0893     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 149) = 434.4474, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 5.5264, p-val = 0.0631
# 
# Model Results:
# 
#                    estimate      se    zval    pval    ci.lb   ci.ub    
# TillageNo_tillage    0.0194  0.0259  0.7482  0.4543  -0.0314  0.0701    
# TillageTillage       0.0348  0.0258  1.3503  0.1769  -0.0157  0.0853    

 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredTillage)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredTillage)
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
Nitrogen_fixation_filteredStraw_retention <- subset(Nitrogen_fixation, Straw_retention %in% c("Retention", "No_retention"))
#
Nitrogen_fixation_filteredStraw_retention$Straw_retention <- droplevels(factor(Nitrogen_fixation_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredStraw_retention %>%
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

overall_model_Nitrogen_fixation_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 55; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   79.7634  -159.5268  -153.5268  -147.6159  -153.0370   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0013  0.0356     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 154.1535, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.7600, p-val = 0.2516
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention   -0.0053  0.0140  -0.3805  0.7036  -0.0327  0.0221    
# Straw_retentionRetention      -0.0183  0.0124  -1.4831  0.1380  -0.0425  0.0059    
  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredStraw_retention)
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
Nitrogen_fixation_filteredPrimer <- subset(Nitrogen_fixation, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
Nitrogen_fixation_filteredPrimer$Primer <- droplevels(factor(Nitrogen_fixation_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Nitrogen_fixation_filteredPrimer %>%
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

overall_model_Nitrogen_fixation_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Nitrogen_fixation_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_fixation_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 265; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  323.9123  -647.8246  -639.8246  -625.5512  -639.6689   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0033  0.0574     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 262) = 795.8341, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 9.9040, p-val = 0.0194
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub     
# PrimerV3-V4    0.0015  0.0115   0.1315  0.8954  -0.0211  0.0241     
# PrimerV4       0.0666  0.0216   3.0785  0.0021   0.0242  0.1090  ** 
# PrimerV4-V5   -0.0138  0.0215  -0.6400  0.5222  -0.0560  0.0284   
# #  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_fixation_filteredPrimer)
vcov_rotation <- vcov(overall_model_Nitrogen_fixation_filteredPrimer)
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
#         "a"         "b"         "a"









#### 9. Linear Mixed Effect Model
# 
Nitrogen_fixation$Wr <- 1 / Nitrogen_fixation$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Nitrogen_fixation)

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
# Model1   Model1 -704.0628 -682.3621 358.0314##################################
# Model2   Model2 -703.9256 -682.2250 357.9628
# Model3   Model3 -704.0603 -682.3596 358.0301
# Model4   Model4 -703.9174 -682.2168 357.9587
# Model5   Model5 -702.6724 -680.9718 357.3362
# Model6   Model6 -702.5403 -680.8397 357.2701
# Model7   Model7 -704.0363 -682.3356 358.0181
# Model8   Model8 -703.9191 -682.2185 357.9596
# Model9   Model9 -664.5987 -628.4310 342.2994
# Model10 Model10 -667.0108 -630.8431 343.5054

##### Model 1 is the best model
summary(Model1)
# Number of obs: 275, groups:  StudyID, 50
anova(Model1) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.2158  0.2158     1 39.527  0.0974 0.7566
# scale(SpeciesRichness) 3.8667  3.8667     1 97.601  1.7454 0.1895
# scale(Duration)        3.1777  3.1777     1 83.447  1.4344 0.2344


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.91232 0.91232     1 48.168  0.4486 0.5062
# scale(SpeciesRichness) 2.17471 2.17471     1 42.578  1.0694 0.3069
# scale(Duration)        0.00746 0.00746     1 38.365  0.0037 0.9520
# scale(pHCK)            0.00009 0.00009     1 30.608  0.0000 0.9949

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 0.0429  0.0429     1 70.687  0.0261 0.87214  
# scale(SpeciesRichness) 2.4545  2.4545     1 77.991  1.4935 0.22536  
# scale(Duration)        7.2539  7.2539     1 68.399  4.4138 0.03934 *
# scale(SOCCK)           4.3364  4.3364     1 19.950  2.6386 0.11999  

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# scale(Rotation_cycles) 1.6312  1.6312     1    47  0.3544 0.5545
# scale(SpeciesRichness) 1.0087  1.0087     1    47  0.2191 0.6419
# scale(Duration)        3.5809  3.5809     1    47  0.7779 0.3823
# scale(TNCK)            0.1860  0.1860     1    47  0.0404 0.8415


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# scale(Rotation_cycles)  2.4072  2.4072     1 11.043  1.5924 0.232991   
# scale(SpeciesRichness)  1.5085  1.5085     1 40.934  0.9979 0.323699   
# scale(Duration)        14.0829 14.0829     1 29.018  9.3158 0.004824 **
# scale(NO3CK)            3.0138  3.0138     1 39.226  1.9936 0.165851   

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 0.1114  0.1114     1  8.196  0.0717 0.79553  
# scale(SpeciesRichness) 1.1896  1.1896     1 37.575  0.7656 0.38713  
# scale(Duration)        7.5527  7.5527     1 21.206  4.8609 0.03865 *
# scale(NH4CK)           0.0070  0.0070     1  5.123  0.0045 0.94916 

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.62451 0.62451     1 58.645  0.3035 0.5838
# scale(SpeciesRichness) 1.23112 1.23112     1 40.420  0.5983 0.4437
# scale(Duration)        0.00628 0.00628     1 54.690  0.0031 0.9561
# scale(APCK)            2.36907 2.36907     1 48.790  1.1513 0.2886

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.1962  0.1962     1  5.486  0.0765 0.7922
# scale(SpeciesRichness) 3.3276  3.3276     1 33.505  1.2976 0.2627
# scale(Duration)        1.3012  1.3012     1  6.175  0.5074 0.5023
# scale(AKCK)            0.0640  0.0640     1 13.396  0.0250 0.8768

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.8933  0.8933     1  9.9138  0.3482 0.5684
# scale(SpeciesRichness) 6.2020  6.2020     1 12.7831  2.4173 0.1444
# scale(Duration)        1.3608  1.3608     1  9.7185  0.5304 0.4836
# scale(ANCK)            1.7229  1.7229     1  6.2860  0.6715 0.4425

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
summary(ModelLatitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.2151  0.2151     1  41.594  0.0979 0.7560
# scale(SpeciesRichness) 3.6088  3.6088     1 116.530  1.6423 0.2026
# scale(Duration)        3.2126  3.2126     1  85.420  1.4621 0.2299
# scale(Latitude)        0.2230  0.2230     1  15.982  0.1015 0.7542

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
summary(ModelLongitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.1329  0.1329     1  39.076  0.0603 0.8074
# scale(SpeciesRichness) 3.2731  3.2731     1 100.594  1.4845 0.2259
# scale(Duration)        3.8771  3.8771     1  84.390  1.7584 0.1884
# scale(Longitude)       0.8074  0.8074     1   8.879  0.3662 0.5602

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
summary(ModelMAPmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.0000  0.0000     1 32.091  0.0000 0.9981
# scale(SpeciesRichness) 2.8179  2.8179     1 83.189  1.2730 0.2625
# scale(Duration)        1.1561  1.1561     1 41.239  0.5222 0.4740
# scale(MAPmean)         4.4481  4.4481     1 15.000  2.0094 0.1768

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Nitrogen_fixation)
summary(ModelMATmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.0874  0.0874     1  39.152  0.0398 0.8429
# scale(SpeciesRichness) 3.5737  3.5737     1 103.910  1.6275 0.2049
# scale(Duration)        2.4602  2.4602     1  61.316  1.1204 0.2940
# scale(MATmean)         0.2838  0.2838     1  13.292  0.1293 0.7249



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Nitrogen_fixation$SpeciesRichness)) ## n = 275
p1 <- ggplot(Nitrogen_fixation, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnNitrogen_fixation275")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Nitrogen_fixation$Duration)) ## n = 275
p2 <- ggplot(Nitrogen_fixation, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnNitrogen_fixation275")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Nitrogen_fixation$Rotation_cycles)) ## n = 275
p3 <- ggplot(Nitrogen_fixation, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnNitrogen_fixation275")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Nitrogen_fixation$Latitude)) ## n = 275
p5 <- ggplot(Nitrogen_fixation, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnNitrogen_fixation275")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Nitrogen_fixation$Longitude)) ## n = 275
p6 <- ggplot(Nitrogen_fixation, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnNitrogen_fixation275")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Nitrogen_fixation$MAPmean)) ## n = 275
p7 <- ggplot(Nitrogen_fixation, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnNitrogen_fixation275")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Nitrogen_fixation$MATmean)) ## n = 275
p8 <- ggplot(Nitrogen_fixation, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnNitrogen_fixation275")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Nitrogen_fixation$pHCK)) ## n = 86
p9 <- ggplot(Nitrogen_fixation, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnNitrogen_fixation86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Nitrogen_fixation$SOCCK)) ## n = 83
p10 <- ggplot(Nitrogen_fixation, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnNitrogen_fixation83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Nitrogen_fixation$TNCK)) ## n = 52
p11 <- ggplot(Nitrogen_fixation, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnNitrogen_fixation52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Nitrogen_fixation$NO3CK)) ## n = 46
p12 <- ggplot(Nitrogen_fixation, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnNitrogen_fixation46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Nitrogen_fixation$NH4CK)) ## n = 45
p13<- ggplot(Nitrogen_fixation, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnNitrogen_fixation45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Nitrogen_fixation$APCK)) ## n = 64
p14 <- ggplot(Nitrogen_fixation, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnNitrogen_fixation64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Nitrogen_fixation$AKCK)) ## n = 39
p15 <- ggplot(Nitrogen_fixation, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnNitrogen_fixation39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Nitrogen_fixation$ANCK)) ## n = 31
p16 <- ggplot(Nitrogen_fixation, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnNitrogen_fixation31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Nitrogen_fixation$RRpH)) ## n = 86
p17 <- ggplot(Nitrogen_fixation, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="RRpH")+
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
sum(!is.na(Nitrogen_fixation$RRSOC)) ## n = 83
p18 <- ggplot(Nitrogen_fixation, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="RRSOC")+
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
sum(!is.na(Nitrogen_fixation$RRTN)) ## n = 52
p19 <- ggplot(Nitrogen_fixation, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="RRTN")+
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
sum(!is.na(Nitrogen_fixation$RRNO3)) ## n = 44
p20 <- ggplot(Nitrogen_fixation, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="RRNO3")+
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
sum(!is.na(Nitrogen_fixation$RRNH4)) ## n = 45
p21 <- ggplot(Nitrogen_fixation, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="RRNH4")+
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
sum(!is.na(Nitrogen_fixation$RRAP)) ## n = 64
p22 <- ggplot(Nitrogen_fixation, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="RRAP")+
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
sum(!is.na(Nitrogen_fixation$RRAK)) ## n = 39
p23 <- ggplot(Nitrogen_fixation, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="RRAK")+
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
sum(!is.na(Nitrogen_fixation$RRAN)) ## n = 31
p24 <- ggplot(Nitrogen_fixation, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_fixation", x="RRAN")+
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
sum(!is.na(Nitrogen_fixation$RRYield)) ## n = 49
p25 <- ggplot(Nitrogen_fixation, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-0.18, 0.42), expand=c(0, 0)) + 
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





