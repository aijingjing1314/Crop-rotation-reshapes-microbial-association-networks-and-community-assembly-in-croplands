
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
ANRA <- read.csv("ANRA.csv", fileEncoding = "latin1")
# Check data
head(ANRA)

# 1. The number of Obversation
total_number <- nrow(ANRA)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 275  

# 2. The number of Study
unique_studyid_number <- length(unique(ANRA$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 


#### 8. Subgroup analysis
### 8.1 Longitude_Sub
ANRA_filteredLongitude_Sub <- subset(ANRA, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
ANRA_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(ANRA_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- ANRA_filteredLongitude_Sub %>%
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

overall_model_ANRA_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = ANRA_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -219.6394   439.2787   445.2787   456.1072   445.3680   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0070  0.0838     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 3665.9038, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.0379, p-val = 0.9812
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0   -0.0022  0.0144  -0.1554  0.8765  -0.0305  0.0260    
# Longitude_SubLongitudeXy0   -0.0034  0.0288  -0.1171  0.9068  -0.0598  0.0530      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_ANRA_filteredLongitude_Sub)
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
ANRA_filteredMAPmean2_Sub <- subset(ANRA, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
ANRA_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(ANRA_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- ANRA_filteredMAPmean2_Sub %>%
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

overall_model_ANRA_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = ANRA_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -220.7095   441.4190   449.4190   463.8422   449.5688   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0070  0.0839     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 3687.2623, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 0.7919, p-val = 0.8514
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean2_SubMAP600Dao1200    0.0039  0.0161   0.2413  0.8093  -0.0277  0.0355    
# MAPmean2_SubMAPDy1200       -0.0014  0.0303  -0.0470  0.9625  -0.0609  0.0581    
# MAPmean2_SubMAPXy600        -0.0085  0.0157  -0.5386  0.5901  -0.0393  0.0223  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_ANRA_filteredMAPmean2_Sub)
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
ANRA_filteredMATmean_Sub <- subset(ANRA, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
ANRA_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(ANRA_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- ANRA_filteredMATmean_Sub %>%
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

overall_model_ANRA_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = ANRA_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -220.9435   441.8870   449.8870   464.3102   450.0368   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0069  0.0832     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 3909.7676, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 0.8229, p-val = 0.8440
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubMAT8Dao15   -0.0041  0.0184  -0.2211  0.8250  -0.0401  0.0319    
# MATmean_SubMATDy15      0.0123  0.0233   0.5263  0.5987  -0.0334  0.0580    
# MATmean_SubMATXy8      -0.0109  0.0160  -0.6828  0.4947  -0.0423  0.0204    
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_ANRA_filteredMATmean_Sub)
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
ANRA_filteredLegumeNonlegume <- subset(ANRA, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
ANRA_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(ANRA_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- ANRA_filteredLegumeNonlegume %>%
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

overall_model_ANRA_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = ANRA_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -79.7582  159.5163  167.5163  181.9395  167.6661   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0071  0.0845     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 2733.9918, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 286.9815, p-val < .0001
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb    ci.ub      
# LegumeNonlegumeLegume to Non-legume        0.0864  0.0146   5.9072  <.0001   0.0577   0.1150  *** 
# LegumeNonlegumeNon-legume to Legume       -0.0315  0.0138  -2.2819  0.0225  -0.0585  -0.0044    * 
# LegumeNonlegumeNon-legume to Non-legume   -0.0281  0.0143  -1.9675  0.0491  -0.0560  -0.0001    * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_ANRA_filteredLegumeNonlegume)
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
ANRA_filteredAMnonAM <- subset(ANRA, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
ANRA_filteredAMnonAM$AMnonAM <- droplevels(factor(ANRA_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- ANRA_filteredAMnonAM %>%
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

overall_model_ANRA_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = ANRA_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -219.7242   439.4485   447.4485   461.8717   447.5983   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0069  0.0832     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 4147.4476, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.2188, p-val = 0.2388
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM      -0.0037  0.0138  -0.2661  0.7902  -0.0307  0.0234    
# AMnonAMAM to nonAM   -0.0241  0.0172  -1.4039  0.1603  -0.0578  0.0095    
# AMnonAMnonAM to AM    0.0292  0.0293   0.9985  0.3180  -0.0282  0.0866  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_ANRA_filteredAMnonAM)
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
ANRA_filteredC3C4 <- subset(ANRA, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
ANRA_filteredC3C4$C3C4 <- droplevels(factor(ANRA_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- ANRA_filteredC3C4 %>%
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

overall_model_ANRA_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = ANRA_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -65.6231  131.2461  139.2461  153.6546  139.3965   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0099  0.0997     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 3239.2805, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 318.7517, p-val < .0001
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb    ci.ub      
# C3C4C3 to C3   -0.0419  0.0176  -2.3779  0.0174  -0.0764  -0.0074    * 
# C3C4C3 to C4    0.0596  0.0162   3.6900  0.0002   0.0279   0.0913  *** 
# C3C4C4 to C3   -0.0768  0.0164  -4.6924  <.0001  -0.1089  -0.0447  *** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredC3C4)
vcov_rotation <- vcov(overall_model_ANRA_filteredC3C4)
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
ANRA_filteredAnnual_Pere <- subset(ANRA, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
ANRA_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(ANRA_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- ANRA_filteredAnnual_Pere %>%
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

overall_model_ANRA_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = ANRA_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -220.9808   441.9616   449.9616   464.3849   450.1115   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0071  0.0840     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 4155.6524, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.5242, p-val = 0.6767
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual       0.0009  0.0137   0.0648  0.9483  -0.0260  0.0278    
# Annual_PereAnnual to Perennial   -0.0281  0.0252  -1.1164  0.2643  -0.0776  0.0213    
# Annual_PerePerennial to Annual   -0.0017  0.0351  -0.0473  0.9623  -0.0704  0.0671 

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_ANRA_filteredAnnual_Pere)
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
ANRA_filteredPlant_stage <- subset(ANRA, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
ANRA_filteredPlant_stage$Plant_stage <- droplevels(factor(ANRA_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- ANRA_filteredPlant_stage %>%
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

overall_model_ANRA_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = ANRA_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredPlant_stage)
# Multivariate Meta-Analysis Model (k = 240; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -249.7769   499.5539   507.5539   521.4261   507.7263   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0063  0.0791     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 237) = 3623.6080, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.1378, p-val = 0.2470
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb   ci.ub    
# Plant_stageHarvest               0.0032  0.0145   0.2205  0.8254  -0.0252  0.0316    
# Plant_stageReproductive stage    0.0082  0.0166   0.4923  0.6225  -0.0243  0.0406    
# Plant_stageVegetative stage     -0.0059  0.0149  -0.3985  0.6903  -0.0351  0.0232 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_ANRA_filteredPlant_stage)
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
ANRA_filteredRotation_cycles2 <- subset(ANRA, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
ANRA_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(ANRA_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- ANRA_filteredRotation_cycles2 %>%
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

overall_model_ANRA_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = ANRA_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -224.4637   448.9274   460.9274   482.5180   461.2468   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0069  0.0833     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 3719.9045, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 3.3856, p-val = 0.6408
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1       0.0061  0.0173   0.3526  0.7244  -0.0278  0.0400    
# Rotation_cycles2D1-3    -0.0007  0.0174  -0.0412  0.9672  -0.0349  0.0334    
# Rotation_cycles2D10     -0.0217  0.0277  -0.7811  0.4348  -0.0760  0.0327    
# Rotation_cycles2D3-5    -0.0123  0.0210  -0.5868  0.5573  -0.0535  0.0288    
# Rotation_cycles2D5-10   -0.0030  0.0203  -0.1481  0.8823  -0.0429  0.0368   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_ANRA_filteredRotation_cycles2)
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
#                   "a"                   "a"                  "a"                  "a"                  "a"                     "a" 


### 8.10 Duration2
ANRA_filteredDuration2 <- subset(ANRA, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
ANRA_filteredDuration2$Duration2 <- droplevels(factor(ANRA_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- ANRA_filteredDuration2 %>%
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

overall_model_ANRA_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = ANRA_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -222.5202   445.0404   455.0404   473.0510   455.2668   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0072  0.0848     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 4067.4929, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 0.6290, p-val = 0.9598
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5      0.0020  0.0182   0.1104  0.9121  -0.0336  0.0376    
# Duration2D11-20   -0.0030  0.0213  -0.1390  0.8895  -0.0448  0.0388    
# Duration2D20-40   -0.0129  0.0225  -0.5748  0.5655  -0.0570  0.0312    
# Duration2D6-10    -0.0079  0.0241  -0.3278  0.7430  -0.0552  0.0393   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredDuration2)
vcov_rotation <- vcov(overall_model_ANRA_filteredDuration2)
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
ANRA_filteredSpeciesRichness2 <- subset(ANRA, SpeciesRichness2 %in% c("R2", "R3"))
#
ANRA_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(ANRA_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- ANRA_filteredSpeciesRichness2 %>%
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

overall_model_ANRA_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = ANRA_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredSpeciesRichness2)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -228.0662   456.1325   462.1325   472.8603   462.2248   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0066  0.0812     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 264) = 4143.3538, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 8.6760, p-val = 0.0131
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb    ci.ub    
# SpeciesRichness2R2   -0.0011  0.0130  -0.0868  0.9309  -0.0266   0.0243    
# SpeciesRichness2R3   -0.0370  0.0169  -2.1893  0.0286  -0.0701  -0.0039  * 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_ANRA_filteredSpeciesRichness2)
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
ANRA_filteredBulk_Rhizosphere <- subset(ANRA, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
ANRA_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(ANRA_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- ANRA_filteredBulk_Rhizosphere %>%
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

overall_model_ANRA_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = ANRA_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -221.0603   442.1206   448.1206   458.9490   448.2098   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0068  0.0824     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 4283.0523, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.0498, p-val = 0.9754
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere   -0.0028  0.0130  -0.2108  0.8330  -0.0283  0.0228    
# Bulk_RhizosphereRhizosphere       -0.0017  0.0140  -0.1201  0.9044  -0.0291  0.0258    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_ANRA_filteredBulk_Rhizosphere)
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
ANRA_filteredSoil_texture <- subset(ANRA, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
ANRA_filteredSoil_texture$Soil_texture <- droplevels(factor(ANRA_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- ANRA_filteredSoil_texture %>%
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

overall_model_ANRA_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = ANRA_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredSoil_texture)
# ivariate Meta-Analysis Model (k = 225; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -244.4571   488.9142   496.9142   510.5249   497.0985   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0070  0.0836     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 222) = 2948.6892, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.5983, p-val = 0.4578
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse    0.0097  0.0326   0.2964  0.7669  -0.0543  0.0736    
# Soil_textureFine      0.0516  0.0351   1.4704  0.1415  -0.0172  0.1205    
# Soil_textureMedium   -0.0126  0.0214  -0.5903  0.5550  -0.0545  0.0293   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_ANRA_filteredSoil_texture)
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
ANRA_filteredTillage <- subset(ANRA, Tillage %in% c("Tillage", "No_tillage"))
#
ANRA_filteredTillage$Tillage <- droplevels(factor(ANRA_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- ANRA_filteredTillage %>%
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

overall_model_ANRA_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = ANRA_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredTillage)
# ultivariate Meta-Analysis Model (k = 151; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  130.0378  -260.0756  -254.0756  -245.0638  -253.9101   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0033  0.0572     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 149) = 386.4748, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.4006, p-val = 0.4964
# 
# Model Results:
# 
#                    estimate      se    zval    pval    ci.lb   ci.ub    
# TillageNo_tillage    0.0058  0.0179  0.3248  0.7453  -0.0294  0.0410    
# TillageTillage       0.0157  0.0177  0.8832  0.3771  -0.0191  0.0505    
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredTillage)
vcov_rotation <- vcov(overall_model_ANRA_filteredTillage)
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
ANRA_filteredStraw_retention <- subset(ANRA, Straw_retention %in% c("Retention", "No_retention"))
#
ANRA_filteredStraw_retention$Straw_retention <- droplevels(factor(ANRA_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- ANRA_filteredStraw_retention %>%
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

overall_model_ANRA_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = ANRA_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredStraw_retention)
# ultivariate Meta-Analysis Model (k = 55; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -358.0010   716.0021   722.0021   727.9130   722.4919   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0122  0.1104     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 2551.7054, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 6.6022, p-val = 0.0368
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention    0.0376  0.0364   1.0310  0.3025  -0.0338  0.1090    
# Straw_retentionRetention      -0.0166  0.0349  -0.4773  0.6332  -0.0850  0.0517    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_ANRA_filteredStraw_retention)
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
#                         "a"                         "b" 

### 8.16 Primer
ANRA_filteredPrimer <- subset(ANRA, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
ANRA_filteredPrimer$Primer <- droplevels(factor(ANRA_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- ANRA_filteredPrimer %>%
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

overall_model_ANRA_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = ANRA_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_ANRA_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 265; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -232.1473   464.2946   472.2946   486.5679   472.4502   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0080  0.0896     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 262) = 3816.2901, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 0.5858, p-val = 0.8997
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4   -0.0044  0.0178  -0.2499  0.8027  -0.0393  0.0304    
# PrimerV4       0.0127  0.0350   0.3629  0.7167  -0.0558  0.0812    
# PrimerV4-V5   -0.0209  0.0333  -0.6258  0.5315  -0.0862  0.0445   
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_ANRA_filteredPrimer)
vcov_rotation <- vcov(overall_model_ANRA_filteredPrimer)
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
ANRA$Wr <- 1 / ANRA$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = ANRA)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = ANRA)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = ANRA)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = ANRA)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = ANRA)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = ANRA)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = ANRA)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = ANRA)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = ANRA)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = ANRA)

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
# Model1   Model1 -275.9930 -254.2924 143.9965
# Model2   Model2 -276.9195 -255.2189 144.4598########################
# Model3   Model3 -276.1213 -254.4206 144.0606
# Model4   Model4 -276.2004 -254.4998 144.1002
# Model5   Model5 -276.1580 -254.4573 144.0790
# Model6   Model6 -276.3039 -254.6032 144.1519
# Model7   Model7 -276.7711 -255.0704 144.3855
# Model8   Model8 -276.3346 -254.6340 144.1673
# Model9   Model9 -240.7266 -204.5589 130.3633
# Model10 Model10 -242.6133 -206.4456 131.3066

##### Model 2 is the best model
summary(Model2)
# Number of obs: 275, groups:  StudyID, 50
anova(Model2) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.00001 0.00001     1 166.66  0.0000 0.9992
# scale(log(SpeciesRichness)) 2.14401 2.14401     1 204.67  0.3571 0.5508
# scale(log(Duration))        0.99873 0.99873     1 194.92  0.1664 0.6838


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(pHCK) + (1 | StudyID), weights = Wr, data = ANRA)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.0114  0.0114     1 35.593  0.0008 0.9777
# scale(log(SpeciesRichness)) 4.5718  4.5718     1 78.442  0.3174 0.5748
# scale(log(Duration))        0.7847  0.7847     1 44.770  0.0545 0.8165
# scale(pHCK)                 3.3561  3.3561     1 54.304  0.2330 0.6313

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = ANRA)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.8335  0.8335     1 33.300  0.0552 0.8157
# scale(log(SpeciesRichness)) 5.5462  5.5462     1 53.110  0.3671 0.5472
# scale(log(Duration))        0.2681  0.2681     1 45.768  0.0177 0.8946
# scale(SOCCK)                6.9226  6.9226     1 16.306  0.4582 0.5080

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(TNCK) + (1 | StudyID), weights = Wr, data = ANRA)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.17886 1.17886     1  3.384  0.5507 0.5062
# scale(log(SpeciesRichness)) 2.08433 2.08433     1 37.119  0.9737 0.3301
# scale(log(Duration))        1.60400 1.60400     1  6.938  0.7493 0.4156
# scale(TNCK)                 0.16729 0.16729     1  4.711  0.0782 0.7917


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = ANRA)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.73768 1.73768     1 11.683  1.0929 0.3170
# scale(log(SpeciesRichness)) 0.01642 0.01642     1 33.038  0.0103 0.9197
# scale(log(Duration))        2.08398 2.08398     1 19.548  1.3107 0.2661
# scale(NO3CK)                0.18031 0.18031     1 15.970  0.1134 0.7407

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = ANRA)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.36434 1.36434     1 13.654  0.8609 0.3696
# scale(log(SpeciesRichness)) 0.02353 0.02353     1 33.636  0.0149 0.9037
# scale(log(Duration))        1.55235 1.55235     1 21.025  0.9795 0.3336
# scale(NH4CK)                0.29071 0.29071     1  9.569  0.1834 0.6779

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(APCK) + (1 | StudyID), weights = Wr, data = ANRA)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# > anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 10.0730 10.0730     1 27.206  0.5388 0.4692
# scale(log(SpeciesRichness))  3.2062  3.2062     1 57.650  0.1715 0.6803
# scale(log(Duration))        11.1036 11.1036     1 30.994  0.5940 0.4467
# scale(APCK)                 25.7328 25.7328     1 44.231  1.3765 0.2470

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(AKCK) + (1 | StudyID), weights = Wr, data = ANRA)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  2.0861  2.0861     1 15.248  0.0666 0.7998
# scale(log(SpeciesRichness))  1.0170  1.0170     1 32.418  0.0325 0.8582
# scale(log(Duration))         6.8176  6.8176     1 18.254  0.2176 0.6464
# scale(AKCK)                 13.8554 13.8554     1 29.250  0.4421 0.5113

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(ANCK) + (1 | StudyID), weights = Wr, data = ANRA)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  88.660  88.660     1  8.4915  2.2491 0.1699
# scale(log(SpeciesRichness))  27.211  27.211     1 25.7662  0.6903 0.4137
# scale(log(Duration))         91.166  91.166     1  7.6229  2.3126 0.1687
# scale(ANCK)                 108.724 108.724     1  6.2628  2.7580 0.1458

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(Latitude) + (1 | StudyID), weights = Wr, data = ANRA)
summary(ModelLatitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.01568 0.01568     1 170.12  0.0026 0.9593
# scale(log(SpeciesRichness)) 2.02193 2.02193     1 203.47  0.3362 0.5627
# scale(log(Duration))        0.70652 0.70652     1 200.04  0.1175 0.7322
# scale(Latitude)             0.95960 0.95960     1 103.16  0.1595 0.6904

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(Longitude) + (1 | StudyID), weights = Wr, data = ANRA)
summary(ModelLongitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.00039 0.00039     1 164.582  0.0001 0.9936
# scale(log(SpeciesRichness)) 1.94977 1.94977     1 198.825  0.3243 0.5697
# scale(log(Duration))        0.59053 0.59053     1 204.750  0.0982 0.7543
# scale(Longitude)            0.42923 0.42923     1  38.125  0.0714 0.7908

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = ANRA)
summary(ModelMAPmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.0403  0.0403     1 167.032  0.0067 0.9348
# scale(log(SpeciesRichness)) 3.0044  3.0044     1 201.747  0.5008 0.4800
# scale(log(Duration))        0.6402  0.6402     1 191.856  0.1067 0.7443
# scale(MAPmean)              5.6272  5.6272     1  78.366  0.9380 0.3358

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(MATmean) + (1 | StudyID), weights = Wr, data = ANRA)
summary(ModelMATmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.2418  0.2418     1 169.001  0.0403 0.8411
# scale(log(SpeciesRichness)) 2.7745  2.7745     1 201.895  0.4624 0.4973
# scale(log(Duration))        0.2230  0.2230     1 195.858  0.0372 0.8473
# scale(MATmean)              6.3948  6.3948     1  60.776  1.0658 0.3060



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(ANRA$SpeciesRichness)) ## n = 275
p1 <- ggplot(ANRA, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnANRA275")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(ANRA$Duration)) ## n = 275
p2 <- ggplot(ANRA, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnANRA275")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(ANRA$Rotation_cycles)) ## n = 275
p3 <- ggplot(ANRA, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnANRA275")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(ANRA$Latitude)) ## n = 275
p5 <- ggplot(ANRA, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnANRA275")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(ANRA$Longitude)) ## n = 275
p6 <- ggplot(ANRA, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnANRA275")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(ANRA$MAPmean)) ## n = 275
p7 <- ggplot(ANRA, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnANRA275")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(ANRA$MATmean)) ## n = 275
p8 <- ggplot(ANRA, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnANRA275")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(ANRA$pHCK)) ## n = 86
p9 <- ggplot(ANRA, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnANRA86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(ANRA$SOCCK)) ## n = 83
p10 <- ggplot(ANRA, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnANRA83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(ANRA$TNCK)) ## n = 52
p11 <- ggplot(ANRA, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnANRA52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(ANRA$NO3CK)) ## n = 46
p12 <- ggplot(ANRA, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnANRA46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(ANRA$NH4CK)) ## n = 45
p13<- ggplot(ANRA, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnANRA45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(ANRA$APCK)) ## n = 64
p14 <- ggplot(ANRA, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnANRA64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(ANRA$AKCK)) ## n = 39
p15 <- ggplot(ANRA, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnANRA39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(ANRA$ANCK)) ## n = 31
p16 <- ggplot(ANRA, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnANRA31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(ANRA$RRpH)) ## n = 86
p17 <- ggplot(ANRA, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="RRpH")+
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
sum(!is.na(ANRA$RRSOC)) ## n = 83
p18 <- ggplot(ANRA, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="RRSOC")+
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
sum(!is.na(ANRA$RRTN)) ## n = 52
p19 <- ggplot(ANRA, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="RRTN")+
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
sum(!is.na(ANRA$RRNO3)) ## n = 44
p20 <- ggplot(ANRA, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="RRNO3")+
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
sum(!is.na(ANRA$RRNH4)) ## n = 45
p21 <- ggplot(ANRA, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="RRNH4")+
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
sum(!is.na(ANRA$RRAP)) ## n = 64
p22 <- ggplot(ANRA, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="RRAP")+
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
sum(!is.na(ANRA$RRAK)) ## n = 39
p23 <- ggplot(ANRA, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="RRAK")+
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
sum(!is.na(ANRA$RRAN)) ## n = 31
p24 <- ggplot(ANRA, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnANRA", x="RRAN")+
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
sum(!is.na(ANRA$RRYield)) ## n = 49
p25 <- ggplot(ANRA, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-0.38, 0.45), expand=c(0, 0)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) 
p25
pdf("RRYield_swapped_limited.pdf", width=8, height=8)
p25
dev.off()


library(patchwork)

# Combine
combined_plot <- (
  p17 | p18 | p19
) / (
  p20 | p21 | p22
) / (
  p23 | p24 | p25
)

combined_plot

# Output
ggsave(
  "Combined_ANRA.pdf",
  combined_plot,
  width = 10,
  height = 10
)











