
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Bacteria_nodes <- read.csv("Bacteria_nodes.csv", fileEncoding = "latin1")
# Check data
head(Bacteria_nodes)

# 1. The number of Obversation
total_number <- nrow(Bacteria_nodes)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 272   

# 2. The number of Study
unique_studyid_number <- length(unique(Bacteria_nodes$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 48  


#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Bacteria_nodes_filteredLongitude_Sub <- subset(Bacteria_nodes, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Bacteria_nodes_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Bacteria_nodes_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredLongitude_Sub %>%
  group_by(Longitude_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Sub Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 LongitudeDa0           111             39
# 2 LongitudeXy0           161              9

overall_model_Bacteria_nodes_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 272; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#   1.4052   -2.8104    3.1896   13.9849    3.2798   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0087  0.0935     48     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 1411.7216, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.0928, p-val = 0.9547
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0   -0.0051  0.0168  -0.3028  0.7620  -0.0380  0.0278    
# Longitude_SubLongitudeXy0    0.0011  0.0322   0.0335  0.9733  -0.0620  0.0642    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredLongitude_Sub)
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
Bacteria_nodes_filteredMAPmean2_Sub <- subset(Bacteria_nodes, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Bacteria_nodes_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Bacteria_nodes_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredMAPmean2_Sub %>%
  group_by(MAPmean2_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
  # MAPmean2_Sub  Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 MAP600Dao1200           66             17
# 2 MAPDy1200              102             11
# 3 MAPXy600               104             21

overall_model_Bacteria_nodes_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 272; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#   4.6237   -9.2473   -1.2473   13.1315   -1.0958   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0077  0.0880     48     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 269) = 1386.6053, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 8.9945, p-val = 0.0294
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean2_SubMAP600Dao1200   -0.0329  0.0185  -1.7783  0.0753  -0.0691  0.0034  . 
# MAPmean2_SubMAPDy1200       -0.0081  0.0346  -0.2332  0.8156  -0.0758  0.0597    
# MAPmean2_SubMAPXy600         0.0231  0.0177   1.3069  0.1912  -0.0116  0.0578  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredMAPmean2_Sub)
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
#                       "a"                      "ab"                       "b" 




### 8.3 MATmean_Sub
Bacteria_nodes_filteredMATmean_Sub <- subset(Bacteria_nodes, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Bacteria_nodes_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Bacteria_nodes_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredMATmean_Sub %>%
  group_by(MATmean_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Sub Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 MAT8Dao15             46             11
# 2 MATDy15              105             15
# 3 MATXy8               121             23

overall_model_Bacteria_nodes_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 272; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#   4.9021   -9.8043   -1.8043   12.5746   -1.6528   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0076  0.0871     48     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 269) = 1405.7713, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 10.0559, p-val = 0.0181
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb    ci.ub    
# MATmean_SubMAT8Dao15   -0.0460  0.0209  -2.1978  0.0280  -0.0869  -0.0050  * 
# MATmean_SubMATDy15     -0.0012  0.0276  -0.0432  0.9655  -0.0552   0.0529    
# MATmean_SubMATXy8       0.0156  0.0174   0.8983  0.3690  -0.0185   0.0497    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredMATmean_Sub)
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
#                  "a"                 "ab"                  "b" 



### 8.5 LegumeNonlegume
Bacteria_nodes_filteredLegumeNonlegume <- subset(Bacteria_nodes, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Bacteria_nodes_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Bacteria_nodes_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredLegumeNonlegume %>%
  group_by(LegumeNonlegume) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   LegumeNonlegume          Observations Unique_StudyID
#   <fct>                           <int>          <int>
# 1 Legume to Non-legume               36             14
# 2 Non-legume to Legume              154             20
# 3 Non-legume to Non-legume           82             23

overall_model_Bacteria_nodes_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 272; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  10.8313  -21.6626  -13.6626    0.7162  -13.5111   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0089  0.0942     48     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 269) = 1393.6298, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 25.7956, p-val < .0001
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume       -0.0242  0.0183  -1.3270  0.1845  -0.0601  0.0116    
# LegumeNonlegumeNon-legume to Legume        0.0175  0.0156   1.1219  0.2619  -0.0131  0.0480    
# LegumeNonlegumeNon-legume to Non-legume   -0.0117  0.0159  -0.7356  0.4619  -0.0428  0.0194    
 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredLegumeNonlegume)
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
    #                                 "a"                                     "b"                                     "a" 




### 8.5 AMnonAM
Bacteria_nodes_filteredAMnonAM <- subset(Bacteria_nodes, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Bacteria_nodes_filteredAMnonAM$AMnonAM <- droplevels(factor(Bacteria_nodes_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             212             37
# 2 AM to nonAM           30              8
# 3 nonAM to AM           30              6

overall_model_Bacteria_nodes_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 272; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#   3.5785   -7.1571    0.8429   15.2218    0.9944   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0092  0.0957     48     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 269) = 1367.3222, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 9.2145, p-val = 0.0266
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM      -0.0155  0.0158  -0.9804  0.3269  -0.0466  0.0155    
# AMnonAMAM to nonAM   -0.0083  0.0188  -0.4392  0.6605  -0.0451  0.0286    
# AMnonAMnonAM to AM    0.0828  0.0325   2.5521  0.0107   0.0192  0.1464  * 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredAMnonAM)
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
               # "a"                "a"                "b" 


### 8.6 C3C4
Bacteria_nodes_filteredC3C4 <- subset(Bacteria_nodes, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Bacteria_nodes_filteredC3C4$C3C4 <- droplevels(factor(Bacteria_nodes_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredC3C4 %>%
  group_by(C3C4) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
  # C3C4     Observations Unique_StudyID
#   <fct>           <int>          <int>
# 1 C3 to C3           73             18
# 2 C3 to C4           63             24
# 3 C4 to C3          135             14

overall_model_Bacteria_nodes_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 271; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#   2.9106   -5.8213    2.1787   16.5427    2.3308   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0094  0.0970     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 268) = 1402.7818, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 9.5008, p-val = 0.0233
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb   ci.ub    
# C3C4C3 to C3    0.0231  0.0184   1.2557  0.2092  -0.0130  0.0592    
# C3C4C3 to C4   -0.0236  0.0171  -1.3785  0.1681  -0.0572  0.0100    
# C3C4C4 to C3    0.0015  0.0179   0.0838  0.9332  -0.0336  0.0366     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredC3C4)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredC3C4)
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
Bacteria_nodes_filteredAnnual_Pere <- subset(Bacteria_nodes, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
Bacteria_nodes_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Bacteria_nodes_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredAnnual_Pere %>%
  group_by(Annual_Pere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
  # Annual_Pere         Observations Unique_StudyID
#   <fct>                      <int>          <int>
# 1 Annual to Annual             238             38
# 2 Annual to Perennial           13              8
# 3 Perennial to Annual           21              7

overall_model_Bacteria_nodes_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 272; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  -1.3413    2.6825   10.6825   25.0614   10.8340   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0086  0.0926     48     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 269) = 1397.1316, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 0.2474, p-val = 0.9696
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual      -0.0042  0.0152  -0.2756  0.7829  -0.0341  0.0257    
# Annual_PereAnnual to Perennial   -0.0075  0.0210  -0.3575  0.7207  -0.0487  0.0337    
# Annual_PerePerennial to Annual    0.0031  0.0292   0.1077  0.9142  -0.0541  0.0604    

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredAnnual_Pere)
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
                           # "a"                            "a"                            "a" 



### 8.8 Plant_stage
Bacteria_nodes_filteredPlant_stage <- subset(Bacteria_nodes, Plant_stage %in% c("Vegetative stage","Reproductive stage","Harvest"))
#
Bacteria_nodes_filteredPlant_stage$Plant_stage <- droplevels(factor(Bacteria_nodes_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredPlant_stage %>%
  group_by(Plant_stage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Plant_stage        Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Harvest                      98             23
# 2 Maturity stage                7              5
# 3 Reproductive stage           61             11
# 4 Vegetative stage             81              7

overall_model_Bacteria_nodes_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredPlant_stage)
# # is Model (k = 240; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -10.8458   21.6917   29.6917   43.5639   29.8641   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0098  0.0990     36     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 237) = 1287.2102, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 25.3809, p-val < .0001
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb   ci.ub    
# Plant_stageHarvest               0.0064  0.0183   0.3521  0.7248  -0.0293  0.0422    
# Plant_stageReproductive stage   -0.0371  0.0208  -1.7857  0.0741  -0.0778  0.0036  . 
# Plant_stageVegetative stage      0.0232  0.0185   1.2506  0.2111  -0.0131  0.0595    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredPlant_stage)
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
           #                "a"                           "b"                           "c" 


### 8.9 Rotation_cycles2
Bacteria_nodes_filteredRotation_cycles2 <- subset(Bacteria_nodes, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10"))
#
Bacteria_nodes_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Bacteria_nodes_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredRotation_cycles2 %>%
  group_by(Rotation_cycles2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Rotation_cycles2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                         83             19
# 2 D1-3                       87             12
# 3 D10                         9              5
# 4 D3-5                       33              7
# 5 D5-10                      60             13
overall_model_Bacteria_nodes_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredRotation_cycles2)
# a-Analysis Model (k = 263; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -15.2723   30.5446   40.5446   58.3287   40.7817   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0094  0.0969     45     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 259) = 1361.1903, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 2.4963, p-val = 0.6453
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1       0.0037  0.0203   0.1831  0.8547  -0.0361  0.0436    
# Rotation_cycles2D1-3     0.0091  0.0204   0.4459  0.6557  -0.0308  0.0490    
# Rotation_cycles2D3-5    -0.0221  0.0267  -0.8286  0.4073  -0.0745  0.0302    
# Rotation_cycles2D5-10   -0.0148  0.0260  -0.5709  0.5681  -0.0657  0.0361     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredRotation_cycles2)
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
# Rotation_cycles2D1  Rotation_cycles2D1-3   Rotation_cycles2D10  Rotation_cycles2D3-5 Rotation_cycles2D5-10 
# "a"                   "a"                   "a"                   "a"                   "a" 

### 8.10 Duration2
Bacteria_nodes_filteredDuration2 <- subset(Bacteria_nodes, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Bacteria_nodes_filteredDuration2$Duration2 <- droplevels(factor(Bacteria_nodes_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredDuration2 %>%
  group_by(Duration2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Duration2 Observations Unique_StudyID
#   <fct>            <int>          <int>
# 1 D1-5               150             26
# 2 D11-20              56              9
# 3 D20-40              35              6
# 4 D6-10               31              9
overall_model_Bacteria_nodes_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 272; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#   1.7922   -3.5844    6.4156   24.3705    6.6446   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0085  0.0925     48     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 268) = 1401.9957, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 6.5769, p-val = 0.1600
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5      0.0105  0.0204   0.5159  0.6059  -0.0294  0.0504    
# Duration2D11-20   -0.0239  0.0236  -1.0142  0.3105  -0.0700  0.0223    
# Duration2D20-40    0.0219  0.0277   0.7904  0.4293  -0.0324  0.0762    
# Duration2D6-10    -0.0417  0.0264  -1.5786  0.1144  -0.0934  0.0101    
# 
# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredDuration2)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredDuration2)
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
  #          "ab"             "a"             "b"             "a" 


### 8.11 SpeciesRichness2
Bacteria_nodes_filteredSpeciesRichness2 <- subset(Bacteria_nodes, SpeciesRichness2 %in% c("R2", "R3"))
#
Bacteria_nodes_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Bacteria_nodes_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredSpeciesRichness2 %>%
  group_by(SpeciesRichness2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   SpeciesRichness2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 R2                        212             42
# 2 R3                         51              8
# 3 R4                          9              3

overall_model_Bacteria_nodes_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredSpeciesRichness2)
# sis Model (k = 263; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  -3.4955    6.9909   12.9909   23.6845   13.0843   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0072  0.0851     45     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 261) = 1350.5277, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.2019, p-val = 0.9040
# 
# Model Results:
# 
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2    0.0063  0.0142  0.4456  0.6559  -0.0215  0.0341    
# SpeciesRichness2R3    0.0065  0.0178  0.3674  0.7133  -0.0283  0.0414    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredSpeciesRichness2)
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
# SpeciesRichness2R2 SpeciesRichness2R3 Spe
#                "a"                "a"        


### 8.12 Bulk_Rhizosphere
Bacteria_nodes_filteredBulk_Rhizosphere <- subset(Bacteria_nodes, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Bacteria_nodes_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Bacteria_nodes_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non_Rhizosphere           213             36
# 2 Rhizosphere                59             22
overall_model_Bacteria_nodes_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredBulk_Rhizosphere)
# ivariate Meta-Analysis Model (k = 272; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#   2.2521   -4.5042    1.4958   12.2911    1.5861   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0082  0.0906     48     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 1410.3002, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 4.9342, p-val = 0.0848
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere   -0.0104  0.0148  -0.7013  0.4831  -0.0394  0.0186    
# Bulk_RhizosphereRhizosphere        0.0094  0.0156   0.6028  0.5466  -0.0212  0.0401    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredBulk_Rhizosphere)
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
Bacteria_nodes_filteredSoil_texture <- subset(Bacteria_nodes, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Bacteria_nodes_filteredSoil_texture$Soil_texture <- droplevels(factor(Bacteria_nodes_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredSoil_texture %>%
  group_by(Soil_texture) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
 # Soil_texture Observations Unique_StudyID
#   <fct>               <int>          <int>
# 1 Coarse                 30              7
# 2 Fine                   20              8
# 3 Medium                178             17

overall_model_Bacteria_nodes_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredSoil_texture)
# ultivariate Meta-Analysis Model (k = 228; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -13.7354   27.4708   35.4708   49.1352   35.6526   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0031  0.0558     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 225) = 1022.5668, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.2759, p-val = 0.3510
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse   -0.0206  0.0232  -0.8890  0.3740  -0.0661  0.0248    
# Soil_textureFine     -0.0462  0.0294  -1.5696  0.1165  -0.1038  0.0115    
# Soil_textureMedium   -0.0022  0.0151  -0.1484  0.8820  -0.0318  0.0273    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredSoil_texture)
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
Bacteria_nodes_filteredTillage <- subset(Bacteria_nodes, Tillage %in% c("Tillage", "No_tillage"))
#
Bacteria_nodes_filteredTillage$Tillage <- droplevels(factor(Bacteria_nodes_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredTillage %>%
  group_by(Tillage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
# Tillage    Observations Unique_StudyID
#   <fct>             <int>          <int>
# 1 No_tillage          123              7
# 2 Tillage              31              9

overall_model_Bacteria_nodes_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredTillage)
# ivariate Meta-Analysis Model (k = 154; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -22.7660   45.5321   51.5321   60.6037   51.6942   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0016  0.0399     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 152) = 667.2296, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 3.0988, p-val = 0.2124
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0252  0.0147  -1.7115  0.0870  -0.0540  0.0037  . 
# TillageTillage      -0.0206  0.0150  -1.3737  0.1695  -0.0499  0.0088     
 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredTillage)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredTillage)
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
Bacteria_nodes_filteredStraw_retention <- subset(Bacteria_nodes, Straw_retention %in% c("Retention", "No_retention"))
#
Bacteria_nodes_filteredStraw_retention$Straw_retention <- droplevels(factor(Bacteria_nodes_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredStraw_retention %>%
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
overall_model_Bacteria_nodes_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredStraw_retention)
# tivariate Meta-Analysis Model (k = 55; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  -3.1402    6.2804   12.2804   18.1913   12.7702   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0018  0.0421     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 219.8583, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 6.4171, p-val = 0.0404
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb    ci.ub    
# Straw_retentionNo_retention   -0.0484  0.0206  -2.3457  0.0190  -0.0888  -0.0080  * 
# Straw_retentionRetention      -0.0302  0.0170  -1.7763  0.0757  -0.0634   0.0031  . 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredStraw_retention)
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
Bacteria_nodes_filteredPrimer <- subset(Bacteria_nodes, Primer %in% c("V3-V4", "V4", "V4-V5", "V5-V7"))
#
Bacteria_nodes_filteredPrimer$Primer <- droplevels(factor(Bacteria_nodes_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Bacteria_nodes_filteredPrimer %>%
  group_by(Primer) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Primer Observations Unique_StudyID
#   <fct>         <int>          <int>
# 1 V3                4              2
# 2 V3-V4           126             28
# 3 V4              103              9
# 4 V4-V5            39              9
overall_model_Bacteria_nodes_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Bacteria_nodes_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Bacteria_nodes_filteredPrimer)
# tivariate Meta-Analysis Model (k = 272; method:odel (k = 268; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#   0.4296   -0.8593    7.1407   21.4597    7.2946   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0084  0.0917     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 265) = 1338.4075, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.5013, p-val = 0.6820
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4    0.0016  0.0184   0.0850  0.9323  -0.0346  0.0377    
# PrimerV4      -0.0267  0.0394  -0.6775  0.4981  -0.1039  0.0505    
# PrimerV4-V5   -0.0348  0.0342  -1.0174  0.3090  -0.1019  0.0323
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_nodes_filteredPrimer)
vcov_rotation <- vcov(overall_model_Bacteria_nodes_filteredPrimer)
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
# Outpuprint(significance_letters)
 # # PrimerV3-V4    PrimerV4 PrimerV4-V5 
 #        "a"        "a"         "a"









#### 9. Linear Mixed Effect Model
# 
Bacteria_nodes$Wr <- 1 / Bacteria_nodes$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Bacteria_nodes)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Bacteria_nodes)

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
# Model1   Model1 -326.2887 -304.6539 169.1444
# Model2   Model2 -325.6771 -304.0423 168.8385
# Model3   Model3 -322.3058 -300.6710 167.1529
# Model4   Model4 -325.5733 -303.9385 168.7867
# Model5   Model5 -325.9711 -304.3363 168.9856
# Model6   Model6 -325.1831 -303.5483 168.5915
# Model7   Model7 -326.5026 -304.8678 169.2513#########################
# Model8   Model8 -321.7244 -300.0896 166.8622
# Model9   Model9 -292.2248 -256.1668 156.1124
# Model10 Model10 -293.4698 -257.4118 156.7349

##### Model 7 is the best model
summary(Model7)
# Number of obs: 272, groups:  StudyID, 48
anova(Model7) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles)) 30.792  30.792     1 210.57  7.8454 0.005569 **
# scale(SpeciesRichness)      21.855  21.855     1 163.77  5.5683 0.019467 * 
# scale(log(Duration))        28.890  28.890     1 187.81  7.3609 0.007286 **


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
summary(ModelpH)
# Number of obs: 80, groups:  StudyID, 30
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.33340 1.33340     1 19.264  0.3246 0.5754
# scale(SpeciesRichness)      0.05640 0.05640     1 68.891  0.0137 0.9071
# scale(log(Duration))        0.03289 0.03289     1 19.408  0.0080 0.9296
# scale(pHCK)                 1.74204 1.74204     1 34.420  0.4241 0.5192

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
summary(ModelSOC)
# Number of obs: 77, groups:  StudyID, 28
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 7.5185  7.5185     1 12.514  1.5334 0.2383
# scale(SpeciesRichness)      4.2336  4.2336     1 54.166  0.8635 0.3569
# scale(log(Duration))        3.4542  3.4542     1 12.778  0.7045 0.4167
# scale(SOCCK)                2.3934  2.3934     1  2.931  0.4881 0.5361

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
summary(ModelTN)
# Number of obs: 46, groups:  StudyID, 17
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.3499  0.3499     1  8.030  0.1839 0.6793
# scale(SpeciesRichness)      0.3761  0.3761     1 40.643  0.1977 0.6589
# scale(log(Duration))        0.6003  0.6003     1 12.103  0.3156 0.5845
# scale(TNCK)                 3.6805  3.6805     1 21.014  1.9350 0.1788


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
summary(ModelNO3)
# Number of obs: 38, groups:  StudyID, 13
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 2.55957 2.55957     1 26.406  3.5044 0.07232 .
# scale(SpeciesRichness)      0.05503 0.05503     1 25.312  0.0753 0.78594  
# scale(log(Duration))        2.66230 2.66230     1 32.989  3.6450 0.06497 .
# scale(NO3CK)                1.25355 1.25355     1 16.486  1.7163 0.20813  

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
summary(ModelNH4)
# Number of obs: 39, groups:  StudyID, 12
anova(ModelNH4) 
# ype III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.44407 1.44407     1 31.635  1.6761 0.2048
# scale(SpeciesRichness)      0.13126 0.13126     1 29.923  0.1524 0.6991
# scale(log(Duration))        1.98769 1.98769     1 33.199  2.3071 0.1382
# scale(NH4CK)                0.34273 0.34273     1  9.626  0.3978 0.5429

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(APCK) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
summary(ModelAP)
# Number of obs: 58, groups:  StudyID, 24
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 4.6165  4.6165     1 23.492  1.1376 0.2970
# scale(SpeciesRichness)      2.2834  2.2834     1 45.722  0.5627 0.4570
# scale(log(Duration))        0.9171  0.9171     1 22.624  0.2260 0.6391
# scale(APCK)                 0.1241  0.1241     1 35.130  0.0306 0.8622

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.0794  0.0794     1  9.1708  0.0155 0.9036
# scale(SpeciesRichness)      6.0379  6.0379     1 29.2416  1.1796 0.2863
# scale(log(Duration))        0.6993  0.6993     1 10.2325  0.1366 0.7192
# scale(AKCK)                 0.7675  0.7675     1 16.9714  0.1499 0.7034

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 10.1627 10.1627     1 19.9311  1.2937 0.2689
# scale(SpeciesRichness)      10.0893 10.0893     1 22.0630  1.2844 0.2693
# scale(log(Duration))         5.8034  5.8034     1 17.4622  0.7388 0.4017
# scale(ANCK)                  0.3575  0.3575     1  5.5411  0.0455 0.8387

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
summary(ModelLatitude)
# Number of obs: 272, groups:  StudyID, 48
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles)) 28.934  28.934     1 209.481  7.3383 0.007308 **
# scale(SpeciesRichness)      23.305  23.305     1 158.370  5.9106 0.016166 * 
# scale(log(Duration))        26.895  26.895     1 189.002  6.8211 0.009732 **
# scale(Latitude)              2.768   2.768     1  78.137  0.7020 0.404661   

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
summary(ModelLongitude)
# Number of obs: 272, groups:  StudyID, 48
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles)) 30.8123 30.8123     1 213.012  7.8723 0.005485 **
# scale(SpeciesRichness)      19.6806 19.6806     1 163.294  5.0282 0.026284 * 
# scale(log(Duration))        29.7698 29.7698     1 202.937  7.6059 0.006349 **
# scale(Longitude)             1.3296  1.3296     1  20.564  0.3397 0.566345   

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
summary(ModelMAPmean)
# Number of obs: 272, groups:  StudyID, 48
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles)) 29.1662 29.1662     1 206.674  7.4261 0.006980 **
# scale(SpeciesRichness)      20.3680 20.3680     1 157.964  5.1860 0.024112 * 
# scale(log(Duration))        27.7558 27.7558     1 181.748  7.0670 0.008552 **
# scale(MAPmean)               3.6414  3.6414     1  44.014  0.9271 0.340866   

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Bacteria_nodes)
summary(ModelMATmean)
# Number of obs: 272, groups:  StudyID, 48
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 26.1749 26.1749     1 209.876  6.6629 0.01053 *
# scale(SpeciesRichness)      19.8638 19.8638     1 161.162  5.0564 0.02589 *
# scale(log(Duration))        23.5169 23.5169     1 187.023  5.9863 0.01534 *
# scale(MATmean)               6.2517  6.2517     1  29.053  1.5914 0.21716  



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Bacteria_nodes$SpeciesRichness)) ## n = 272
p1 <- ggplot(Bacteria_nodes, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnBacteria_nodes272")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Bacteria_nodes$Duration)) ## n = 272
p2 <- ggplot(Bacteria_nodes, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnBacteria_nodes272")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Bacteria_nodes$Rotation_cycles)) ## n = 272
p3 <- ggplot(Bacteria_nodes, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnBacteria_nodes272")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Bacteria_nodes$Latitude)) ## n = 272
p5 <- ggplot(Bacteria_nodes, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnBacteria_nodes272")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Bacteria_nodes$Longitude)) ## n = 272
p6 <- ggplot(Bacteria_nodes, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnBacteria_nodes272")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Bacteria_nodes$MAPmean)) ## n = 272
p7 <- ggplot(Bacteria_nodes, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnBacteria_nodes272")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Bacteria_nodes$MATmean)) ## n = 272
p8 <- ggplot(Bacteria_nodes, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnBacteria_nodes272")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Bacteria_nodes$pHCK)) ## n = 80
p9 <- ggplot(Bacteria_nodes, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnBacteria_nodes80")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Bacteria_nodes$SOCCK)) ## n = 77
p10 <- ggplot(Bacteria_nodes, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnBacteria_nodes77")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Bacteria_nodes$TNCK)) ## n = 46
p11 <- ggplot(Bacteria_nodes, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnBacteria_nodes46")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Bacteria_nodes$NO3CK)) ## n = 38
p12 <- ggplot(Bacteria_nodes, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnBacteria_nodes38")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Bacteria_nodes$NH4CK)) ## n = 39
p13<- ggplot(Bacteria_nodes, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnBacteria_nodes39")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Bacteria_nodes$APCK)) ## n = 58
p14 <- ggplot(Bacteria_nodes, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnBacteria_nodes58")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Bacteria_nodes$AKCK)) ## n = 39
p15 <- ggplot(Bacteria_nodes, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnBacteria_nodes39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Bacteria_nodes$ANCK)) ## n = 31
p16 <- ggplot(Bacteria_nodes, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnBacteria_nodes31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Bacteria_nodes$RRpH)) ## n = 80
p17 <- ggplot(Bacteria_nodes, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="RRpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +  scale_y_continuous(
      limits = c(-1, 0.75),
      expand = c(0, 0)
    ) +
  labs(x="RRpH" , y="RR80")
p17
pdf("RRpH.pdf",width=8,height=8)
p17
dev.off() 

## RRSOC
sum(!is.na(Bacteria_nodes$RRSOC)) ## n = 77
p18 <- ggplot(Bacteria_nodes, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="RRSOC")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) + scale_y_continuous(
      limits = c(-1, 0.75),
      expand = c(0, 0)
    ) +
  labs(x="RRSOC" , y="RR77")
p18
pdf("RRSOC.pdf",width=8,height=8)
p18
dev.off() 

## RRTN
sum(!is.na(Bacteria_nodes$RRTN)) ## n = 46
p19 <- ggplot(Bacteria_nodes, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="RRTN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) + scale_y_continuous(
      limits = c(-1, 0.75),
      expand = c(0, 0)
    ) +
  labs(x="RRTN" , y="RR46")
p19
pdf("RRTN.pdf",width=8,height=8)
p19
dev.off() 

## RRNO3
sum(!is.na(Bacteria_nodes$RRNO3)) ## n = 38
p20 <- ggplot(Bacteria_nodes, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="RRNO3")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) + scale_y_continuous(
      limits = c(-1, 0.75),
      expand = c(0, 0)
    ) +
  labs(x="RRNO3" , y="RR38")
p20
pdf("RRNO3.pdf",width=8,height=8)
p20
dev.off() 

## RRNH4
sum(!is.na(Bacteria_nodes$RRNH4)) ## n = 39
p21 <- ggplot(Bacteria_nodes, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="RRNH4")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) + scale_y_continuous(
      limits = c(-1, 0.75),
      expand = c(0, 0)
    ) +
  labs(x="RRNH4" , y="RR39")
p21
pdf("RRNH4.pdf",width=8,height=8)
p21
dev.off() 

## RRAP
sum(!is.na(Bacteria_nodes$RRAP)) ## n = 58
p22 <- ggplot(Bacteria_nodes, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="RRAP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) + scale_y_continuous(
      limits = c(-1, 0.75),
      expand = c(0, 0)
    ) +
  labs(x="RRAP" , y="RR58")
p22
pdf("RRAP.pdf",width=8,height=8)
p22
dev.off() 

## RRAK
sum(!is.na(Bacteria_nodes$RRAK)) ## n = 39
p23 <- ggplot(Bacteria_nodes, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="RRAK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) + scale_y_continuous(
      limits = c(-1, 0.75),
      expand = c(0, 0)
    ) +
  labs(x="RRAK" , y="RR39")
p23
pdf("RRAK.pdf",width=8,height=8)
p23
dev.off() 

## RRAN
sum(!is.na(Bacteria_nodes$RRAN)) ## n = 31
p24 <- ggplot(Bacteria_nodes, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_nodes", x="RRAN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) + scale_y_continuous(
      limits = c(-1, 0.75),
      expand = c(0, 0)
    ) +
  labs(x="RRAN" , y="RR31")
p24
pdf("RRAN.pdf",width=8,height=8)
p24
dev.off() 

## RRYield
sum(!is.na(Bacteria_nodes$RRYield)) ## n = 49
p25 <- ggplot(Bacteria_nodes, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-0.85, 0.3), expand=c(0, 0)) + 
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

