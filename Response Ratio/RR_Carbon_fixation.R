
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Carbon_fixation <- read.csv("Carbon_fixation.csv", fileEncoding = "latin1")
# Check data
head(Carbon_fixation)

# 1. The number of Obversation
total_number <- nrow(Carbon_fixation)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 275  

# 2. The number of Study
unique_studyid_number <- length(unique(Carbon_fixation$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 



#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Carbon_fixation_filteredLongitude_Sub <- subset(Carbon_fixation, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Carbon_fixation_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Carbon_fixation_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredLongitude_Sub %>%
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

overall_model_Carbon_fixation_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Carbon_fixation_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  329.9951  -659.9903  -653.9903  -643.1619  -653.9011   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0024  0.0487     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 1189.8532, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.5070, p-val = 0.7761
# 
# Model Results:
# 
#                            estimate      se    zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0    0.0053  0.0082  0.6422  0.5208  -0.0108  0.0213    
# Longitude_SubLongitudeXy0    0.0051  0.0166  0.3075  0.7584  -0.0275  0.0377       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredLongitude_Sub)
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
Carbon_fixation_filteredMAPmean2_Sub <- subset(Carbon_fixation, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Carbon_fixation_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Carbon_fixation_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredMAPmean2_Sub %>%
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

overall_model_Carbon_fixation_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Carbon_fixation_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  335.4171  -670.8342  -662.8342  -648.4110  -662.6844   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0027  0.0524     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1191.9736, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 15.1428, p-val = 0.0017
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub     
# MAPmean2_SubMAP600Dao1200    0.0269  0.0104   2.5888  0.0096   0.0065  0.0473  ** 
# MAPmean2_SubMAPDy1200        0.0062  0.0181   0.3450  0.7301  -0.0292  0.0417     
# MAPmean2_SubMAPXy600        -0.0136  0.0100  -1.3632  0.1728  -0.0331  0.0059     


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredMAPmean2_Sub)
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
#                       "a"                       "ab"                       "b" 




### 8.3 MATmean_Sub
Carbon_fixation_filteredMATmean_Sub <- subset(Carbon_fixation, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Carbon_fixation_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Carbon_fixation_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredMATmean_Sub %>%
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

overall_model_Carbon_fixation_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Carbon_fixation_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  334.8300  -669.6600  -661.6600  -647.2368  -661.5102   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0028  0.0526     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1198.3831, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 14.4440, p-val = 0.0024
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubMAT8Dao15    0.0284  0.0122   2.3186  0.0204   0.0044  0.0524  * 
# MATmean_SubMATDy15      0.0165  0.0142   1.1654  0.2439  -0.0113  0.0444    
# MATmean_SubMATXy8      -0.0120  0.0101  -1.1919  0.2333  -0.0317  0.0077      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredMATmean_Sub)
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
#                "a"                 "ab"                  "b" 




### 8.4 LegumeNonlegume
Carbon_fixation_filteredLegumeNonlegume <- subset(Carbon_fixation, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Carbon_fixation_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Carbon_fixation_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredLegumeNonlegume %>%
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

overall_model_Carbon_fixation_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Carbon_fixation_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  330.7557  -661.5114  -653.5114  -639.0882  -653.3616   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0022  0.0474     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1200.1756, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 10.4682, p-val = 0.0150
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume        0.0104  0.0082   1.2711  0.2037  -0.0056  0.0264    
# LegumeNonlegumeNon-legume to Legume        0.0125  0.0076   1.6563  0.0977  -0.0023  0.0273  . 
# LegumeNonlegumeNon-legume to Non-legume   -0.0031  0.0077  -0.3974  0.6911  -0.0182  0.0121    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredLegumeNonlegume)
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
# LegumeNonlegumeLegume to Non-legume     LegumeNonlegumeNon-legume to Legume 
#                                     "a"                                     "a" 
# LegumeNonlegumeNon-legume to Non-legume 
                                    # "b" 




### 8.5 AMnonAM
Carbon_fixation_filteredAMnonAM <- subset(Carbon_fixation, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Carbon_fixation_filteredAMnonAM$AMnonAM <- droplevels(factor(Carbon_fixation_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredAMnonAM %>%
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

overall_model_Carbon_fixation_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Carbon_fixation_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredAMnonAM)
# ltivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  336.0665  -672.1330  -664.1330  -649.7098  -663.9832   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0025  0.0503     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1170.4200, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 17.5020, p-val = 0.0006
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb    ci.ub     
# AMnonAMAM to AM       0.0119  0.0082   1.4472  0.1478  -0.0042   0.0281     
# AMnonAMAM to nonAM   -0.0300  0.0115  -2.6120  0.0090  -0.0525  -0.0075  ** 
# AMnonAMnonAM to AM    0.0079  0.0183   0.4306  0.6668  -0.0280   0.0437     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredAMnonAM)
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
   #      "a"                "b"               "ab" 


### 8.6 C3C4
Carbon_fixation_filteredC3C4 <- subset(Carbon_fixation, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Carbon_fixation_filteredC3C4$C3C4 <- droplevels(factor(Carbon_fixation_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredC3C4 %>%
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

overall_model_Carbon_fixation_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Carbon_fixation_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredC3C4)
# ltivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  333.8877  -667.7754  -659.7754  -645.3669  -659.6250   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0020  0.0448     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1158.7493, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 12.9361, p-val = 0.0048
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb   ci.ub     
# C3C4C3 to C3   -0.0020  0.0080  -0.2547  0.7990  -0.0177  0.0136     
# C3C4C3 to C4    0.0080  0.0073   1.0900  0.2757  -0.0064  0.0223     
# C3C4C4 to C3    0.0212  0.0080   2.6611  0.0078   0.0056  0.0368  ** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredC3C4)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredC3C4)
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
Carbon_fixation_filteredAnnual_Pere <- subset(Carbon_fixation, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
Carbon_fixation_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Carbon_fixation_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredAnnual_Pere %>%
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

overall_model_Carbon_fixation_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Carbon_fixation_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredAnnual_Pere)
# ltivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  328.1028  -656.2055  -648.2055  -633.7823  -648.0557   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0023  0.0481     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1153.1391, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.9646, p-val = 0.3971
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual       0.0087  0.0077   1.1250  0.2606  -0.0065  0.0238    
# Annual_PereAnnual to Perennial    0.0028  0.0117   0.2410  0.8095  -0.0201  0.0258    
# Annual_PerePerennial to Annual   -0.0177  0.0164  -1.0835  0.2786  -0.0498  0.0143    

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredAnnual_Pere)
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
Carbon_fixation_filteredPlant_stage <- subset(Carbon_fixation, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
Carbon_fixation_filteredPlant_stage$Plant_stage <- droplevels(factor(Carbon_fixation_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredPlant_stage %>%
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

overall_model_Carbon_fixation_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Carbon_fixation_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredPlant_stage)
# tivariate Meta-Analysis Model (k = 240; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  278.6945  -557.3890  -549.3890  -535.5168  -549.2166   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0378     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 237) = 994.7572, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 11.7848, p-val = 0.0082
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb    ci.ub    
# Plant_stageHarvest              -0.0012  0.0071  -0.1741  0.8618  -0.0152   0.0127    
# Plant_stageReproductive stage   -0.0187  0.0088  -2.1156  0.0344  -0.0360  -0.0014  * 
# Plant_stageVegetative stage      0.0067  0.0076   0.8911  0.3729  -0.0081   0.0216    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredPlant_stage)
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
#       "a"                           "b"                           "c" 


### 8.9 Rotation_cycles2
Carbon_fixation_filteredRotation_cycles2 <- subset(Carbon_fixation, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
Carbon_fixation_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Carbon_fixation_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredRotation_cycles2 %>%
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

overall_model_Carbon_fixation_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Carbon_fixation_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredRotation_cycles2)
# ltivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  324.3072  -648.6144  -636.6144  -615.0239  -636.2950   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0025  0.0500     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 1195.4532, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 4.9171, p-val = 0.4261
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1       0.0049  0.0099   0.4921  0.6226  -0.0146  0.0243    
# Rotation_cycles2D1-3     0.0067  0.0100   0.6649  0.5061  -0.0130  0.0263    
# Rotation_cycles2D10      0.0151  0.0181   0.8371  0.4025  -0.0203  0.0506    
# Rotation_cycles2D3-5     0.0081  0.0127   0.6373  0.5239  -0.0168  0.0329    
# Rotation_cycles2D5-10   -0.0018  0.0124  -0.1451  0.8847  -0.0261  0.0225   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredRotation_cycles2)
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
#                "a"                   "a"                   "a"                   "a"                   "a" 


### 8.10 Duration2
Carbon_fixation_filteredDuration2 <- subset(Carbon_fixation, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Carbon_fixation_filteredDuration2$Duration2 <- droplevels(factor(Carbon_fixation_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredDuration2 %>%
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

overall_model_Carbon_fixation_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Carbon_fixation_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredDuration2)
# ltivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  330.9309  -661.8617  -651.8617  -633.8512  -651.6353   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0026  0.0511     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1194.4011, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 11.0550, p-val = 0.0260
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5      0.0067  0.0107   0.6277  0.5302  -0.0143  0.0277    
# Duration2D11-20    0.0185  0.0122   1.5078  0.1316  -0.0055  0.0424    
# Duration2D20-40   -0.0161  0.0140  -1.1476  0.2511  -0.0435  0.0114    
# Duration2D6-10     0.0027  0.0128   0.2140  0.8306  -0.0223  0.0278    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredDuration2)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredDuration2)
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
#   "ab"             "a"             "b"            "ab" 


### 8.11 SpeciesRichness2
Carbon_fixation_filteredSpeciesRichness2 <- subset(Carbon_fixation, SpeciesRichness2 %in% c("R2", "R3"))
#
Carbon_fixation_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Carbon_fixation_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredSpeciesRichness2 %>%
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

overall_model_Carbon_fixation_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Carbon_fixation_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredSpeciesRichness2)
# ultivariate Meta-Analysis Model (k = 266; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  316.6109  -633.2218  -627.2218  -616.4939  -627.1295   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0024  0.0486     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 264) = 1175.7048, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.4776, p-val = 0.7876
# 
# Model Results:
# 
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2    0.0051  0.0076  0.6748  0.4998  -0.0097  0.0199    
# SpeciesRichness2R3    0.0039  0.0092  0.4221  0.6729  -0.0142  0.0219   


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredSpeciesRichness2)
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
Carbon_fixation_filteredBulk_Rhizosphere <- subset(Carbon_fixation, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Carbon_fixation_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Carbon_fixation_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredBulk_Rhizosphere %>%
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

overall_model_Carbon_fixation_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Carbon_fixation_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredBulk_Rhizosphere)
# ltivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  328.6237  -657.2474  -651.2474  -640.4190  -651.1582   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0023  0.0483     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 1176.6571, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.7329, p-val = 0.6932
# 
# Model Results:
# 
#                                  estimate      se    zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere    0.0059  0.0074  0.7984  0.4246  -0.0086  0.0205    
# Bulk_RhizosphereRhizosphere        0.0037  0.0080  0.4572  0.6476  -0.0120  0.0193      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredBulk_Rhizosphere)
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
Carbon_fixation_filteredSoil_texture <- subset(Carbon_fixation, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Carbon_fixation_filteredSoil_texture$Soil_texture <- droplevels(factor(Carbon_fixation_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredSoil_texture %>%
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

overall_model_Carbon_fixation_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Carbon_fixation_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredSoil_texture)
# ultivariate Meta-Analysis Model (k = 225; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  273.6212  -547.2425  -539.2425  -525.6318  -539.0581   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0014  0.0379     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 222) = 821.8878, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.7808, p-val = 0.2861
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse    0.0224  0.0148   1.5158  0.1296  -0.0066  0.0513    
# Soil_textureFine     -0.0158  0.0168  -0.9432  0.3456  -0.0486  0.0170    
# Soil_textureMedium   -0.0076  0.0099  -0.7704  0.4411  -0.0269  0.0117     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredSoil_texture)
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
Carbon_fixation_filteredTillage <- subset(Carbon_fixation, Tillage %in% c("Tillage", "No_tillage"))
#
Carbon_fixation_filteredTillage$Tillage <- droplevels(factor(Carbon_fixation_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredTillage %>%
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

overall_model_Carbon_fixation_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Carbon_fixation_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredTillage)
# ltivariate Meta-Analysis Model (k = 151; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  208.6885  -417.3771  -411.3771  -402.3652  -411.2116   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0009  0.0292     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 149) = 345.8303, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 3.1158, p-val = 0.2106
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0161  0.0092  -1.7530  0.0796  -0.0342  0.0019  . 
# TillageTillage      -0.0143  0.0090  -1.5913  0.1115  -0.0319  0.0033      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredTillage)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredTillage)
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
Carbon_fixation_filteredStraw_retention <- subset(Carbon_fixation, Straw_retention %in% c("Retention", "No_retention"))
#
Carbon_fixation_filteredStraw_retention$Straw_retention <- droplevels(factor(Carbon_fixation_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredStraw_retention %>%
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

overall_model_Carbon_fixation_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Carbon_fixation_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredStraw_retention)
# ultivariate Meta-Analysis Model (k = 55; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   98.6280  -197.2560  -191.2560  -185.3451  -190.7662   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0015  0.0383     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 283.7269, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.4323, p-val = 0.8056
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention   -0.0025  0.0129  -0.1938  0.8464  -0.0277  0.0227    
# Straw_retentionRetention      -0.0064  0.0123  -0.5171  0.6051  -0.0305  0.0178    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredStraw_retention)
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
Carbon_fixation_filteredPrimer <- subset(Carbon_fixation, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
Carbon_fixation_filteredPrimer$Primer <- droplevels(factor(Carbon_fixation_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Carbon_fixation_filteredPrimer %>%
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

overall_model_Carbon_fixation_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Carbon_fixation_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Carbon_fixation_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 265; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  313.6057  -627.2114  -619.2114  -604.9380  -619.0558   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0025  0.0499     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 262) = 1132.2221, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.5700, p-val = 0.3118
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4    0.0145  0.0099   1.4657  0.1427  -0.0049  0.0339    
# PrimerV4       0.0081  0.0181   0.4459  0.6556  -0.0273  0.0434    
# PrimerV4-V5   -0.0201  0.0182  -1.1058  0.2688  -0.0557  0.0155   
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Carbon_fixation_filteredPrimer)
vcov_rotation <- vcov(overall_model_Carbon_fixation_filteredPrimer)
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
Carbon_fixation$Wr <- 1 / Carbon_fixation$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Carbon_fixation)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Carbon_fixation)

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
# Model1   Model1 -736.8524 -715.1517 374.4262
# Model2   Model2 -735.4200 -713.7194 373.7100
# Model3   Model3 -737.0645 -715.3639 374.5323############################
# Model4   Model4 -736.8001 -715.0995 374.4001
# Model5   Model5 -735.1201 -713.4195 373.5601
# Model6   Model6 -735.0641 -713.3635 373.5321
# Model7   Model7 -735.4766 -713.7760 373.7383
# Model8   Model8 -737.0120 -715.3114 374.5060
# Model9   Model9 -701.9631 -665.7954 360.9816
# Model10 Model10 -702.7425 -666.5748 361.3712

##### Model 3 is the best model
summary(Model3)
# Number of obs: 275, groups:  StudyID, 50
anova(Model3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.2199  1.2199     1 129.38  0.4976 0.4818
# scale(SpeciesRichness)      0.2628  0.2628     1 232.07  0.1072 0.7436
# scale(Duration)             4.8201  4.8201     1 131.74  1.9662 0.1632


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.88156 0.88156     1 41.872  0.6753 0.4159
# scale(SpeciesRichness)      0.36940 0.36940     1 49.011  0.2830 0.5972
# scale(Duration)             0.66039 0.66039     1 70.262  0.5059 0.4793
# scale(pHCK)                 0.32194 0.32194     1 49.920  0.2466 0.6216

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.76835 0.76835     1 38.476  0.4370 0.5125
# scale(SpeciesRichness)      1.13169 1.13169     1 75.058  0.6437 0.4249
# scale(Duration)             0.85249 0.85249     1 77.957  0.4849 0.4883
# scale(SOCCK)                0.09036 0.09036     1 24.603  0.0514 0.8225

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.27052 0.27052     1 24.803  0.1490 0.7027
# scale(SpeciesRichness)      0.22060 0.22060     1 29.793  0.1215 0.7298
# scale(Duration)             0.86974 0.86974     1 46.071  0.4792 0.4923
# scale(TNCK)                 0.12437 0.12437     1 15.488  0.0685 0.7970


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 3.8808  3.8808     1  8.904  2.4974 0.1489
# scale(SpeciesRichness)      0.1531  0.1531     1 29.951  0.0985 0.7558
# scale(Duration)             3.3197  3.3197     1 36.957  2.1363 0.1523
# scale(NO3CK)                0.0002  0.0002     1 26.640  0.0001 0.9906

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 3.2351  3.2351     1 11.526  2.1527 0.1691
# scale(SpeciesRichness)      0.1916  0.1916     1 30.841  0.1275 0.7234
# scale(Duration)             2.9666  2.9666     1 38.937  1.9741 0.1679
# scale(NH4CK)                0.1411  0.1411     1 10.048  0.0939 0.7655

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# > anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.00519 0.00519     1 31.291  0.0036 0.9528
# scale(SpeciesRichness)      0.05651 0.05651     1 35.251  0.0388 0.8450
# scale(Duration)             0.06144 0.06144     1 56.705  0.0422 0.8380
# scale(APCK)                 1.75483 1.75483     1 49.535  1.2051 0.2776

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.10745 0.10745     1 13.394  0.0614 0.8081
# scale(SpeciesRichness)      0.17659 0.17659     1 18.326  0.1009 0.7544
# scale(Duration)             0.03155 0.03155     1 18.677  0.0180 0.8946
# scale(AKCK)                 0.28597 0.28597     1 11.847  0.1634 0.6933

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.04187 0.04187     1  9.2653  0.0302 0.8657
# scale(SpeciesRichness)      0.11863 0.11863     1 11.3336  0.0856 0.7751
# scale(Duration)             0.01028 0.01028     1  7.5981  0.0074 0.9336
# scale(ANCK)                 1.49131 1.49131     1  6.3849  1.0761 0.3372

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
summary(ModelLatitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.1676  1.1676     1 130.903  0.4819 0.4888
# scale(SpeciesRichness)      0.5067  0.5067     1 237.651  0.2091 0.6479
# scale(Duration)             5.1874  5.1874     1 134.171  2.1408 0.1458
# scale(Latitude)             6.3727  6.3727     1  44.965  2.6300 0.1119

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
summary(ModelLongitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.2306  1.2306     1 132.940  0.5037 0.4791
# scale(SpeciesRichness)      0.2029  0.2029     1 235.722  0.0830 0.7735
# scale(Duration)             5.0184  5.0184     1 131.036  2.0541 0.1542
# scale(Longitude)            0.2246  0.2246     1  23.469  0.0919 0.7644

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
summary(ModelMAPmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.7441  0.7441     1 135.713  0.3112 0.5779
# scale(SpeciesRichness)      0.0705  0.0705     1 235.831  0.0295 0.8638
# scale(Duration)             3.3125  3.3125     1 114.771  1.3853 0.2416
# scale(MAPmean)              6.1565  6.1565     1  39.406  2.5747 0.1166

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Carbon_fixation)
summary(ModelMATmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)  
# scale(log(Rotation_cycles))  0.5227  0.5227     1 127.83  0.2184 0.6411  
# scale(SpeciesRichness)       0.1425  0.1425     1 230.74  0.0595 0.8074  
# scale(Duration)              3.0678  3.0678     1 121.41  1.2817 0.2598  
# scale(MATmean)              12.1834 12.1834     1  29.65  5.0900 0.0316 *



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Carbon_fixation$SpeciesRichness)) ## n = 275
p1 <- ggplot(Carbon_fixation, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnCarbon_fixation275")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Carbon_fixation$Duration)) ## n = 275
p2 <- ggplot(Carbon_fixation, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnCarbon_fixation275")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Carbon_fixation$Rotation_cycles)) ## n = 275
p3 <- ggplot(Carbon_fixation, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnCarbon_fixation275")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Carbon_fixation$Latitude)) ## n = 275
p5 <- ggplot(Carbon_fixation, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnCarbon_fixation275")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Carbon_fixation$Longitude)) ## n = 275
p6 <- ggplot(Carbon_fixation, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnCarbon_fixation275")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Carbon_fixation$MAPmean)) ## n = 275
p7 <- ggplot(Carbon_fixation, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnCarbon_fixation275")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Carbon_fixation$MATmean)) ## n = 275
p8 <- ggplot(Carbon_fixation, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnCarbon_fixation275")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Carbon_fixation$pHCK)) ## n = 86
p9 <- ggplot(Carbon_fixation, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnCarbon_fixation86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Carbon_fixation$SOCCK)) ## n = 83
p10 <- ggplot(Carbon_fixation, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnCarbon_fixation83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Carbon_fixation$TNCK)) ## n = 52
p11 <- ggplot(Carbon_fixation, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnCarbon_fixation52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Carbon_fixation$NO3CK)) ## n = 46
p12 <- ggplot(Carbon_fixation, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnCarbon_fixation46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Carbon_fixation$NH4CK)) ## n = 45
p13<- ggplot(Carbon_fixation, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnCarbon_fixation45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Carbon_fixation$APCK)) ## n = 64
p14 <- ggplot(Carbon_fixation, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnCarbon_fixation64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Carbon_fixation$AKCK)) ## n = 39
p15 <- ggplot(Carbon_fixation, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnCarbon_fixation39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Carbon_fixation$ANCK)) ## n = 31
p16 <- ggplot(Carbon_fixation, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnCarbon_fixation31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Carbon_fixation$RRpH)) ## n = 86
p17 <- ggplot(Carbon_fixation, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="RRpH")+
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
sum(!is.na(Carbon_fixation$RRSOC)) ## n = 83
p18 <- ggplot(Carbon_fixation, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="RRSOC")+
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
sum(!is.na(Carbon_fixation$RRTN)) ## n = 52
p19 <- ggplot(Carbon_fixation, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="RRTN")+
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
sum(!is.na(Carbon_fixation$RRNO3)) ## n = 44
p20 <- ggplot(Carbon_fixation, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="RRNO3")+
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
sum(!is.na(Carbon_fixation$RRNH4)) ## n = 45
p21 <- ggplot(Carbon_fixation, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="RRNH4")+
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
sum(!is.na(Carbon_fixation$RRAP)) ## n = 64
p22 <- ggplot(Carbon_fixation, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="RRAP")+
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
sum(!is.na(Carbon_fixation$RRAK)) ## n = 39
p23 <- ggplot(Carbon_fixation, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="RRAK")+
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
sum(!is.na(Carbon_fixation$RRAN)) ## n = 31
p24 <- ggplot(Carbon_fixation, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnCarbon_fixation", x="RRAN")+
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
sum(!is.na(Carbon_fixation$RRYield)) ## n = 49
p25 <- ggplot(Carbon_fixation, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-0.185, 0.29), expand=c(0, 0)) + 
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
