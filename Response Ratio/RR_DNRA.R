
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
DNRA <- read.csv("DNRA.csv", fileEncoding = "latin1")
# Check data
head(DNRA)

# 1. The number of Obversation
total_number <- nrow(DNRA)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 274  

# 2. The number of Study
unique_studyid_number <- length(unique(DNRA$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 





#### 8. Subgroup analysis
### 8.1 Longitude_Sub
DNRA_filteredLongitude_Sub <- subset(DNRA, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
DNRA_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(DNRA_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- DNRA_filteredLongitude_Sub %>%
  group_by(Longitude_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Sub Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 LongitudeDa0           117             41
# 2 LongitudeXy0           157              9

overall_model_DNRA_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = DNRA_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -173.4971   346.9942   352.9942   363.8116   353.0838   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0301  0.1735     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 2100.6386, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.5349, p-val = 0.7653
# 
# Model Results:
# 
#                            estimate      se    zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0    0.0054  0.0290  0.1871  0.8516  -0.0514  0.0622    
# Longitude_SubLongitudeXy0    0.0416  0.0589  0.7070  0.4795  -0.0737  0.1570     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_DNRA_filteredLongitude_Sub)
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
DNRA_filteredMAPmean2_Sub <- subset(DNRA, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
DNRA_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(DNRA_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- DNRA_filteredMAPmean2_Sub %>%
  group_by(MAPmean2_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MAPmean2_Sub  Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 MAP600Dao1200           72             19
# 2 MAPDy1200               98             11
# 3 MAPXy600               104             21

overall_model_DNRA_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = DNRA_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -170.1325   340.2650   348.2650   362.6735   348.4154   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0337  0.1836     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 2130.5320, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 8.4301, p-val = 0.0379
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean2_SubMAP600Dao1200   -0.0487  0.0364  -1.3380  0.1809  -0.1200  0.0226    
# MAPmean2_SubMAPDy1200        0.0329  0.0635   0.5174  0.6049  -0.0916  0.1574    
# MAPmean2_SubMAPXy600         0.0570  0.0348   1.6371  0.1016  -0.0112  0.1251    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_DNRA_filteredMAPmean2_Sub)
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
#                   "a"                      "ab"                       "b" 




### 8.3 MATmean_Sub
DNRA_filteredMATmean_Sub <- subset(DNRA, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
DNRA_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(DNRA_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- DNRA_filteredMATmean_Sub %>%
  group_by(MATmean_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Sub Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 MAT8Dao15             46             11
# 2 MATDy15              107             17
# 3 MATXy8               121             23

overall_model_DNRA_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = DNRA_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -168.4223   336.8446   344.8446   359.2530   344.9949   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0313  0.1769     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 2152.0310, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 12.1732, p-val = 0.0068
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb    ci.ub    
# MATmean_SubMAT8Dao15   -0.0908  0.0418  -2.1696  0.0300  -0.1728  -0.0088  * 
# MATmean_SubMATDy15      0.0374  0.0477   0.7829  0.4337  -0.0562   0.1309    
# MATmean_SubMATXy8       0.0414  0.0339   1.2223  0.2216  -0.0250   0.1077    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_DNRA_filteredMATmean_Sub)
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
#               "a"                  "b"                  "b" 




### 8.4 LegumeNonlegume
DNRA_filteredLegumeNonlegume <- subset(DNRA, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
DNRA_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(DNRA_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- DNRA_filteredLegumeNonlegume %>%
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
# 3 Non-legume to Non-legume           87             25

overall_model_DNRA_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = DNRA_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -168.2508   336.5015   344.5015   358.9100   344.6519   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0327  0.1808     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1866.6508, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 16.6658, p-val = 0.0008
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume       -0.0301  0.0308  -0.9781  0.3280  -0.0905  0.0302    
# LegumeNonlegumeNon-legume to Legume        0.0295  0.0290   1.0165  0.3094  -0.0274  0.0864    
# LegumeNonlegumeNon-legume to Non-legume    0.0222  0.0299   0.7445  0.4566  -0.0363  0.0808  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_DNRA_filteredLegumeNonlegume)
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
DNRA_filteredAMnonAM <- subset(DNRA, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
DNRA_filteredAMnonAM$AMnonAM <- droplevels(factor(DNRA_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- DNRA_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             208             37
# 2 AM to nonAM           30              8
# 3 nonAM to AM           36              8

overall_model_DNRA_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = DNRA_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -171.5589   343.1177   351.1177   365.5262   351.2681   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0298  0.1725     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 2102.6737, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 9.5392, p-val = 0.0229
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM       0.0152  0.0267   0.5674  0.5705  -0.0372  0.0676    
# AMnonAMAM to nonAM   -0.0294  0.0295  -0.9982  0.3182  -0.0872  0.0283    
# AMnonAMnonAM to AM    0.0347  0.0389   0.8924  0.3722  -0.0415  0.1110   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_DNRA_filteredAMnonAM)
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
   #             "a"                "b"                "a" 


### 8.6 C3C4
DNRA_filteredC3C4 <- subset(DNRA, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
DNRA_filteredC3C4$C3C4 <- droplevels(factor(DNRA_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- DNRA_filteredC3C4 %>%
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
# 3 C4 to C3          131             14

overall_model_DNRA_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = DNRA_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 273; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -171.3863   342.7727   350.7727   365.1663   350.9236   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0308  0.1754     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 2006.7586, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 11.6913, p-val = 0.0085
# 
# Model Results:
# 
#               estimate      se    zval    pval    ci.lb   ci.ub    
# C3C4C3 to C3    0.0063  0.0299  0.2116  0.8324  -0.0523  0.0650    
# C3C4C3 to C4    0.0029  0.0277  0.1054  0.9160  -0.0514  0.0573    
# C3C4C4 to C3    0.0498  0.0292  1.7087  0.0875  -0.0073  0.1070  . 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredC3C4)
vcov_rotation <- vcov(overall_model_DNRA_filteredC3C4)
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
#         "ab"          "a"          "b" 


### 8.7 Annual_Pere
DNRA_filteredAnnual_Pere <- subset(DNRA, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
DNRA_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(DNRA_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- DNRA_filteredAnnual_Pere %>%
  group_by(Annual_Pere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Annual_Pere         Observations Unique_StudyID
#   <fct>                      <int>          <int>
# 1 Annual to Annual             240             40
# 2 Annual to Perennial           13              8
# 3 Perennial to Annual           21              7

overall_model_DNRA_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = DNRA_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -155.2813   310.5626   318.5626   332.9710   318.7129   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0348  0.1867     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 2063.9341, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 43.7160, p-val < .0001
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb    ci.ub     
# Annual_PereAnnual to Annual       0.0384  0.0283   1.3575  0.1746  -0.0170   0.0938     
# Annual_PereAnnual to Perennial   -0.1023  0.0329  -3.1083  0.0019  -0.1667  -0.0378  ** 
# Annual_PerePerennial to Annual   -0.0768  0.0398  -1.9291  0.0537  -0.1547   0.0012   . 

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_DNRA_filteredAnnual_Pere)
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
   #                 "a"                            "b"                            "b" 



### 8.8 Plant_stage
DNRA_filteredPlant_stage <- subset(DNRA, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
DNRA_filteredPlant_stage$Plant_stage <- droplevels(factor(DNRA_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- DNRA_filteredPlant_stage %>%
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
# 3 Vegetative stage             80              8

overall_model_DNRA_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = DNRA_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredPlant_stage)
# Multivariate Meta-Analysis Model (k = 239; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -134.7803   269.5607   277.5607   291.4160   277.7338   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0294  0.1714     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 236) = 1661.3944, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 51.7334, p-val < .0001
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb    ci.ub      
# Plant_stageHarvest               0.0690  0.0308   2.2414  0.0250   0.0087   0.1294    * 
# Plant_stageReproductive stage   -0.0938  0.0360  -2.6016  0.0093  -0.1644  -0.0231   ** 
# Plant_stageVegetative stage      0.1058  0.0317   3.3339  0.0009   0.0436   0.1680  *** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_DNRA_filteredPlant_stage)
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
# "      "a"                           "b"                           "c" 

### 8.9 Rotation_cycles2
DNRA_filteredRotation_cycles2 <- subset(DNRA, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
DNRA_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(DNRA_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- DNRA_filteredRotation_cycles2 %>%
  group_by(Rotation_cycles2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Rotation_cycles2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                         79             19
# 2 D1-3                       87             12
# 3 D10                        11              7
# 4 D3-5                       35              9
# 5 D5-10                      62             15

overall_model_DNRA_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = DNRA_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -173.3414   346.6828   358.6828   380.2511   359.0034   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0297  0.1725     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 269) = 2101.5067, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 10.5039, p-val = 0.0622
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1       0.0304  0.0345   0.8828  0.3774  -0.0371  0.0980    
# Rotation_cycles2D1-3     0.0306  0.0348   0.8805  0.3786  -0.0376  0.0989    
# Rotation_cycles2D10     -0.0860  0.0462  -1.8615  0.0627  -0.1766  0.0046  . 
# Rotation_cycles2D3-5     0.0056  0.0417   0.1337  0.8937  -0.0762  0.0874    
# Rotation_cycles2D5-10    0.0106  0.0409   0.2605  0.7945  -0.0695  0.0908 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_DNRA_filteredRotation_cycles2)
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
#                  "a"                   "a"                   "b"                   "a"                   "a"               "a" 


### 8.10 Duration2
DNRA_filteredDuration2 <- subset(DNRA, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
DNRA_filteredDuration2$Duration2 <- droplevels(factor(DNRA_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- DNRA_filteredDuration2 %>%
  group_by(Duration2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Duration2 Observations Unique_StudyID
#   <fct>            <int>          <int>
# 1 D1-5               146             26
# 2 D11-20              58             11
# 3 D20-40              37              8
# 4 D6-10               33             11

overall_model_DNRA_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = DNRA_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredDuration2)
# ultivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -167.6596   335.3193   345.3193   363.3114   345.5466   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0308  0.1755     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 2069.5279, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 16.4046, p-val = 0.0025
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5     -0.0039  0.0364  -0.1079  0.9141  -0.0754  0.0675    
# Duration2D11-20   -0.0448  0.0430  -1.0424  0.2972  -0.1290  0.0394    
# Duration2D20-40    0.0240  0.0430   0.5576  0.5771  -0.0603  0.1082    
# Duration2D6-10     0.1021  0.0429   2.3811  0.0173   0.0181  0.1861  * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredDuration2)
vcov_rotation <- vcov(overall_model_DNRA_filteredDuration2)
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
#    "abc"             "a"             "b"             "c" 


### 8.11 SpeciesRichness2
DNRA_filteredSpeciesRichness2 <- subset(DNRA, SpeciesRichness2 %in% c("R2", "R3"))
#
DNRA_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(DNRA_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- DNRA_filteredSpeciesRichness2 %>%
  group_by(SpeciesRichness2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   SpeciesRichness2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 R2                        214             44
# 2 R3                         51              8

overall_model_DNRA_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = DNRA_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredSpeciesRichness2)
# ultivariate Meta-Analysis Model (k = 265; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -171.8607   343.7214   349.7214   360.4378   349.8140   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0314  0.1771     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 263) = 2122.1126, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 8.0810, p-val = 0.0176
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2    0.0138  0.0273   0.5054  0.6132  -0.0397  0.0673    
# SpeciesRichness2R3   -0.0270  0.0301  -0.8954  0.3706  -0.0860  0.0321    



# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_DNRA_filteredSpeciesRichness2)
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
DNRA_filteredBulk_Rhizosphere <- subset(DNRA, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
DNRA_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(DNRA_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- DNRA_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
 # Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non_Rhizosphere           209             36
# 2 Rhizosphere                65             24

overall_model_DNRA_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = DNRA_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -169.6372   339.2745   345.2745   356.0919   345.3640   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0302  0.1739     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 2090.2630, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 11.8937, p-val = 0.0026
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere    0.0265  0.0264   1.0057  0.3146  -0.0252  0.0782    
# Bulk_RhizosphereRhizosphere       -0.0158  0.0273  -0.5769  0.5640  -0.0693  0.0378   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_DNRA_filteredBulk_Rhizosphere)
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
DNRA_filteredSoil_texture <- subset(DNRA, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
DNRA_filteredSoil_texture$Soil_texture <- droplevels(factor(DNRA_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- DNRA_filteredSoil_texture %>%
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
# 3 Medium                174             17

overall_model_DNRA_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = DNRA_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredSoil_texture)
# ltivariate Meta-Analysis Model (k = 224; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -172.4184   344.8368   352.8368   366.4295   353.0220   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0289  0.1700     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 221) = 1764.0442, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.1772, p-val = 0.5365
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse   -0.0002  0.0661  -0.0025  0.9980  -0.1298  0.1294    
# Soil_textureFine      0.0974  0.0680   1.4320  0.1521  -0.0359  0.2307    
# Soil_textureMedium    0.0152  0.0428   0.3557  0.7221  -0.0686  0.0991   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_DNRA_filteredSoil_texture)
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
DNRA_filteredTillage <- subset(DNRA, Tillage %in% c("Tillage", "No_tillage"))
#
DNRA_filteredTillage$Tillage <- droplevels(factor(DNRA_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- DNRA_filteredTillage %>%
  group_by(Tillage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Tillage    Observations Unique_StudyID
#   <fct>             <int>          <int>
# 1 No_tillage          119              7
# 2 Tillage              31              9

overall_model_DNRA_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = DNRA_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredTillage)
# Multivariate Meta-Analysis Model (k = 150; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -159.0706   318.1412   324.1412   333.1328   324.3079   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0067  0.0818     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 148) = 1052.0272, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 5.3600, p-val = 0.0686
# 
# Model Results:
# 
#                    estimate      se    zval    pval    ci.lb   ci.ub    
# TillageNo_tillage    0.0249  0.0261  0.9516  0.3413  -0.0263  0.0761    
# TillageTillage       0.0522  0.0259  2.0164  0.0438   0.0015  0.1029  * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredTillage)
vcov_rotation <- vcov(overall_model_DNRA_filteredTillage)
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
DNRA_filteredStraw_retention <- subset(DNRA, Straw_retention %in% c("Retention", "No_retention"))
#
DNRA_filteredStraw_retention$Straw_retention <- droplevels(factor(DNRA_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- DNRA_filteredStraw_retention %>%
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

overall_model_DNRA_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = DNRA_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredStraw_retention)
# ultivariate Meta-Analysis Model (k = 55; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -77.6256  155.2513  161.2513  167.1622  161.7411   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0254  0.1595     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 863.2963, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.8751, p-val = 0.3916
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention    0.0259  0.0544   0.4760  0.6340  -0.0807  0.1325    
# Straw_retentionRetention      -0.0211  0.0514  -0.4114  0.6808  -0.1218  0.0795      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_DNRA_filteredStraw_retention)
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
DNRA_filteredPrimer <- subset(DNRA, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
DNRA_filteredPrimer$Primer <- droplevels(factor(DNRA_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- DNRA_filteredPrimer %>%
  group_by(Primer) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Primer Observations Unique_StudyID
#   <fct>         <int>          <int>
# 1 V3-V4           126             28
# 2 V4               99              9
# 3 V4-V5            39              9

overall_model_DNRA_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = DNRA_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_DNRA_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 264; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -174.1995   348.3990   356.3990   370.6571   356.5553   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0304  0.1744     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 261) = 1963.9274, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 7.0429, p-val = 0.0705
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4   -0.0195  0.0343  -0.5686  0.5696  -0.0867  0.0477    
# PrimerV4       0.1694  0.0673   2.5151  0.0119   0.0374  0.3014  * 
# PrimerV4-V5   -0.0384  0.0611  -0.6276  0.5303  -0.1582  0.0814     
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_DNRA_filteredPrimer)
vcov_rotation <- vcov(overall_model_DNRA_filteredPrimer)
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
DNRA$Wr <- 1 / DNRA$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = DNRA)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = DNRA)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = DNRA)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = DNRA)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = DNRA)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = DNRA)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = DNRA)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = DNRA)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = DNRA)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = DNRA)

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
#           Model       AIC        BIC   logLik
# Model1   Model1 -56.94920 -35.270434 34.47460
# Model2   Model2 -56.10785 -34.429078 34.05392
# Model3   Model3 -55.56951 -33.890744 33.78476
# Model4   Model4 -57.06148 -35.382716 34.53074
# Model5   Model5 -57.02483 -35.346057 34.51241
# Model6   Model6 -57.11496 -35.436192 34.55748$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Model7   Model7 -56.06957 -34.390802 34.03479
# Model8   Model8 -55.59588 -33.917111 33.79794
# Model9   Model9 -34.41738   1.713901 27.20869
# Model10 Model10 -41.09832  -4.967038 30.54916

##### Model 6 is the best model
summary(Model6)
# Number of obs: 274, groups:  StudyID, 50
anova(Model6)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)      8.1304  8.1304     1 131.188  1.8105 0.1808
# scale(log(SpeciesRichness)) 6.7366  6.7366     1 268.846  1.5001 0.2217
# scale(log(Duration))        9.5176  9.5176     1  98.508  2.1194 0.1486


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(pHCK) + (1 | StudyID), weights = Wr, data = DNRA)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles)       6.1705  6.1705     1 38.037  0.9210 0.34326  
# scale(log(SpeciesRichness)) 15.7827 15.7827     1 60.524  2.3558 0.13003  
# scale(log(Duration))         2.9280  2.9280     1 25.716  0.4370 0.51443  
# scale(pHCK)                 27.4008 27.4008     1 29.183  4.0899 0.05239 .

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = DNRA)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# scale(Rotation_cycles)       6.986   6.986     1 49.664  1.3092 0.258029   
# scale(log(SpeciesRichness)) 52.211  52.211     1 66.714  9.7842 0.002609 **
# scale(log(Duration))         4.497   4.497     1 29.712  0.8427 0.366004   
# scale(SOCCK)                 0.752   0.752     1 19.190  0.1409 0.711553   

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(TNCK) + (1 | StudyID), weights = Wr, data = DNRA)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles)      9.5282  9.5282     1 46.355  3.3385 0.07412 .
# scale(log(SpeciesRichness)) 1.2487  1.2487     1 32.749  0.4375 0.51295  
# scale(log(Duration))        3.6670  3.6670     1 27.659  1.2848 0.26673  
# scale(TNCK)                 0.8363  0.8363     1 13.246  0.2930 0.59729  


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = DNRA)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)      2.8454  2.8454     1 30.516  0.9829 0.3293
# scale(log(SpeciesRichness)) 3.3849  3.3849     1 35.012  1.1693 0.2869
# scale(log(Duration))        0.6146  0.6146     1 13.201  0.2123 0.6525
# scale(NO3CK)                1.6480  1.6480     1 36.976  0.5693 0.4553

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = DNRA)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)      3.9615  3.9615     1 30.515  1.4210 0.2424
# scale(log(SpeciesRichness)) 3.4144  3.4144     1 36.513  1.2248 0.2757
# scale(log(Duration))        1.9408  1.9408     1 12.782  0.6962 0.4194
# scale(NH4CK)                3.6595  3.6595     1 10.558  1.3127 0.2772

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(APCK) + (1 | StudyID), weights = Wr, data = DNRA)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# > anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)      11.3126 11.3126     1 58.769  2.1725 0.1458
# scale(log(SpeciesRichness))  1.5683  1.5683     1 42.337  0.3012 0.5860
# scale(log(Duration))         5.1544  5.1544     1 35.586  0.9899 0.3265
# scale(APCK)                 14.4583 14.4583     1 57.504  2.7766 0.1011

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(AKCK) + (1 | StudyID), weights = Wr, data = DNRA)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)      0.6402  0.6402     1  9.1446  0.0768 0.7878
# scale(log(SpeciesRichness)) 1.3656  1.3656     1 26.0288  0.1639 0.6889
# scale(log(Duration))        0.8668  0.8668     1  9.3316  0.1040 0.7541
# scale(AKCK)                 7.9747  7.9747     1 24.6979  0.9572 0.3374

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(ANCK) + (1 | StudyID), weights = Wr, data = DNRA)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(Rotation_cycles)      43.111  43.111     1  7.2842  3.8527 0.088832 . 
# scale(log(SpeciesRichness)) 99.683  99.683     1 20.4325  8.9084 0.007208 **
# scale(log(Duration))        42.701  42.701     1  6.1699  3.8161 0.097256 . 
# scale(ANCK)                 25.744  25.744     1  4.2790  2.3007 0.199347   

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(Latitude) + (1 | StudyID), weights = Wr, data = DNRA)
summary(ModelLatitude)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles)       5.0946  5.0946     1 129.030  1.1387 0.28791  
# scale(log(SpeciesRichness))  7.5043  7.5043     1 266.843  1.6773 0.19640  
# scale(log(Duration))         7.3158  7.3158     1  97.835  1.6352 0.20401  
# scale(Latitude)             19.4280 19.4280     1  84.124  4.3425 0.04021 *

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(Longitude) + (1 | StudyID), weights = Wr, data = DNRA)
summary(ModelLongitude)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)      8.4328  8.4328     1 131.243  1.8846 0.1722
# scale(log(SpeciesRichness)) 7.9707  7.9707     1 268.422  1.7813 0.1831
# scale(log(Duration))        4.9843  4.9843     1 100.654  1.1139 0.2938
# scale(Longitude)            5.0337  5.0337     1  29.968  1.1250 0.2973

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = DNRA)
summary(ModelMAPmean)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)      7.8861  7.8861     1 134.216  1.7559 0.1874
# scale(log(SpeciesRichness)) 6.6316  6.6316     1 268.143  1.4766 0.2254
# scale(log(Duration))        9.2861  9.2861     1  97.921  2.0676 0.1536
# scale(MAPmean)              0.0548  0.0548     1  61.265  0.0122 0.9124

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + scale(MATmean) + (1 | StudyID), weights = Wr, data = DNRA)
summary(ModelMATmean)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles)      6.3035  6.3035     1 142.806  1.4083 0.2373
# scale(log(SpeciesRichness)) 6.3504  6.3504     1 268.094  1.4188 0.2347
# scale(log(Duration))        7.7410  7.7410     1 103.428  1.7294 0.1914
# scale(MATmean)              1.8252  1.8252     1  40.594  0.4078 0.5267


############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(DNRA$SpeciesRichness)) ## n = 274
p1 <- ggplot(DNRA, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnDNRA274")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(DNRA$Duration)) ## n = 274
p2 <- ggplot(DNRA, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnDNRA274")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(DNRA$Rotation_cycles)) ## n = 274
p3 <- ggplot(DNRA, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnDNRA274")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(DNRA$Latitude)) ## n = 274
p5 <- ggplot(DNRA, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnDNRA274")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(DNRA$Longitude)) ## n = 274
p6 <- ggplot(DNRA, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnDNRA274")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(DNRA$MAPmean)) ## n = 274
p7 <- ggplot(DNRA, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnDNRA274")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(DNRA$MATmean)) ## n = 274
p8 <- ggplot(DNRA, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnDNRA274")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(DNRA$pHCK)) ## n = 86
p9 <- ggplot(DNRA, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnDNRA86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(DNRA$SOCCK)) ## n = 83
p10 <- ggplot(DNRA, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnDNRA83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(DNRA$TNCK)) ## n = 52
p11 <- ggplot(DNRA, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnDNRA52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(DNRA$NO3CK)) ## n = 46
p12 <- ggplot(DNRA, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnDNRA46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(DNRA$NH4CK)) ## n = 45
p13<- ggplot(DNRA, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnDNRA45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(DNRA$APCK)) ## n = 64
p14 <- ggplot(DNRA, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnDNRA64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(DNRA$AKCK)) ## n = 39
p15 <- ggplot(DNRA, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnDNRA39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(DNRA$ANCK)) ## n = 31
p16 <- ggplot(DNRA, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnDNRA31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(DNRA$RRpH)) ## n = 86
p17 <- ggplot(DNRA, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="RRpH")+
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
sum(!is.na(DNRA$RRSOC)) ## n = 83
p18 <- ggplot(DNRA, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="RRSOC")+
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
sum(!is.na(DNRA$RRTN)) ## n = 52
p19 <- ggplot(DNRA, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="RRTN")+
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
sum(!is.na(DNRA$RRNO3)) ## n = 44
p20 <- ggplot(DNRA, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="RRNO3")+
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
sum(!is.na(DNRA$RRNH4)) ## n = 45
p21 <- ggplot(DNRA, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="RRNH4")+
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
sum(!is.na(DNRA$RRAP)) ## n = 64
p22 <- ggplot(DNRA, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="RRAP")+
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
sum(!is.na(DNRA$RRAK)) ## n = 39
p23 <- ggplot(DNRA, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="RRAK")+
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
sum(!is.na(DNRA$RRAN)) ## n = 31
p24 <- ggplot(DNRA, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnDNRA", x="RRAN")+
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
sum(!is.na(DNRA$RRYield)) ## n = 49
p25 <- ggplot(DNRA, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-0.62, 0.85), expand=c(0, 0)) + 
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






