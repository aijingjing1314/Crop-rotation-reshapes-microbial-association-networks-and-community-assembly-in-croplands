
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Methane_metabolism <- read.csv("Methane_metabolism.csv", fileEncoding = "latin1")
# Check data
head(Methane_metabolism)

# 1. The number of Obversation
total_number <- nrow(Methane_metabolism)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 274  

# 2. The number of Study
unique_studyid_number <- length(unique(Methane_metabolism$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 







#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Methane_metabolism_filteredLongitude_Sub <- subset(Methane_metabolism, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Methane_metabolism_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Methane_metabolism_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredLongitude_Sub %>%
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

overall_model_Methane_metabolism_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Methane_metabolism_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -104.5762   209.1523   215.1523   225.9697   215.2419   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0282  0.1678     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1304.3065, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.3397, p-val = 0.5118
# 
# Model Results:
# 
#                            estimate      se    zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0    0.0316  0.0288  1.0956  0.2733  -0.0249  0.0880    
# Longitude_SubLongitudeXy0    0.0213  0.0570  0.3734  0.7088  -0.0904  0.1330        

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredLongitude_Sub)
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
Methane_metabolism_filteredMAPmean2_Sub <- subset(Methane_metabolism, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Methane_metabolism_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Methane_metabolism_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredMAPmean2_Sub %>%
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

overall_model_Methane_metabolism_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Methane_metabolism_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -90.6824  181.3648  189.3648  203.7732  189.5152   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0470  0.2168     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1298.4027, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 33.4276, p-val < .0001
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub      
# MAPmean2_SubMAP600Dao1200    0.1491  0.0411   3.6271  0.0003   0.0685  0.2296  *** 
# MAPmean2_SubMAPDy1200        0.0148  0.0749   0.1981  0.8430  -0.1320  0.1616      
# MAPmean2_SubMAPXy600        -0.0625  0.0399  -1.5656  0.1174  -0.1408  0.0157     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredMAPmean2_Sub)
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
Methane_metabolism_filteredMATmean_Sub <- subset(Methane_metabolism, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Methane_metabolism_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Methane_metabolism_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredMATmean_Sub %>%
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

overall_model_Methane_metabolism_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Methane_metabolism_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -91.9106  183.8212  191.8212  206.2297  191.9716   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0498  0.2232     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1302.8020, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 32.2196, p-val < .0001
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb   ci.ub     
# MATmean_SubMAT8Dao15    0.1577  0.0481   3.2786  0.0010   0.0634  0.2520  ** 
# MATmean_SubMATDy15      0.0812  0.0601   1.3508  0.1768  -0.0366  0.1989     
# MATmean_SubMATXy8      -0.0558  0.0418  -1.3347  0.1820  -0.1377  0.0261  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredMATmean_Sub)
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
#               "a"                 "ab"                  "b" 




### 8.4 LegumeNonlegume
Methane_metabolism_filteredLegumeNonlegume <- subset(Methane_metabolism, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Methane_metabolism_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Methane_metabolism_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredLegumeNonlegume %>%
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

overall_model_Methane_metabolism_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Methane_metabolism_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -103.1830   206.3660   214.3660   228.7744   214.5163   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0258  0.1606     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1275.5730, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 8.3404, p-val = 0.0395
# 
# Model Results:
# 
#                                          estimate      se    zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume        0.0805  0.0328  2.4552  0.0141   0.0162  0.1448  * 
# LegumeNonlegumeNon-legume to Legume        0.0183  0.0280  0.6524  0.5141  -0.0366  0.0732    
# LegumeNonlegumeNon-legume to Non-legume    0.0146  0.0292  0.5006  0.6167  -0.0426  0.0718     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredLegumeNonlegume)
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
    #          "a"                                     "b"                                    "ab" 




### 8.5 AMnonAM
Methane_metabolism_filteredAMnonAM <- subset(Methane_metabolism, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Methane_metabolism_filteredAMnonAM$AMnonAM <- droplevels(factor(Methane_metabolism_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredAMnonAM %>%
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

overall_model_Methane_metabolism_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Methane_metabolism_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredAMnonAM)
# ltivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -104.3752   208.7504   216.7504   231.1588   216.9007   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0286  0.1692     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1302.2209, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.1745, p-val = 0.1595
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM       0.0434  0.0274   1.5860  0.1127  -0.0102  0.0971    
# AMnonAMAM to nonAM   -0.0109  0.0341  -0.3185  0.7501  -0.0778  0.0560    
# AMnonAMnonAM to AM   -0.0042  0.0502  -0.0826  0.9341  -0.1026  0.0943    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredAMnonAM)
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
Methane_metabolism_filteredC3C4 <- subset(Methane_metabolism, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Methane_metabolism_filteredC3C4$C3C4 <- droplevels(factor(Methane_metabolism_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredC3C4 %>%
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

overall_model_Methane_metabolism_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Methane_metabolism_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 273; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -89.3813  178.7627  186.7627  201.1564  186.9136   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0372  0.1930     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 1281.5965, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 37.1731, p-val < .0001
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb   ci.ub      
# C3C4C3 to C3   -0.0519  0.0336  -1.5464  0.1220  -0.1177  0.0139      
# C3C4C3 to C4    0.0572  0.0316   1.8108  0.0702  -0.0047  0.1191    . 
# C3C4C4 to C3    0.1158  0.0340   3.4080  0.0007   0.0492  0.1824  *** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredC3C4)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredC3C4)
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
#    "a"          "b"          "c" 


### 8.7 Annual_Pere
Methane_metabolism_filteredAnnual_Pere <- subset(Methane_metabolism, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
Methane_metabolism_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Methane_metabolism_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredAnnual_Pere %>%
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

overall_model_Methane_metabolism_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Methane_metabolism_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -106.0125   212.0249   220.0249   234.4334   220.1753   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0272  0.1651     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1297.6966, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.8389, p-val = 0.4171
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual       0.0324  0.0263   1.2336  0.2173  -0.0191  0.0839    
# Annual_PereAnnual to Perennial    0.0419  0.0380   1.1033  0.2699  -0.0325  0.1164    
# Annual_PerePerennial to Annual   -0.0140  0.0525  -0.2667  0.7897  -0.1168  0.0888    

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredAnnual_Pere)
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
   #           "a"                            "a"                            "a" 



### 8.8 Plant_stage
Methane_metabolism_filteredPlant_stage <- subset(Methane_metabolism, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
Methane_metabolism_filteredPlant_stage$Plant_stage <- droplevels(factor(Methane_metabolism_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredPlant_stage %>%
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

overall_model_Methane_metabolism_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Methane_metabolism_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredPlant_stage)
# Multivariate Meta-Analysis Model (k = 239; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -117.5995   235.1990   243.1990   257.0544   243.3722   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0314  0.1773     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 236) = 1116.1677, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.2085, p-val = 0.5303
# 
# Model Results:
# 
#                                estimate      se    zval    pval    ci.lb   ci.ub    
# Plant_stageHarvest               0.0174  0.0329  0.5278  0.5976  -0.0471  0.0819    
# Plant_stageReproductive stage    0.0509  0.0397  1.2827  0.1996  -0.0269  0.1288    
# Plant_stageVegetative stage      0.0068  0.0357  0.1917  0.8480  -0.0631  0.0768 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredPlant_stage)
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
# "      "a"                           "a"                           "a" 

### 8.9 Rotation_cycles2
Methane_metabolism_filteredRotation_cycles2 <- subset(Methane_metabolism, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
Methane_metabolism_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Methane_metabolism_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredRotation_cycles2 %>%
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

overall_model_Methane_metabolism_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Methane_metabolism_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -92.9130  185.8259  197.8259  219.3942  198.1466   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0311  0.1764     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 269) = 1227.0964, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 32.2392, p-val < .0001
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1       0.0527  0.0373   1.4120  0.1579  -0.0205  0.1259    
# Rotation_cycles2D1-3     0.0228  0.0384   0.5938  0.5526  -0.0525  0.0982    
# Rotation_cycles2D10      0.0261  0.0504   0.5173  0.6049  -0.0727  0.1249    
# Rotation_cycles2D3-5     0.0870  0.0434   2.0061  0.0449   0.0020  0.1720  * 
# Rotation_cycles2D5-10   -0.0304  0.0421  -0.7207  0.4711  -0.1129  0.0522    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredRotation_cycles2)
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
#           "ab"                  "ab"                  "ab"                   "a"                   "b"         


### 8.10 Duration2
Methane_metabolism_filteredDuration2 <- subset(Methane_metabolism, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Methane_metabolism_filteredDuration2$Duration2 <- droplevels(factor(Methane_metabolism_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredDuration2 %>%
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

overall_model_Methane_metabolism_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Methane_metabolism_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -92.2728  184.5455  194.5455  212.5376  194.7728   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0279  0.1670     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 1246.9745, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 29.6399, p-val < .0001
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub     
# Duration2D1-5     -0.0006  0.0358  -0.0179  0.9857  -0.0708  0.0695     
# Duration2D11-20    0.1175  0.0424   2.7730  0.0056   0.0345  0.2006  ** 
# Duration2D20-40   -0.0523  0.0435  -1.2028  0.2291  -0.1376  0.0329     
# Duration2D6-10     0.0935  0.0439   2.1305  0.0331   0.0075  0.1796   * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredDuration2)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredDuration2)
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
#    "ab"             "c"             "a"            "bc" 


### 8.11 SpeciesRichness2
Methane_metabolism_filteredSpeciesRichness2 <- subset(Methane_metabolism, SpeciesRichness2 %in% c("R2", "R3"))
#
Methane_metabolism_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Methane_metabolism_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredSpeciesRichness2 %>%
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

overall_model_Methane_metabolism_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Methane_metabolism_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredSpeciesRichness2)
# ultivariate Meta-Analysis Model (k = 265; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -109.2329   218.4658   224.4658   235.1823   224.5585   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0282  0.1680     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 263) = 1287.5248, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.7343, p-val = 0.2548
# 
# Model Results:
# 
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2    0.0218  0.0266  0.8211  0.4116  -0.0303  0.0740    
# SpeciesRichness2R3    0.0562  0.0350  1.6066  0.1081  -0.0124  0.1248      



# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredSpeciesRichness2)
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
Methane_metabolism_filteredBulk_Rhizosphere <- subset(Methane_metabolism, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Methane_metabolism_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Methane_metabolism_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredBulk_Rhizosphere %>%
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

overall_model_Methane_metabolism_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Methane_metabolism_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -105.8231   211.6461   217.6461   228.4635   217.7357   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0274  0.1656     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1306.0227, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.6961, p-val = 0.4283
# 
# Model Results:
# 
#                                  estimate      se    zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere    0.0333  0.0262  1.2680  0.2048  -0.0182  0.0847    
# Bulk_RhizosphereRhizosphere        0.0225  0.0279  0.8075  0.4194  -0.0321  0.0771  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredBulk_Rhizosphere)
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
Methane_metabolism_filteredSoil_texture <- subset(Methane_metabolism, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Methane_metabolism_filteredSoil_texture$Soil_texture <- droplevels(factor(Methane_metabolism_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredSoil_texture %>%
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

overall_model_Methane_metabolism_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Methane_metabolism_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 224; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -122.4785   244.9571   252.9571   266.5497   253.1423   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0279  0.1670     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 221) = 1043.6156, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.1927, p-val = 0.3629
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse    0.0913  0.0673   1.3570  0.1748  -0.0406  0.2232    
# Soil_textureFine     -0.0095  0.0667  -0.1417  0.8873  -0.1402  0.1213    
# Soil_textureMedium    0.0493  0.0427   1.1538  0.2486  -0.0344  0.1330    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredSoil_texture)
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
#            "a"                "a"                "a" 


### 8.14 Tillage
Methane_metabolism_filteredTillage <- subset(Methane_metabolism, Tillage %in% c("Tillage", "No_tillage"))
#
Methane_metabolism_filteredTillage$Tillage <- droplevels(factor(Methane_metabolism_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredTillage %>%
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

overall_model_Methane_metabolism_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Methane_metabolism_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredTillage)
# ltivariate Meta-Analysis Model (k = 150; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -118.7870   237.5739   243.5739   252.5655   243.7406   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0121  0.1100     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 148) = 715.6579, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 8.7841, p-val = 0.0124
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0353  0.0352  -1.0025  0.3161  -0.1042  0.0337    
# TillageTillage       0.0345  0.0344   1.0024  0.3162  -0.0329  0.1018    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredTillage)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredTillage)
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
Methane_metabolism_filteredStraw_retention <- subset(Methane_metabolism, Straw_retention %in% c("Retention", "No_retention"))
#
Methane_metabolism_filteredStraw_retention$Straw_retention <- droplevels(factor(Methane_metabolism_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredStraw_retention %>%
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

overall_model_Methane_metabolism_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Methane_metabolism_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 55; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -13.0385   26.0770   32.0770   37.9879   32.5668   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0439  0.2095     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 372.6301, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.0265, p-val = 0.3630
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention   -0.0031  0.0708  -0.0432  0.9655  -0.1418  0.1356    
# Straw_retentionRetention       0.0497  0.0672   0.7395  0.4596  -0.0820  0.1814    
  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredStraw_retention)
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
Methane_metabolism_filteredPrimer <- subset(Methane_metabolism, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
Methane_metabolism_filteredPrimer$Primer <- droplevels(factor(Methane_metabolism_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Methane_metabolism_filteredPrimer %>%
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

overall_model_Methane_metabolism_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Methane_metabolism_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Methane_metabolism_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 264; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -107.8060   215.6121   223.6121   237.8701   223.7683   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0284  0.1685     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 261) = 1184.0684, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.8078, p-val = 0.1864
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4    0.0481  0.0343   1.4030  0.1606  -0.0191  0.1153    
# PrimerV4       0.0756  0.0640   1.1809  0.2377  -0.0499  0.2011    
# PrimerV4-V5   -0.0728  0.0606  -1.2020  0.2294  -0.1916  0.0459  
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Methane_metabolism_filteredPrimer)
vcov_rotation <- vcov(overall_model_Methane_metabolism_filteredPrimer)
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
Methane_metabolism$Wr <- 1 / Methane_metabolism$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Methane_metabolism)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Methane_metabolism)

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
#           Model       AIC      BIC      logLik
# Model1   Model1  9.338806 31.01757  1.33059695
# Model2   Model2 12.125310 33.80408 -0.06265501
# Model3   Model3  8.980908 30.65968  1.50954599
# Model4   Model4  9.317230 30.99600  1.34138498
# Model5   Model5 12.395136 34.07390 -0.19756791
# Model6   Model6 12.422644 34.10141 -0.21132204
# Model7   Model7 12.095076 33.77384 -0.04753792
# Model8   Model8  8.961836 30.64061  1.51908181#########################################
# Model9   Model9 31.651517 67.78280 -5.82575847
# Model10 Model10 30.004853 66.13613 -5.00242659

##### Model 8 is the best model
summary(Model8)
# Number of obs: 274, groups:  StudyID, 50
anova(Model8)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  1.7861  1.7861     1  81.199  0.5631 0.45518  
# scale(log(SpeciesRichness))  5.1671  5.1671     1 198.327  1.6291 0.20332  
# scale(Duration)             12.6373 12.6373     1 149.469  3.9842 0.04774 *


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.1504  0.1504     1 37.311  0.0537 0.8179
# scale(log(SpeciesRichness)) 4.0672  4.0672     1 58.622  1.4539 0.2327
# scale(Duration)             0.0032  0.0032     1 64.370  0.0012 0.9730
# scale(pHCK)                 0.3164  0.3164     1 32.119  0.1131 0.7388.

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.01721 0.01721     1 32.065  0.0060 0.9389
# scale(log(SpeciesRichness)) 1.32878 1.32878     1 75.536  0.4615 0.4990
# scale(Duration)             0.64744 0.64744     1 63.959  0.2248 0.6370
# scale(SOCCK)                0.05679 0.05679     1 18.765  0.0197 0.8898

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.22088 0.22088     1 19.985  0.1414 0.7108
# scale(log(SpeciesRichness)) 1.79409 1.79409     1 45.770  1.1489 0.2894
# scale(Duration)             1.09411 1.09411     1 38.815  0.7007 0.4077
# scale(TNCK)                 1.05519 1.05519     1 12.067  0.6757 0.4270


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.28771 0.28771     1 18.918  0.2746 0.6064
# scale(log(SpeciesRichness)) 2.91744 2.91744     1 28.240  2.7844 0.1062
# scale(Duration)             1.29724 1.29724     1 38.620  1.2381 0.2727
# scale(NO3CK)                1.06324 1.06324     1  9.101  1.0148 0.3398

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.04142 0.04142     1 12.295  0.0360 0.8525
# scale(log(SpeciesRichness)) 2.57939 2.57939     1 29.789  2.2445 0.1446
# scale(Duration)             0.32439 0.32439     1 34.559  0.2823 0.5986
# scale(NH4CK)                0.01215 0.01215     1  6.770  0.0106 0.9211

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.12251 0.12251     1 32.865  0.0612 0.8062
# scale(log(SpeciesRichness)) 3.13905 3.13905     1 41.921  1.5675 0.2175
# scale(Duration)             0.07238 0.07238     1 58.513  0.0361 0.8499
# scale(APCK)                 0.02407 0.02407     1 52.272  0.0120 0.9131

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.95952 1.95952     1 13.033  0.9395 0.3500
# scale(log(SpeciesRichness)) 2.71893 2.71893     1 26.533  1.3036 0.2638
# scale(Duration)             0.70150 0.70150     1 13.754  0.3363 0.5713
# scale(AKCK)                 0.63518 0.63518     1 14.010  0.3045 0.5897
# 
# #### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 25.9240 25.9240     1 16.6746  6.2243 0.02342 *
# scale(log(SpeciesRichness))  1.3050  1.3050     1 25.8603  0.3133 0.58047  
# scale(Duration)             22.5172 22.5172     1 20.4006  5.4063 0.03048 *
# scale(ANCK)                  3.3426  3.3426     1  4.6018  0.8025 0.41473  

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
summary(ModelLatitude)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  1.6870  1.6870     1  80.422  0.5351 0.46658  
# scale(log(SpeciesRichness))  5.7276  5.7276     1 212.313  1.8168 0.17913  
# scale(Duration)             13.0366 13.0366     1 156.503  4.1354 0.04369 *
# scale(Latitude)              1.5026  1.5026     1  48.028  0.4766 0.49327  

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
summary(ModelLongitude)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  1.6504  1.6504     1  79.641  0.5215 0.47232  
# scale(log(SpeciesRichness))  4.7497  4.7497     1 205.300  1.5008 0.22196  
# scale(Duration)             12.2630 12.2630     1 157.146  3.8748 0.05078 .
# scale(Longitude)             0.1022  0.1022     1  24.902  0.0323 0.85885  

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
summary(ModelMAPmean)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelMAPmean) 
# ype III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  1.5623  1.5623     1  73.027  0.4937 0.48451  
# scale(log(SpeciesRichness))  5.0059  5.0059     1 190.536  1.5820 0.21001  
# scale(Duration)             11.7621 11.7621     1 121.503  3.7171 0.05619 .
# scale(MAPmean)               0.0690  0.0690     1  41.055  0.0218 0.88330  

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Methane_metabolism)
summary(ModelMATmean)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  2.1851  2.1851     1  70.527  0.6853 0.41054  
# scale(log(SpeciesRichness))  5.3299  5.3299     1 187.673  1.6717 0.19761  
# scale(Duration)             13.4926 13.4926     1 126.506  4.2320 0.04172 *
# scale(MATmean)               1.3739  1.3739     1  28.913  0.4309 0.51673  


############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Methane_metabolism$SpeciesRichness)) ## n = 274
p1 <- ggplot(Methane_metabolism, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnMethane_metabolism274")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Methane_metabolism$Duration)) ## n = 274
p2 <- ggplot(Methane_metabolism, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnMethane_metabolism274")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Methane_metabolism$Rotation_cycles)) ## n = 274
p3 <- ggplot(Methane_metabolism, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnMethane_metabolism274")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Methane_metabolism$Latitude)) ## n = 274
p5 <- ggplot(Methane_metabolism, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnMethane_metabolism274")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Methane_metabolism$Longitude)) ## n = 274
p6 <- ggplot(Methane_metabolism, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnMethane_metabolism274")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Methane_metabolism$MAPmean)) ## n = 274
p7 <- ggplot(Methane_metabolism, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnMethane_metabolism274")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Methane_metabolism$MATmean)) ## n = 274
p8 <- ggplot(Methane_metabolism, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnMethane_metabolism274")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Methane_metabolism$pHCK)) ## n = 86
p9 <- ggplot(Methane_metabolism, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnMethane_metabolism86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Methane_metabolism$SOCCK)) ## n = 83
p10 <- ggplot(Methane_metabolism, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnMethane_metabolism83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Methane_metabolism$TNCK)) ## n = 52
p11 <- ggplot(Methane_metabolism, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnMethane_metabolism52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Methane_metabolism$NO3CK)) ## n = 46
p12 <- ggplot(Methane_metabolism, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnMethane_metabolism46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Methane_metabolism$NH4CK)) ## n = 45
p13<- ggplot(Methane_metabolism, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnMethane_metabolism45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Methane_metabolism$APCK)) ## n = 64
p14 <- ggplot(Methane_metabolism, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnMethane_metabolism64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Methane_metabolism$AKCK)) ## n = 39
p15 <- ggplot(Methane_metabolism, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnMethane_metabolism39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Methane_metabolism$ANCK)) ## n = 31
p16 <- ggplot(Methane_metabolism, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnMethane_metabolism31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Methane_metabolism$RRpH)) ## n = 86
p17 <- ggplot(Methane_metabolism, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="RRpH")+
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
sum(!is.na(Methane_metabolism$RRSOC)) ## n = 83
p18 <- ggplot(Methane_metabolism, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="RRSOC")+
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
sum(!is.na(Methane_metabolism$RRTN)) ## n = 52
p19 <- ggplot(Methane_metabolism, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="RRTN")+
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
sum(!is.na(Methane_metabolism$RRNO3)) ## n = 44
p20 <- ggplot(Methane_metabolism, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="RRNO3")+
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
sum(!is.na(Methane_metabolism$RRNH4)) ## n = 45
p21 <- ggplot(Methane_metabolism, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="RRNH4")+
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
sum(!is.na(Methane_metabolism$RRAP)) ## n = 64
p22 <- ggplot(Methane_metabolism, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="RRAP")+
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
sum(!is.na(Methane_metabolism$RRAK)) ## n = 39
p23 <- ggplot(Methane_metabolism, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="RRAK")+
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
sum(!is.na(Methane_metabolism$RRAN)) ## n = 31
p24 <- ggplot(Methane_metabolism, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMethane_metabolism", x="RRAN")+
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
sum(!is.na(Methane_metabolism$RRYield)) ## n = 49
p25 <- ggplot(Methane_metabolism, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-0.32, 0.7), expand=c(0, 0)) + 
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





