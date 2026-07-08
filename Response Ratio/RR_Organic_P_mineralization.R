
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Organic_P_mineralization <- read.csv("Organic_P_mineralization.csv", fileEncoding = "latin1")
# Check data
head(Organic_P_mineralization)

# 1. The number of Obversation
total_number <- nrow(Organic_P_mineralization)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 275  

# 2. The number of Study
unique_studyid_number <- length(unique(Organic_P_mineralization$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 






#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Organic_P_mineralization_filteredLongitude_Sub <- subset(Organic_P_mineralization, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Organic_P_mineralization_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Organic_P_mineralization_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredLongitude_Sub %>%
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

overall_model_Organic_P_mineralization_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  260.5193  -521.0386  -515.0386  -504.2102  -514.9494   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0028  0.0532     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 1246.3828, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.2626, p-val = 0.8769
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0   -0.0040  0.0092  -0.4372  0.6620  -0.0220  0.0140    
# Longitude_SubLongitudeXy0    0.0049  0.0183   0.2674  0.7892  -0.0309  0.0407    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredLongitude_Sub)
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
Organic_P_mineralization_filteredMAPmean2_Sub <- subset(Organic_P_mineralization, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Organic_P_mineralization_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Organic_P_mineralization_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredMAPmean2_Sub %>%
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

overall_model_Organic_P_mineralization_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  260.7390  -521.4781  -513.4781  -499.0549  -513.3283   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0031  0.0557     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1244.7327, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.1289, p-val = 0.2479
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean2_SubMAP600Dao1200    0.0115  0.0112   1.0260  0.3049  -0.0105  0.0334    
# MAPmean2_SubMAPDy1200       -0.0118  0.0201  -0.5901  0.5551  -0.0512  0.0275    
# MAPmean2_SubMAPXy600        -0.0107  0.0109  -0.9863  0.3240  -0.0320  0.0106 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredMAPmean2_Sub)
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
Organic_P_mineralization_filteredMATmean_Sub <- subset(Organic_P_mineralization, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Organic_P_mineralization_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Organic_P_mineralization_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredMATmean_Sub %>%
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

overall_model_Organic_P_mineralization_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  259.6127  -519.2254  -511.2254  -496.8022  -511.0756   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0031  0.0552     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1232.9360, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.3151, p-val = 0.5096
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubMAT8Dao15    0.0042  0.0131   0.3193  0.7495  -0.0214  0.0298    
# MATmean_SubMATDy15      0.0073  0.0154   0.4724  0.6367  -0.0229  0.0374    
# MATmean_SubMATXy8      -0.0114  0.0108  -1.0506  0.2934  -0.0326  0.0098    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredMATmean_Sub)
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
Organic_P_mineralization_filteredLegumeNonlegume <- subset(Organic_P_mineralization, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Organic_P_mineralization_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Organic_P_mineralization_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredLegumeNonlegume %>%
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

overall_model_Organic_P_mineralization_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  257.6722  -515.3444  -507.3444  -492.9212  -507.1946   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0030  0.0545     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1175.8573, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.2100, p-val = 0.3604
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume       -0.0097  0.0094  -1.0353  0.3005  -0.0281  0.0087    
# LegumeNonlegumeNon-legume to Legume       -0.0026  0.0087  -0.2955  0.7676  -0.0196  0.0144    
# LegumeNonlegumeNon-legume to Non-legume    0.0018  0.0088   0.2021  0.8399  -0.0155  0.0191    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredLegumeNonlegume)
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
    #                   "a"                                     "a"                                     "a" 




### 8.5 AMnonAM
Organic_P_mineralization_filteredAMnonAM <- subset(Organic_P_mineralization, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Organic_P_mineralization_filteredAMnonAM$AMnonAM <- droplevels(factor(Organic_P_mineralization_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredAMnonAM %>%
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

overall_model_Organic_P_mineralization_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  260.7220  -521.4441  -513.4441  -499.0209  -513.2943   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0028  0.0531     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1253.8438, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 6.8059, p-val = 0.0783
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM       0.0057  0.0087   0.6507  0.5152  -0.0114  0.0228    
# AMnonAMAM to nonAM   -0.0184  0.0113  -1.6344  0.1022  -0.0405  0.0037    
# AMnonAMnonAM to AM   -0.0254  0.0135  -1.8842  0.0595  -0.0518  0.0010  .   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredAMnonAM)
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
   #         "a"                "b"                "b" 


### 8.6 C3C4
Organic_P_mineralization_filteredC3C4 <- subset(Organic_P_mineralization, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Organic_P_mineralization_filteredC3C4$C3C4 <- droplevels(factor(Organic_P_mineralization_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredC3C4 %>%
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

overall_model_Organic_P_mineralization_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  260.8591  -521.7183  -513.7183  -499.3098  -513.5679   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0027  0.0523     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1207.7119, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 10.0577, p-val = 0.0181
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb   ci.ub    
# C3C4C3 to C3   -0.0110  0.0099  -1.1177  0.2637  -0.0303  0.0083    
# C3C4C3 to C4   -0.0007  0.0087  -0.0759  0.9395  -0.0178  0.0164    
# C3C4C4 to C3    0.0118  0.0091   1.2974  0.1945  -0.0060  0.0296 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredC3C4)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredC3C4)
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
Organic_P_mineralization_filteredAnnual_Pere <- subset(Organic_P_mineralization, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
Organic_P_mineralization_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Organic_P_mineralization_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredAnnual_Pere %>%
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

overall_model_Organic_P_mineralization_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  262.9130  -525.8259  -517.8259  -503.4027  -517.6761   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0025  0.0502     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1096.3389, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 11.9746, p-val = 0.0075
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb    ci.ub     
# Annual_PereAnnual to Annual       0.0051  0.0083   0.6168  0.5373  -0.0111   0.0213     
# Annual_PereAnnual to Perennial   -0.0155  0.0135  -1.1448  0.2523  -0.0420   0.0110     
# Annual_PerePerennial to Annual   -0.0441  0.0150  -2.9353  0.0033  -0.0736  -0.0147  ** 

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredAnnual_Pere)
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
   #                         "a"                           "a"                            "b" 



### 8.8 Plant_stage
Organic_P_mineralization_filteredPlant_stage <- subset(Organic_P_mineralization, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
Organic_P_mineralization_filteredPlant_stage$Plant_stage <- droplevels(factor(Organic_P_mineralization_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredPlant_stage %>%
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

overall_model_Organic_P_mineralization_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredPlant_stage)
# Multivariate Meta-Analysis Model (k = 240; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  220.2116  -440.4232  -432.4232  -418.5510  -432.2508   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0020  0.0451     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 237) = 1030.1846, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.6107, p-val = 0.3067
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb   ci.ub    
# Plant_stageHarvest              -0.0016  0.0087  -0.1863  0.8522  -0.0187  0.0154    
# Plant_stageReproductive stage   -0.0063  0.0111  -0.5708  0.5681  -0.0281  0.0154    
# Plant_stageVegetative stage     -0.0078  0.0090  -0.8631  0.3881  -0.0254  0.0099   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredPlant_stage)
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
Organic_P_mineralization_filteredRotation_cycles2 <- subset(Organic_P_mineralization, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
Organic_P_mineralization_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Organic_P_mineralization_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredRotation_cycles2 %>%
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

overall_model_Organic_P_mineralization_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  256.7767  -513.5534  -501.5534  -479.9628  -501.2340   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0029  0.0538     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 1277.3880, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 8.4217, p-val = 0.1345
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1      -0.0013  0.0112  -0.1196  0.9048  -0.0232  0.0205    
# Rotation_cycles2D1-3    -0.0095  0.0112  -0.8420  0.3998  -0.0315  0.0126    
# Rotation_cycles2D10     -0.0025  0.0175  -0.1403  0.8884  -0.0369  0.0319    
# Rotation_cycles2D3-5     0.0083  0.0134   0.6160  0.5379  -0.0180  0.0345    
# Rotation_cycles2D5-10   -0.0018  0.0131  -0.1394  0.8892  -0.0276  0.0239 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredRotation_cycles2)
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
#           "a"                   "b"                  "ab"                  "ab"                  "ab" 


### 8.10 Duration2
Organic_P_mineralization_filteredDuration2 <- subset(Organic_P_mineralization, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Organic_P_mineralization_filteredDuration2$Duration2 <- droplevels(factor(Organic_P_mineralization_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredDuration2 %>%
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

overall_model_Organic_P_mineralization_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  257.6086  -515.2172  -505.2172  -487.2066  -504.9908   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0031  0.0552     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1272.1334, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 2.5433, p-val = 0.6369
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5     -0.0041  0.0119  -0.3493  0.7269  -0.0274  0.0191    
# Duration2D11-20    0.0063  0.0140   0.4508  0.6521  -0.0211  0.0337    
# Duration2D20-40   -0.0113  0.0147  -0.7698  0.4414  -0.0402  0.0175    
# Duration2D6-10     0.0003  0.0146   0.0186  0.9852  -0.0283  0.0289   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredDuration2)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredDuration2)
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
Organic_P_mineralization_filteredSpeciesRichness2 <- subset(Organic_P_mineralization, SpeciesRichness2 %in% c("R2", "R3"))
#
Organic_P_mineralization_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Organic_P_mineralization_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredSpeciesRichness2 %>%
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

overall_model_Organic_P_mineralization_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredSpeciesRichness2)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  251.3901  -502.7802  -496.7802  -486.0524  -496.6879   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0027  0.0518     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 264) = 1261.8404, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.0934, p-val = 0.5788
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2   -0.0037  0.0083  -0.4443  0.6568  -0.0198  0.0125    
# SpeciesRichness2R3   -0.0108  0.0109  -0.9943  0.3201  -0.0322  0.0105


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredSpeciesRichness2)
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
Organic_P_mineralization_filteredBulk_Rhizosphere <- subset(Organic_P_mineralization, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Organic_P_mineralization_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Organic_P_mineralization_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredBulk_Rhizosphere %>%
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

overall_model_Organic_P_mineralization_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  270.9092  -541.8184  -535.8184  -524.9900  -535.7292   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0030  0.0547     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 1227.2413, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 23.7560, p-val < .0001
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb    ci.ub    
# Bulk_RhizosphereNon_Rhizosphere    0.0083  0.0087   0.9615  0.3363  -0.0087   0.0254    
# Bulk_RhizosphereRhizosphere       -0.0228  0.0094  -2.4232  0.0154  -0.0412  -0.0044  *  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredBulk_Rhizosphere)
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
Organic_P_mineralization_filteredSoil_texture <- subset(Organic_P_mineralization, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Organic_P_mineralization_filteredSoil_texture$Soil_texture <- droplevels(factor(Organic_P_mineralization_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredSoil_texture %>%
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

overall_model_Organic_P_mineralization_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 225; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  195.0479  -390.0959  -382.0959  -368.4851  -381.9115   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0016  0.0397     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 222) = 946.6731, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.9623, p-val = 0.5803
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse   -0.0015  0.0157  -0.0971  0.9227  -0.0323  0.0292    
# Soil_textureFine      0.0142  0.0199   0.7129  0.4759  -0.0248  0.0532    
# Soil_textureMedium   -0.0125  0.0104  -1.2020  0.2294  -0.0328  0.0079 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredSoil_texture)
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
Organic_P_mineralization_filteredTillage <- subset(Organic_P_mineralization, Tillage %in% c("Tillage", "No_tillage"))
#
Organic_P_mineralization_filteredTillage$Tillage <- droplevels(factor(Organic_P_mineralization_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredTillage %>%
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

overall_model_Organic_P_mineralization_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredTillage)
# Multivariate Meta-Analysis Model (k = 151; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  155.7083  -311.4166  -305.4166  -296.4048  -305.2511   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0019  0.0436     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 149) = 457.1475, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.8955, p-val = 0.6391
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0047  0.0137  -0.3434  0.7313  -0.0315  0.0221    
# TillageTillage       0.0015  0.0135   0.1120  0.9108  -0.0250  0.0281    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredTillage)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredTillage)
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
Organic_P_mineralization_filteredStraw_retention <- subset(Organic_P_mineralization, Straw_retention %in% c("Retention", "No_retention"))
#
Organic_P_mineralization_filteredStraw_retention$Straw_retention <- droplevels(factor(Organic_P_mineralization_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredStraw_retention %>%
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

overall_model_Organic_P_mineralization_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 55; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  42.4402  -84.8804  -78.8804  -72.9695  -78.3906   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0027  0.0517     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 411.3277, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.0988, p-val = 0.5773
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention   -0.0018  0.0177  -0.1011  0.9195  -0.0365  0.0329    
# Straw_retentionRetention      -0.0123  0.0167  -0.7358  0.4619  -0.0450  0.0204   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredStraw_retention)
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
Organic_P_mineralization_filteredPrimer <- subset(Organic_P_mineralization, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
Organic_P_mineralization_filteredPrimer$Primer <- droplevels(factor(Organic_P_mineralization_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Organic_P_mineralization_filteredPrimer %>%
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

overall_model_Organic_P_mineralization_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Organic_P_mineralization_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Organic_P_mineralization_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 265; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  296.1397  -592.2795  -584.2795  -570.0061  -584.1239   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0035  0.0588     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 262) = 1123.8281, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.2362, p-val = 0.7443
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4    0.0020  0.0118   0.1717  0.8637  -0.0210  0.0250    
# PrimerV4      -0.0212  0.0228  -0.9312  0.3517  -0.0658  0.0234    
# PrimerV4-V5   -0.0125  0.0215  -0.5827  0.5601  -0.0547  0.0296     
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Organic_P_mineralization_filteredPrimer)
vcov_rotation <- vcov(overall_model_Organic_P_mineralization_filteredPrimer)
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
Organic_P_mineralization$Wr <- 1 / Organic_P_mineralization$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Organic_P_mineralization)

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
# Model1   Model1 -667.8213 -646.1207 339.9107
# Model2   Model2 -668.2151 -646.5145 340.1076
# Model3   Model3 -668.2919 -646.5912 340.1459###############################################################
# Model4   Model4 -667.7438 -646.0432 339.8719
# Model5   Model5 -667.8126 -646.1120 339.9063
# Model6   Model6 -667.6910 -645.9904 339.8455
# Model7   Model7 -668.2899 -646.5893 340.1449
# Model8   Model8 -668.2461 -646.5455 340.1231
# Model9   Model9 -628.9008 -592.7331 324.4504
# Model10 Model10 -630.1799 -594.0122 325.0900

##### Model 3 is the best model
summary(Model3)
# Number of obs: 275, groups:  StudyID, 50
anova(Model3)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.66368 1.66368     1 126.67  0.5594 0.4559
# scale(SpeciesRichness)      0.00399 0.00399     1 236.17  0.0013 0.9708
# scale(Duration)             1.10734 1.10734     1 154.50  0.3723 0.5426


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.55902 0.55902     1 26.706  0.1156 0.7365
# scale(SpeciesRichness)      0.84822 0.84822     1 67.416  0.1754 0.6767
# scale(Duration)             2.13643 2.13643     1 51.391  0.4418 0.5092
# scale(pHCK)                 0.80433 0.80433     1 34.483  0.1663 0.6859

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.33482 0.33482     1 27.085  0.0842 0.7739
# scale(SpeciesRichness)      0.55108 0.55108     1 70.329  0.1385 0.7109
# scale(Duration)             0.04114 0.04114     1 63.301  0.0103 0.9193
# scale(SOCCK)                0.34886 0.34886     1 13.724  0.0877 0.7716

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 3.2031  3.2031     1  8.556  0.8127 0.3920
# scale(SpeciesRichness)      1.4945  1.4945     1 42.818  0.3792 0.5413
# scale(Duration)             4.1836  4.1836     1 40.263  1.0614 0.3090
# scale(TNCK)                 5.8916  5.8916     1 26.598  1.4947 0.2322


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 3.9363  3.9363     1 17.764  4.8683 0.04077 *
# scale(SpeciesRichness)      0.5870  0.5870     1 28.970  0.7260 0.40116  
# scale(Duration)             5.4061  5.4061     1 38.519  6.6863 0.01362 *
# scale(NO3CK)                0.0985  0.0985     1 11.851  0.1218 0.73321  

# #### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 4.6532  4.6532     1 18.213  5.7537 0.02736 *
# scale(SpeciesRichness)      0.6005  0.6005     1 29.712  0.7425 0.39575  
# scale(Duration)             5.7310  5.7310     1 37.516  7.0864 0.01137 *
# scale(NH4CK)                0.3799  0.3799     1  9.369  0.4697 0.50972  

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  2.3196  2.3196     1 14.928  0.4631 0.50660  
# scale(SpeciesRichness)       0.5717  0.5717     1 51.576  0.1141 0.73686  
# scale(Duration)              0.7625  0.7625     1 35.626  0.1522 0.69874  
# scale(APCK)                 23.3240 23.3240     1 24.508  4.6563 0.04095 *

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.40651 0.40651     1 14.400  0.0639 0.8040
# scale(SpeciesRichness)      0.41534 0.41534     1 29.065  0.0653 0.8001
# scale(Duration)             2.04764 2.04764     1 18.180  0.3219 0.5774
# scale(AKCK)                 1.04418 1.04418     1 24.119  0.1642 0.6889

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  0.0180  0.0180     1 18.2507  0.0016 0.9682
# scale(SpeciesRichness)       0.8838  0.8838     1 25.4509  0.0802 0.7793
# scale(Duration)              0.2688  0.2688     1 19.3464  0.0244 0.8775
# scale(ANCK)                 20.6098 20.6098     1  4.0269  1.8706 0.2428

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
summary(ModelLatitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.80039 1.80039     1 125.014  0.6053 0.4380
# scale(SpeciesRichness)      0.02006 0.02006     1 238.213  0.0067 0.9346
# scale(Duration)             1.16339 1.16339     1 152.199  0.3911 0.5326
# scale(Latitude)             1.32505 1.32505     1  78.534  0.4455 0.5064

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
summary(ModelLongitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLongitude)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.75440 1.75440     1 126.970  0.5913 0.4434
# scale(SpeciesRichness)      0.00012 0.00012     1 237.204  0.0000 0.9949
# scale(Duration)             0.72993 0.72993     1 152.763  0.2460 0.6206
# scale(Longitude)            0.16610 0.16610     1  24.503  0.0560 0.8149

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
summary(ModelMAPmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 2.9319  2.9319     1 121.067  0.9881 0.3222
# scale(SpeciesRichness)      0.0829  0.0829     1 229.814  0.0279 0.8674
# scale(Duration)             2.2269  2.2269     1 135.441  0.7505 0.3878
# scale(MAPmean)              7.9621  7.9621     1  48.193  2.6833 0.1079

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Organic_P_mineralization)
summary(ModelMATmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 2.3822  2.3822     1 118.245  0.8015 0.3725
# scale(SpeciesRichness)      0.0002  0.0002     1 230.995  0.0001 0.9933
# scale(Duration)             1.7806  1.7806     1 140.817  0.5991 0.4402
# scale(MATmean)              4.8624  4.8624     1  34.829  1.6360 0.2093



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Organic_P_mineralization$SpeciesRichness)) ## n = 275
p1 <- ggplot(Organic_P_mineralization, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnOrganic_P_mineralization275")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Organic_P_mineralization$Duration)) ## n = 275
p2 <- ggplot(Organic_P_mineralization, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnOrganic_P_mineralization275")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Organic_P_mineralization$Rotation_cycles)) ## n = 275
p3 <- ggplot(Organic_P_mineralization, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnOrganic_P_mineralization275")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Organic_P_mineralization$Latitude)) ## n = 275
p5 <- ggplot(Organic_P_mineralization, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnOrganic_P_mineralization275")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Organic_P_mineralization$Longitude)) ## n = 275
p6 <- ggplot(Organic_P_mineralization, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnOrganic_P_mineralization275")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Organic_P_mineralization$MAPmean)) ## n = 275
p7 <- ggplot(Organic_P_mineralization, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnOrganic_P_mineralization275")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Organic_P_mineralization$MATmean)) ## n = 275
p8 <- ggplot(Organic_P_mineralization, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnOrganic_P_mineralization275")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Organic_P_mineralization$pHCK)) ## n = 86
p9 <- ggplot(Organic_P_mineralization, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnOrganic_P_mineralization86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Organic_P_mineralization$SOCCK)) ## n = 83
p10 <- ggplot(Organic_P_mineralization, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnOrganic_P_mineralization83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Organic_P_mineralization$TNCK)) ## n = 52
p11 <- ggplot(Organic_P_mineralization, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnOrganic_P_mineralization52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Organic_P_mineralization$NO3CK)) ## n = 46
p12 <- ggplot(Organic_P_mineralization, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnOrganic_P_mineralization46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Organic_P_mineralization$NH4CK)) ## n = 45
p13<- ggplot(Organic_P_mineralization, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnOrganic_P_mineralization45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Organic_P_mineralization$APCK)) ## n = 64
p14 <- ggplot(Organic_P_mineralization, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnOrganic_P_mineralization64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Organic_P_mineralization$AKCK)) ## n = 39
p15 <- ggplot(Organic_P_mineralization, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnOrganic_P_mineralization39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Organic_P_mineralization$ANCK)) ## n = 31
p16 <- ggplot(Organic_P_mineralization, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnOrganic_P_mineralization31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Organic_P_mineralization$RRpH)) ## n = 86
p17 <- ggplot(Organic_P_mineralization, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="RRpH")+
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
sum(!is.na(Organic_P_mineralization$RRSOC)) ## n = 83
p18 <- ggplot(Organic_P_mineralization, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="RRSOC")+
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
sum(!is.na(Organic_P_mineralization$RRTN)) ## n = 52
p19 <- ggplot(Organic_P_mineralization, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="RRTN")+
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
sum(!is.na(Organic_P_mineralization$RRNO3)) ## n = 44
p20 <- ggplot(Organic_P_mineralization, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="RRNO3")+
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
sum(!is.na(Organic_P_mineralization$RRNH4)) ## n = 45
p21 <- ggplot(Organic_P_mineralization, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="RRNH4")+
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
sum(!is.na(Organic_P_mineralization$RRAP)) ## n = 64
p22 <- ggplot(Organic_P_mineralization, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="RRAP")+
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
sum(!is.na(Organic_P_mineralization$RRAK)) ## n = 39
p23 <- ggplot(Organic_P_mineralization, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="RRAK")+
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
sum(!is.na(Organic_P_mineralization$RRAN)) ## n = 31
p24 <- ggplot(Organic_P_mineralization, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnOrganic_P_mineralization", x="RRAN")+
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
sum(!is.na(Organic_P_mineralization$RRYield)) ## n = 49
p25 <- ggplot(Organic_P_mineralization, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-0.18, 0.12), expand=c(0, 0)) + 
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





