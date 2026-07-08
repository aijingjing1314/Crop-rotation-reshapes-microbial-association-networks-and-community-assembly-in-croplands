
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
P_uptake_and_transport <- read.csv("P_uptake_and_transport.csv", fileEncoding = "latin1")
# Check data
head(P_uptake_and_transport)

# 1. The number of Obversation
total_number <- nrow(P_uptake_and_transport)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 275  

# 2. The number of Study
unique_studyid_number <- length(unique(P_uptake_and_transport$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 






#### 8. Subgroup analysis
### 8.1 Longitude_Sub
P_uptake_and_transport_filteredLongitude_Sub <- subset(P_uptake_and_transport, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
P_uptake_and_transport_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(P_uptake_and_transport_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredLongitude_Sub %>%
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

overall_model_P_uptake_and_transport_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  387.8522  -775.7044  -769.7044  -758.8760  -769.6152   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0311     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 914.5996, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.7180, p-val = 0.4236
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0   -0.0072  0.0055  -1.3101  0.1902  -0.0179  0.0035    
# Longitude_SubLongitudeXy0   -0.0004  0.0106  -0.0398  0.9683  -0.0213  0.0204        

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredLongitude_Sub)
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
P_uptake_and_transport_filteredMAPmean2_Sub <- subset(P_uptake_and_transport, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
P_uptake_and_transport_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(P_uptake_and_transport_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredMAPmean2_Sub %>%
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

overall_model_P_uptake_and_transport_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  389.3429  -778.6859  -770.6859  -756.2627  -770.5361   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0309     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 898.3612, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 8.6212, p-val = 0.0348
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb    ci.ub     
# MAPmean2_SubMAP600Dao1200    0.0053  0.0071   0.7433  0.4573  -0.0086   0.0191     
# MAPmean2_SubMAPDy1200       -0.0000  0.0114  -0.0043  0.9966  -0.0224   0.0223     
# MAPmean2_SubMAPXy600        -0.0169  0.0065  -2.5881  0.0097  -0.0297  -0.0041  ** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredMAPmean2_Sub)
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
#           "a"                      "ab"                       "b" 




### 8.3 MATmean_Sub
P_uptake_and_transport_filteredMATmean_Sub <- subset(P_uptake_and_transport, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
P_uptake_and_transport_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(P_uptake_and_transport_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredMATmean_Sub %>%
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

overall_model_P_uptake_and_transport_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  388.5363  -777.0726  -769.0726  -754.6494  -768.9228   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0313     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 905.2882, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 7.3307, p-val = 0.0621
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb    ci.ub    
# MATmean_SubMAT8Dao15    0.0042  0.0089   0.4673  0.6403  -0.0133   0.0216    
# MATmean_SubMATDy15      0.0036  0.0089   0.4097  0.6821  -0.0138   0.0211    
# MATmean_SubMATXy8      -0.0154  0.0064  -2.3969  0.0165  -0.0280  -0.0028  *   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredMATmean_Sub)
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
#                 "a"                 "ab"                  "b" 




### 8.4 LegumeNonlegume
P_uptake_and_transport_filteredLegumeNonlegume <- subset(P_uptake_and_transport, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
P_uptake_and_transport_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(P_uptake_and_transport_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredLegumeNonlegume %>%
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

overall_model_P_uptake_and_transport_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  384.8211  -769.6423  -761.6423  -747.2191  -761.4925   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0309     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 902.7731, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.5922, p-val = 0.2042
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume       -0.0097  0.0056  -1.7394  0.0820  -0.0206  0.0012  . 
# LegumeNonlegumeNon-legume to Legume       -0.0044  0.0051  -0.8584  0.3907  -0.0144  0.0056    
# LegumeNonlegumeNon-legume to Non-legume   -0.0049  0.0053  -0.9283  0.3533  -0.0153  0.0055 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredLegumeNonlegume)
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
    #                  "a"                                     "a"                                     "a" 




### 8.5 AMnonAM
P_uptake_and_transport_filteredAMnonAM <- subset(P_uptake_and_transport, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
P_uptake_and_transport_filteredAMnonAM$AMnonAM <- droplevels(factor(P_uptake_and_transport_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredAMnonAM %>%
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

overall_model_P_uptake_and_transport_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  388.6778  -777.3557  -769.3557  -754.9325  -769.2059   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0310     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 915.8551, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 9.2976, p-val = 0.0256
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb    ci.ub     
# AMnonAMAM to AM      -0.0028  0.0053  -0.5209  0.6025  -0.0133   0.0077     
# AMnonAMAM to nonAM   -0.0231  0.0082  -2.8138  0.0049  -0.0391  -0.0070  ** 
# AMnonAMnonAM to AM   -0.0048  0.0095  -0.5052  0.6134  -0.0235   0.0139     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredAMnonAM)
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
P_uptake_and_transport_filteredC3C4 <- subset(P_uptake_and_transport, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
P_uptake_and_transport_filteredC3C4$C3C4 <- droplevels(factor(P_uptake_and_transport_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredC3C4 %>%
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

overall_model_P_uptake_and_transport_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  392.1442  -784.2883  -776.2883  -761.8799  -776.1380   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0007  0.0267     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 879.3653, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 11.5283, p-val = 0.0092
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb    ci.ub    
# C3C4C3 to C3   -0.0111  0.0056  -1.9927  0.0463  -0.0221  -0.0002  * 
# C3C4C3 to C4   -0.0034  0.0047  -0.7220  0.4703  -0.0126   0.0058    
# C3C4C4 to C3    0.0041  0.0050   0.8241  0.4099  -0.0057   0.0139    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredC3C4)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredC3C4)
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
#        "a"          "a"          "b" 


### 8.7 Annual_Pere
P_uptake_and_transport_filteredAnnual_Pere <- subset(P_uptake_and_transport, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
P_uptake_and_transport_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(P_uptake_and_transport_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredAnnual_Pere %>%
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

overall_model_P_uptake_and_transport_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  389.8706  -779.7412  -771.7412  -757.3180  -771.5914   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0314     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 879.4577, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 13.5728, p-val = 0.0035
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb    ci.ub     
# Annual_PereAnnual to Annual      -0.0037  0.0052  -0.7158  0.4741  -0.0140   0.0065     
# Annual_PereAnnual to Perennial   -0.0032  0.0084  -0.3878  0.6981  -0.0196   0.0132     
# Annual_PerePerennial to Annual   -0.0241  0.0093  -2.6032  0.0092  -0.0423  -0.0060  ** 

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredAnnual_Pere)
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
P_uptake_and_transport_filteredPlant_stage <- subset(P_uptake_and_transport, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
P_uptake_and_transport_filteredPlant_stage$Plant_stage <- droplevels(factor(P_uptake_and_transport_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredPlant_stage %>%
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

overall_model_P_uptake_and_transport_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredPlant_stage)
# Multivariate Meta-Analysis Model (k = 240; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  335.1541  -670.3082  -662.3082  -648.4360  -662.1358   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0011  0.0338     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 237) = 715.6071, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.2971, p-val = 0.1513
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb   ci.ub    
# Plant_stageHarvest              -0.0053  0.0065  -0.8208  0.4118  -0.0181  0.0074    
# Plant_stageReproductive stage   -0.0162  0.0084  -1.9380  0.0526  -0.0326  0.0002  . 
# Plant_stageVegetative stage     -0.0076  0.0066  -1.1446  0.2524  -0.0206  0.0054    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredPlant_stage)
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
P_uptake_and_transport_filteredRotation_cycles2 <- subset(P_uptake_and_transport, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
P_uptake_and_transport_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(P_uptake_and_transport_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredRotation_cycles2 %>%
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

overall_model_P_uptake_and_transport_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  379.3245  -758.6490  -746.6490  -725.0585  -746.3296   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0310     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 907.0038, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 3.2377, p-val = 0.6634
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1      -0.0105  0.0064  -1.6411  0.1008  -0.0231  0.0020    
# Rotation_cycles2D1-3    -0.0105  0.0065  -1.6237  0.1044  -0.0231  0.0022    
# Rotation_cycles2D10     -0.0027  0.0110  -0.2419  0.8088  -0.0242  0.0189    
# Rotation_cycles2D3-5     0.0032  0.0082   0.3947  0.6931  -0.0129  0.0194    
# Rotation_cycles2D5-10    0.0009  0.0079   0.1163  0.9074  -0.0146  0.0164   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredRotation_cycles2)
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
P_uptake_and_transport_filteredDuration2 <- subset(P_uptake_and_transport, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
P_uptake_and_transport_filteredDuration2$Duration2 <- droplevels(factor(P_uptake_and_transport_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredDuration2 %>%
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

overall_model_P_uptake_and_transport_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  384.6451  -769.2903  -759.2903  -741.2797  -759.0639   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0011  0.0325     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 912.5290, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 3.7931, p-val = 0.4347
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5     -0.0114  0.0070  -1.6197  0.1053  -0.0252  0.0024    
# Duration2D11-20    0.0025  0.0093   0.2693  0.7877  -0.0157  0.0207    
# Duration2D20-40   -0.0078  0.0105  -0.7450  0.4563  -0.0284  0.0128    
# Duration2D6-10     0.0028  0.0102   0.2766  0.7821  -0.0172  0.0228    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredDuration2)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredDuration2)
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
P_uptake_and_transport_filteredSpeciesRichness2 <- subset(P_uptake_and_transport, SpeciesRichness2 %in% c("R2", "R3"))
#
P_uptake_and_transport_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(P_uptake_and_transport_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredSpeciesRichness2 %>%
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

overall_model_P_uptake_and_transport_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredSpeciesRichness2)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  379.5165  -759.0331  -753.0331  -742.3052  -752.9408   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0311     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 264) = 890.4397, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.1768, p-val = 0.3368
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2   -0.0068  0.0050  -1.3522  0.1763  -0.0166  0.0030    
# SpeciesRichness2R3   -0.0088  0.0064  -1.3803  0.1675  -0.0214  0.0037  


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredSpeciesRichness2)
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
P_uptake_and_transport_filteredBulk_Rhizosphere <- subset(P_uptake_and_transport, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
P_uptake_and_transport_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(P_uptake_and_transport_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredBulk_Rhizosphere %>%
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

overall_model_P_uptake_and_transport_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 275; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  387.7197  -775.4395  -769.4395  -758.6110  -769.3502   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0010  0.0315     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 273) = 913.4757, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 3.9411, p-val = 0.1394
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere   -0.0037  0.0051  -0.7220  0.4703  -0.0136  0.0063    
# Bulk_RhizosphereRhizosphere       -0.0101  0.0056  -1.7991  0.0720  -0.0211  0.0009  .   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredBulk_Rhizosphere)
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
P_uptake_and_transport_filteredSoil_texture <- subset(P_uptake_and_transport, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
P_uptake_and_transport_filteredSoil_texture$Soil_texture <- droplevels(factor(P_uptake_and_transport_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredSoil_texture %>%
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

overall_model_P_uptake_and_transport_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 225; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  305.8759  -611.7519  -603.7519  -590.1412  -603.5675   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0013  0.0358     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 222) = 751.4526, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.1106, p-val = 0.1639
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse   -0.0049  0.0139  -0.3507  0.7258  -0.0322  0.0224    
# Soil_textureFine     -0.0304  0.0158  -1.9302  0.0536  -0.0613  0.0005  . 
# Soil_textureMedium   -0.0103  0.0092  -1.1233  0.2613  -0.0284  0.0077    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredSoil_texture)
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
P_uptake_and_transport_filteredTillage <- subset(P_uptake_and_transport, Tillage %in% c("Tillage", "No_tillage"))
#
P_uptake_and_transport_filteredTillage$Tillage <- droplevels(factor(P_uptake_and_transport_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredTillage %>%
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

overall_model_P_uptake_and_transport_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredTillage)
# Multivariate Meta-Analysis Model (k = 151; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  235.1014  -470.2027  -464.2027  -455.1909  -464.0372   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0035  0.0595     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 149) = 390.6578, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 3.2365, p-val = 0.1982
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0194  0.0171  -1.1337  0.2569  -0.0529  0.0141    
# TillageTillage      -0.0244  0.0170  -1.4332  0.1518  -0.0578  0.0090    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredTillage)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredTillage)
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
P_uptake_and_transport_filteredStraw_retention <- subset(P_uptake_and_transport, Straw_retention %in% c("Retention", "No_retention"))
#
P_uptake_and_transport_filteredStraw_retention$Straw_retention <- droplevels(factor(P_uptake_and_transport_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredStraw_retention %>%
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

overall_model_P_uptake_and_transport_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 55; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   96.6969  -193.3938  -187.3938  -181.4829  -186.9040   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0002  0.0128     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 174.7688, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 6.7640, p-val = 0.0340
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb    ci.ub    
# Straw_retentionNo_retention   -0.0035  0.0056  -0.6220  0.5339  -0.0146   0.0075    
# Straw_retentionRetention      -0.0119  0.0049  -2.4378  0.0148  -0.0215  -0.0023  * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredStraw_retention)
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
P_uptake_and_transport_filteredPrimer <- subset(P_uptake_and_transport, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
P_uptake_and_transport_filteredPrimer$Primer <- droplevels(factor(P_uptake_and_transport_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- P_uptake_and_transport_filteredPrimer %>%
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

overall_model_P_uptake_and_transport_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = P_uptake_and_transport_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_P_uptake_and_transport_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 265; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  392.6423  -785.2846  -777.2846  -763.0112  -777.1290   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0012  0.0339     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 262) = 820.0101, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.0620, p-val = 0.5596
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4   -0.0028  0.0069  -0.4120  0.6803  -0.0164  0.0107    
# PrimerV4      -0.0174  0.0128  -1.3613  0.1734  -0.0425  0.0077    
# PrimerV4-V5   -0.0025  0.0129  -0.1977  0.8433  -0.0278  0.0227 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_P_uptake_and_transport_filteredPrimer)
vcov_rotation <- vcov(overall_model_P_uptake_and_transport_filteredPrimer)
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
P_uptake_and_transport$Wr <- 1 / P_uptake_and_transport$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = P_uptake_and_transport)

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
# Model1   Model1 -896.2029 -874.5022 454.1014
# Model2   Model2 -899.8324 -878.1318 455.9162
# Model3   Model3 -896.5331 -874.8325 454.2666
# Model4   Model4 -896.1085 -874.4079 454.0543
# Model5   Model5 -897.3607 -875.6600 454.6803
# Model6   Model6 -897.2298 -875.5292 454.6149
# Model7   Model7 -900.0031 -878.3025 456.0015##############################
# Model8   Model8 -896.4363 -874.7356 454.2181
# Model9   Model9 -852.4098 -816.2421 436.2049
# Model10 Model10 -854.9074 -818.7397 437.4537

##### Model 7 is the best model
summary(Model7)
# Number of obs: 275, groups:  StudyID, 50
anova(Model7) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 10.8972 10.8972     1 42.961  3.7104 0.06071 .
# scale(SpeciesRichness)       3.9416  3.9416     1 92.850  1.3421 0.24964  
# scale(log(Duration))        11.2685 11.2685     1 49.900  3.8368 0.05574 .

#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(pHCK) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.08588 0.08588     1  7.812  0.0270 0.8738
# scale(SpeciesRichness)      0.21481 0.21481     1 41.560  0.0674 0.7964
# scale(log(Duration))        0.10021 0.10021     1 10.081  0.0315 0.8627
# scale(pHCK)                 0.86340 0.86340     1  8.930  0.2710 0.6153

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.33765 1.33765     1  5.1437  0.4269 0.5416
# scale(SpeciesRichness)      0.98826 0.98826     1 15.1276  0.3154 0.5826
# scale(log(Duration))        1.43694 1.43694     1  6.4537  0.4586 0.5218
# scale(SOCCK)                0.00229 0.00229     1  2.0576  0.0007 0.9808

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(TNCK) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 4.1312  4.1312     1  4.455  1.0659 0.3546
# scale(SpeciesRichness)      2.8509  2.8509     1 43.619  0.7356 0.3958
# scale(log(Duration))        0.2654  0.2654     1 18.210  0.0685 0.7965
# scale(TNCK)                 0.3342  0.3342     1 46.090  0.0862 0.7703


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.36169 0.36169     1  3.3961  0.3128 0.6108
# scale(SpeciesRichness)      0.13815 0.13815     1 22.9004  0.1195 0.7328
# scale(log(Duration))        0.22697 0.22697     1  5.1566  0.1963 0.6757
# scale(NO3CK)                0.35745 0.35745     1 11.1849  0.3091 0.5892

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 1.29668 1.29668     1  3.6911  1.2513 0.3308
# scale(SpeciesRichness)      0.14979 0.14979     1 21.6384  0.1446 0.7075
# scale(log(Duration))        0.95355 0.95355     1  5.5086  0.9202 0.3776
# scale(NH4CK)                0.03623 0.03623     1  4.1012  0.0350 0.8606

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(APCK) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.6650  0.6650     1  5.616  0.1924 0.6773
# scale(SpeciesRichness)      1.9522  1.9522     1 22.905  0.5647 0.4600
# scale(log(Duration))        0.0003  0.0003     1  5.112  0.0001 0.9935
# scale(APCK)                 6.0074  6.0074     1 48.984  1.7376 0.1936

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(AKCK) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 2.95187 2.95187     1  4.9267  0.7322 0.4318
# scale(SpeciesRichness)      1.12194 1.12194     1 27.7401  0.2783 0.6020
# scale(log(Duration))        3.05561 3.05561     1  6.3268  0.7579 0.4158
# scale(AKCK)                 0.02571 0.02571     1 18.6582  0.0064 0.9372

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(ANCK) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 6.8062  6.8062     1    26  1.3546 0.2550
# scale(SpeciesRichness)      7.3427  7.3427     1    26  1.4614 0.2376
# scale(log(Duration))        9.5823  9.5823     1    26  1.9071 0.1790
# scale(ANCK)                 0.7954  0.7954     1    26  0.1583 0.6940

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(Latitude) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
summary(ModelLatitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLatitude)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 11.6270 11.6270     1 40.638  3.9562 0.05346 .
# scale(SpeciesRichness)       4.8447  4.8447     1 97.274  1.6485 0.20222  
# scale(log(Duration))        12.4708 12.4708     1 46.877  4.2433 0.04498 *
# scale(Latitude)              1.8450  1.8450     1 34.648  0.6278 0.43356  

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(Longitude) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
summary(ModelLongitude)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 12.3205 12.3205     1 40.214  4.1869 0.04731 *
# scale(SpeciesRichness)       3.7249  3.7249     1 78.079  1.2658 0.26399  
# scale(log(Duration))        15.6732 15.6732     1 55.040  5.3263 0.02479 *
# scale(Longitude)             4.5737  4.5737     1 11.464  1.5543 0.23737  

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
summary(ModelMAPmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 12.2767 12.2767     1 42.621  4.1946 0.04675 *
# scale(SpeciesRichness)       4.4335  4.4335     1 93.215  1.5148 0.22151  
# scale(log(Duration))        13.8620 13.8620     1 55.787  4.7362 0.03378 *
# scale(MAPmean)               4.3154  4.3154     1 17.785  1.4744 0.24052  

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(MATmean) + (1 | StudyID), weights = Wr, data = P_uptake_and_transport)
summary(ModelMATmean)
# Number of obs: 275, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 11.2033 11.2033     1 31.633  3.8406 0.05889 .
# scale(SpeciesRichness)       5.7633  5.7633     1 72.478  1.9757 0.16412  
# scale(log(Duration))        11.6272 11.6272     1 38.739  3.9859 0.05295 .
# scale(MATmean)              16.4890 16.4890     1 17.381  5.6526 0.02915 *



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(P_uptake_and_transport$SpeciesRichness)) ## n = 275
p1 <- ggplot(P_uptake_and_transport, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnP_uptake_and_transport275")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(P_uptake_and_transport$Duration)) ## n = 275
p2 <- ggplot(P_uptake_and_transport, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnP_uptake_and_transport275")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(P_uptake_and_transport$Rotation_cycles)) ## n = 275
p3 <- ggplot(P_uptake_and_transport, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnP_uptake_and_transport275")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(P_uptake_and_transport$Latitude)) ## n = 275
p5 <- ggplot(P_uptake_and_transport, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnP_uptake_and_transport275")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(P_uptake_and_transport$Longitude)) ## n = 275
p6 <- ggplot(P_uptake_and_transport, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnP_uptake_and_transport275")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(P_uptake_and_transport$MAPmean)) ## n = 275
p7 <- ggplot(P_uptake_and_transport, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnP_uptake_and_transport275")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(P_uptake_and_transport$MATmean)) ## n = 275
p8 <- ggplot(P_uptake_and_transport, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnP_uptake_and_transport275")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(P_uptake_and_transport$pHCK)) ## n = 86
p9 <- ggplot(P_uptake_and_transport, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnP_uptake_and_transport86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(P_uptake_and_transport$SOCCK)) ## n = 83
p10 <- ggplot(P_uptake_and_transport, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnP_uptake_and_transport83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(P_uptake_and_transport$TNCK)) ## n = 52
p11 <- ggplot(P_uptake_and_transport, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnP_uptake_and_transport52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(P_uptake_and_transport$NO3CK)) ## n = 46
p12 <- ggplot(P_uptake_and_transport, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnP_uptake_and_transport46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(P_uptake_and_transport$NH4CK)) ## n = 45
p13<- ggplot(P_uptake_and_transport, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnP_uptake_and_transport45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(P_uptake_and_transport$APCK)) ## n = 64
p14 <- ggplot(P_uptake_and_transport, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnP_uptake_and_transport64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(P_uptake_and_transport$AKCK)) ## n = 39
p15 <- ggplot(P_uptake_and_transport, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnP_uptake_and_transport39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(P_uptake_and_transport$ANCK)) ## n = 31
p16 <- ggplot(P_uptake_and_transport, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnP_uptake_and_transport31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(P_uptake_and_transport$RRpH)) ## n = 86
p17 <- ggplot(P_uptake_and_transport, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="RRpH")+
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
sum(!is.na(P_uptake_and_transport$RRSOC)) ## n = 83
p18 <- ggplot(P_uptake_and_transport, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="RRSOC")+
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
sum(!is.na(P_uptake_and_transport$RRTN)) ## n = 52
p19 <- ggplot(P_uptake_and_transport, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="RRTN")+
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
sum(!is.na(P_uptake_and_transport$RRNO3)) ## n = 44
p20 <- ggplot(P_uptake_and_transport, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="RRNO3")+
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
sum(!is.na(P_uptake_and_transport$RRNH4)) ## n = 45
p21 <- ggplot(P_uptake_and_transport, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="RRNH4")+
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
sum(!is.na(P_uptake_and_transport$RRAP)) ## n = 64
p22 <- ggplot(P_uptake_and_transport, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="RRAP")+
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
sum(!is.na(P_uptake_and_transport$RRAK)) ## n = 39
p23 <- ggplot(P_uptake_and_transport, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="RRAK")+
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
sum(!is.na(P_uptake_and_transport$RRAN)) ## n = 31
p24 <- ggplot(P_uptake_and_transport, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnP_uptake_and_transport", x="RRAN")+
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
sum(!is.na(P_uptake_and_transport$RRYield)) ## n = 49
p25 <- ggplot(P_uptake_and_transport, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-0.21, 0.13), expand=c(0, 0)) + 
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






