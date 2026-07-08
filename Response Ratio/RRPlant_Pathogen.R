
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Plant_Pathogen <- read.csv("Plant_Pathogen.csv", fileEncoding = "latin1")
# Check data
head(Plant_Pathogen)

# 1. The number of Obversation
total_number <- nrow(Plant_Pathogen)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 312  

# 2. The number of Study
unique_studyid_number <- length(unique(Plant_Pathogen$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 45 





#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Plant_Pathogen_filteredLongitude_Sub <- subset(Plant_Pathogen, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Plant_Pathogen_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Plant_Pathogen_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredLongitude_Sub %>%
  group_by(Longitude_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Sub Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 LongitudeDa0           139             39
# 2 LongitudeXy0           173              6

overall_model_Plant_Pathogen_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 312; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -823.0969  1646.1937  1652.1937  1663.4035  1652.2722   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3326  0.5767     45     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 310) = 3265.2183, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 4.1755, p-val = 0.1240
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0   -0.1748  0.0973  -1.7968  0.0724  -0.3655  0.0159  . 
# Longitude_SubLongitudeXy0    0.2367  0.2432   0.9732  0.3305  -0.2400  0.7134    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredLongitude_Sub)
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
Plant_Pathogen_filteredMAPmean2_Sub <- subset(Plant_Pathogen, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Plant_Pathogen_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Plant_Pathogen_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredMAPmean2_Sub %>%
  group_by(MAPmean2_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MAPmean2_Sub  Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 MAP600Dao1200          122             17
# 2 MAPDy1200              106             10
# 3 MAPXy600                84             19

overall_model_Plant_Pathogen_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 312; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -821.0916  1642.1832  1650.1832  1665.1166  1650.3148   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3833  0.6191     45     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 309) = 3238.9771, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 6.3323, p-val = 0.0965
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean2_SubMAP600Dao1200   -0.2527  0.1395  -1.8112  0.0701  -0.5262  0.0208  . 
# MAPmean2_SubMAPDy1200       -0.2682  0.2061  -1.3012  0.1932  -0.6722  0.1358    
# MAPmean2_SubMAPXy600         0.0846  0.1344   0.6298  0.5288  -0.1787  0.3479     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredMAPmean2_Sub)
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
Plant_Pathogen_filteredMATmean_Sub <- subset(Plant_Pathogen, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Plant_Pathogen_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Plant_Pathogen_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredMATmean_Sub %>%
  group_by(MATmean_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Sub Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 MAT8Dao15             40             11
# 2 MATDy15              115             15
# 3 MATXy8               157             20

overall_model_Plant_Pathogen_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 312; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -821.1472  1642.2943  1650.2943  1665.2277  1650.4259   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3872  0.6223     45     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 309) = 3280.1345, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 6.7719, p-val = 0.0795
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubMAT8Dao15   -0.3161  0.1650  -1.9163  0.0553  -0.6394  0.0072  . 
# MATmean_SubMATDy15     -0.2182  0.1673  -1.3044  0.1921  -0.5461  0.1097    
# MATmean_SubMATXy8       0.0703  0.1344   0.5230  0.6010  -0.1932  0.3338

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredMATmean_Sub)
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
#        "a"                 "ab"                  "b"




### 8.4 LegumeNonlegume
Plant_Pathogen_filteredLegumeNonlegume <- subset(Plant_Pathogen, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Plant_Pathogen_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Plant_Pathogen_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredLegumeNonlegume %>%
  group_by(LegumeNonlegume) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   LegumeNonlegume          Observations Unique_StudyID
#   <fct>                           <int>          <int>
# 1 Legume to Non-legume               55             10
# 2 Non-legume to Legume              152             15
# 3 Non-legume to Non-legume          105             29

overall_model_Plant_Pathogen_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 312; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -768.3850  1536.7699  1544.7699  1559.7033  1544.9015   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3942  0.6278     45     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 309) = 3357.1513, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 117.7805, p-val < .0001
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb    ci.ub      
# LegumeNonlegumeLegume to Non-legume       -0.1541  0.1080  -1.4272  0.1535  -0.3657   0.0575      
# LegumeNonlegumeNon-legume to Legume       -0.4166  0.1016  -4.1006  <.0001  -0.6158  -0.2175  *** 
# LegumeNonlegumeNon-legume to Non-legume    0.0209  0.0997   0.2096  0.8340  -0.1745   0.2163      
 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredLegumeNonlegume)
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
    #             "a"                                     "b"                                     "c" 




### 8.5 AMnonAM
Plant_Pathogen_filteredAMnonAM <- subset(Plant_Pathogen, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Plant_Pathogen_filteredAMnonAM$AMnonAM <- droplevels(factor(Plant_Pathogen_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             271             37
# 2 AM to nonAM           25              7
# 3 nonAM to AM           16              6

overall_model_Plant_Pathogen_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 312; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -804.2058  1608.4116  1616.4116  1631.3449  1616.5432   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3582  0.5985     45     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 309) = 3390.7249, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 45.7116, p-val < .0001
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb    ci.ub      
# AMnonAMAM to AM      -0.2238  0.0948  -2.3600  0.0183  -0.4097  -0.0379    * 
# AMnonAMAM to nonAM    0.0179  0.1080   0.1660  0.8681  -0.1937   0.2295      
# AMnonAMnonAM to AM    0.4697  0.1402   3.3507  0.0008   0.1949   0.7444  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredAMnonAM)
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
   #       "a"                "b"                "c" 


### 8.6 C3C4
Plant_Pathogen_filteredC3C4 <- subset(Plant_Pathogen, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Plant_Pathogen_filteredC3C4$C3C4 <- droplevels(factor(Plant_Pathogen_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredC3C4 %>%
  group_by(C3C4) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   C3C4     Observations Unique_StudyID
#   <fct>           <int>          <int>
# 1 C3 to C3           68             17
# 2 C3 to C4           91             22
# 3 C4 to C3          148             13

overall_model_Plant_Pathogen_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 307; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -749.1572  1498.3144  1506.3144  1521.1825  1506.4481   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.8485  0.9211     43     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 304) = 3264.4569, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 152.9558, p-val < .0001
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb    ci.ub      
# C3C4C3 to C3    0.7902  0.1671   4.7299  <.0001   0.4627   1.1176  *** 
# C3C4C3 to C4   -0.5168  0.1513  -3.4166  0.0006  -0.8133  -0.2203  *** 
# C3C4C4 to C3   -0.8101  0.1538  -5.2682  <.0001  -1.1115  -0.5087  *** 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredC3C4)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredC3C4)
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
#      "a"          "b"          "c" 


### 8.7 Annual_Pere
Plant_Pathogen_filteredAnnual_Pere <- subset(Plant_Pathogen, Annual_Pere %in% c("Annual to Annual",  "Perennial to Annual"))
#
Plant_Pathogen_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Plant_Pathogen_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredAnnual_Pere %>%
  group_by(Annual_Pere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Annual_Pere            Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 Annual to Annual                263             33
# 2 Annual to Perennial               9              6
# 3 Perennial to Annual              34              9
# 4 Perennial to Perennial            6              3

overall_model_Plant_Pathogen_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredAnnual_Pere)
# ultivariate Meta-Analysis Model (k = 297; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -665.8089  1331.6177  1337.6177  1348.6787  1337.7002   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1825  0.4272     42     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 295) = 2542.0411, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 12.0035, p-val = 0.0025
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb    ci.ub      
# Annual_PereAnnual to Annual       0.0078  0.0812   0.0958  0.9237  -0.1514   0.1669      
# Annual_PerePerennial to Annual   -0.5423  0.1566  -3.4633  0.0005  -0.8492  -0.2354  *** 

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredAnnual_Pere)
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
     # Annual_PereAnnual to Annual Annual_PerePerennial to Annual 
     #                       "a"                            "b" 

### 8.8 Plant_stage
Plant_Pathogen_filteredPlant_stage <- subset(Plant_Pathogen, Plant_stage %in% c("Vegetative stage","Reproductive stage","Harvest"))
#
Plant_Pathogen_filteredPlant_stage$Plant_stage <- droplevels(factor(Plant_Pathogen_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredPlant_stage %>%
  group_by(Plant_stage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
  # Plant_stage        Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Harvest                     117             21
# 2 Maturity stage                9              6
# 3 Reproductive stage           54             10
# 4 Vegetative stage             92              9

overall_model_Plant_Pathogen_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredPlant_stage)
# tivariate Meta-Analysis Model (k = 263; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -763.6792  1527.3584  1535.3584  1549.6011  1535.5152   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4130  0.6426     31     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 260) = 2838.5857, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 12.9714, p-val = 0.0047
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb   ci.ub    
# Plant_stageHarvest              -0.1307  0.1215  -1.0757  0.2821  -0.3687  0.1074    
# Plant_stageReproductive stage   -0.1183  0.1268  -0.9325  0.3511  -0.3668  0.1303    
# Plant_stageVegetative stage      0.0191  0.1249   0.1528  0.8786  -0.2257  0.2638    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredPlant_stage)
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
            #               "a"                           "a"                           "b" 

### 8.9 Rotation_cycles2
Plant_Pathogen_filteredRotation_cycles2 <- subset(Plant_Pathogen, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
Plant_Pathogen_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Plant_Pathogen_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredRotation_cycles2 %>%
  group_by(Rotation_cycles2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Rotation_cycles2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                         79             18
# 2 D1-3                       99             15
# 3 D10                        20              7
# 4 D3-5                       80              7
# 5 D5-10                      34             11

overall_model_Plant_Pathogen_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 312; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -806.9863  1613.9726  1625.9726  1648.3337  1626.2526   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3667  0.6056     45     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 307) = 3302.3838, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 40.6805, p-val < .0001
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb    ci.ub    
# Rotation_cycles2D1      -0.1328  0.1073  -1.2373  0.2160  -0.3432   0.0776    
# Rotation_cycles2D1-3    -0.0502  0.1066  -0.4711  0.6376  -0.2591   0.1587    
# Rotation_cycles2D10     -0.1993  0.1350  -1.4761  0.1399  -0.4638   0.0653    
# Rotation_cycles2D3-5    -0.2965  0.1259  -2.3546  0.0185  -0.5434  -0.0497  * 
# Rotation_cycles2D5-10   -0.0113  0.1229  -0.0919  0.9268  -0.2521   0.2295      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredRotation_cycles2)
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
                 # "ab"                  "cd"                  "ac"                  "ac"                  "bd" 


### 8.10 Duration2
Plant_Pathogen_filteredDuration2 <- subset(Plant_Pathogen, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Plant_Pathogen_filteredDuration2$Duration2 <- droplevels(factor(Plant_Pathogen_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredDuration2 %>%
  group_by(Duration2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Duration2 Observations Unique_StudyID
#   <fct>            <int>          <int>
# 1 D1-5               156             25
# 2 D11-20              35             11
# 3 D20-40              87              6
# 4 D6-10               34             10
overall_model_Plant_Pathogen_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 312; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -820.7045  1641.4091  1651.4091  1670.0596  1651.6078   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3530  0.5941     45     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 308) = 3315.2470, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 10.5330, p-val = 0.0323
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5     -0.1243  0.1086  -1.1444  0.2524  -0.3373  0.0886    
# Duration2D11-20   -0.2042  0.1147  -1.7813  0.0749  -0.4290  0.0205  . 
# Duration2D20-40    0.0852  0.1511   0.5640  0.5728  -0.2109  0.3813    
# Duration2D6-10    -0.1031  0.1191  -0.8663  0.3863  -0.3365  0.1302    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredDuration2)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredDuration2)
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
 #       "ab"             "a"             "b"             "b" 

### 8.11 SpeciesRichness2
Plant_Pathogen_filteredSpeciesRichness2 <- subset(Plant_Pathogen, SpeciesRichness2 %in% c("R2", "R3"))
#
Plant_Pathogen_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Plant_Pathogen_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredSpeciesRichness2 %>%
  group_by(SpeciesRichness2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   SpeciesRichness2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 R2                        278             42
# 2 R3                         24              4
# 3 R4                         10              2

overall_model_Plant_Pathogen_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredSpeciesRichness2)
# riate Meta-Analysis Model (k = 302; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -820.2116  1640.4232  1646.4232  1657.5346  1646.5043   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.3577  0.5981     43     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 300) = 3410.3916, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 7.5117, p-val = 0.0234
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb    ci.ub    
# SpeciesRichness2R2   -0.1080  0.0956  -1.1302  0.2584  -0.2954   0.0793    
# SpeciesRichness2R3   -0.2669  0.1136  -2.3484  0.0189  -0.4896  -0.0441  * 



# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredSpeciesRichness2)
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
#      "a"                "b"               



### 8.12 Bulk_Rhizosphere
Plant_Pathogen_filteredBulk_Rhizosphere <- subset(Plant_Pathogen, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Plant_Pathogen_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Plant_Pathogen_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non_Rhizosphere           250             33
# 2 Rhizosphere                62             20

overall_model_Plant_Pathogen_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 312; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -771.8234  1543.6469  1549.6469  1560.8566  1549.7253   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4592  0.6777     45     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 310) = 3437.5959, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 111.2983, p-val < .0001
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb    ci.ub      
# Bulk_RhizosphereNon_Rhizosphere    0.0724  0.1064   0.6806  0.4961  -0.1362   0.2810      
# Bulk_RhizosphereRhizosphere       -0.4808  0.1105  -4.3524  <.0001  -0.6974  -0.2643  *** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredBulk_Rhizosphere)
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
#               "a"                             "b" 


### 8.13 Soil_texture
Plant_Pathogen_filteredSoil_texture <- subset(Plant_Pathogen, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Plant_Pathogen_filteredSoil_texture$Soil_texture <- droplevels(factor(Plant_Pathogen_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredSoil_texture %>%
  group_by(Soil_texture) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
  # Soil_texture Observations Unique_StudyID
#   <fct>               <int>          <int>
# 1 Coarse                 30              5
# 2 Fine                   97             10
# 3 Medium                129             11

overall_model_Plant_Pathogen_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 256; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -722.1296  1444.2591  1452.2591  1466.3927  1452.4204   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4326  0.6577     26     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 253) = 2727.2424, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.0195, p-val = 0.1704
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse   -0.2036  0.2966  -0.6863  0.4925  -0.7850  0.3778    
# Soil_textureFine     -0.4182  0.2181  -1.9180  0.0551  -0.8456  0.0092  . 
# Soil_textureMedium    0.1914  0.2052   0.9326  0.3510  -0.2108  0.5936


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredSoil_texture)
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
#      "ab"                "a"                "b"


### 8.14 Tillage
Plant_Pathogen_filteredTillage <- subset(Plant_Pathogen, Tillage %in% c("Tillage", "No_tillage"))
#
Plant_Pathogen_filteredTillage$Tillage <- droplevels(factor(Plant_Pathogen_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredTillage %>%
  group_by(Tillage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Tillage    Observations Unique_StudyID
#   <fct>             <int>          <int>
# 1 No_tillage          103              5
# 2 Tillage              91              8

overall_model_Plant_Pathogen_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredTillage)
# Multivariate Meta-Analysis Model (k = 194; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -304.3440   608.6881   614.6881   624.4605   614.8157   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.4734  0.6881     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 192) = 1036.1175, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 70.6838, p-val < .0001
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.4246  0.2220  -1.9124  0.0558  -0.8598  0.0106  . 
# TillageTillage       0.2767  0.2191   1.2630  0.2066  -0.1527  0.7061   
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredTillage)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredTillage)
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
Plant_Pathogen_filteredStraw_retention <- subset(Plant_Pathogen, Straw_retention %in% c("Retention", "No_retention"))
#
Plant_Pathogen_filteredStraw_retention$Straw_retention <- droplevels(factor(Plant_Pathogen_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredStraw_retention %>%
  group_by(Straw_retention) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Straw_retention Observations Unique_StudyID
#   <fct>                  <int>          <int>
# 1 No_retention              11              5
# 2 Retention                 34              8

overall_model_Plant_Pathogen_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 45; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -133.1594   266.3188   272.3188   277.6024   272.9342   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.2691  0.5188     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 43) = 435.2584, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.5125, p-val = 0.7739
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention   -0.1577  0.2210  -0.7135  0.4756  -0.5908  0.2754    
# Straw_retentionRetention      -0.0379  0.1840  -0.2058  0.8370  -0.3985  0.3228 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredStraw_retention)
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
Plant_Pathogen_filteredPrimer <- subset(Plant_Pathogen, Primer %in% c("ITS1", "ITS2"))
#
Plant_Pathogen_filteredPrimer$Primer <- droplevels(factor(Plant_Pathogen_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Plant_Pathogen_filteredPrimer %>%
  group_by(Primer) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Primer             Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 ITS1                        226             40
# 2 ITS1 + 5.8S + ITS2            3              1
# 3 ITS2                         83              4

overall_model_Plant_Pathogen_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Plant_Pathogen_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Plant_Pathogen_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 309; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -703.5432  1407.0865  1413.0865  1424.2670  1413.1657   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1671  0.4088     44     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 307) = 2295.1557, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.9929, p-val = 0.3692
# 
# Model Results:
# 
#             estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerITS1   -0.0814  0.0708  -1.1494  0.2504  -0.2201  0.0574    
# PrimerITS2    0.1766  0.2154   0.8196  0.4124  -0.2457  0.5988         
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Plant_Pathogen_filteredPrimer)
vcov_rotation <- vcov(overall_model_Plant_Pathogen_filteredPrimer)
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
  # PrimerITS1 PrimerITS2 
  #      "a"        "a" 









#### 9. Linear Mixed Effect Model
# 
Plant_Pathogen$Wr <- 1 / Plant_Pathogen$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Plant_Pathogen)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Plant_Pathogen)

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
#           Model      AIC      BIC    logLik
# Model1   Model1 700.0860 722.5440 -344.0430
# Model2   Model2 696.2983 718.7564 -342.1492
# Model3   Model3 698.0972 720.5553 -343.0486
# Model4   Model4 700.1599 722.6180 -344.0800
# Model5   Model5 700.3603 722.8183 -344.1801
# Model6   Model6 700.4380 722.8960 -344.2190
# Model7   Model7 696.1872 718.6452 -342.0936######################
# Model8   Model8 698.1914 720.6494 -343.0957
# Model9   Model9 707.5260 744.9560 -343.7630
# Model10 Model10 703.8537 741.2838 -341.9269

##### Model 7 is the best model
summary(Model7)
# Number of obs: 312, groups:  StudyID, 45
anova(Model7)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 32.144  32.144     1 304.34  5.0900 0.02477 *
# scale(SpeciesRichness)       0.041   0.041     1 224.00  0.0065 0.93568  
# scale(log(Duration))        17.532  17.532     1 189.15  2.7763 0.09733 .


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
summary(ModelpH)
# Number of obs: 92, groups:  StudyID, 29
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  0.8754  0.8754     1 78.905  0.1182 0.73188  
# scale(SpeciesRichness)      27.7097 27.7097     1 72.295  3.7423 0.05697 .
# scale(log(Duration))         0.0323  0.0323     1 53.408  0.0044 0.94757  
# scale(pHCK)                 19.0553 19.0553     1 58.749  2.5735 0.11403  

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
summary(ModelSOC)
# Number of obs: 102, groups:  StudyID, 29
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  0.5330  0.5330     1 83.517  0.0765 0.78274  
# scale(SpeciesRichness)      20.5623 20.5623     1 96.008  2.9524 0.08898 .
# scale(log(Duration))         1.1654  1.1654     1 50.260  0.1673 0.68423  
# scale(SOCCK)                 7.3569  7.3569     1 29.981  1.0563 0.31228  

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
summary(ModelTN)
# Number of obs: 64, groups:  StudyID, 18
anova(ModelTN) 
# > anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.63663 0.63663     1 58.693  0.3948 0.5322
# scale(SpeciesRichness)      2.12473 2.12473     1 14.473  1.3176 0.2696
# scale(log(Duration))        0.39337 0.39337     1 26.464  0.2439 0.6254
# scale(TNCK)                 0.12618 0.12618     1  6.129  0.0783 0.7889


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
summary(ModelNO3)
# Number of obs: 29, groups:  StudyID, 12
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 22.7250 22.7250     1 10.5746  6.4218 0.02853 *
# scale(SpeciesRichness)      12.1796 12.1796     1  9.2324  3.4418 0.09572 .
# scale(log(Duration))        18.6825 18.6825     1 10.3069  5.2794 0.04370 *
# scale(NO3CK)                 0.9767  0.9767     1 23.8258  0.2760 0.60419  

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
summary(ModelNH4)
# Number of obs: 29, groups:  StudyID, 12
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 19.0629 19.0629     1    24  5.6891 0.02532 *
# scale(SpeciesRichness)      12.0737 12.0737     1    24  3.6032 0.06976 .
# scale(log(Duration))         8.5889  8.5889     1    24  2.5632 0.12246  
# scale(NH4CK)                 9.5676  9.5676     1    24  2.8553 0.10402  

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(APCK) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
summary(ModelAP)
# Number of obs: 57, groups:  StudyID, 24
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 11.8985 11.8985     1 25.046  1.3203 0.2614
# scale(SpeciesRichness)       2.9053  2.9053     1 50.678  0.3224 0.5727
# scale(log(Duration))        13.4312 13.4312     1 24.439  1.4904 0.2338
# scale(APCK)                 10.8597 10.8597     1 21.347  1.2051 0.2845

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
summary(ModelAK)
# Number of obs: 43, groups:  StudyID, 19
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 17.7358 17.7358     1 16.928  1.3898 0.2548
# scale(SpeciesRichness)       4.4579  4.4579     1 37.281  0.3493 0.5581
# scale(log(Duration))        20.9669 20.9669     1 15.669  1.6429 0.2186
# scale(AKCK)                  2.8785  2.8785     1 26.371  0.2256 0.6387

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
summary(ModelAN)
# Number of obs: 42, groups:  StudyID, 13
anova(ModelAN)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
# scale(log(Rotation_cycles)) 468.45  468.45     1    37 27.1782 7.286e-06 ***
# scale(SpeciesRichness)      167.00  167.00     1    37  9.6887  0.003566 ** 
# scale(log(Duration))        482.21  482.21     1    37 27.9764 5.754e-06 ***
# scale(ANCK)                  15.14   15.14     1    37  0.8784  0.354710

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
summary(ModelLatitude)
# Number of obs: 312, groups:  StudyID, 45
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 33.379  33.379     1 302.868  5.2996 0.02201 *
# scale(SpeciesRichness)       0.025   0.025     1 220.692  0.0040 0.94945  
# scale(log(Duration))        22.318  22.318     1 196.549  3.5434 0.06126 .
# scale(Latitude)             13.508  13.508     1  37.533  2.1446 0.15140  

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
summary(ModelLongitude)
# Number of obs: 312, groups:  StudyID, 45
anova(ModelLongitude)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 36.167  36.167     1 304.036  5.7533 0.01706 *
# scale(SpeciesRichness)       0.046   0.046     1 219.766  0.0073 0.93186  
# scale(log(Duration))        25.545  25.545     1 207.444  4.0637 0.04511 *
# scale(Longitude)            18.586  18.586     1  33.657  2.9567 0.09471 .

# ### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
summary(ModelMAPmean)
# Number of obs: 312, groups:  StudyID, 45
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 33.033  33.033     1 303.175  5.2397 0.02276 *
# scale(SpeciesRichness)       0.000   0.000     1 226.010  0.0000 0.99642  
# scale(log(Duration))        19.914  19.914     1 190.026  3.1588 0.07712 .
# scale(MAPmean)               7.907   7.907     1  40.821  1.2541 0.26931 

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Plant_Pathogen)
summary(ModelMATmean)
# Number of obs: 312, groups:  StudyID, 45
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 33.897  33.897     1 303.13  5.3783 0.02105 *
# scale(SpeciesRichness)       0.005   0.005     1 226.75  0.0008 0.97692  
# scale(log(Duration))        20.923  20.923     1 193.84  3.3198 0.06999 .
# scale(MATmean)              10.690  10.690     1  36.35  1.6961 0.20099  
# ---


############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Plant_Pathogen$SpeciesRichness)) ## n = 312
p1 <- ggplot(Plant_Pathogen, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnPlant_Pathogen312")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Plant_Pathogen$Duration)) ## n = 312
p2 <- ggplot(Plant_Pathogen, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnPlant_Pathogen312")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Plant_Pathogen$Rotation_cycles)) ## n = 312
p3 <- ggplot(Plant_Pathogen, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnPlant_Pathogen312")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Plant_Pathogen$Latitude)) ## n = 312
p5 <- ggplot(Plant_Pathogen, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnPlant_Pathogen312")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Plant_Pathogen$Longitude)) ## n = 312
p6 <- ggplot(Plant_Pathogen, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnPlant_Pathogen312")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Plant_Pathogen$MAPmean)) ## n = 312
p7 <- ggplot(Plant_Pathogen, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnPlant_Pathogen312")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Plant_Pathogen$MATmean)) ## n = 312
p8 <- ggplot(Plant_Pathogen, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnPlant_Pathogen312")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Plant_Pathogen$pHCK)) ## n = 92
p9 <- ggplot(Plant_Pathogen, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnPlant_Pathogen92")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Plant_Pathogen$SOCCK)) ## n = 102
p10 <- ggplot(Plant_Pathogen, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnPlant_Pathogen102")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Plant_Pathogen$TNCK)) ## n = 64
p11 <- ggplot(Plant_Pathogen, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnPlant_Pathogen64")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Plant_Pathogen$NO3CK)) ## n = 29
p12 <- ggplot(Plant_Pathogen, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnPlant_Pathogen29")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Plant_Pathogen$NH4CK)) ## n = 29
p13<- ggplot(Plant_Pathogen, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnPlant_Pathogen29")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Plant_Pathogen$APCK)) ## n = 57
p14 <- ggplot(Plant_Pathogen, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnPlant_Pathogen57")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Plant_Pathogen$AKCK)) ## n = 43
p15 <- ggplot(Plant_Pathogen, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnPlant_Pathogen43")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Plant_Pathogen$ANCK)) ## n = 42
p16 <- ggplot(Plant_Pathogen, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnPlant_Pathogen42")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Plant_Pathogen$RRpH)) ## n = 92
p17 <- ggplot(Plant_Pathogen, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="RRpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRpH" , y="RR92")
p17
pdf("RRpH.pdf",width=8,height=8)
p17
dev.off() 

## RRSOC
sum(!is.na(Plant_Pathogen$RRSOC)) ## n = 102
p18 <- ggplot(Plant_Pathogen, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="RRSOC")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRSOC" , y="RR102")
p18
pdf("RRSOC.pdf",width=8,height=8)
p18
dev.off() 

## RRTN
sum(!is.na(Plant_Pathogen$RRTN)) ## n = 64
p19 <- ggplot(Plant_Pathogen, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="RRTN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRTN" , y="RR64")
p19
pdf("RRTN.pdf",width=8,height=8)
p19
dev.off() 

## RRNO3
sum(!is.na(Plant_Pathogen$RRNO3)) ## n = 29
p20 <- ggplot(Plant_Pathogen, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="RRNO3")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNO3" , y="RR29")
p20
pdf("RRNO3.pdf",width=8,height=8)
p20
dev.off() 

## RRNH4
sum(!is.na(Plant_Pathogen$RRNH4)) ## n = 29
p21 <- ggplot(Plant_Pathogen, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="RRNH4")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNH4" , y="RR29")
p21
pdf("RRNH4.pdf",width=8,height=8)
p21
dev.off() 

## RRAP
sum(!is.na(Plant_Pathogen$RRAP)) ## n = 57
p22 <- ggplot(Plant_Pathogen, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="RRAP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAP" , y="RR57")
p22
pdf("RRAP.pdf",width=8,height=8)
p22
dev.off() 

## RRAK
sum(!is.na(Plant_Pathogen$RRAK)) ## n = 43
p23 <- ggplot(Plant_Pathogen, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="RRAK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAK" , y="RR43")
p23
pdf("RRAK.pdf",width=8,height=8)
p23
dev.off() 

## RRAN
sum(!is.na(Plant_Pathogen$RRAN)) ## n = 42
p24 <- ggplot(Plant_Pathogen, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnPlant_Pathogen", x="RRAN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAN" , y="RR42")
p24
pdf("RRAN.pdf",width=8,height=8)
p24
dev.off() 

## RRYield
sum(!is.na(Plant_Pathogen$RRYield)) ## n = 50
p25 <- ggplot(Plant_Pathogen, aes(x=RR, y=RRYield)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm, color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") +
  theme_bw() +
  theme(text = element_text(family = "serif", size=20)) +
  geom_vline(aes(xintercept=0), colour="black", linewidth=0.5, linetype="dashed") +
  labs(x="RR", y="RRYield50") +
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
  scale_x_continuous(limits=c(-2.22, 1.25), expand=c(0, 0)) + 
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



