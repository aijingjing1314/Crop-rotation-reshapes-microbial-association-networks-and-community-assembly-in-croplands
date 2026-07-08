
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Nitrification <- read.csv("Nitrification.csv", fileEncoding = "latin1")
# Check data
head(Nitrification)

# 1. The number of Obversation
total_number <- nrow(Nitrification)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 274  

# 2. The number of Study
unique_studyid_number <- length(unique(Nitrification$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 






#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Nitrification_filteredLongitude_Sub <- subset(Nitrification, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Nitrification_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Nitrification_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredLongitude_Sub %>%
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

overall_model_Nitrification_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Nitrification_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -277.4716   554.9432   560.9432   571.7606   561.0328   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1242  0.3524     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1613.1687, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.9473, p-val = 0.6227
# 
# Model Results:
# 
#                            estimate      se    zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0    0.0497  0.0604  0.8224  0.4109  -0.0687  0.1680    
# Longitude_SubLongitudeXy0    0.0617  0.1185  0.5206  0.6027  -0.1706  0.2940       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Nitrification_filteredLongitude_Sub)
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
Nitrification_filteredMAPmean2_Sub <- subset(Nitrification, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Nitrification_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Nitrification_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredMAPmean2_Sub %>%
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

overall_model_Nitrification_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Nitrification_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -254.4900   508.9801   516.9801   531.3886   517.1305   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.2096  0.4578     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1515.0979, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 49.2604, p-val < .0001
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb    ci.ub     
# MAPmean2_SubMAP600Dao1200   -0.2884  0.0877  -3.2890  0.0010  -0.4602  -0.1165  ** 
# MAPmean2_SubMAPDy1200        0.2960  0.1540   1.9224  0.0546  -0.0058   0.5979   . 
# MAPmean2_SubMAPXy600         0.2534  0.0849   2.9862  0.0028   0.0871   0.4197  **    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Nitrification_filteredMAPmean2_Sub)
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
#             "a"                       "b"                       "b" 




### 8.3 MATmean_Sub
Nitrification_filteredMATmean_Sub <- subset(Nitrification, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Nitrification_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Nitrification_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredMATmean_Sub %>%
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

overall_model_Nitrification_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Nitrification_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -250.9071   501.8142   509.8142   524.2227   509.9646   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1823  0.4270     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1540.9993, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 55.5649, p-val < .0001
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb    ci.ub      
# MATmean_SubMAT8Dao15   -0.4213  0.0959  -4.3949  <.0001  -0.6091  -0.2334  *** 
# MATmean_SubMATDy15      0.2158  0.1145   1.8846  0.0595  -0.0086   0.4402    . 
# MATmean_SubMATXy8       0.1709  0.0815   2.0971  0.0360   0.0112   0.3306    * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Nitrification_filteredMATmean_Sub)
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
#         "a"                  "b"                  "b" 




### 8.4 LegumeNonlegume
Nitrification_filteredLegumeNonlegume <- subset(Nitrification, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Nitrification_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Nitrification_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredLegumeNonlegume %>%
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

overall_model_Nitrification_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Nitrification_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -216.4654   432.9308   440.9308   455.3393   441.0812   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1298  0.3602     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1526.0803, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 127.5985, p-val < .0001
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb   ci.ub     
# LegumeNonlegumeLegume to Non-legume        0.0635  0.0658   0.9651  0.3345  -0.0655  0.1925     
# LegumeNonlegumeNon-legume to Legume       -0.1010  0.0568  -1.7776  0.0755  -0.2124  0.0104   . 
# LegumeNonlegumeNon-legume to Non-legume    0.1552  0.0570   2.7201  0.0065   0.0434  0.2670  ** 
 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Nitrification_filteredLegumeNonlegume)
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
    #        "a"                                     "b"                                     "a" 




### 8.5 AMnonAM
Nitrification_filteredAMnonAM <- subset(Nitrification, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Nitrification_filteredAMnonAM$AMnonAM <- droplevels(factor(Nitrification_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredAMnonAM %>%
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

overall_model_Nitrification_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Nitrification_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -252.8799   505.7598   513.7598   528.1683   513.9102   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1232  0.3509     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1474.2508, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 51.5672, p-val < .0001
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub      
# AMnonAMAM to AM      -0.0162  0.0572  -0.2824  0.7777  -0.1283  0.0960      
# AMnonAMAM to nonAM    0.3593  0.0704   5.1037  <.0001   0.2213  0.4972  *** 
# AMnonAMnonAM to AM    0.0785  0.1343   0.5850  0.5585  -0.1846  0.3417      
  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Nitrification_filteredAMnonAM)
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
   #        "a"                "b"               "ab"


### 8.6 C3C4
Nitrification_filteredC3C4 <- subset(Nitrification, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Nitrification_filteredC3C4$C3C4 <- droplevels(factor(Nitrification_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredC3C4 %>%
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

overall_model_Nitrification_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Nitrification_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 273; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -257.5426   515.0851   523.0851   537.4788   523.2361   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0915  0.3025     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 1525.3162, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 38.3629, p-val < .0001
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb    ci.ub    
# C3C4C3 to C3    0.1362  0.0568   2.3975  0.0165   0.0249   0.2476  * 
# C3C4C3 to C4    0.1135  0.0516   2.1996  0.0278   0.0124   0.2147  * 
# C3C4C4 to C3   -0.1267  0.0572  -2.2125  0.0269  -0.2388  -0.0145  * 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredC3C4)
vcov_rotation <- vcov(overall_model_Nitrification_filteredC3C4)
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
#   "a"          "a"          "b" 


### 8.7 Annual_Pere
Nitrification_filteredAnnual_Pere <- subset(Nitrification, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
Nitrification_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Nitrification_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredAnnual_Pere %>%
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

overall_model_Nitrification_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Nitrification_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -276.1774   552.3548   560.3548   574.7633   560.5052   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1111  0.3333     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 1613.0774, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 6.0785, p-val = 0.1079
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual       0.0261  0.0532   0.4904  0.6239  -0.0782  0.1305    
# Annual_PereAnnual to Perennial   -0.0006  0.0716  -0.0090  0.9928  -0.1410  0.1397    
# Annual_PerePerennial to Annual    0.3355  0.1431   2.3442  0.0191   0.0550  0.6159  * 

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Nitrification_filteredAnnual_Pere)
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
   #              "a"                            "a"                            "b" 



### 8.8 Plant_stage
Nitrification_filteredPlant_stage <- subset(Nitrification, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
Nitrification_filteredPlant_stage$Plant_stage <- droplevels(factor(Nitrification_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredPlant_stage %>%
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

overall_model_Nitrification_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Nitrification_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredPlant_stage)
# Multivariate Meta-Analysis Model (k = 239; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -225.2993   450.5986   458.5986   472.4539   458.7718   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1105  0.3324     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 236) = 1344.7225, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 30.7270, p-val < .0001
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb    ci.ub     
# Plant_stageHarvest               0.1298  0.0618   2.0992  0.0358   0.0086   0.2510   * 
# Plant_stageReproductive stage   -0.2018  0.0751  -2.6859  0.0072  -0.3491  -0.0545  ** 
# Plant_stageVegetative stage      0.0505  0.0629   0.8039  0.4215  -0.0727   0.1738     



# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Nitrification_filteredPlant_stage)
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
# "       "a"                           "b"                           "c" 

### 8.9 Rotation_cycles2
Nitrification_filteredRotation_cycles2 <- subset(Nitrification, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
Nitrification_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Nitrification_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredRotation_cycles2 %>%
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

overall_model_Nitrification_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Nitrification_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -245.2581   490.5163   502.5163   524.0846   502.8369   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1153  0.3395     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 269) = 1524.9415, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 70.7918, p-val < .0001
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub     
# Rotation_cycles2D1       0.1907  0.0706   2.7013  0.0069   0.0523  0.3291  ** 
# Rotation_cycles2D1-3     0.0161  0.0712   0.2266  0.8208  -0.1233  0.1556     
# Rotation_cycles2D10     -0.0570  0.0915  -0.6227  0.5335  -0.2364  0.1224     
# Rotation_cycles2D3-5    -0.1093  0.0836  -1.3078  0.1910  -0.2731  0.0545     
# Rotation_cycles2D5-10    0.0220  0.0801   0.2752  0.7832  -0.1349  0.1790    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Nitrification_filteredRotation_cycles2)
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
#            "a"                  "bc"                  "bc"                   "b"                  "ac"      


### 8.10 Duration2
Nitrification_filteredDuration2 <- subset(Nitrification, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Nitrification_filteredDuration2$Duration2 <- droplevels(factor(Nitrification_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredDuration2 %>%
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

overall_model_Nitrification_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Nitrification_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -252.0598   504.1196   514.1196   532.1117   514.3468   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1265  0.3557     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 1543.6873, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 53.0349, p-val < .0001
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb    ci.ub     
# Duration2D1-5      0.0360  0.0766   0.4697  0.6386  -0.1142   0.1862     
# Duration2D11-20   -0.1915  0.0852  -2.2467  0.0247  -0.3586  -0.0244   * 
# Duration2D20-40    0.2392  0.0901   2.6548  0.0079   0.0626   0.4159  ** 
# Duration2D6-10     0.1781  0.0870   2.0475  0.0406   0.0076   0.3485   * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredDuration2)
vcov_rotation <- vcov(overall_model_Nitrification_filteredDuration2)
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
#    "a"             "b"             "a"             "a" 


### 8.11 SpeciesRichness2
Nitrification_filteredSpeciesRichness2 <- subset(Nitrification, SpeciesRichness2 %in% c("R2", "R3"))
#
Nitrification_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Nitrification_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredSpeciesRichness2 %>%
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

overall_model_Nitrification_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Nitrification_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredSpeciesRichness2)
# Multivariate Meta-Analysis Model (k = 265; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -274.5503   549.1007   555.1007   565.8172   555.1934   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1306  0.3613     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 263) = 1583.9324, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 3.5699, p-val = 0.1678
# 
# Model Results:
# 
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2    0.0431  0.0570  0.7558  0.4498  -0.0686  0.1548    
# SpeciesRichness2R3    0.1212  0.0708  1.7125  0.0868  -0.0175  0.2599  . 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Nitrification_filteredSpeciesRichness2)
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
Nitrification_filteredBulk_Rhizosphere <- subset(Nitrification, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Nitrification_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Nitrification_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredBulk_Rhizosphere %>%
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

overall_model_Nitrification_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Nitrification_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -261.2556   522.5113   528.5113   539.3287   528.6008   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1305  0.3612     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 1604.8941, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 36.5721, p-val < .0001
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb   ci.ub     
# Bulk_RhizosphereNon_Rhizosphere   -0.0178  0.0563  -0.3158  0.7521  -0.1281  0.0925     
# Bulk_RhizosphereRhizosphere        0.1752  0.0587   2.9838  0.0028   0.0601  0.2902  ** 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Nitrification_filteredBulk_Rhizosphere)
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
Nitrification_filteredSoil_texture <- subset(Nitrification, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Nitrification_filteredSoil_texture$Soil_texture <- droplevels(factor(Nitrification_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredSoil_texture %>%
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

overall_model_Nitrification_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Nitrification_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 224; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -242.9862   485.9725   493.9725   507.5651   494.1577   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1229  0.3506     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 221) = 1256.4788, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.5469, p-val = 0.6715
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse    0.0310  0.1464   0.2117  0.8323  -0.2560  0.3180    
# Soil_textureFine     -0.0769  0.1398  -0.5503  0.5821  -0.3510  0.1971    
# Soil_textureMedium    0.0976  0.0891   1.0951  0.2735  -0.0771  0.2723       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Nitrification_filteredSoil_texture)
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
Nitrification_filteredTillage <- subset(Nitrification, Tillage %in% c("Tillage", "No_tillage"))
#
Nitrification_filteredTillage$Tillage <- droplevels(factor(Nitrification_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredTillage %>%
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

overall_model_Nitrification_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Nitrification_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredTillage)
# Multivariate Meta-Analysis Model (k = 150; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -160.0894   320.1788   326.1788   335.1704   326.3454   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0817  0.2858     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 148) = 740.5980, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 3.1678, p-val = 0.2052
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0374  0.0864  -0.4335  0.6647  -0.2068  0.1319    
# TillageTillage       0.0433  0.0854   0.5070  0.6122  -0.1241  0.2107     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredTillage)
vcov_rotation <- vcov(overall_model_Nitrification_filteredTillage)
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
Nitrification_filteredStraw_retention <- subset(Nitrification, Straw_retention %in% c("Retention", "No_retention"))
#
Nitrification_filteredStraw_retention$Straw_retention <- droplevels(factor(Nitrification_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredStraw_retention %>%
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

overall_model_Nitrification_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Nitrification_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 55; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -42.8726   85.7452   91.7452   97.6561   92.2350   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1308  0.3617     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 456.5720, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 5.5120, p-val = 0.0635
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention   -0.0900  0.1267  -0.7103  0.4775  -0.3382  0.1583    
# Straw_retentionRetention       0.1000  0.1198   0.8348  0.4038  -0.1348  0.3348      
  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Nitrification_filteredStraw_retention)
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
Nitrification_filteredPrimer <- subset(Nitrification, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
Nitrification_filteredPrimer$Primer <- droplevels(factor(Nitrification_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Nitrification_filteredPrimer %>%
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

overall_model_Nitrification_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Nitrification_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Nitrification_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 264; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -244.8668   489.7335   497.7335   511.9916   497.8898   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1325  0.3640     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 261) = 1359.3744, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.4321, p-val = 0.2184
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4    0.0303  0.0739   0.4093  0.6823  -0.1146  0.1752    
# PrimerV4       0.2804  0.1360   2.0611  0.0393   0.0138  0.5470  * 
# PrimerV4-V5   -0.0167  0.1298  -0.1286  0.8977  -0.2712  0.2378    
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrification_filteredPrimer)
vcov_rotation <- vcov(overall_model_Nitrification_filteredPrimer)
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
Nitrification$Wr <- 1 / Nitrification$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Nitrification)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Nitrification)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Nitrification)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Nitrification)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Nitrification)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Nitrification)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Nitrification)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Nitrification)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Nitrification)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Nitrification)

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
# Model1   Model1 295.1695 316.8483 -141.5848
# Model2   Model2 292.2347 313.9135 -140.1174
# Model3   Model3 286.1946 307.8734 -137.0973######################################
# Model4   Model4 295.1480 316.8268 -141.5740
# Model5   Model5 293.6548 315.3335 -140.8274
# Model6   Model6 293.5025 315.1812 -140.7512
# Model7   Model7 292.3332 314.0120 -140.1666
# Model8   Model8 286.2044 307.8832 -137.1022
# Model9   Model9 308.9290 345.0603 -144.4645
# Model10 Model10 300.9175 337.0487 -140.4587

##### Model 3 is the best model
summary(Model3)
# Number of obs: 274, groups:  StudyID, 50
anova(Model3)
# Type III Analysis of Variance Table with Satterthwaite's method
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)    
# scale(log(Rotation_cycles)) 42.589  42.589     1 112.61 12.4817 0.000597 ***
# scale(SpeciesRichness)       0.488   0.488     1 170.68  0.1430 0.705779    
# scale(Duration)             23.188  23.188     1 111.36  6.7958 0.010388 *  


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Nitrification)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.1866  0.1866     1 34.114  0.0400 0.8426
# scale(SpeciesRichness)      4.6233  4.6233     1 70.645  0.9916 0.3228
# scale(Duration)             0.9089  0.9089     1 45.680  0.1949 0.6609
# scale(pHCK)                 4.5394  4.5394     1 26.097  0.9736 0.3329

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Nitrification)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.26744 0.26744     1 33.421  0.0693 0.7940
# scale(SpeciesRichness)      0.01100 0.01100     1 62.852  0.0028 0.9576
# scale(Duration)             1.71939 1.71939     1 58.342  0.4456 0.5071
# scale(SOCCK)                0.86984 0.86984     1 22.320  0.2254 0.6396

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Nitrification)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles))  0.0123  0.0123     1  9.907  0.0026 0.9606
# scale(SpeciesRichness)       0.6403  0.6403     1 34.918  0.1333 0.7172
# scale(Duration)              0.0132  0.0132     1 18.876  0.0028 0.9587
# scale(TNCK)                 10.5614 10.5614     1  8.990  2.1987 0.1723


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Nitrification)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 2.92313 2.92313     1 40.279  4.0162 0.05182 .
# scale(SpeciesRichness)      0.58261 0.58261     1 32.009  0.8005 0.37763  
# scale(Duration)             0.94772 0.94772     1 40.769  1.3021 0.26049  
# scale(NO3CK)                0.01774 0.01774     1 18.757  0.0244 0.87761  

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Nitrification)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 3.8259  3.8259     1 39.424  5.3942 0.02545 *
# scale(SpeciesRichness)      0.4177  0.4177     1 31.858  0.5889 0.44850  
# scale(Duration)             1.3854  1.3854     1 39.169  1.9534 0.17009  
# scale(NH4CK)                0.2305  0.2305     1 11.330  0.3250 0.57973 

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = Nitrification)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.8156  0.8156     1 30.359  0.1850 0.6701
# scale(SpeciesRichness)      2.0451  2.0451     1 58.721  0.4639 0.4985
# scale(Duration)             3.4170  3.4170     1 41.933  0.7751 0.3836
# scale(APCK)                 0.2697  0.2697     1 47.594  0.0612 0.8057

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Nitrification)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.0074  0.0074     1 13.705  0.0016 0.9691
# scale(SpeciesRichness)      2.1236  2.1236     1 31.971  0.4436 0.5102
# scale(Duration)             7.5063  7.5063     1 17.657  1.5679 0.2268
# scale(AKCK)                 6.0448  6.0448     1  8.725  1.2627 0.2911
# 
# #### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Nitrification)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 8.0614  8.0614     1 20.3583  1.1243 0.3014
# scale(SpeciesRichness)      2.0942  2.0942     1 24.0492  0.2921 0.5939
# scale(Duration)             6.4525  6.4525     1 21.0920  0.8999 0.3535
# scale(ANCK)                 0.0094  0.0094     1  5.8285  0.0013 0.9723  

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Nitrification)
summary(ModelLatitude)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# scale(log(Rotation_cycles)) 42.663  42.663     1 114.370 12.5167 0.0005839 ***
# scale(SpeciesRichness)       0.443   0.443     1 182.871  0.1300 0.7188091    
# scale(Duration)             22.640  22.640     1 127.670  6.6424 0.0110939 *  
# scale(Latitude)              0.001   0.001     1  30.886  0.0002 0.9891034  

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Nitrification)
summary(ModelLongitude)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# scale(log(Rotation_cycles)) 42.712  42.712     1 111.590 12.5428 0.0005816 ***
# scale(SpeciesRichness)       0.353   0.353     1 167.847  0.1036 0.7479262    
# scale(Duration)             17.491  17.491     1 126.219  5.1365 0.0251306 *  
# scale(Longitude)             0.553   0.553     1  27.371  0.1623 0.6901839   

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Nitrification)
summary(ModelMAPmean)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# scale(log(Rotation_cycles)) 40.457  40.457     1 103.272 11.8634 0.0008285 ***
# scale(SpeciesRichness)       0.545   0.545     1 162.511  0.1597 0.6899184    
# scale(Duration)             21.620  21.620     1  98.672  6.3398 0.0134148 *  
# scale(MAPmean)               1.625   1.625     1  34.949  0.4764 0.4946025    

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Nitrification)
summary(ModelMATmean)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# scale(log(Rotation_cycles)) 40.618  40.618     1 100.138 11.8549 0.0008409 ***
# scale(SpeciesRichness)       0.517   0.517     1 162.382  0.1510 0.6980735    
# scale(Duration)             22.401  22.401     1  99.723  6.5381 0.0120659 *  
# scale(MATmean)               1.923   1.923     1  28.136  0.5611 0.4600213  


############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Nitrification$SpeciesRichness)) ## n = 274
p1 <- ggplot(Nitrification, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnNitrification274")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Nitrification$Duration)) ## n = 274
p2 <- ggplot(Nitrification, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnNitrification274")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Nitrification$Rotation_cycles)) ## n = 274
p3 <- ggplot(Nitrification, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnNitrification274")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Nitrification$Latitude)) ## n = 274
p5 <- ggplot(Nitrification, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnNitrification274")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Nitrification$Longitude)) ## n = 274
p6 <- ggplot(Nitrification, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnNitrification274")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Nitrification$MAPmean)) ## n = 274
p7 <- ggplot(Nitrification, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnNitrification274")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Nitrification$MATmean)) ## n = 274
p8 <- ggplot(Nitrification, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnNitrification274")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Nitrification$pHCK)) ## n = 86
p9 <- ggplot(Nitrification, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnNitrification86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Nitrification$SOCCK)) ## n = 83
p10 <- ggplot(Nitrification, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnNitrification83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Nitrification$TNCK)) ## n = 52
p11 <- ggplot(Nitrification, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnNitrification52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Nitrification$NO3CK)) ## n = 46
p12 <- ggplot(Nitrification, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnNitrification46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Nitrification$NH4CK)) ## n = 45
p13<- ggplot(Nitrification, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnNitrification45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Nitrification$APCK)) ## n = 64
p14 <- ggplot(Nitrification, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnNitrification64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Nitrification$AKCK)) ## n = 39
p15 <- ggplot(Nitrification, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnNitrification39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Nitrification$ANCK)) ## n = 31
p16 <- ggplot(Nitrification, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnNitrification31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Nitrification$RRpH)) ## n = 86
p17 <- ggplot(Nitrification, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="RRpH")+
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
sum(!is.na(Nitrification$RRSOC)) ## n = 83
p18 <- ggplot(Nitrification, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="RRSOC")+
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
sum(!is.na(Nitrification$RRTN)) ## n = 52
p19 <- ggplot(Nitrification, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="RRTN")+
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
sum(!is.na(Nitrification$RRNO3)) ## n = 44
p20 <- ggplot(Nitrification, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="RRNO3")+
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
sum(!is.na(Nitrification$RRNH4)) ## n = 45
p21 <- ggplot(Nitrification, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="RRNH4")+
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
sum(!is.na(Nitrification$RRAP)) ## n = 64
p22 <- ggplot(Nitrification, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="RRAP")+
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
sum(!is.na(Nitrification$RRAK)) ## n = 39
p23 <- ggplot(Nitrification, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="RRAK")+
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
sum(!is.na(Nitrification$RRAN)) ## n = 31
p24 <- ggplot(Nitrification, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrification", x="RRAN")+
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
sum(!is.na(Nitrification$RRYield)) ## n = 49
p25 <- ggplot(Nitrification, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-1.32, 1), expand=c(0, 0)) + 
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




