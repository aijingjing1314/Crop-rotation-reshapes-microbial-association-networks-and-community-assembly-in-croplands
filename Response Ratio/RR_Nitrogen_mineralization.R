
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Nitrogen_mineralization <- read.csv("Nitrogen_mineralization.csv", fileEncoding = "latin1")
# Check data
head(Nitrogen_mineralization)

# 1. The number of Obversation
total_number <- nrow(Nitrogen_mineralization)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 274  

# 2. The number of Study
unique_studyid_number <- length(unique(Nitrogen_mineralization$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 50 






#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Nitrogen_mineralization_filteredLongitude_Sub <- subset(Nitrogen_mineralization, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Nitrogen_mineralization_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Nitrogen_mineralization_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredLongitude_Sub %>%
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

overall_model_Nitrogen_mineralization_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredLongitude_Sub)
# ultivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  265.5735  -531.1470  -525.1470  -514.3296  -525.0575   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0036  0.0604     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 960.6213, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.5745, p-val = 0.2760
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0   -0.0169  0.0105  -1.6043  0.1087  -0.0375  0.0037    
# Longitude_SubLongitudeXy0   -0.0006  0.0204  -0.0276  0.9780  -0.0406  0.0395       

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredLongitude_Sub)
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
Nitrogen_mineralization_filteredMAPmean2_Sub <- subset(Nitrogen_mineralization, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Nitrogen_mineralization_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Nitrogen_mineralization_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredMAPmean2_Sub %>%
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

overall_model_Nitrogen_mineralization_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  265.5668  -531.1337  -523.1337  -508.7252  -522.9833   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0034  0.0585     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 949.0853, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 6.0776, p-val = 0.1079
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean2_SubMAP600Dao1200    0.0003  0.0114   0.0292  0.9767  -0.0220  0.0227    
# MAPmean2_SubMAPDy1200       -0.0385  0.0224  -1.7203  0.0854  -0.0824  0.0054  . 
# MAPmean2_SubMAPXy600        -0.0161  0.0112  -1.4475  0.1478  -0.0380  0.0057     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredMAPmean2_Sub)
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
#              "a"                       "a"                       "a" 




### 8.3 MATmean_Sub
Nitrogen_mineralization_filteredMATmean_Sub <- subset(Nitrogen_mineralization, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Nitrogen_mineralization_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Nitrogen_mineralization_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredMATmean_Sub %>%
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

overall_model_Nitrogen_mineralization_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  264.1957  -528.3915  -520.3915  -505.9830  -520.2411   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0035  0.0593     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 931.2230, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.7929, p-val = 0.2847
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubMAT8Dao15   -0.0153  0.0132  -1.1569  0.2473  -0.0411  0.0106    
# MATmean_SubMATDy15      0.0031  0.0172   0.1820  0.8556  -0.0306  0.0368    
# MATmean_SubMATXy8      -0.0222  0.0115  -1.9388  0.0525  -0.0447  0.0002  . 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredMATmean_Sub)
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
#               "a"                  "a"                  "a" 




### 8.4 LegumeNonlegume
Nitrogen_mineralization_filteredLegumeNonlegume <- subset(Nitrogen_mineralization, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Nitrogen_mineralization_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Nitrogen_mineralization_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredLegumeNonlegume %>%
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

overall_model_Nitrogen_mineralization_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  261.4761  -522.9522  -514.9522  -500.5437  -514.8018   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0036  0.0596     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 930.2633, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.2773, p-val = 0.5169
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume       -0.0144  0.0105  -1.3674  0.1715  -0.0349  0.0062    
# LegumeNonlegumeNon-legume to Legume       -0.0123  0.0097  -1.2665  0.2053  -0.0314  0.0067    
# LegumeNonlegumeNon-legume to Non-legume   -0.0138  0.0100  -1.3757  0.1689  -0.0334  0.0059    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredLegumeNonlegume)
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
    #        "a"                                     "a"                                     "a" 



### 8.5 AMnonAM
Nitrogen_mineralization_filteredAMnonAM <- subset(Nitrogen_mineralization, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Nitrogen_mineralization_filteredAMnonAM$AMnonAM <- droplevels(factor(Nitrogen_mineralization_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredAMnonAM %>%
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

overall_model_Nitrogen_mineralization_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  264.5122  -529.0243  -521.0243  -506.6159  -520.8740   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0035  0.0594     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 964.0369, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.9250, p-val = 0.1153
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb    ci.ub    
# AMnonAMAM to AM      -0.0146  0.0098  -1.4926  0.1355  -0.0338   0.0046    
# AMnonAMAM to nonAM   -0.0252  0.0120  -2.0957  0.0361  -0.0488  -0.0016  * 
# AMnonAMnonAM to AM    0.0072  0.0183   0.3946  0.6931  -0.0287   0.0432    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredAMnonAM)
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
Nitrogen_mineralization_filteredC3C4 <- subset(Nitrogen_mineralization, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Nitrogen_mineralization_filteredC3C4$C3C4 <- droplevels(factor(Nitrogen_mineralization_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredC3C4 %>%
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

overall_model_Nitrogen_mineralization_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 273; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  262.7805  -525.5610  -517.5610  -503.1673  -517.4100   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0035  0.0588     49     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 936.2026, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.4085, p-val = 0.2206
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb   ci.ub    
# C3C4C3 to C3   -0.0160  0.0110  -1.4556  0.1455  -0.0374  0.0055    
# C3C4C3 to C4   -0.0130  0.0100  -1.3029  0.1926  -0.0325  0.0065    
# C3C4C4 to C3   -0.0048  0.0102  -0.4659  0.6413  -0.0247  0.0152  


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredC3C4)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredC3C4)
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
#    "a"          "a"          "a"


### 8.7 Annual_Pere
Nitrogen_mineralization_filteredAnnual_Pere <- subset(Nitrogen_mineralization, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
Nitrogen_mineralization_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Nitrogen_mineralization_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredAnnual_Pere %>%
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

overall_model_Nitrogen_mineralization_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  262.8515  -525.7029  -517.7029  -503.2944  -517.5525   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0036  0.0597     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 271) = 896.5755, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.9733, p-val = 0.3958
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual      -0.0120  0.0097  -1.2335  0.2174  -0.0310  0.0070    
# Annual_PereAnnual to Perennial   -0.0238  0.0151  -1.5769  0.1148  -0.0535  0.0058    
# Annual_PerePerennial to Annual   -0.0128  0.0210  -0.6098  0.5420  -0.0539  0.0283  

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredAnnual_Pere)
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
Nitrogen_mineralization_filteredPlant_stage <- subset(Nitrogen_mineralization, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
Nitrogen_mineralization_filteredPlant_stage$Plant_stage <- droplevels(factor(Nitrogen_mineralization_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredPlant_stage %>%
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

overall_model_Nitrogen_mineralization_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredPlant_stage)
# Multivariate Meta-Analysis Model (k = 239; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  223.8802  -447.7603  -439.7603  -425.9050  -439.5871   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0047  0.0689     37     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 236) = 821.1406, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.8674, p-val = 0.6004
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb   ci.ub    
# Plant_stageHarvest              -0.0159  0.0125  -1.2767  0.2017  -0.0403  0.0085    
# Plant_stageReproductive stage   -0.0098  0.0147  -0.6646  0.5063  -0.0386  0.0190    
# Plant_stageVegetative stage     -0.0170  0.0127  -1.3364  0.1814  -0.0420  0.0079    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredPlant_stage)
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
Nitrogen_mineralization_filteredRotation_cycles2 <- subset(Nitrogen_mineralization, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
Nitrogen_mineralization_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Nitrogen_mineralization_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredRotation_cycles2 %>%
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

overall_model_Nitrogen_mineralization_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  258.0618  -516.1237  -504.1237  -482.5554  -503.8031   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0037  0.0608     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 269) = 951.6400, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 2.1782, p-val = 0.8240
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1      -0.0156  0.0128  -1.2206  0.2222  -0.0407  0.0095    
# Rotation_cycles2D1-3    -0.0151  0.0129  -1.1690  0.2424  -0.0403  0.0102    
# Rotation_cycles2D10     -0.0124  0.0204  -0.6066  0.5441  -0.0524  0.0276    
# Rotation_cycles2D3-5    -0.0118  0.0153  -0.7739  0.4390  -0.0417  0.0181    
# Rotation_cycles2D5-10   -0.0103  0.0148  -0.6980  0.4852  -0.0393  0.0187     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredRotation_cycles2)
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
#            "a"                   "a"                   "a"                   "a"                   "a" 


### 8.10 Duration2
Nitrogen_mineralization_filteredDuration2 <- subset(Nitrogen_mineralization, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Nitrogen_mineralization_filteredDuration2$Duration2 <- droplevels(factor(Nitrogen_mineralization_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredDuration2 %>%
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

overall_model_Nitrogen_mineralization_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  262.4373  -524.8747  -514.8747  -496.8826  -514.6474   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0036  0.0601     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 270) = 975.1423, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 3.8518, p-val = 0.4264
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5     -0.0244  0.0130  -1.8677  0.0618  -0.0499  0.0012  . 
# Duration2D11-20    0.0014  0.0151   0.0945  0.9247  -0.0282  0.0310    
# Duration2D20-40   -0.0045  0.0164  -0.2769  0.7819  -0.0366  0.0276    
# Duration2D6-10    -0.0042  0.0172  -0.2418  0.8089  -0.0380  0.0296 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredDuration2)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredDuration2)
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
#     "a"             "a"             "a"             "a" 


### 8.11 SpeciesRichness2
Nitrogen_mineralization_filteredSpeciesRichness2 <- subset(Nitrogen_mineralization, SpeciesRichness2 %in% c("R2", "R3"))
#
Nitrogen_mineralization_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Nitrogen_mineralization_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredSpeciesRichness2 %>%
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

overall_model_Nitrogen_mineralization_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredSpeciesRichness2)
# Multivariate Meta-Analysis Model (k = 265; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  253.6221  -507.2441  -501.2441  -490.5277  -501.1515   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0038  0.0613     47     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 263) = 964.8526, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 4.4135, p-val = 0.1101
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb    ci.ub    
# SpeciesRichness2R2   -0.0149  0.0098  -1.5156  0.1296  -0.0341   0.0044    
# SpeciesRichness2R3   -0.0256  0.0122  -2.0991  0.0358  -0.0494  -0.0017  *     



# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredSpeciesRichness2)
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
Nitrogen_mineralization_filteredBulk_Rhizosphere <- subset(Nitrogen_mineralization, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Nitrogen_mineralization_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Nitrogen_mineralization_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredBulk_Rhizosphere %>%
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

overall_model_Nitrogen_mineralization_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 274; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  265.8992  -531.7984  -525.7984  -514.9810  -525.7088   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0037  0.0610     50     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 272) = 974.0560, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 5.9109, p-val = 0.0521
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb    ci.ub    
# Bulk_RhizosphereNon_Rhizosphere   -0.0094  0.0097  -0.9665  0.3338  -0.0283   0.0096    
# Bulk_RhizosphereRhizosphere       -0.0231  0.0106  -2.1724  0.0298  -0.0439  -0.0023  * 
 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredBulk_Rhizosphere)
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
Nitrogen_mineralization_filteredSoil_texture <- subset(Nitrogen_mineralization, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Nitrogen_mineralization_filteredSoil_texture$Soil_texture <- droplevels(factor(Nitrogen_mineralization_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredSoil_texture %>%
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

overall_model_Nitrogen_mineralization_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 224; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  197.1863  -394.3725  -386.3725  -372.7799  -386.1874   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0058  0.0761     32     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 221) = 863.9732, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.2476, p-val = 0.5226
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse   -0.0055  0.0293  -0.1868  0.8518  -0.0629  0.0520    
# Soil_textureFine     -0.0363  0.0314  -1.1564  0.2475  -0.0979  0.0252    
# Soil_textureMedium   -0.0182  0.0194  -0.9356  0.3495  -0.0562  0.0199     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredSoil_texture)
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
Nitrogen_mineralization_filteredTillage <- subset(Nitrogen_mineralization, Tillage %in% c("Tillage", "No_tillage"))
#
Nitrogen_mineralization_filteredTillage$Tillage <- droplevels(factor(Nitrogen_mineralization_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredTillage %>%
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

overall_model_Nitrogen_mineralization_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredTillage)
# Multivariate Meta-Analysis Model (k = 150; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  182.4261  -364.8522  -358.8522  -349.8606  -358.6856   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0079  0.0890     13     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 148) = 349.1826, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.7584, p-val = 0.6844
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0206  0.0258  -0.7979  0.4249  -0.0712  0.0300    
# TillageTillage      -0.0222  0.0257  -0.8649  0.3871  -0.0725  0.0281      

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredTillage)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredTillage)
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
Nitrogen_mineralization_filteredStraw_retention <- subset(Nitrogen_mineralization, Straw_retention %in% c("Retention", "No_retention"))
#
Nitrogen_mineralization_filteredStraw_retention$Straw_retention <- droplevels(factor(Nitrogen_mineralization_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredStraw_retention %>%
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

overall_model_Nitrogen_mineralization_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 55; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  36.3084  -72.6168  -66.6168  -60.7059  -66.1270   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0037  0.0605     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 279.7669, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 3.7920, p-val = 0.1502
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention    0.0166  0.0209   0.7928  0.4279  -0.0244  0.0576    
# Straw_retentionRetention      -0.0113  0.0196  -0.5781  0.5632  -0.0497  0.0271    
  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredStraw_retention)
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
Nitrogen_mineralization_filteredPrimer <- subset(Nitrogen_mineralization, Primer %in% c("V3-V4", "V4", "V4-V5"))
#
Nitrogen_mineralization_filteredPrimer$Primer <- droplevels(factor(Nitrogen_mineralization_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Nitrogen_mineralization_filteredPrimer %>%
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

overall_model_Nitrogen_mineralization_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Nitrogen_mineralization_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Nitrogen_mineralization_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 264; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  286.3986  -572.7972  -564.7972  -550.5391  -564.6410   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0043  0.0658     46     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 261) = 863.9494, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.9587, p-val = 0.2660
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4   -0.0063  0.0132  -0.4822  0.6297  -0.0321  0.0194    
# PrimerV4      -0.0460  0.0260  -1.7725  0.0763  -0.0969  0.0049  . 
# PrimerV4-V5   -0.0188  0.0246  -0.7645  0.4446  -0.0670  0.0294    

#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Nitrogen_mineralization_filteredPrimer)
vcov_rotation <- vcov(overall_model_Nitrogen_mineralization_filteredPrimer)
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
Nitrogen_mineralization$Wr <- 1 / Nitrogen_mineralization$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)

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
# Model1   Model1 -631.2236 -609.5448 321.6118##############################
# Model2   Model2 -631.0593 -609.3806 321.5297
# Model3   Model3 -630.6369 -608.9581 321.3184
# Model4   Model4 -631.1953 -609.5166 321.5977
# Model5   Model5 -631.1546 -609.4758 321.5773
# Model6   Model6 -631.1070 -609.4282 321.5535
# Model7   Model7 -631.0360 -609.3572 321.5180
# Model8   Model8 -630.6523 -608.9735 321.3261
# Model9   Model9 -590.1608 -554.0295 305.0804
# Model10 Model10 -592.6474 -556.5162 306.3237

##### Model 1 is the best model
summary(Model1)
# Number of obs: 274, groups:  StudyID, 50
anova(Model1)
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles) 2.13205 2.13205     1  15.511  0.7543 0.3984
# scale(SpeciesRichness) 0.00792 0.00792     1 115.014  0.0028 0.9579
# scale(Duration)        0.00623 0.00623     1  27.269  0.0022 0.9629


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
summary(ModelpH)
# Number of obs: 86, groups:  StudyID, 32
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles) 3.6240  3.6240     1  2.7134  0.7916 0.4454
# scale(SpeciesRichness) 0.2262  0.2262     1 31.0711  0.0494 0.8255
# scale(Duration)        0.0489  0.0489     1  4.3386  0.0107 0.9223
# scale(pHCK)            0.7379  0.7379     1  4.0491  0.1612 0.7084

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
summary(ModelSOC)
# Number of obs: 83, groups:  StudyID, 30
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 6.3501  6.3501     1 16.658  1.6570 0.2156
# scale(SpeciesRichness) 0.0510  0.0510     1 48.182  0.0133 0.9086
# scale(Duration)        1.5794  1.5794     1 45.219  0.4121 0.5241
# scale(SOCCK)           1.5613  1.5613     1  6.374  0.4074 0.5455

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
summary(ModelTN)
# Number of obs: 52, groups:  StudyID, 19
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 7.5198  7.5198     1 10.469  2.0449 0.1819
# scale(SpeciesRichness) 1.2800  1.2800     1 25.118  0.3481 0.5605
# scale(Duration)        3.4097  3.4097     1 22.775  0.9272 0.3457
# scale(TNCK)            0.5659  0.5659     1 38.198  0.1539 0.6970


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
summary(ModelNO3)
# Number of obs: 46, groups:  StudyID, 15
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.82625 0.82625     1 35.206  0.6900 0.4118
# scale(SpeciesRichness) 0.10930 0.10930     1 27.953  0.0913 0.7648
# scale(Duration)        2.74722 2.74722     1 40.621  2.2941 0.1376
# scale(NO3CK)           0.00079 0.00079     1 26.109  0.0007 0.9796

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
summary(ModelNH4)
# Number of obs: 45, groups:  StudyID, 14
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.41736 0.41736     1 39.868  0.3597 0.5521
# scale(SpeciesRichness) 0.13988 0.13988     1 28.857  0.1205 0.7310
# scale(Duration)        2.06036 2.06036     1 37.534  1.7755 0.1907
# scale(NH4CK)           0.00302 0.00302     1 10.872  0.0026 0.9602

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
summary(ModelAP)
# Number of obs: 64, groups:  StudyID, 26
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 30.9679 30.9679     1 27.295  5.0019 0.03368 *
# scale(SpeciesRichness)  1.7288  1.7288     1 35.134  0.2792 0.60052  
# scale(Duration)         8.3259  8.3259     1 32.691  1.3448 0.25459  
# scale(APCK)            10.0906 10.0906     1 27.580  1.6298 0.21237

# #### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 31.4647 31.4647     1 13.387  4.0842 0.06375 .
# scale(SpeciesRichness)  0.7013  0.7013     1 33.362  0.0910 0.76474  
# scale(Duration)        21.4475 21.4475     1 17.017  2.7840 0.11351  
# scale(AKCK)             0.2520  0.2520     1 16.866  0.0327 0.85862
# # 
# #### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles) 11.8035 11.8035     1  5.9202  1.5065 0.2662
# scale(SpeciesRichness)  4.1469  4.1469     1 25.2758  0.5293 0.4736
# scale(Duration)        15.3901 15.3901     1  5.8586  1.9643 0.2117
# scale(ANCK)            18.3852 18.3852     1  3.8044  2.3465 0.2039

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
summary(ModelLatitude)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles) 1.65582 1.65582     1  15.661  0.5918 0.4532
# scale(SpeciesRichness) 0.02241 0.02241     1 131.103  0.0080 0.9288
# scale(Duration)        0.00379 0.00379     1  30.445  0.0014 0.9709
# scale(Latitude)        0.64845 0.64845     1  11.791  0.2317 0.6390

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
summary(ModelLongitude)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles) 1.63676 1.63676     1  18.227  0.5877 0.4531
# scale(SpeciesRichness) 0.00266 0.00266     1 123.492  0.0010 0.9754
# scale(Duration)        0.01129 0.01129     1  41.051  0.0041 0.9495
# scale(Longitude)       0.00374 0.00374     1   7.799  0.0013 0.9717

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
summary(ModelMAPmean)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles) 1.26458 1.26458     1  14.420  0.4505 0.5127
# scale(SpeciesRichness) 0.00533 0.00533     1 104.008  0.0019 0.9653
# scale(Duration)        0.08578 0.08578     1  22.031  0.0306 0.8628
# scale(MAPmean)         1.24459 1.24459     1  14.775  0.4434 0.5157  

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Nitrogen_mineralization)
summary(ModelMATmean)
# Number of obs: 274, groups:  StudyID, 50
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.6875  0.6875     1 11.270  0.2447 0.6304
# scale(SpeciesRichness) 0.0032  0.0032     1 90.449  0.0011 0.9733
# scale(Duration)        0.2943  0.2943     1 18.700  0.1047 0.7498
# scale(MATmean)         6.0006  6.0006     1  7.211  2.1355 0.1861 


############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Nitrogen_mineralization$SpeciesRichness)) ## n = 274
p1 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnNitrogen_mineralization274")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Nitrogen_mineralization$Duration)) ## n = 274
p2 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnNitrogen_mineralization274")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Nitrogen_mineralization$Rotation_cycles)) ## n = 274
p3 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnNitrogen_mineralization274")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Nitrogen_mineralization$Latitude)) ## n = 274
p5 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnNitrogen_mineralization274")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Nitrogen_mineralization$Longitude)) ## n = 274
p6 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnNitrogen_mineralization274")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Nitrogen_mineralization$MAPmean)) ## n = 274
p7 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnNitrogen_mineralization274")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Nitrogen_mineralization$MATmean)) ## n = 274
p8 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnNitrogen_mineralization274")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Nitrogen_mineralization$pHCK)) ## n = 86
p9 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnNitrogen_mineralization86")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Nitrogen_mineralization$SOCCK)) ## n = 83
p10 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnNitrogen_mineralization83")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Nitrogen_mineralization$TNCK)) ## n = 52
p11 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnNitrogen_mineralization52")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Nitrogen_mineralization$NO3CK)) ## n = 46
p12 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnNitrogen_mineralization46")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Nitrogen_mineralization$NH4CK)) ## n = 45
p13<- ggplot(Nitrogen_mineralization, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnNitrogen_mineralization45")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Nitrogen_mineralization$APCK)) ## n = 64
p14 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnNitrogen_mineralization64")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Nitrogen_mineralization$AKCK)) ## n = 39
p15 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnNitrogen_mineralization39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Nitrogen_mineralization$ANCK)) ## n = 31
p16 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnNitrogen_mineralization31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Nitrogen_mineralization$RRpH)) ## n = 86
p17 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="RRpH")+
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
sum(!is.na(Nitrogen_mineralization$RRSOC)) ## n = 83
p18 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="RRSOC")+
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
sum(!is.na(Nitrogen_mineralization$RRTN)) ## n = 52
p19 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="RRTN")+
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
sum(!is.na(Nitrogen_mineralization$RRNO3)) ## n = 44
p20 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="RRNO3")+
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
sum(!is.na(Nitrogen_mineralization$RRNH4)) ## n = 45
p21 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="RRNH4")+
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
sum(!is.na(Nitrogen_mineralization$RRAP)) ## n = 64
p22 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="RRAP")+
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
sum(!is.na(Nitrogen_mineralization$RRAK)) ## n = 39
p23 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="RRAK")+
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
sum(!is.na(Nitrogen_mineralization$RRAN)) ## n = 31
p24 <- ggplot(Nitrogen_mineralization, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnNitrogen_mineralization", x="RRAN")+
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
sum(!is.na(Nitrogen_mineralization$RRYield)) ## n = 49
p25 <- ggplot(Nitrogen_mineralization, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-0.36, 0.15), expand=c(0, 0)) + 
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





