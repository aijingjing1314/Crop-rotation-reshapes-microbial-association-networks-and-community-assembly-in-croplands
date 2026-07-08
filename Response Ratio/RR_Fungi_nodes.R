
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Fungi_nodes <- read.csv("Fungi_nodes.csv", fileEncoding = "latin1")
# Check data
head(Fungi_nodes)

# 1. The number of Obversation
total_number <- nrow(Fungi_nodes)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 302   

# 2. The number of Study
unique_studyid_number <- length(unique(Fungi_nodes$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 41 





#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Fungi_nodes_filteredLongitude_Sub <- subset(Fungi_nodes, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Fungi_nodes_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Fungi_nodes_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredLongitude_Sub %>%
  group_by(Longitude_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Sub Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 LongitudeDa0           130             36
# 2 LongitudeXy0           172              5
overall_model_Fungi_nodes_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Fungi_nodes_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 302; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   97.3351  -194.6703  -188.6703  -177.5589  -188.5892   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0449  0.2118     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 300) = 1815.6314, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.8308, p-val = 0.6601
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0    0.0207  0.0367   0.5646  0.5724  -0.0513  0.0928    
# Longitude_SubLongitudeXy0   -0.0689  0.0964  -0.7156  0.4743  -0.2578  0.1199   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredLongitude_Sub)
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
Fungi_nodes_filteredMAPmean2_Sub <- subset(Fungi_nodes, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Fungi_nodes_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Fungi_nodes_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredMAPmean2_Sub %>%
  group_by(MAPmean2_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MAPmean2_Sub  Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 MAP600Dao1200          111             14
# 2 MAPDy1200              107             10
# 3 MAPXy600                84             18

overall_model_Fungi_nodes_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Fungi_nodes_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 302; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   99.0294  -198.0587  -190.0587  -175.2570  -189.9227   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0384  0.1959     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 299) = 1486.1372, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.0391, p-val = 0.1690
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean2_SubMAP600Dao1200   -0.0006  0.0518  -0.0115  0.9908  -0.1022  0.1010    
# MAPmean2_SubMAPDy1200       -0.0994  0.0659  -1.5084  0.1315  -0.2285  0.0298    
# MAPmean2_SubMAPXy600         0.0760  0.0461   1.6486  0.0992  -0.0144  0.1664  .

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredMAPmean2_Sub)
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
#           "ab"                       "a"                       "b" 




### 8.3 MATmean_Sub
Fungi_nodes_filteredMATmean_Sub <- subset(Fungi_nodes, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Fungi_nodes_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Fungi_nodes_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredMATmean_Sub %>%
  group_by(MATmean_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Sub Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 MAT8Dao15             39             10
# 2 MATDy15              115             14
# 3 MATXy8               148             18

overall_model_Fungi_nodes_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Fungi_nodes_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 302; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   99.1101  -198.2202  -190.2202  -175.4184  -190.0841   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0387  0.1967     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 299) = 1369.4556, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.5955, p-val = 0.1330
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubMAT8Dao15    0.0294  0.0611   0.4807  0.6307  -0.0903  0.1490    
# MATmean_SubMATDy15     -0.0918  0.0552  -1.6610  0.0967  -0.2000  0.0165  . 
# MATmean_SubMATXy8       0.0772  0.0463   1.6665  0.0956  -0.0136  0.1680  . 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredMATmean_Sub)
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
#       "ab"                  "a"                  "b" 




### 8.4 LegumeNonlegume
Fungi_nodes_filteredLegumeNonlegume <- subset(Fungi_nodes, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Fungi_nodes_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Fungi_nodes_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredLegumeNonlegume %>%
  group_by(LegumeNonlegume) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   LegumeNonlegume          Observations Unique_StudyID
#   <fct>                           <int>          <int>
# 1 Legume to Non-legume               55             10
# 2 Non-legume to Legume              151             13
# 3 Non-legume to Non-legume           96             27

overall_model_Fungi_nodes_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Fungi_nodes_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 302; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   95.6059  -191.2117  -183.2117  -168.4100  -183.0757   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0441  0.2101     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 299) = 1390.5859, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.9746, p-val = 0.2642
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume        0.0235  0.0372   0.6307  0.5282  -0.0495  0.0964    
# LegumeNonlegumeNon-legume to Legume       -0.0023  0.0361  -0.0641  0.9489  -0.0730  0.0684    
# LegumeNonlegumeNon-legume to Non-legume    0.0094  0.0350   0.2686  0.7883  -0.0591  0.0779    
 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredLegumeNonlegume)
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
    #              "a"                                     "a"                                     "a"




### 8.5 AMnonAM
Fungi_nodes_filteredAMnonAM <- subset(Fungi_nodes, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Fungi_nodes_filteredAMnonAM$AMnonAM <- droplevels(factor(Fungi_nodes_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
  # AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             262             34
# 2 AM to nonAM           25              7
# 3 nonAM to AM           15              5

overall_model_Fungi_nodes_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Fungi_nodes_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredAMnonAM)
# ultivariate Meta-Analysis Model (k = 302; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   96.7883  -193.5766  -185.5766  -170.7748  -185.4405   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0419  0.2046     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 299) = 1750.8695, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 6.1078, p-val = 0.1065
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM       0.0115  0.0337   0.3409  0.7332  -0.0545  0.0774    
# AMnonAMAM to nonAM   -0.0319  0.0379  -0.8415  0.4000  -0.1063  0.0424    
# AMnonAMnonAM to AM    0.0440  0.0488   0.9010  0.3676  -0.0517  0.1397 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredAMnonAM)
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
   #       "a"                "b"               "ab" 


### 8.6 C3C4
Fungi_nodes_filteredC3C4 <- subset(Fungi_nodes, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Fungi_nodes_filteredC3C4$C3C4 <- droplevels(factor(Fungi_nodes_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredC3C4 %>%
  group_by(C3C4) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   C3C4     Observations Unique_StudyID
#   <fct>           <int>          <int>
# 1 C3 to C3           58             14
# 2 C3 to C4           91             22
# 3 C4 to C3          148             12

overall_model_Fungi_nodes_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Fungi_nodes_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredC3C4)
# ultivariate Meta-Analysis Model (k = 297; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  115.4228  -230.8457  -222.8457  -208.1114  -222.7073   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0353  0.1878     39     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 294) = 1196.3484, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 14.2755, p-val = 0.0026
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb    ci.ub    
# C3C4C3 to C3   -0.1172  0.0489  -2.3968  0.0165  -0.2131  -0.0214  * 
# C3C4C3 to C4    0.0656  0.0357   1.8366  0.0663  -0.0044   0.1355  . 
# C3C4C4 to C3    0.0381  0.0366   1.0416  0.2976  -0.0336   0.1097    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredC3C4)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredC3C4)
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
#   "a"          "b"          "c" 


### 8.7 Annual_Pere
Fungi_nodes_filteredAnnual_Pere <- subset(Fungi_nodes, Annual_Pere %in% c("Annual to Annual", "Perennial to Annual"))
#
Fungi_nodes_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Fungi_nodes_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredAnnual_Pere %>%
  group_by(Annual_Pere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Annual_Pere            Observations Unique_StudyID
#   <fct>                         <int>          <int>
# 1 Annual to Annual                253             30
# 2 Annual to Perennial               9              6
# 3 Perennial to Annual              35              9
# 4 Perennial to Perennial            5              2

overall_model_Fungi_nodes_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Fungi_nodes_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 302; method: d: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  121.2558  -242.5116  -236.5116  -225.5437  -236.4265   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0415  0.2038     39     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 286) = 1706.9007, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 4.5236, p-val = 0.1042
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual       0.0482  0.0388   1.2426  0.2140  -0.0278  0.1242    
# Annual_PerePerennial to Annual   -0.1262  0.0731  -1.7262  0.0843  -0.2695  0.0171  . 
 

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredAnnual_Pere)
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
# "a"                            "b"

### 8.8 Plant_stage
Fungi_nodes_filteredPlant_stage <- subset(Fungi_nodes, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
Fungi_nodes_filteredPlant_stage$Plant_stage <- droplevels(factor(Fungi_nodes_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredPlant_stage %>%
  group_by(Plant_stage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Plant_stage        Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Harvest                     118             21
# 2 Maturity stage                9              5
# 3 Reproductive stage           52              9
# 4 Vegetative stage             93              9

overall_model_Fungi_nodes_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Fungi_nodes_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredPlant_stage)
# # Multivari(k = 263; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  103.2913  -206.5826  -198.5826  -184.3398  -198.4257   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0403  0.2007     30     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 260) = 1453.6441, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.9935, p-val = 0.2622
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb   ci.ub    
# Plant_stageHarvest              -0.0113  0.0383  -0.2960  0.7673  -0.0865  0.0638    
# Plant_stageReproductive stage   -0.0311  0.0396  -0.7856  0.4321  -0.1088  0.0465    
# Plant_stageVegetative stage     -0.0269  0.0387  -0.6946  0.4873  -0.1028  0.0490    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredPlant_stage)
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
Fungi_nodes_filteredRotation_cycles2 <- subset(Fungi_nodes, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
Fungi_nodes_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Fungi_nodes_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredRotation_cycles2 %>%
  group_by(Rotation_cycles2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Rotation_cycles2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                         68             15
# 2 D1-3                      100             15
# 3 D10                        18              6
# 4 D3-5                       80              7
# 5 D5-10                      36             11
overall_model_Fungi_nodes_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Fungi_nodes_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 302; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   98.3095  -196.6191  -184.6191  -162.4567  -184.3294   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0624  0.2498     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 297) = 1639.3706, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 15.7204, p-val = 0.0077
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub     
# Rotation_cycles2D1      -0.0633  0.0470  -1.3468  0.1780  -0.1554  0.0288     
# Rotation_cycles2D1-3    -0.0711  0.0468  -1.5195  0.1286  -0.1629  0.0206     
# Rotation_cycles2D10      0.1597  0.0580   2.7553  0.0059   0.0461  0.2734  ** 
# Rotation_cycles2D3-5     0.1494  0.0567   2.6361  0.0084   0.0383  0.2604  ** 
# Rotation_cycles2D5-10    0.1103  0.0560   1.9695  0.0489   0.0005  0.2200   *   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredRotation_cycles2)
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
#         "a"                   "a"                   "b"                   "b"                   "c" 


### 8.10 Duration2
Fungi_nodes_filteredDuration2 <- subset(Fungi_nodes, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Fungi_nodes_filteredDuration2$Duration2 <- droplevels(factor(Fungi_nodes_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredDuration2 %>%
  group_by(Duration2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Duration2 Observations Unique_StudyID
#   <fct>            <int>          <int>
# 1 D1-5               146             22
# 2 D11-20              35             10
# 3 D20-40              87              6
# 4 D6-10               34             10
overall_model_Fungi_nodes_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Fungi_nodes_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 302; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  100.2347  -200.4695  -190.4695  -171.9840  -190.2640   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0588  0.2425     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 298) = 1682.4920, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 12.0391, p-val = 0.0171
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5     -0.0768  0.0486  -1.5815  0.1138  -0.1720  0.0184    
# Duration2D11-20    0.1211  0.0517   2.3413  0.0192   0.0197  0.2225  * 
# Duration2D20-40    0.0456  0.0646   0.7063  0.4800  -0.0809  0.1722    
# Duration2D6-10     0.1188  0.0531   2.2393  0.0251   0.0148  0.2228  * 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredDuration2)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredDuration2)
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
 #       "a"             "b"            "ab"             "b"

### 8.11 SpeciesRichness2
Fungi_nodes_filteredSpeciesRichness2 <- subset(Fungi_nodes, SpeciesRichness2 %in% c("R2", "R3"))
#
Fungi_nodes_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Fungi_nodes_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredSpeciesRichness2 %>%
  group_by(SpeciesRichness2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   SpeciesRichness2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 R2                        277             39
# 2 R3                         24              4
# 3 R4                          1              1

overall_model_Fungi_nodes_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Fungi_nodes_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredSpeciesRichness2)
# # tivariate Meta-Analysis Model (k = 302; meth 301; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   97.5199  -195.0397  -189.0397  -177.9384  -188.9584   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0409  0.2023     40     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 299) = 1793.6733, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 3.1890, p-val = 0.2030
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2    0.0208  0.0333   0.6266  0.5309  -0.0444  0.0860    
# SpeciesRichness2R3   -0.0151  0.0388  -0.3899  0.6966  -0.0912  0.0609    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredSpeciesRichness2)
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
#      a                    a



### 8.12 Bulk_Rhizosphere
Fungi_nodes_filteredBulk_Rhizosphere <- subset(Fungi_nodes, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Fungi_nodes_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Fungi_nodes_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non_Rhizosphere           249             30
# 2 Rhizosphere                53             19

overall_model_Fungi_nodes_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Fungi_nodes_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 302; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  110.4680  -220.9360  -214.9360  -203.8247  -214.8549   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0453  0.2129     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 300) = 1851.6284, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 31.3807, p-val < .0001
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere   -0.0208  0.0349  -0.5946  0.5521  -0.0892  0.0477    
# Bulk_RhizosphereRhizosphere        0.0676  0.0360   1.8765  0.0606  -0.0030  0.1382  . 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredBulk_Rhizosphere)
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
Fungi_nodes_filteredSoil_texture <- subset(Fungi_nodes, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Fungi_nodes_filteredSoil_texture$Soil_texture <- droplevels(factor(Fungi_nodes_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredSoil_texture %>%
  group_by(Soil_texture) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Soil_texture Observations Unique_StudyID
#   <fct>               <int>          <int>
# 1 Coarse                 20              4
# 2 Fine                   98              9
# 3 Medium                129             10

overall_model_Fungi_nodes_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Fungi_nodes_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 247; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#   73.8691  -147.7381  -139.7381  -125.7495  -139.5708   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0358  0.1892     23     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 244) = 1259.6450, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 5.2984, p-val = 0.1512
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb    ci.ub    
# Soil_textureCoarse   -0.2181  0.1006  -2.1667  0.0303  -0.4153  -0.0208  * 
# Soil_textureFine      0.0477  0.0659   0.7244  0.4688  -0.0814   0.1768    
# Soil_textureMedium    0.0173  0.0616   0.2807  0.7789  -0.1035   0.1381    


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredSoil_texture)
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
#         "a"                "b"                "b" 


### 8.14 Tillage
Fungi_nodes_filteredTillage <- subset(Fungi_nodes, Tillage %in% c("Tillage", "No_tillage"))
#
Fungi_nodes_filteredTillage$Tillage <- droplevels(factor(Fungi_nodes_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredTillage %>%
  group_by(Tillage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Tillage    Observations Unique_StudyID
#   <fct>             <int>          <int>
# 1 No_tillage          103              5
# 2 Tillage              90              7

overall_model_Fungi_nodes_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Fungi_nodes_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredTillage)
# riate Meta-Analysis Model (k = 193; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  167.8894  -335.7789  -329.7789  -320.0221  -329.6505   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0567  0.2380     10     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 191) = 456.7277, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.9174, p-val = 0.2325
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.0755  0.0792  -0.9525  0.3409  -0.2307  0.0798    
# TillageTillage      -0.0315  0.0789  -0.3995  0.6895  -0.1862  0.1231    

#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredTillage)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredTillage)
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
Fungi_nodes_filteredStraw_retention <- subset(Fungi_nodes, Straw_retention %in% c("Retention", "No_retention"))
#
Fungi_nodes_filteredStraw_retention$Straw_retention <- droplevels(factor(Fungi_nodes_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredStraw_retention %>%
  group_by(Straw_retention) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Straw_retention Observations Unique_StudyID
#   <fct>                  <int>          <int>
# 1 No_retention              11              5
# 2 Retention                 32              7

overall_model_Fungi_nodes_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Fungi_nodes_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredStraw_retention)
# # Miate Meta-Analysis Model (k = 43; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  -7.6062   15.2124   21.2124   26.3531   21.8611   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0165  0.1285     10     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 41) = 200.1313, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.3462, p-val = 0.8411
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention    0.0056  0.0485   0.1146  0.9088  -0.0895  0.1006    
# Straw_retentionRetention      -0.0159  0.0462  -0.3450  0.7301  -0.1064  0.0746    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredStraw_retention)
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
Fungi_nodes_filteredPrimer <- subset(Fungi_nodes, Primer %in% c("ITS1", "ITS2" ))
#
Fungi_nodes_filteredPrimer$Primer <- droplevels(factor(Fungi_nodes_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Fungi_nodes_filteredPrimer %>%
  group_by(Primer) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Primer             Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 ITS1                        218             38
# 2 ITS1 + 5.8S + ITS2            4              1
# 3 ITS2                         80              2
overall_model_Fungi_nodes_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Fungi_nodes_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Fungi_nodes_filteredPrimer)
# alysis Model (k = 298; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
#  101.9074  -203.8148  -197.8148  -186.7437  -197.7326   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0453  0.2129     40     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 296) = 1445.9561, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.2260, p-val = 0.8932
# 
# Model Results:
# 
#             estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerITS1    0.0170  0.0359   0.4737  0.6357  -0.0534  0.0874    
# PrimerITS2   -0.0060  0.1509  -0.0399  0.9682  -0.3018  0.2898       

#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_nodes_filteredPrimer)
vcov_rotation <- vcov(overall_model_Fungi_nodes_filteredPrimer)
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
              # PrimerITS1    PrimerITS2 
              #        "a"        "a" 









#### 9. Linear Mixed Effect Model
# 
Fungi_nodes$Wr <- 1 / Fungi_nodes$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Fungi_nodes)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Fungi_nodes)

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
# Model1   Model1 -289.4569 -267.1943 150.7284
# Model2   Model2 -289.5651 -267.3026 150.7826
# Model3   Model3 -290.4118 -268.1492 151.2059##############################
# Model4   Model4 -289.2632 -267.0006 150.6316
# Model5   Model5 -288.7382 -266.4757 150.3691
# Model6   Model6 -288.5428 -266.2802 150.2714
# Model7   Model7 -289.7570 -267.4945 150.8785
# Model8   Model8 -290.2220 -267.9594 151.1110
# Model9   Model9 -274.4576 -237.3534 147.2288
# Model10 Model10 -272.1547 -235.0504 146.0773

##### Model 3 is the best model
summary(Model3)
# Number of obs: 302, groups:  StudyID, 41
anova(Model3)
# Type III Ane III Analysis of Variance Table with Satterthwaite's method
# Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(SpeciesRichness)      5.8109  5.8109     1 289.39  2.4274 0.1203
# scale(Duration)             0.0403  0.0403     1 125.34  0.0168 0.8970
# scale(log(Rotation_cycles)) 0.3824  0.3824     1 287.66  0.1597 0.6897

#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
summary(ModelpH)
# Number of obs: 90, groups:  StudyID, 27
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles))  0.0960  0.0960     1 80.063  0.0405 0.841097   
# scale(SpeciesRichness)       1.1713  1.1713     1 63.917  0.4939 0.484746   
# scale(Duration)              0.0050  0.0050     1 71.304  0.0021 0.963404   
# scale(pHCK)                 25.8883 25.8883     1 84.298 10.9159 0.001401 **  

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
summary(ModelSOC)
# Number of obs: 101, groups:  StudyID, 28
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 0.1425  0.1425     1 91.146  0.0603 0.80663  
# scale(SpeciesRichness)      7.3457  7.3457     1 80.235  3.1070 0.08177 .
# scale(Duration)             0.0777  0.0777     1 73.788  0.0329 0.85663  
# scale(SOCCK)                0.3839  0.3839     1 42.295  0.1624 0.68902 

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
summary(ModelTN)
# Number of obs: 63, groups:  StudyID, 17
anova(ModelTN) 
# > anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.11476 0.11476     1 48.416  0.0773 0.7822
# scale(SpeciesRichness)      0.30707 0.30707     1 54.087  0.2068 0.6511
# scale(Duration)             0.51987 0.51987     1 57.576  0.3500 0.5564
# scale(TNCK)                 0.00096 0.00096     1 13.214  0.0006 0.9801


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
summary(ModelNO3)
# Number of obs: 29, groups:  StudyID, 12
anova(ModelNO3)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  0.0623  0.0623     1 10.3961  0.0289 0.86834  
# scale(SpeciesRichness)      16.0277 16.0277     1 12.6970  7.4284 0.01764 *
# scale(Duration)              0.0381  0.0381     1 23.1931  0.0176 0.89549  
# scale(NO3CK)                 0.0614  0.0614     1  7.5286  0.0285 0.87044  

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
summary(ModelNH4)
# Number of obs: 29, groups:  StudyID, 12
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  0.1253  0.1253     1  9.8863  0.0578 0.81492  
# scale(SpeciesRichness)      16.1275 16.1275     1 12.0927  7.4415 0.01823 *
# scale(Duration)              0.0091  0.0091     1 23.3961  0.0042 0.94893  
# scale(NH4CK)                 0.0030  0.0030     1  7.3658  0.0014 0.97145  

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
summary(ModelAP)
# Number of obs: 57, groups:  StudyID, 23
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 10.5533 10.5533     1 30.772  3.8356 0.05929 .
# scale(SpeciesRichness)       2.3742  2.3742     1 47.320  0.8629 0.35764  
# scale(Duration)              4.1529  4.1529     1 51.922  1.5094 0.22477  
# scale(APCK)                  5.5888  5.5888     1 45.186  2.0313 0.16096 

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
summary(ModelAK)
# Number of obs: 43, groups:  StudyID, 18
anova(ModelAK)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 6.3010  6.3010     1 17.893  2.0306 0.1714
# scale(SpeciesRichness)      3.5153  3.5153     1 37.508  1.1329 0.2940
# scale(Duration)             2.0714  2.0714     1 15.084  0.6675 0.4266
# scale(AKCK)                 0.6524  0.6524     1 33.973  0.2102 0.6495

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
summary(ModelAN)
# Number of obs: 42, groups:  StudyID, 12
anova(ModelAN)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 12.1373 12.1373     1 14.0774  2.6443 0.12609  
# scale(SpeciesRichness)      15.8343 15.8343     1 17.6519  3.4497 0.08003 .
# scale(Duration)              5.5142  5.5142     1 11.2808  1.2013 0.29591  
# scale(ANCK)                  0.0008  0.0008     1  5.1234  0.0002 0.99020  

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
summary(ModelLatitude)
# Number of obs: 302, groups:  StudyID, 41
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles))  0.3470  0.3470     1 290.776  0.1435 0.70506  
# scale(SpeciesRichness)       7.0997  7.0997     1 290.547  2.9369 0.08764 .
# scale(Duration)              0.4573  0.4573     1  98.356  0.1892 0.66454  
# scale(Latitude)             12.8526 12.8526     1  29.651  5.3168 0.02830 *

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
summary(ModelLongitude)
# Number of obs: 302, groups:  StudyID, 41
anova(ModelLongitude)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.1798  0.1798     1 282.98  0.0754 0.7839
# scale(SpeciesRichness)      5.6462  5.6462     1 289.06  2.3665 0.1251
# scale(Duration)             0.0698  0.0698     1 140.34  0.0292 0.8645
# scale(Longitude)            3.1398  3.1398     1  30.47  1.3159 0.2602

# ### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
summary(ModelMAPmean)
# Number of obs: 302, groups:  StudyID, 41
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 0.3468  0.3468     1 289.14  0.1439 0.70475  
# scale(SpeciesRichness)      6.5869  6.5869     1 289.56  2.7322 0.09943 .
# scale(Duration)             0.2899  0.2899     1 109.22  0.1202 0.72945  
# scale(MAPmean)              5.5060  5.5060     1  32.60  2.2838 0.14036  

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Fungi_nodes)
summary(ModelMATmean)
# Number of obs: 302, groups:  StudyID, 41
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles))  0.2667  0.2667     1 292.305  0.1102 0.740099   
# scale(SpeciesRichness)       7.6294  7.6294     1 291.286  3.1535 0.076808 . 
# scale(Duration)              0.4720  0.4720     1  94.468  0.1951 0.659732   
# scale(MATmean)              18.8456 18.8456     1  29.159  7.7895 0.009173 **
# # ---


############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Fungi_nodes$SpeciesRichness)) ## n = 302
p1 <- ggplot(Fungi_nodes, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnFungi_nodes302")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Fungi_nodes$Duration)) ## n = 302
p2 <- ggplot(Fungi_nodes, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnFungi_nodes302")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Fungi_nodes$Rotation_cycles)) ## n = 302
p3 <- ggplot(Fungi_nodes, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnFungi_nodes302")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Fungi_nodes$Latitude)) ## n = 302
p5 <- ggplot(Fungi_nodes, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnFungi_nodes302")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Fungi_nodes$Longitude)) ## n = 302
p6 <- ggplot(Fungi_nodes, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnFungi_nodes302")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Fungi_nodes$MAPmean)) ## n = 302
p7 <- ggplot(Fungi_nodes, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnFungi_nodes302")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Fungi_nodes$MATmean)) ## n = 302
p8 <- ggplot(Fungi_nodes, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnFungi_nodes302")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Fungi_nodes$pHCK)) ## n = 90
p9 <- ggplot(Fungi_nodes, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnFungi_nodes90")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Fungi_nodes$SOCCK)) ## n = 101
p10 <- ggplot(Fungi_nodes, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnFungi_nodes101")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Fungi_nodes$TNCK)) ## n = 63
p11 <- ggplot(Fungi_nodes, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnFungi_nodes63")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Fungi_nodes$NO3CK)) ## n = 29
p12 <- ggplot(Fungi_nodes, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnFungi_nodes29")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Fungi_nodes$NH4CK)) ## n = 29
p13<- ggplot(Fungi_nodes, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnFungi_nodes29")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Fungi_nodes$APCK)) ## n = 57
p14 <- ggplot(Fungi_nodes, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnFungi_nodes57")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Fungi_nodes$AKCK)) ## n = 43
p15 <- ggplot(Fungi_nodes, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnFungi_nodes43")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Fungi_nodes$ANCK)) ## n = 42
p16 <- ggplot(Fungi_nodes, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnFungi_nodes42")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Fungi_nodes$RRpH)) ## n = 90
p17 <- ggplot(Fungi_nodes, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="RRpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRpH" , y="RR90")
p17
pdf("RRpH.pdf",width=8,height=8)
p17
dev.off() 

## RRSOC
sum(!is.na(Fungi_nodes$RRSOC)) ## n = 101
p18 <- ggplot(Fungi_nodes, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="RRSOC")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRSOC" , y="RR101")
p18
pdf("RRSOC.pdf",width=8,height=8)
p18
dev.off() 

## RRTN
sum(!is.na(Fungi_nodes$RRTN)) ## n = 63
p19 <- ggplot(Fungi_nodes, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="RRTN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRTN" , y="RR63")
p19
pdf("RRTN.pdf",width=8,height=8)
p19
dev.off() 

## RRNO3
sum(!is.na(Fungi_nodes$RRNO3)) ## n = 29
p20 <- ggplot(Fungi_nodes, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="RRNO3")+
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
sum(!is.na(Fungi_nodes$RRNH4)) ## n = 29
p21 <- ggplot(Fungi_nodes, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="RRNH4")+
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
sum(!is.na(Fungi_nodes$RRAP)) ## n = 57
p22 <- ggplot(Fungi_nodes, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="RRAP")+
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
sum(!is.na(Fungi_nodes$RRAK)) ## n = 43
p23 <- ggplot(Fungi_nodes, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="RRAK")+
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
sum(!is.na(Fungi_nodes$RRAN)) ## n = 42
p24 <- ggplot(Fungi_nodes, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_nodes", x="RRAN")+
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
sum(!is.na(Fungi_nodes$RRYield)) ## n = 51
p25 <- ggplot(Fungi_nodes, aes(x=RR, y=RRYield)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm, color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") +
  theme_bw() +
  theme(text = element_text(family = "serif", size=20)) +
  geom_vline(aes(xintercept=0), colour="black", linewidth=0.5, linetype="dashed") +
  labs(x="RR", y="RRYield51") +
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
  scale_x_continuous(limits=c(-0.55, 0.55), expand=c(0, 0)) + 
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








