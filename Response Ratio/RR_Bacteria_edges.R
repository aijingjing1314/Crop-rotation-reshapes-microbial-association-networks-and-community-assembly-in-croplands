
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Bacteria_edges <- read.csv("Bacteria_edges.csv", fileEncoding = "latin1")
# Check data
head(Bacteria_edges)

# 1. The number of Obversation
total_number <- nrow(Bacteria_edges)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 266  

# 2. The number of Study
unique_studyid_number <- length(unique(Bacteria_edges$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 47 


#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Bacteria_edges_filteredLongitude_Sub <- subset(Bacteria_edges, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Bacteria_edges_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Bacteria_edges_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredLongitude_Sub %>%
  group_by(Longitude_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Sub Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 LongitudeDa0           109             38
# 2 LongitudeXy0           157              9

overall_model_Bacteria_edges_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Bacteria_edges_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -51.9342  103.8684  109.8684  120.5962  109.9607   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0199  0.1410     47     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 264) = 1259.3854, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.0532, p-val = 0.9737
# Model Results:
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0    0.0021  0.0266   0.0806  0.9358  -0.0499  0.0542    
# Longitude_SubLongitudeXy0   -0.0107  0.0497  -0.2161  0.8289  -0.1081  0.0866  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredLongitude_Sub)
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
Bacteria_edges_filteredMAPmean2_Sub <- subset(Bacteria_edges, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Bacteria_edges_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Bacteria_edges_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredMAPmean2_Sub %>%
  group_by(MAPmean2_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MAPmean2_Sub  Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 MAP600Dao1200           66             17
# 2 MAPDy1200               98             11
# 3 MAPXy600               102             20

overall_model_Bacteria_edges_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Bacteria_edges_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -50.5527  101.1054  109.1054  123.3940  109.2605   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0181  0.1346     47     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 263) = 1234.0073, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 4.0189, p-val = 0.2594
# Model Results:
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean2_SubMAP600Dao1200   -0.0392  0.0298  -1.3150  0.1885  -0.0976  0.0192    
# MAPmean2_SubMAPDy1200        0.0283  0.0590   0.4797  0.6314  -0.0873  0.1439    
# MAPmean2_SubMAPXy600         0.0216  0.0282   0.7640  0.4449  -0.0338  0.0769 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredMAPmean2_Sub)
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
Bacteria_edges_filteredMATmean_Sub <- subset(Bacteria_edges, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Bacteria_edges_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Bacteria_edges_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredMATmean_Sub %>%
  group_by(MATmean_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Sub Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 MAT8Dao15             46             11
# 2 MATDy15              101             15
# 3 MATXy8               119             22

overall_model_Bacteria_edges_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Bacteria_edges_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -51.2739  102.5477  110.5477  124.8363  110.7027   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0187  0.1369     47     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 263) = 1251.8557, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.9706, p-val = 0.3962
# Model Results:
#                       estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubMAT8Dao15   -0.0415  0.0350  -1.1853  0.2359  -0.1100  0.0271    
# MATmean_SubMATDy15      0.0014  0.0462   0.0295  0.9765  -0.0892  0.0919    
# MATmean_SubMATXy8       0.0159  0.0282   0.5654  0.5718  -0.0393  0.0711 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredMATmean_Sub)
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
Bacteria_edges_filteredLegumeNonlegume <- subset(Bacteria_edges, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Bacteria_edges_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Bacteria_edges_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredLegumeNonlegume %>%
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
# 3 Non-legume to Non-legume           79             22

overall_model_Bacteria_edges_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Bacteria_edges_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -46.2012   92.4024  100.4024  114.6910  100.5574   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0210  0.1450     47     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 263) = 1252.1522, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 17.8413, p-val = 0.0005
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume       -0.0377  0.0308  -1.2252  0.2205  -0.0981  0.0226    
# LegumeNonlegumeNon-legume to Legume        0.0247  0.0249   0.9896  0.3224  -0.0242  0.0735    
# LegumeNonlegumeNon-legume to Non-legume   -0.0029  0.0252  -0.1152  0.9083  -0.0522  0.0464   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredLegumeNonlegume)
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
Bacteria_edges_filteredAMnonAM <- subset(Bacteria_edges, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Bacteria_edges_filteredAMnonAM$AMnonAM <- droplevels(factor(Bacteria_edges_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
  # AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             206             36
# 2 AM to nonAM           30              8
# 3 nonAM to AM           30              6

overall_model_Bacteria_edges_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Bacteria_edges_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -51.5933  103.1867  111.1867  125.4753  111.3417   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0196  0.1399     47     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 263) = 1162.0365, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.3330, p-val = 0.3431
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM      -0.0165  0.0249  -0.6614  0.5084  -0.0653  0.0324    
# AMnonAMAM to nonAM   -0.0037  0.0331  -0.1107  0.9119  -0.0684  0.0611    
# AMnonAMnonAM to AM    0.1085  0.0656   1.6533  0.0983  -0.0201  0.2372  . 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredAMnonAM)
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
Bacteria_edges_filteredC3C4 <- subset(Bacteria_edges, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Bacteria_edges_filteredC3C4$C3C4 <- droplevels(factor(Bacteria_edges_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredC3C4 %>%
  group_by(C3C4) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
  # C3C4     Observations Unique_StudyID
#   <fct>           <int>          <int>
# 1 C3 to C3           72             17
# 2 C3 to C4           62             23
# 3 C4 to C3          131             14

overall_model_Bacteria_edges_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Bacteria_edges_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredC3C4)
# Multivariate Meta-Analysis Model (k = 265; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -47.6482   95.2963  103.2963  117.5697  103.4520   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0209  0.1446     46     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 262) = 1249.8961, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 12.9358, p-val = 0.0048
# Model Results:
#               estimate      se     zval    pval    ci.lb   ci.ub    
# C3C4C3 to C3    0.0356  0.0273   1.3028  0.1926  -0.0180  0.0892    
# C3C4C3 to C4   -0.0330  0.0263  -1.2545  0.2097  -0.0845  0.0185    
# C3C4C4 to C3    0.0161  0.0287   0.5619  0.5742  -0.0401  0.0723   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredC3C4)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredC3C4)
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
Bacteria_edges_filteredAnnual_Pere <- subset(Bacteria_edges, Annual_Pere %in% c("Annual to Annual", "Perennial to Perennial", "Annual to Perennial", "Perennial to Annual"))
#
Bacteria_edges_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Bacteria_edges_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredAnnual_Pere %>%
  group_by(Annual_Pere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Annual_Pere         Observations Unique_StudyID
#   <fct>                      <int>          <int>
# 1 Annual to Annual             233             37
# 2 Annual to Perennial           12              7
# 3 Perennial to Annual           21              7

overall_model_Bacteria_edges_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Bacteria_edges_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredAnnual_Pere)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -52.2057  104.4114  112.4114  126.7000  112.5664   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0217  0.1472     47     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 263) = 1224.6750, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.9942, p-val = 0.2621
# Model Results:
#                                 estimate      se     zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual       0.0106  0.0252   0.4191  0.6751  -0.0389  0.0600    
# Annual_PereAnnual to Perennial    0.0017  0.0349   0.0494  0.9606  -0.0666  0.0700    
# Annual_PerePerennial to Annual   -0.0896  0.0515  -1.7389  0.0820  -0.1906  0.0114  . 

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredAnnual_Pere)
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
Bacteria_edges_filteredPlant_stage <- subset(Bacteria_edges, Plant_stage %in% c("Vegetative stage","Reproductive stage","Harvest"))
#
Bacteria_edges_filteredPlant_stage$Plant_stage <- droplevels(factor(Bacteria_edges_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredPlant_stage %>%
  group_by(Plant_stage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Plant_stage        Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Harvest                      96             22
# 2 Maturity stage                7              5
# 3 Reproductive stage           61             11
# 4 Vegetative stage             77              7

overall_model_Bacteria_edges_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Bacteria_edges_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredPlant_stage)
# ltivariate Meta-Analysis Model (k = 234; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -43.4565   86.9130   94.9130  108.6827   95.0900   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0335  0.1830     35     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 231) = 1165.5998, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 55.0876, p-val < .0001
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb    ci.ub     
# Plant_stageHarvest               0.0419  0.0340   1.2325  0.2178  -0.0247   0.1086     
# Plant_stageReproductive stage   -0.1233  0.0375  -3.2855  0.0010  -0.1969  -0.0498  ** 
# Plant_stageVegetative stage      0.0545  0.0342   1.5943  0.1109  -0.0125   0.1215     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredPlant_stage)
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
Bacteria_edges_filteredRotation_cycles2 <- subset(Bacteria_edges, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10"))
#
Bacteria_edges_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Bacteria_edges_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredRotation_cycles2 %>%
  group_by(Rotation_cycles2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Rotation_cycles2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                         77             18
# 2 D1-3                       87             12
# 3 D10                         9              5
# 4 D3-5                       33              7
# 5 D5-10                      60             13
overall_model_Bacteria_edges_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Bacteria_edges_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredRotation_cycles2)
# tivariate Meta-Analysis Model (k = 257; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -55.5441  111.0883  121.0883  138.7552  121.3312   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0215  0.1468     44     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 253) = 1209.6663, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 10.3375, p-val = 0.0351
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1       0.0017  0.0320   0.0547  0.9564  -0.0610  0.0645    
# Rotation_cycles2D1-3     0.0173  0.0320   0.5415  0.5882  -0.0454  0.0800    
# Rotation_cycles2D3-5     0.0025  0.0436   0.0567  0.9548  -0.0830  0.0880    
# Rotation_cycles2D5-10   -0.0088  0.0414  -0.2115  0.8325  -0.0900  0.0725     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredRotation_cycles2)
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
# Rotation_cycles2D1  Rotation_cycles2D1-3  Rotation_cycles2D3-5 Rotation_cycles2D5-10 
# "a"                   "b"                  "ab"                  "ab" 

### 8.10 Duration2
Bacteria_edges_filteredDuration2 <- subset(Bacteria_edges, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Bacteria_edges_filteredDuration2$Duration2 <- droplevels(factor(Bacteria_edges_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredDuration2 %>%
  group_by(Duration2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Duration2 Observations Unique_StudyID
#   <fct>            <int>          <int>
# 1 D1-5               144             25
# 2 D11-20              56              9
# 3 D20-40              35              6
# 4 D6-10               31              9
overall_model_Bacteria_edges_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Bacteria_edges_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -51.2796  102.5591  112.5591  130.4009  112.7935   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0193  0.1388     47     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 262) = 1235.8534, p-val < .0001
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 4.2397, p-val = 0.3745
# Model Results:
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5      0.0288  0.0321   0.8959  0.3703  -0.0342  0.0918    
# Duration2D11-20   -0.0177  0.0397  -0.4450  0.6563  -0.0956  0.0602    
# Duration2D20-40    0.0116  0.0451   0.2574  0.7969  -0.0767  0.0999    
# Duration2D6-10    -0.0858  0.0499  -1.7172  0.0859  -0.1836  0.0121  . 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredDuration2)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredDuration2)
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
Bacteria_edges_filteredSpeciesRichness2 <- subset(Bacteria_edges, SpeciesRichness2 %in% c("R2", "R3"))
#
Bacteria_edges_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Bacteria_edges_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredSpeciesRichness2 %>%
  group_by(SpeciesRichness2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   SpeciesRichness2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 R2                        206             41
# 2 R3                         51              8
# 3 R4                          9              3

overall_model_Bacteria_edges_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Bacteria_edges_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredSpeciesRichness2)
# is Model (k = 257; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -52.9470  105.8941  111.8941  122.5179  111.9897   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0183  0.1354     44     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 255) = 1210.2775, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 3.0472, p-val = 0.2179
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2    0.0133  0.0234   0.5673  0.5705  -0.0326  0.0591    
# SpeciesRichness2R3   -0.0145  0.0273  -0.5319  0.5948  -0.0681  0.0391  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredSpeciesRichness2)
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
# SpeciesRichness2R2 SpeciesRichness2R3 Speci
#                "a"               "a"             


### 8.12 Bulk_Rhizosphere
Bacteria_edges_filteredBulk_Rhizosphere <- subset(Bacteria_edges, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Bacteria_edges_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Bacteria_edges_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non_Rhizosphere           207             35
# 2 Rhizosphere                59             22

overall_model_Bacteria_edges_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Bacteria_edges_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 266; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -52.3128  104.6255  110.6255  121.3534  110.7178   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0193  0.1389     47     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 264) = 1256.4359, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 2.5304, p-val = 0.2822
# Model Results:
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere   -0.0079  0.0235  -0.3368  0.7363  -0.0541  0.0382    
# Bulk_RhizosphereRhizosphere        0.0129  0.0247   0.5245  0.5999  -0.0354  0.0613    

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredBulk_Rhizosphere)
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
Bacteria_edges_filteredSoil_texture <- subset(Bacteria_edges, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Bacteria_edges_filteredSoil_texture$Soil_texture <- droplevels(factor(Bacteria_edges_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredSoil_texture %>%
  group_by(Soil_texture) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#  Soil_texture Observations Unique_StudyID
#   <fct>               <int>          <int>
# 1 Coarse                 30              7
# 2 Fine                   18              7
# 3 Medium                174             17

overall_model_Bacteria_edges_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Bacteria_edges_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 222; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
# -33.6711   67.3422   75.3422   88.8985   75.5291   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0118  0.1087     31     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 219) = 982.3855, p-val < .0001
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 2.2626, p-val = 0.5197
# Model Results:
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse   -0.0387  0.0462  -0.8384  0.4018  -0.1292  0.0518    
# Soil_textureFine     -0.0670  0.0555  -1.2060  0.2278  -0.1759  0.0419    
# Soil_textureMedium   -0.0095  0.0293  -0.3245  0.7456  -0.0669  0.0479   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredSoil_texture)
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
Bacteria_edges_filteredTillage <- subset(Bacteria_edges, Tillage %in% c("Tillage", "No_tillage"))
#
Bacteria_edges_filteredTillage$Tillage <- droplevels(factor(Bacteria_edges_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredTillage %>%
  group_by(Tillage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Tillage    Observations Unique_StudyID
#   <fct>             <int>          <int>
# 1 No_tillage          119              7
# 2 Tillage              29              8

overall_model_Bacteria_edges_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Bacteria_edges_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredTillage)
# Multivariate Meta-Analysis Model (k = 148; method: REML)
#   logLik  Deviance       AIC       BIC      AICc   
#  35.9174  -71.8349  -65.8349  -56.8841  -65.6659   
# Variance Components:
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0058  0.0763     12     no  StudyID 
# Test for Residual Heterogeneity:
# QE(df = 146) = 534.9127, p-val < .0001
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 4.7988, p-val = 0.0908
# Model Results:
#                    estimate      se     zval    pval    ci.lb    ci.ub    
# TillageNo_tillage   -0.0599  0.0277  -2.1648  0.0304  -0.1141  -0.0057  * 
# TillageTillage      -0.0275  0.0287  -0.9582  0.3380  -0.0836   0.0287      
 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredTillage)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredTillage)
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
Bacteria_edges_filteredStraw_retention <- subset(Bacteria_edges, Straw_retention %in% c("Retention", "No_retention"))
#
Bacteria_edges_filteredStraw_retention$Straw_retention <- droplevels(factor(Bacteria_edges_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredStraw_retention %>%
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
overall_model_Bacteria_edges_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Bacteria_edges_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 55; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -34.0091   68.0182   74.0182   79.9291   74.5080   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0078  0.0881     11     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 53) = 257.7742, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 6.6680, p-val = 0.0357
# 
# Model Results:
#                              estimate      se     zval    pval    ci.lb    ci.ub    
# Straw_retentionNo_retention   -0.0849  0.0359  -2.3641  0.0181  -0.1553  -0.0145  * 
# Straw_retentionRetention      -0.0774  0.0325  -2.3847  0.0171  -0.1410  -0.0138  * 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredStraw_retention)
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
Bacteria_edges_filteredPrimer <- subset(Bacteria_edges, Primer %in% c("V3-V4", "V4", "V4-V5", "V5-V7"))
#
Bacteria_edges_filteredPrimer$Primer <- droplevels(factor(Bacteria_edges_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Bacteria_edges_filteredPrimer %>%
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
# 3 V4               97              8
# 4 V4-V5            39              9
overall_model_Bacteria_edges_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Bacteria_edges_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Bacteria_edges_filteredPrimer)
# # Multivariate Meta-Analysis Model (k = 266; method: REML) Model (k = 262; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -47.3909   94.7818  102.7818  117.0091  102.9393   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.0194  0.1394     45     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 259) = 1195.1326, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.9966, p-val = 0.5731
# 
# Model Results:
# 
#              estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerV3-V4    0.0146  0.0295   0.4949  0.6207  -0.0432  0.0723    
# PrimerV4      -0.0752  0.0646  -1.1642  0.2443  -0.2019  0.0514    
# PrimerV4-V5   -0.0327  0.0519  -0.6296  0.5290  -0.1344  0.0691     
 
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Bacteria_edges_filteredPrimer)
vcov_rotation <- vcov(overall_model_Bacteria_edges_filteredPrimer)
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
   # PrimerV3 PrimerV3-V4    PrimerV4 PrimerV4-V5 
   #      "a"         "a"         "a"         "a" 









#### 9. Linear Mixed Effect Model
# 
Bacteria_edges$Wr <- 1 / Bacteria_edges$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Bacteria_edges)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Bacteria_edges)

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
# 
model_comparison <- data.frame(
  Model = names(models),
  AIC = sapply(models, AIC),
  BIC = sapply(models, BIC),
  logLik = sapply(models, logLik)
)
# 查看结果
print(model_comparison)
#           Model       AIC        BIC   logLik
# Model1   Model1 -147.7632 -126.26220 79.88159
# Model2   Model2 -144.6985 -123.19754 78.34926
# Model3   Model3 -143.4818 -121.98086 77.74092
# Model4   Model4 -147.2075 -125.70650 79.60374
# Model5   Model5 -150.0062 -128.50518 81.00308############################
# Model6   Model6 -149.3947 -127.89377 80.69737
# Model7   Model7 -145.1910 -123.69001 78.59549
# Model8   Model8 -143.0845 -121.58351 77.54225
# Model9   Model9 -120.8890  -85.05401 70.44449
# Model10 Model10 -119.3528  -83.51786 69.67641

##### Model 5 is the best model
summary(Model5)
# Number of obs: 266, groups:  StudyID, 47
anova(Model5) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 18.372  18.372     1  90.696  5.3029 0.02358 *
# scale(SpeciesRichness) 20.704  20.704     1 203.564  5.9759 0.01535 *
# scale(log(Duration))   22.251  22.251     1 230.375  6.4225 0.01193 *


#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
summary(ModelpH)
# Number of obs: 78, groups:  StudyID, 29
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.0685  0.0685     1 16.305  0.0154 0.9026
# scale(SpeciesRichness) 7.1606  7.1606     1 69.871  1.6145 0.2081
# scale(log(Duration))   4.1129  4.1129     1 26.723  0.9273 0.3442
# scale(pHCK)            0.8528  0.8528     1 14.855  0.1923 0.6673

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
summary(ModelSOC)
# Number of obs: 75, groups:  StudyID, 27
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles)  1.1366  1.1366     1 20.443  0.2134 0.6490
# scale(SpeciesRichness)  2.8038  2.8038     1 51.525  0.5265 0.4714
# scale(log(Duration))    0.1597  0.1597     1 32.786  0.0300 0.8636
# scale(SOCCK)           10.5892 10.5892     1  7.582  1.9884 0.1982

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
summary(ModelTN)
# Number of obs: 44, groups:  StudyID, 16
anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)  
# scale(Rotation_cycles) 7.2758  7.2758     1  6.025  2.5419 0.1618  
# scale(SpeciesRichness) 9.1704  9.1704     1 38.172  3.2038 0.0814 .
# scale(log(Duration))   6.5773  6.5773     1 18.528  2.2979 0.1464  
# scale(TNCK)            2.7876  2.7876     1 35.286  0.9739 0.3304  


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
summary(ModelNO3)
# Number of obs: 38, groups:  StudyID, 13
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 6.7253  6.7253     1 31.453  6.9118 0.01313 *
# scale(SpeciesRichness) 3.8999  3.8999     1 21.025  4.0080 0.05835 .
# scale(log(Duration))   0.5931  0.5931     1 29.092  0.6095 0.44128  
# scale(NO3CK)           0.0484  0.0484     1 16.916  0.0498 0.82614  

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
summary(ModelNH4)
# Number of obs: 39, groups:  StudyID, 12
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 4.5860  4.5860     1 33.803  4.2366 0.04734 *
# scale(SpeciesRichness) 4.0350  4.0350     1 23.738  3.7276 0.06555 .
# scale(log(Duration))   0.0414  0.0414     1 30.419  0.0383 0.84617  
# scale(NH4CK)           0.0229  0.0229     1  8.618  0.0212 0.88763  

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + scale(APCK) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
summary(ModelAP)
# Number of obs: 56, groups:  StudyID, 23
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Rotation_cycles) 2.5625  2.5625     1 18.498  0.6106 0.4445
# scale(SpeciesRichness) 4.8174  4.8174     1 44.145  1.1479 0.2898
# scale(log(Duration))   0.0089  0.0089     1 29.735  0.0021 0.9635
# scale(APCK)            0.9348  0.9348     1 47.596  0.2227 0.6391

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
summary(ModelAK)
# Number of obs: 39, groups:  StudyID, 18
anova(ModelAK) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# scale(Rotation_cycles) 1.2153  1.2153     1    34  0.2070 0.6520
# scale(SpeciesRichness) 9.1261  9.1261     1    34  1.5547 0.2210
# scale(log(Duration))   1.2388  1.2388     1    34  0.2110 0.6489
# scale(AKCK)            2.4558  2.4558     1    34  0.4184 0.5221

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
summary(ModelAN)
# Number of obs: 31, groups:  StudyID, 12
anova(ModelAN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Rotation_cycles) 0.28590 0.28590     1 10.7840  0.0380 0.8491
# scale(SpeciesRichness) 1.27034 1.27034     1 24.5714  0.1687 0.6848
# scale(log(Duration))   1.88882 1.88882     1  9.4448  0.2509 0.6279
# scale(ANCK)            0.28427 0.28427     1  4.4966  0.0378 0.8544

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
summary(ModelLatitude)
# Number of obs: 266, groups:  StudyID, 47
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 16.7280 16.7280     1  87.963  4.8047 0.03102 *
# scale(SpeciesRichness) 21.7796 21.7796     1 202.370  6.2556 0.01317 *
# scale(log(Duration))   20.0678 20.0678     1 233.608  5.7639 0.01714 *
# scale(Latitude)         1.2596  1.2596     1  85.076  0.3618 0.54911  

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
summary(ModelLongitude)
# Number of obs: 266, groups:  StudyID, 47
anova(ModelLongitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 16.4741 16.4741     1  81.082  4.7641 0.03195 *
# scale(SpeciesRichness) 18.4841 18.4841     1 204.465  5.3454 0.02177 *
# scale(log(Duration))   22.5788 22.5788     1 243.119  6.5295 0.01122 *
# scale(Longitude)        0.4012  0.4012     1  21.710  0.1160 0.73665  

### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
summary(ModelMAPmean)
# Number of obs: 266, groups:  StudyID, 47
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 17.3499 17.3499     1  89.325  4.9960 0.02790 *
# scale(SpeciesRichness) 20.1562 20.1562     1 198.664  5.8041 0.01690 *
# scale(log(Duration))   21.5361 21.5361     1 226.313  6.2014 0.01348 *
# scale(MAPmean)          1.1818  1.1818     1  47.467  0.3403 0.56241  

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Bacteria_edges)
summary(ModelMATmean)
# Number of obs: 266, groups:  StudyID, 47
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                         Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(Rotation_cycles) 16.7831 16.7831     1  89.836  4.8373 0.03042 *
# scale(SpeciesRichness) 20.4199 20.4199     1 200.583  5.8856 0.01615 *
# scale(log(Duration))   20.3323 20.3323     1 230.921  5.8603 0.01626 *
# scale(MATmean)          0.8786  0.8786     1  30.207  0.2532 0.61845  



############# 12. Plot
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggpmisc)

## SpeciesRichness
sum(!is.na(Bacteria_edges$SpeciesRichness)) ## n = 266
p1 <- ggplot(Bacteria_edges, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnBacteria_edges266")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Bacteria_edges$Duration)) ## n = 266
p2 <- ggplot(Bacteria_edges, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnBacteria_edges266")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Bacteria_edges$Rotation_cycles)) ## n = 266
p3 <- ggplot(Bacteria_edges, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnBacteria_edges266")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Bacteria_edges$Latitude)) ## n = 266
p5 <- ggplot(Bacteria_edges, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnBacteria_edges266")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Bacteria_edges$Longitude)) ## n = 266
p6 <- ggplot(Bacteria_edges, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnBacteria_edges266")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Bacteria_edges$MAPmean)) ## n = 266
p7 <- ggplot(Bacteria_edges, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnBacteria_edges266")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Bacteria_edges$MATmean)) ## n = 266
p8 <- ggplot(Bacteria_edges, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnBacteria_edges266")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Bacteria_edges$pHCK)) ## n = 78
p9 <- ggplot(Bacteria_edges, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnBacteria_edges78")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Bacteria_edges$SOCCK)) ## n = 75
p10 <- ggplot(Bacteria_edges, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnBacteria_edges75")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Bacteria_edges$TNCK)) ## n = 44
p11 <- ggplot(Bacteria_edges, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnBacteria_edges44")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Bacteria_edges$NO3CK)) ## n = 38
p12 <- ggplot(Bacteria_edges, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnBacteria_edges38")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Bacteria_edges$NH4CK)) ## n = 39
p13<- ggplot(Bacteria_edges, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnBacteria_edges39")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Bacteria_edges$APCK)) ## n = 56
p14 <- ggplot(Bacteria_edges, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnBacteria_edges56")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Bacteria_edges$AKCK)) ## n = 39
p15 <- ggplot(Bacteria_edges, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnBacteria_edges39")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Bacteria_edges$ANCK)) ## n = 31
p16 <- ggplot(Bacteria_edges, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnBacteria_edges31")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Bacteria_edges$RRpH)) ## n = 78
p17 <- ggplot(Bacteria_edges, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="RRpH")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRpH" , y="RR78")
p17
pdf("RRpH.pdf",width=8,height=8)
p17
dev.off() 

## RRSOC
sum(!is.na(Bacteria_edges$RRSOC)) ## n = 75
p18 <- ggplot(Bacteria_edges, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="RRSOC")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRSOC" , y="RR75")
p18
pdf("RRSOC.pdf",width=8,height=8)
p18
dev.off() 

## RRTN
sum(!is.na(Bacteria_edges$RRTN)) ## n = 44
p19 <- ggplot(Bacteria_edges, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="RRTN")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRTN" , y="RR44")
p19
pdf("RRTN.pdf",width=8,height=8)
p19
dev.off() 

## RRNO3
sum(!is.na(Bacteria_edges$RRNO3)) ## n = 38
p20 <- ggplot(Bacteria_edges, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="RRNO3")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNO3" , y="RR38")
p20
pdf("RRNO3.pdf",width=8,height=8)
p20
dev.off() 

## RRNH4
sum(!is.na(Bacteria_edges$RRNH4)) ## n = 39
p21 <- ggplot(Bacteria_edges, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="RRNH4")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRNH4" , y="RR39")
p21
pdf("RRNH4.pdf",width=8,height=8)
p21
dev.off() 

## RRAP
sum(!is.na(Bacteria_edges$RRAP)) ## n = 56
p22 <- ggplot(Bacteria_edges, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="RRAP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="RRAP" , y="RR56")
p22
pdf("RRAP.pdf",width=8,height=8)
p22
dev.off() 

## RRAK
sum(!is.na(Bacteria_edges$RRAK)) ## n = 39
p23 <- ggplot(Bacteria_edges, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="RRAK")+
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
sum(!is.na(Bacteria_edges$RRAN)) ## n = 31
p24 <- ggplot(Bacteria_edges, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnBacteria_edges", x="RRAN")+
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
sum(!is.na(Bacteria_edges$RRYield)) ## n = 47
p25 <- ggplot(Bacteria_edges, aes(x=RR, y=RRYield)) +
  geom_point(color="gray", size=10, shape=21) +
  geom_smooth(method=lm, color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") +
  theme_bw() +
  theme(text = element_text(family = "serif", size=20)) +
  geom_vline(aes(xintercept=0), colour="black", linewidth=0.5, linetype="dashed") +
  labs(x="RR", y="RRYield47") +
  theme(panel.grid=element_blank()) + 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = x ~ y,  
    parse = TRUE,
    color="black",
    size = 5, # 公式字体大小
    label.x = 0.05,  # 公式位置
    label.y = 0.85
  ) +
  stat_cor(method = "spearman", size = 5) +
  scale_x_continuous(limits=c(-0.8, 0.6), expand=c(0, 0)) + 
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



