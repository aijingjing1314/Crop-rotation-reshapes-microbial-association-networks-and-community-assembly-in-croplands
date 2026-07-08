
library(metafor)
library(boot)
library(parallel)
library(dplyr)
library(multcompView)
library(lme4)
library(MuMIn)
library(lmerTest)
Fungi_edges <- read.csv("Fungi_edges.csv", fileEncoding = "latin1")
# Check data
head(Fungi_edges)

# 1. The number of Obversation
total_number <- nrow(Fungi_edges)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 301   

# 2. The number of Study
unique_studyid_number <- length(unique(Fungi_edges$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 41 






#### 8. Subgroup analysis
### 8.1 Longitude_Sub
Fungi_edges_filteredLongitude_Sub <- subset(Fungi_edges, Longitude_Sub %in% c("LongitudeDa0", "LongitudeXy0"))
#
Fungi_edges_filteredLongitude_Sub$Longitude_Sub <- droplevels(factor(Fungi_edges_filteredLongitude_Sub$Longitude_Sub))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredLongitude_Sub %>%
  group_by(Longitude_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Longitude_Sub Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 LongitudeDa0           129             36
# 2 LongitudeXy0           172              5
overall_model_Fungi_edges_filteredLongitude_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Longitude_Sub, random = ~ 1 | StudyID, data = Fungi_edges_filteredLongitude_Sub, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredLongitude_Sub)
# Multivariate Meta-Analysis Model (k = 301; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -134.1089   268.2178   274.2178   285.3192   274.2992   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1672  0.4089     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 299) = 1665.1984, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.6943, p-val = 0.7067
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# Longitude_SubLongitudeDa0    0.0448  0.0711   0.6302  0.5286  -0.0946  0.1843    
# Longitude_SubLongitudeXy0   -0.1023  0.1877  -0.5452  0.5857  -0.4703  0.2656     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredLongitude_Sub)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredLongitude_Sub)
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
Fungi_edges_filteredMAPmean2_Sub <- subset(Fungi_edges, MAPmean2_Sub %in% c("MAPXy600", "MAP600Dao1200", "MAPDy1200"))
#
Fungi_edges_filteredMAPmean2_Sub$MAPmean2_Sub <- droplevels(factor(Fungi_edges_filteredMAPmean2_Sub$MAPmean2_Sub))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredMAPmean2_Sub %>%
  group_by(MAPmean2_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MAPmean2_Sub  Observations Unique_StudyID
#   <fct>                <int>          <int>
# 1 MAP600Dao1200          111             14
# 2 MAPDy1200              106             10
# 3 MAPXy600                84             18

overall_model_Fungi_edges_filteredMAPmean2_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MAPmean2_Sub, random = ~ 1 | StudyID, data = Fungi_edges_filteredMAPmean2_Sub, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredMAPmean2_Sub)
# Multivariate Meta-Analysis Model (k = 301; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -133.0842   266.1684   274.1684   288.9568   274.3049   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1586  0.3982     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 298) = 1435.3598, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.7222, p-val = 0.6320
# 
# Model Results:
# 
#                            estimate      se     zval    pval    ci.lb   ci.ub    
# MAPmean2_SubMAP600Dao1200   -0.0593  0.1065  -0.5569  0.5776  -0.2680  0.1494    
# MAPmean2_SubMAPDy1200        0.0021  0.1337   0.0155  0.9876  -0.2600  0.2642    
# MAPmean2_SubMAPXy600         0.1036  0.0932   1.1113  0.2664  -0.0791  0.2863  

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredMAPmean2_Sub)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredMAPmean2_Sub)
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
#           "a"                       "a"                       "a" 




### 8.3 MATmean_Sub
Fungi_edges_filteredMATmean_Sub <- subset(Fungi_edges, MATmean_Sub %in% c("MATXy8", "MAT8Dao15", "MATDy15"))
#
Fungi_edges_filteredMATmean_Sub$MATmean_Sub <- droplevels(factor(Fungi_edges_filteredMATmean_Sub$MATmean_Sub))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredMATmean_Sub %>%
  group_by(MATmean_Sub) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   MATmean_Sub Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 MAT8Dao15             39             10
# 2 MATDy15              114             14
# 3 MATXy8               148             18

overall_model_Fungi_edges_filteredMATmean_Sub <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + MATmean_Sub, random = ~ 1 | StudyID, data = Fungi_edges_filteredMATmean_Sub, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredMATmean_Sub)
# Multivariate Meta-Analysis Model (k = 301; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -133.1986   266.3972   274.3972   289.1855   274.5337   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1639  0.4049     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 298) = 1391.8467, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.8619, p-val = 0.6016
# 
# Model Results:
# 
#                       estimate      se     zval    pval    ci.lb   ci.ub    
# MATmean_SubMAT8Dao15    0.0543  0.1259   0.4313  0.6663  -0.1925  0.3011    
# MATmean_SubMATDy15     -0.0904  0.1137  -0.7953  0.4265  -0.3132  0.1324    
# MATmean_SubMATXy8       0.1020  0.0950   1.0737  0.2830  -0.0842  0.2881  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredMATmean_Sub)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredMATmean_Sub)
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
#       "a"                  "a"                  "a" 




### 8.4 LegumeNonlegume
Fungi_edges_filteredLegumeNonlegume <- subset(Fungi_edges, LegumeNonlegume %in% c("Non-legume to Non-legume", "Non-legume to Legume", "Legume to Non-legume"))
#
Fungi_edges_filteredLegumeNonlegume$LegumeNonlegume <- droplevels(factor(Fungi_edges_filteredLegumeNonlegume$LegumeNonlegume))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredLegumeNonlegume %>%
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
# 3 Non-legume to Non-legume           95             27

overall_model_Fungi_edges_filteredLegumeNonlegume <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + LegumeNonlegume, random = ~ 1 | StudyID, data = Fungi_edges_filteredLegumeNonlegume, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredLegumeNonlegume)
# Multivariate Meta-Analysis Model (k = 301; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -135.5758   271.1515   279.1515   293.9399   279.2881   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1647  0.4058     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 298) = 1466.4795, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 3.8270, p-val = 0.2808
# 
# Model Results:
# 
#                                          estimate      se     zval    pval    ci.lb   ci.ub    
# LegumeNonlegumeLegume to Non-legume       -0.0173  0.0705  -0.2452  0.8063  -0.1555  0.1209    
# LegumeNonlegumeNon-legume to Legume        0.0014  0.0683   0.0211  0.9832  -0.1324  0.1353    
# LegumeNonlegumeNon-legume to Non-legume    0.0496  0.0672   0.7381  0.4605  -0.0821  0.1812  
 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredLegumeNonlegume)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredLegumeNonlegume)
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
Fungi_edges_filteredAMnonAM <- subset(Fungi_edges, AMnonAM %in% c("AM to AM", "AM to nonAM", "nonAM to AM"))
#
Fungi_edges_filteredAMnonAM$AMnonAM <- droplevels(factor(Fungi_edges_filteredAMnonAM$AMnonAM))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredAMnonAM %>%
  group_by(AMnonAM) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   AMnonAM     Observations Unique_StudyID
#   <fct>              <int>          <int>
# 1 AM to AM             261             34
# 2 AM to nonAM           25              7
# 3 nonAM to AM           15              5

overall_model_Fungi_edges_filteredAMnonAM <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + AMnonAM, random = ~ 1 | StudyID, data = Fungi_edges_filteredAMnonAM, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredAMnonAM)
# Multivariate Meta-Analysis Model (k = 301; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -133.4241   266.8481   274.8481   289.6365   274.9847   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1491  0.3862     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 298) = 1636.9788, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 6.3469, p-val = 0.0959
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# AMnonAMAM to AM       0.0288  0.0644   0.4464  0.6553  -0.0975  0.1550    
# AMnonAMAM to nonAM   -0.0711  0.0782  -0.9090  0.3634  -0.2244  0.0822    
# AMnonAMnonAM to AM    0.1255  0.1040   1.2066  0.2276  -0.0784  0.3294
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredAMnonAM)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredAMnonAM)
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
   #       "ab"                "a"                "b" 


### 8.6 C3C4
Fungi_edges_filteredC3C4 <- subset(Fungi_edges, C3C4 %in% c("C3 to C3", "C3 to C4", "C4 to C3"))
#
Fungi_edges_filteredC3C4$C3C4 <- droplevels(factor(Fungi_edges_filteredC3C4$C3C4))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredC3C4 %>%
  group_by(C3C4) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   C3C4     Observations Unique_StudyID
#   <fct>           <int>          <int>
# 1 C3 to C3           57             14
# 2 C3 to C4           91             22
# 3 C4 to C3          148             12

overall_model_Fungi_edges_filteredC3C4 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + C3C4, random = ~ 1 | StudyID, data = Fungi_edges_filteredC3C4, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredC3C4)
# #Multivariate Meta-Analysis Model (k = 296; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -113.6338   227.2677   235.2677   249.9883   235.4065   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1471  0.3835     39     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 293) = 1343.5202, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 9.7041, p-val = 0.0213
# 
# Model Results:
# 
#               estimate      se     zval    pval    ci.lb    ci.ub    
# C3C4C3 to C3   -0.2219  0.1025  -2.1644  0.0304  -0.4228  -0.0210  * 
# C3C4C3 to C4    0.1174  0.0735   1.5984  0.1100  -0.0266   0.2615    
# C3C4C4 to C3    0.1330  0.0744   1.7874  0.0739  -0.0128   0.2788  . 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredC3C4)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredC3C4)
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
#      "a"          "b"          "b" 


### 8.7 Annual_Pere
Fungi_edges_filteredAnnual_Pere <- subset(Fungi_edges, Annual_Pere %in% c("Annual to Annual",  "Perennial to Annual"))
#
Fungi_edges_filteredAnnual_Pere$Annual_Pere <- droplevels(factor(Fungi_edges_filteredAnnual_Pere$Annual_Pere))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredAnnual_Pere %>%
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
# 3 Perennial to Annual              34              9
# 4 Perennial to Perennial            5              2

overall_model_Fungi_edges_filteredAnnual_Pere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Annual_Pere, random = ~ 1 | StudyID, data = Fungi_edges_filteredAnnual_Pere, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredAnnual_Pere)
# ate Meta-Analysis Model (k = 287; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -112.4544   224.9088   230.9088   241.8663   230.9942   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1750  0.4183     39     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 285) = 1578.8307, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.2852, p-val = 0.8671
# 
# Model Results:
# 
#                                 estimate      se     zval    pval    ci.lb   ci.ub    
# Annual_PereAnnual to Annual       0.0412  0.0797   0.5175  0.6048  -0.1150  0.1974    
# Annual_PerePerennial to Annual   -0.0197  0.1490  -0.1319  0.8950  -0.3118  0.2724  
 

# # Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredAnnual_Pere)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredAnnual_Pere)
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
      # Annual_PereAnnual to Annual    Annual_PereAnnual to Perennial    Annual_PerePerennial to Annual Annual_PerePerennial to Perennial 
      #                   "a"                               "a"                               "a"                               "a"" 


### 8.8 Plant_stage
Fungi_edges_filteredPlant_stage <- subset(Fungi_edges, Plant_stage %in% c("Vegetative stage","Reproductive stage", "Harvest"))
#
Fungi_edges_filteredPlant_stage$Plant_stage <- droplevels(factor(Fungi_edges_filteredPlant_stage$Plant_stage))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredPlant_stage %>%
  group_by(Plant_stage) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Plant_stage        Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 Harvest                     117             21
# 2 Maturity stage                9              5
# 3 Reproductive stage           52              9
# 4 Vegetative stage             93              9

overall_model_Fungi_edges_filteredPlant_stage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Plant_stage, random = ~ 1 | StudyID, data = Fungi_edges_filteredPlant_stage, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredPlant_stage)
# # Multel (k = 262; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -82.1756  164.3512  172.3512  186.5785  172.5087   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1558  0.3948     30     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 259) = 1283.8542, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.5854, p-val = 0.6627
# 
# Model Results:
# 
#                                estimate      se     zval    pval    ci.lb   ci.ub    
# Plant_stageHarvest              -0.0186  0.0756  -0.2454  0.8061  -0.1668  0.1297    
# Plant_stageReproductive stage   -0.0437  0.0774  -0.5645  0.5724  -0.1954  0.1080    
# Plant_stageVegetative stage     -0.0323  0.0762  -0.4243  0.6714  -0.1817  0.1170     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredPlant_stage)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredPlant_stage)
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
Fungi_edges_filteredRotation_cycles2 <- subset(Fungi_edges, Rotation_cycles2 %in% c("D1", "D1-3", "D3-5", "D5-10", "D10"))
#
Fungi_edges_filteredRotation_cycles2$Rotation_cycles2 <- droplevels(factor(Fungi_edges_filteredRotation_cycles2$Rotation_cycles2))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredRotation_cycles2 %>%
  group_by(Rotation_cycles2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Rotation_cycles2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 D1                         68             15
# 2 D1-3                       99             15
# 3 D10                        18              6
# 4 D3-5                       80              7
# 5 D5-10                      36             11
overall_model_Fungi_edges_filteredRotation_cycles2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Rotation_cycles2, random = ~ 1 | StudyID, data = Fungi_edges_filteredRotation_cycles2, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredRotation_cycles2)
# Multivariate Meta-Analysis Model (k = 301; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -130.8103   261.6205   273.6205   295.7627   273.9112   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.2305  0.4801     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 296) = 1503.2615, p-val < .0001
# 
# Test of Moderators (coefficients 1:5):
# QM(df = 5) = 16.4570, p-val = 0.0057
# 
# Model Results:
# 
#                        estimate      se     zval    pval    ci.lb   ci.ub    
# Rotation_cycles2D1      -0.0762  0.0924  -0.8243  0.4098  -0.2573  0.1049    
# Rotation_cycles2D1-3    -0.0967  0.0922  -1.0485  0.2944  -0.2775  0.0841    
# Rotation_cycles2D10      0.2601  0.1161   2.2403  0.0251   0.0325  0.4877  * 
# Rotation_cycles2D3-5     0.2602  0.1144   2.2738  0.0230   0.0359  0.4845  * 
# Rotation_cycles2D5-10    0.1426  0.1135   1.2564  0.2090  -0.0798  0.3650     

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredRotation_cycles2)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredRotation_cycles2)
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
#         "a"                   "a"                   "b"                   "b"                   "a" 


### 8.10 Duration2
Fungi_edges_filteredDuration2 <- subset(Fungi_edges, Duration2 %in% c("D1-5", "D6-10", "D11-20", "D20-40"))
#
Fungi_edges_filteredDuration2$Duration2 <- droplevels(factor(Fungi_edges_filteredDuration2$Duration2))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredDuration2 %>%
  group_by(Duration2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Duration2 Observations Unique_StudyID
#   <fct>            <int>          <int>
# 1 D1-5               145             22
# 2 D11-20              35             10
# 3 D20-40              87              6
# 4 D6-10               34             10
overall_model_Fungi_edges_filteredDuration2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Duration2, random = ~ 1 | StudyID, data = Fungi_edges_filteredDuration2, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredDuration2)
# Multivariate Meta-Analysis Model (k = 301; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -130.1547   260.3093   270.3093   288.7780   270.5155   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.2245  0.4738     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 297) = 1580.5678, p-val < .0001
# 
# Test of Moderators (coefficients 1:4):
# QM(df = 4) = 11.3625, p-val = 0.0228
# 
# Model Results:
# 
#                  estimate      se     zval    pval    ci.lb   ci.ub    
# Duration2D1-5     -0.1132  0.0978  -1.1578  0.2470  -0.3049  0.0784    
# Duration2D11-20    0.1498  0.1051   1.4246  0.1543  -0.0563  0.3558    
# Duration2D20-40    0.0734  0.1266   0.5797  0.5621  -0.1748  0.3216    
# Duration2D6-10     0.2552  0.1067   2.3929  0.0167   0.0462  0.4643  * 


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredDuration2)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredDuration2)
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
 #       "a"             "b"           "abc"             "c" 

### 8.11 SpeciesRichness2
Fungi_edges_filteredSpeciesRichness2 <- subset(Fungi_edges, SpeciesRichness2 %in% c("R2", "R3"))
#
Fungi_edges_filteredSpeciesRichness2$SpeciesRichness2 <- droplevels(factor(Fungi_edges_filteredSpeciesRichness2$SpeciesRichness2))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredSpeciesRichness2 %>%
  group_by(SpeciesRichness2) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   SpeciesRichness2 Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 R2                        276             39
# 2 R3                         24              4
# 3 R4                          1              1

overall_model_Fungi_edges_filteredSpeciesRichness2 <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + SpeciesRichness2, random = ~ 1 | StudyID, data = Fungi_edges_filteredSpeciesRichness2, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredSpeciesRichness2)
# Analysis Model (k = 300; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -132.9017   265.8034   271.8034   282.8947   271.8851   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1467  0.3831     40     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 298) = 1661.4655, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 1.0013, p-val = 0.6061
# 
# Model Results:
# 
#                     estimate      se    zval    pval    ci.lb   ci.ub    
# SpeciesRichness2R2    0.0472  0.0633  0.7450  0.4563  -0.0769  0.1712    
# SpeciesRichness2R3    0.0820  0.0820  1.0006  0.3170  -0.0787  0.2427   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredSpeciesRichness2)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredSpeciesRichness2)
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
#       "a"                "a"             



### 8.12 Bulk_Rhizosphere
Fungi_edges_filteredBulk_Rhizosphere <- subset(Fungi_edges, Bulk_Rhizosphere %in% c("Non_Rhizosphere", "Rhizosphere"))
#
Fungi_edges_filteredBulk_Rhizosphere$Bulk_Rhizosphere <- droplevels(factor(Fungi_edges_filteredBulk_Rhizosphere$Bulk_Rhizosphere))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredBulk_Rhizosphere %>%
  group_by(Bulk_Rhizosphere) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Bulk_Rhizosphere Observations Unique_StudyID
#   <fct>                   <int>          <int>
# 1 Non_Rhizosphere           249             30
# 2 Rhizosphere                52             19

overall_model_Fungi_edges_filteredBulk_Rhizosphere <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Bulk_Rhizosphere, random = ~ 1 | StudyID, data = Fungi_edges_filteredBulk_Rhizosphere, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredBulk_Rhizosphere)
# Multivariate Meta-Analysis Model (k = 301; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -128.6845   257.3690   263.3690   274.4703   263.4504   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1658  0.4072     41     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 299) = 1678.1348, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 15.8479, p-val = 0.0004
# 
# Model Results:
# 
#                                  estimate      se     zval    pval    ci.lb   ci.ub    
# Bulk_RhizosphereNon_Rhizosphere   -0.0151  0.0671  -0.2257  0.8215  -0.1466  0.1163    
# Bulk_RhizosphereRhizosphere        0.1058  0.0692   1.5279  0.1265  -0.0299  0.2414 

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredBulk_Rhizosphere)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredBulk_Rhizosphere)
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
Fungi_edges_filteredSoil_texture <- subset(Fungi_edges, Soil_texture %in% c("Fine", "Medium", "Coarse"))
#
Fungi_edges_filteredSoil_texture$Soil_texture <- droplevels(factor(Fungi_edges_filteredSoil_texture$Soil_texture))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredSoil_texture %>%
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
# 3 Medium                128             10

overall_model_Fungi_edges_filteredSoil_texture <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Soil_texture, random = ~ 1 | StudyID, data = Fungi_edges_filteredSoil_texture, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredSoil_texture)
# Multivariate Meta-Analysis Model (k = 246; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -119.6259   239.2519   247.2519   261.2241   247.4199   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1846  0.4296     23     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 243) = 1382.1146, p-val < .0001
# 
# Test of Moderators (coefficients 1:3):
# QM(df = 3) = 1.9110, p-val = 0.5911
# 
# Model Results:
# 
#                     estimate      se     zval    pval    ci.lb   ci.ub    
# Soil_textureCoarse   -0.2792  0.2198  -1.2706  0.2039  -0.7099  0.1515    
# Soil_textureFine      0.0649  0.1499   0.4327  0.6652  -0.2290  0.3587    
# Soil_textureMedium    0.0459  0.1388   0.3308  0.7408  -0.2261  0.3179  


# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredSoil_texture)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredSoil_texture)
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
#         "a"                "a"                "a" 


### 8.14 Tillage
Fungi_edges_filteredTillage <- subset(Fungi_edges, Tillage %in% c("Tillage", "No_tillage"))
#
Fungi_edges_filteredTillage$Tillage <- droplevels(factor(Fungi_edges_filteredTillage$Tillage))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredTillage %>%
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

overall_model_Fungi_edges_filteredTillage <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Tillage, random = ~ 1 | StudyID, data = Fungi_edges_filteredTillage, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredTillage)
# Multivariate Meta-Analysis Model (k = 193; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  23.1431  -46.2862  -40.2862  -30.5294  -40.1579   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1860  0.4313     10     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 191) = 521.5903, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.9972, p-val = 0.6074
# 
# Model Results:
# 
#                    estimate      se     zval    pval    ci.lb   ci.ub    
# TillageNo_tillage   -0.1449  0.1470  -0.9861  0.3241  -0.4329  0.1431    
# TillageTillage      -0.1267  0.1462  -0.8670  0.3859  -0.4132  0.1597    
#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredTillage)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredTillage)
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
Fungi_edges_filteredStraw_retention <- subset(Fungi_edges, Straw_retention %in% c("Retention", "No_retention"))
#
Fungi_edges_filteredStraw_retention$Straw_retention <- droplevels(factor(Fungi_edges_filteredStraw_retention$Straw_retention))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredStraw_retention %>%
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

overall_model_Fungi_edges_filteredStraw_retention <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Straw_retention, random = ~ 1 | StudyID, data = Fungi_edges_filteredStraw_retention, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredStraw_retention)
# Multivariate Meta-Analysis Model (k = 43; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -59.3360  118.6720  124.6720  129.8127  125.3207   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1758  0.4193     10     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 41) = 280.3471, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.7555, p-val = 0.6854
# 
# Model Results:
# 
#                              estimate      se     zval    pval    ci.lb   ci.ub    
# Straw_retentionNo_retention   -0.0766  0.1405  -0.5454  0.5855  -0.3519  0.1987    
# Straw_retentionRetention      -0.1089  0.1381  -0.7885  0.4304  -0.3796  0.1618   

# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredStraw_retention)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredStraw_retention)
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
Fungi_edges_filteredPrimer <- subset(Fungi_edges, Primer %in% c("ITS1", "ITS2"))
#
Fungi_edges_filteredPrimer$Primer <- droplevels(factor(Fungi_edges_filteredPrimer$Primer))
# The number of Observations and StudyID
group_summary <- Fungi_edges_filteredPrimer %>%
  group_by(Primer) %>%
  summarise(
    Observations = n(),                   
    Unique_StudyID = n_distinct(StudyID)  
  )
print(group_summary)
#   Primer             Observations Unique_StudyID
#   <fct>                     <int>          <int>
# 1 ITS1                        217             38
# 2 ITS1 + 5.8S + ITS2            4              1
# 3 ITS2                         80              2
overall_model_Fungi_edges_filteredPrimer <- rma.mv(yi = RR, V = Vi, mods = ~ 0 + Primer, random = ~ 1 | StudyID, data = Fungi_edges_filteredPrimer, method = "REML")
# QM and p value
summary(overall_model_Fungi_edges_filteredPrimer)
# Multivariate Meta-Analysis Model (k = 297; method: REML)
# 
#    logLik   Deviance        AIC        BIC       AICc   
# -131.7762   263.5525   269.5525   280.6134   269.6349   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed   factor 
# sigma^2    0.1748  0.4181     40     no  StudyID 
# 
# Test for Residual Heterogeneity:
# QE(df = 295) = 1377.9545, p-val < .0001
# 
# Test of Moderators (coefficients 1:2):
# QM(df = 2) = 0.2081, p-val = 0.9012
# 
# Model Results:
# 
#             estimate      se     zval    pval    ci.lb   ci.ub    
# PrimerITS1    0.0319  0.0707   0.4511  0.6519  -0.1066  0.1704    
# PrimerITS2   -0.0202  0.2964  -0.0680  0.9458  -0.6012  0.5608    

#  
# Extract model coefficients and covariance matrix
coef_rotation <- coef(overall_model_Fungi_edges_filteredPrimer)
vcov_rotation <- vcov(overall_model_Fungi_edges_filteredPrimer)
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
              # PrimerITS1               PrimerITS2 
              #        "a"                  "a" 









#### 9. Linear Mixed Effect Model
# 
Fungi_edges$Wr <- 1 / Fungi_edges$Vi
# Model selection
Model1 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Fungi_edges)
Model2 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Fungi_edges)
Model3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + (1 | StudyID), weights = Wr, data = Fungi_edges)
Model4 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Fungi_edges)
Model5 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Fungi_edges)
Model6 <- lmer(RR ~ scale(Rotation_cycles) + scale(log(SpeciesRichness)) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Fungi_edges)
Model7 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(log(Duration)) + (1 | StudyID), weights = Wr, data = Fungi_edges)
Model8 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(Duration) + (1 | StudyID), weights = Wr, data = Fungi_edges)
Model9 <- lmer(RR ~ scale(Rotation_cycles) + scale(SpeciesRichness) + scale(Duration) + 
                 scale(Rotation_cycles) * scale(SpeciesRichness) * scale(Duration) + 
                 (1 | StudyID), weights = Wr, data = Fungi_edges)
Model10 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(log(SpeciesRichness)) + scale(log(Duration)) + 
                  scale(log(Rotation_cycles)) * scale(log(SpeciesRichness)) * scale(log(Duration)) + 
                  (1 | StudyID), weights = Wr, data = Fungi_edges)

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
#           Model       AIC      BIC    logLik
# Model1   Model1  98.31353 120.5562 -43.15676
# Model2   Model2  98.47983 120.7225 -43.23991
# Model3   Model3  97.29299 119.5357 -42.64650###############################################
# Model4   Model4  98.45146 120.6941 -43.22573
# Model5   Model5  99.39894 121.6416 -43.69947
# Model6   Model6  99.53365 121.7763 -43.76683
# Model7   Model7  98.34173 120.5844 -43.17086
# Model8   Model8  97.43526 119.6779 -42.71763
# Model9   Model9 102.92846 139.9996 -41.46423
# Model10 Model10 108.53809 145.6092 -44.26904

##### Model 3 is the best model
summary(Model3)
# Number of obs: 301, groups:  StudyID, 41
anova(Model3)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.3671  0.3671     1 287.618  0.1255 0.7234
# scale(SpeciesRichness)      1.5468  1.5468     1 296.988  0.5288 0.4677
# scale(Duration)             3.8280  3.8280     1  68.083  1.3086 0.2567

#### 10.1. ModelpH
ModelpH <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(pHCK) + (1 | StudyID), weights = Wr, data = Fungi_edges)
summary(ModelpH)
# Number of obs: 90, groups:  StudyID, 27
anova(ModelpH) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.1558  0.1558     1 70.999  0.0984 0.7547
# scale(SpeciesRichness)      0.4642  0.4642     1 58.920  0.2931 0.5903
# scale(Duration)             0.0238  0.0238     1 79.535  0.0150 0.9028
# scale(pHCK)                 3.7714  3.7714     1 84.980  2.3811 0.1265  

#### 10.2. ModelSOC
ModelSOC <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(SOCCK) + (1 | StudyID), weights = Wr, data = Fungi_edges)
summary(ModelSOC)
# Number of obs: 101, groups:  StudyID, 28
anova(ModelSOC) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.10226 0.10226     1 77.705  0.0716 0.7897
# scale(SpeciesRichness)      0.07298 0.07298     1 78.308  0.0511 0.8217
# scale(Duration)             0.48845 0.48845     1 94.004  0.3420 0.5601
# scale(SOCCK)                0.29497 0.29497     1 65.895  0.2065 0.6510  

#### 10.3. ModelTN
ModelTN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(TNCK) + (1 | StudyID), weights = Wr, data = Fungi_edges)
summary(ModelTN)
# Number of obs: 63, groups:  StudyID, 17
anova(ModelTN) 
# > anova(ModelTN) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                               Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.026663 0.026663     1 44.446  0.0284 0.8670
# scale(SpeciesRichness)      0.058129 0.058129     1 57.846  0.0618 0.8045
# scale(Duration)             0.303228 0.303228     1 52.160  0.3226 0.5725
# scale(TNCK)                 0.250924 0.250924     1 15.254  0.2670 0.6128


#### 10.4. ModelNO3
ModelNO3 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(NO3CK) + (1 | StudyID), weights = Wr, data = Fungi_edges)
summary(ModelNO3)
# Number of obs: 29, groups:  StudyID, 12
anova(ModelNO3) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles))  0.6999  0.6999     1 10.6191  0.4471 0.517967   
# scale(SpeciesRichness)      18.6636 18.6636     1 14.4359 11.9238 0.003729 **
# scale(Duration)              0.0048  0.0048     1 22.8233  0.0031 0.956401   
# scale(NO3CK)                 0.1198  0.1198     1  7.8449  0.0765 0.789236   

#### 10.5. ModelNH4
ModelNH4 <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(NH4CK) + (1 | StudyID), weights = Wr, data = Fungi_edges)
summary(ModelNH4)
# Number of obs: 29, groups:  StudyID, 12
anova(ModelNH4) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(log(Rotation_cycles))  0.7615  0.7615     1 10.1406  0.4820 0.503128   
# scale(SpeciesRichness)      18.8099 18.8099     1 13.9778 11.9058 0.003909 **
# scale(Duration)              0.0090  0.0090     1 23.0032  0.0057 0.940599   
# scale(NH4CK)                 0.1901  0.1901     1  6.8753  0.1203 0.739067 

#### 10.6. ModelAP
ModelAP <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(APCK) + (1 | StudyID), weights = Wr, data = Fungi_edges)
summary(ModelAP)
# Number of obs: 57, groups:  StudyID, 23
anova(ModelAP) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
# scale(log(Rotation_cycles)) 9.7750  9.7750     1 25.327  5.2182 0.03099 *
# scale(SpeciesRichness)      1.4900  1.4900     1 45.765  0.7954 0.37714  
# scale(Duration)             3.8641  3.8641     1 51.705  2.0628 0.15696  
# scale(APCK)                 2.0029  2.0029     1 39.693  1.0692 0.30739  

#### 10.7. ModelAK
ModelAK <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(AKCK) + (1 | StudyID), weights = Wr, data = Fungi_edges)
summary(ModelAK)
# Number of obs: 43, groups:  StudyID, 18
anova(ModelAK)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 3.2540  3.2540     1 17.175  1.6979 0.2098
# scale(SpeciesRichness)      1.6951  1.6951     1 37.546  0.8844 0.3530
# scale(Duration)             0.2640  0.2640     1 15.223  0.1377 0.7157
# scale(AKCK)                 0.3806  0.3806     1 37.493  0.1986 0.6584

#### 10.8. ModelAN
ModelAN <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(ANCK) + (1 | StudyID), weights = Wr, data = Fungi_edges)
summary(ModelAN)
# Number of obs: 42, groups:  StudyID, 12
anova(ModelAN)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 5.7575  5.7575     1 12.4852  2.7804 0.1203
# scale(SpeciesRichness)      3.6307  3.6307     1 13.5901  1.7533 0.2073
# scale(Duration)             1.4182  1.4182     1 10.7398  0.6849 0.4259
# scale(ANCK)                 0.0020  0.0020     1  6.8797  0.0010 0.9762

#### 11. Latitude, Longitude
### 11.1. Latitude
ModelLatitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(Latitude) + (1 | StudyID), weights = Wr, data = Fungi_edges)
summary(ModelLatitude)
# Number of obs: 301, groups:  StudyID, 41
anova(ModelLatitude) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.3803  0.3803     1 287.902  0.1298 0.7189
# scale(SpeciesRichness)      1.7969  1.7969     1 295.937  0.6135 0.4341
# scale(Duration)             4.6818  4.6818     1  64.258  1.5984 0.2107
# scale(Latitude)             3.2616  3.2616     1  22.932  1.1136 0.3023

### 11.2. Longitude
ModelLongitude <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(Longitude) + (1 | StudyID), weights = Wr, data = Fungi_edges)
summary(ModelLongitude)
# Number of obs: 301, groups:  StudyID, 41
anova(ModelLongitude)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.48107 0.48107     1 281.576  0.1652 0.6847
# scale(SpeciesRichness)      1.28650 1.28650     1 295.989  0.4417 0.5068
# scale(Duration)             2.25468 2.25468     1  83.780  0.7742 0.3814
# scale(Longitude)            1.31167 1.31167     1  23.176  0.4504 0.5088

# ### 11.3. MAPmean
ModelMAPmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(MAPmean) + (1 | StudyID), weights = Wr, data = Fungi_edges)
summary(ModelMAPmean)
# Number of obs: 301, groups:  StudyID, 41
anova(ModelMAPmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.3382  0.3382     1 286.656  0.1157 0.7340
# scale(SpeciesRichness)      1.5159  1.5159     1 296.000  0.5187 0.4720
# scale(Duration)             3.9321  3.9321     1  67.177  1.3454 0.2502
# scale(MAPmean)              0.2815  0.2815     1  26.546  0.0963 0.7587

### 11.4. MATmean
ModelMATmean <- lmer(RR ~ scale(log(Rotation_cycles)) + scale(SpeciesRichness) + scale(Duration) + scale(MATmean) + (1 | StudyID), weights = Wr, data = Fungi_edges)
summary(ModelMATmean)
# Number of obs: 301, groups:  StudyID, 41
anova(ModelMATmean) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(log(Rotation_cycles)) 0.4074  0.4074     1 287.865  0.1392 0.7094
# scale(SpeciesRichness)      1.8309  1.8309     1 295.957  0.6254 0.4297
# scale(Duration)             4.5664  4.5664     1  65.068  1.5599 0.2162
# scale(MATmean)              3.6929  3.6929     1  21.985  1.2615 0.2735
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
sum(!is.na(Fungi_edges$SpeciesRichness)) ## n = 301
p1 <- ggplot(Fungi_edges, aes(y=RR, x=SpeciesRichness)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="SpeciesRichness")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SpeciesRichness" , y="lnFungi_edges301")
p1
pdf("SpeciesRichness.pdf",width=8,height=8)
p1
dev.off() 

## Duration
sum(!is.na(Fungi_edges$Duration)) ## n = 301
p2 <- ggplot(Fungi_edges, aes(y=RR, x=Duration)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="Duration")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Duration" , y="lnFungi_edges301")
p2
pdf("Duration.pdf",width=8,height=8)
p2
dev.off() 

## Rotation_cycles
sum(!is.na(Fungi_edges$Rotation_cycles)) ## n = 301
p3 <- ggplot(Fungi_edges, aes(y=RR, x=Rotation_cycles)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="Rotation_cycles")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Rotation_cycles" , y="lnFungi_edges301")
p3
pdf("Rotation_cycles.pdf",width=8,height=8)
p3
dev.off() 

## Latitude
sum(!is.na(Fungi_edges$Latitude)) ## n = 301
p5 <- ggplot(Fungi_edges, aes(y=RR, x=Latitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="Latitude")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Latitude" , y="lnFungi_edges301")
p5
pdf("Latitude.pdf",width=8,height=8)
p5
dev.off() 

## Longitude
sum(!is.na(Fungi_edges$Longitude)) ## n = 301
p6 <- ggplot(Fungi_edges, aes(y=RR, x=Longitude)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="Longitude")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="Longitude" , y="lnFungi_edges301")
p6
pdf("Longitude.pdf",width=8,height=8)
p6
dev.off() 


## MAPmean
sum(!is.na(Fungi_edges$MAPmean)) ## n = 301
p7 <- ggplot(Fungi_edges, aes(y=RR, x=MAPmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="MAPmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MAPmean" , y="lnFungi_edges301")
p7
pdf("MAPmean.pdf",width=8,height=8)
p7
dev.off() 

## MATmean
sum(!is.na(Fungi_edges$MATmean)) ## n = 301
p8 <- ggplot(Fungi_edges, aes(y=RR, x=MATmean)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="MATmean")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="MATmean" , y="lnFungi_edges301")
p8
pdf("MATmean.pdf",width=8,height=8)
p8
dev.off() 


## pHCK
sum(!is.na(Fungi_edges$pHCK)) ## n = 90
p9 <- ggplot(Fungi_edges, aes(y=RR, x=pHCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="pHCK" , y="lnFungi_edges90")
p9
pdf("pHCK.pdf",width=8,height=8)
p9
dev.off() 

## SOCCK
sum(!is.na(Fungi_edges$SOCCK)) ## n = 101
p10 <- ggplot(Fungi_edges, aes(y=RR, x=SOCCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="SOCCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="SOCCK" , y="lnFungi_edges101")
p10
pdf("SOCCK.pdf",width=8,height=8)
p10
dev.off() 

## TNCK
sum(!is.na(Fungi_edges$TNCK)) ## n = 63
p11 <- ggplot(Fungi_edges, aes(y=RR, x=TNCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="TNCK" , y="lnFungi_edges63")
p11
pdf("TNCK.pdf",width=8,height=8)
p11
dev.off() 

## NO3CK
sum(!is.na(Fungi_edges$NO3CK)) ## n = 29
p12 <- ggplot(Fungi_edges, aes(y=RR, x=NO3CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="NO3CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NO3CK" , y="lnFungi_edges29")
p12
pdf("NO3CK.pdf",width=8,height=8)
p12
dev.off() 

## NH4CK
sum(!is.na(Fungi_edges$NH4CK)) ## n = 29
p13<- ggplot(Fungi_edges, aes(y=RR, x=NH4CK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="NH4CK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="NH4CK" , y="lnFungi_edges29")
p13
pdf("NH4CK.pdf",width=8,height=8)
p13
dev.off() 

## APCK
sum(!is.na(Fungi_edges$APCK)) ## n = 57
p14 <- ggplot(Fungi_edges, aes(y=RR, x=APCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="APCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="APCK" , y="lnFungi_edges57")
p14
pdf("APCK.pdf",width=8,height=8)
p14
dev.off() 

## AKCK
sum(!is.na(Fungi_edges$AKCK)) ## n = 43
p15 <- ggplot(Fungi_edges, aes(y=RR, x=AKCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="AKCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="AKCK" , y="lnFungi_edges43")
p15
pdf("AKCK.pdf",width=8,height=8)
p15
dev.off() 

## ANCK
sum(!is.na(Fungi_edges$ANCK)) ## n = 42
p16 <- ggplot(Fungi_edges, aes(y=RR, x=ANCK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="ANCK")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "spearman", size = 5) +
  labs(x="ANCK" , y="lnFungi_edges42")
p16
pdf("ANCK.pdf",width=8,height=8)
p16
dev.off() 

## RRpH
sum(!is.na(Fungi_edges$RRpH)) ## n = 90
p17 <- ggplot(Fungi_edges, aes(y=RR, x=RRpH)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="RRpH")+
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
sum(!is.na(Fungi_edges$RRSOC)) ## n = 101
p18 <- ggplot(Fungi_edges, aes(y=RR, x=RRSOC)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="RRSOC")+
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
sum(!is.na(Fungi_edges$RRTN)) ## n = 63
p19 <- ggplot(Fungi_edges, aes(y=RR, x=RRTN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="RRTN")+
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
sum(!is.na(Fungi_edges$RRNO3)) ## n = 29
p20 <- ggplot(Fungi_edges, aes(y=RR, x=RRNO3)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="RRNO3")+
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
sum(!is.na(Fungi_edges$RRNH4)) ## n = 29
p21 <- ggplot(Fungi_edges, aes(y=RR, x=RRNH4)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="RRNH4")+
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
sum(!is.na(Fungi_edges$RRAP)) ## n = 57
p22 <- ggplot(Fungi_edges, aes(y=RR, x=RRAP)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="RRAP")+
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
sum(!is.na(Fungi_edges$RRAK)) ## n = 43
p23 <- ggplot(Fungi_edges, aes(y=RR, x=RRAK)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="RRAK")+
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
sum(!is.na(Fungi_edges$RRAN)) ## n = 42
p24 <- ggplot(Fungi_edges, aes(y=RR, x=RRAN)) +
  geom_point(color="gray", size=10, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnFungi_edges", x="RRAN")+
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
sum(!is.na(Fungi_edges$RRYield)) ## n = 51
p25 <- ggplot(Fungi_edges, aes(x=RR, y=RRYield)) +
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
  scale_x_continuous(limits=c(-1.42, 1.25), expand=c(0, 0)) + 
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






