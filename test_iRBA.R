install.packages("randomForestSRC", "ranger", "feather","parallel")
## After VSLRelief - AML acute GVHD
library(randomForestSRC)
library(ranger)
library(feather)
library(parallel)
options(rf.cores=detectCores(), mc.cores=detectCores())

G = read.csv('AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75.csv', row.names = 1, stringsAsFactors = F)

G400 <- G[, 1:401]

rfsrc_obj400 <- rfsrc(agvhi24 ~ ., data = G400, ntree=300, na.action='na.impute', nimpute = 50,
                      importance=TRUE, membership=TRUE)
print(rfsrc_obj400)
#                          Sample size: 331
#            Frequency of class labels: 146, 185
#                      Number of trees: 300
#            Forest terminal node size: 1
#        Average no. of terminal nodes: 51.48
# No. of variables tried at each split: 20
#               Total no. of variables: 400
#        Resampling used to grow trees: swr
#     Resample size used to grow trees: 331
#                             Analysis: RF-C
#                               Family: class
#                       Splitting rule: gini *random*
#        Number of random split points: 10
#               Normalized brier score: 71.01
#                                  AUC: 96.38
#                           Error rate: 0.16, 0.34, 0.01
#
# Confusion matrix:
#
#           predicted
#   observed aGVHD Non-GVHD class.error
#   aGVHD       96       50      0.3425
#   Non-GVHD     3      182      0.0162
#
#         Overall error rate: 15.71%

save(rfsrc_obj400, file='AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rfSRC_model_importance_member_top400.RData')

####

BMT_cases <-read_feather('AML_patient_donor_autosome_SNP_mat_TXoutcome.feather')
BMT_cases$BMT_case <- as.integer(BMT_cases$BMT_case)
HLA_table <- read.csv('../GWAS_cleaned_BMT_HLA_typing.csv')
meta_info <- read.csv('../GWAS_cleaned_BMT_cases_info.csv')

convert_ind <- sapply(1:dim(BMT_cases)[1], function(x) which(HLA_table$bmt_case == BMT_cases$BMT_case[x]))

G400_hla <- cbind(G400, HLA_table[convert_ind, -c(1,2,3)])
rfsrc_obj400_hla <- rfsrc(agvhi24 ~ ., data = G400_hla, ntree=300, na.action='na.impute', nimpute = 50,
                          importance=TRUE, membership=TRUE)
print(rfsrc_obj400_hla)
save(rfsrc_obj400_hla, file='AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rfSRC_model_importance_member_wHLAtyping_top400.RData')

rfsrc_oo_sbsample400 <- subsample(rfsrc_obj400, joint=TRUE, bootstrap=TRUE)
save(rfsrc_oo_sbsample400, file = 'AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rfSRC_subsample_importance_dobule_bootstraped_top400.RData')

rfsrc_oo_sbsample400_2 <- subsample(rfsrc_obj400, joint=TRUE)
save(rfsrc_oo_sbsample400_2, file = 'AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rfSRC_subsample_importance_subsample_top400.RData')

obj400_subtree <- max.subtree(rfsrc_obj400, max.order =2, sub.order=T)
print(round(obj400_subtree$order, 3))
save(obj400_subtree, file='AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rfSRC_maxSubtree_top400.RData')


rfsrc_oo_sbsample400_wHLA <- subsample(rfsrc_obj400_hla, joint=TRUE, bootstrap=TRUE)
save(rfsrc_oo_sbsample400_wHLA, file = 'AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rfSRC_subsample_importance_dobule_bootstraped_wHLA_top400.RData')

rfsrc_oo_sbsample400_2_wHLA <- subsample(rfsrc_obj400_hla, joint=TRUE)
save(rfsrc_oo_sbsample400_2_wHLA, file = 'AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rfSRC_subsample_importance_subsample_wHLA_top400.RData')

obj400_subtree_wHLA <- max.subtree(rfsrc_obj400_hla, max.order =2, sub.order=T)
print(round(obj400_subtree_wHLA$order, 3))
save(obj400_subtree_wHLA, file='AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rfSRC_maxSubtree_wHLA_top400.RData')


# with sex
convert_ind2 <- sapply(1:dim(BMT_cases)[1], function(x) which(meta_info$bmt_case == BMT_cases$BMT_case[x]))

G400_hla_sex <- cbind(G400_hla, sexmatch=meta_info$sexmatch[convert_ind2])
rfsrc_obj400_hla_sex_random <- rfsrc(agvhi24 ~ ., data = G400_hla_sex, ntree=300, na.action='na.impute', nimpute = 50,
                                     importance='random', membership=TRUE)
print(rfsrc_obj400_hla_sex_random)
save(rfsrc_obj400_hla_sex_random, file='AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rfSRC_model_importance_member_wHLAtyping_wSexmm_top400_random.RData')

rfsrc_oo_sbsample400_wHLA_wSex <- subsample(rfsrc_obj400_hla_sex_random, joint=TRUE, bootstrap=TRUE)
save(rfsrc_oo_sbsample400_wHLA_wSex, file = 'AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rfSRC_subsample_importance_dobule_bootstraped_wHLA_wSexmm_top400.RData')

rfsrc_oo_sbsample400_2_wHLA_wSex <- subsample(rfsrc_obj400_hla_sex_random, joint=TRUE)
save(rfsrc_oo_sbsample400_2_wHLA_wSex, file = 'AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rfSRC_subsample_importance_subsample_wHLA_wSexmm_top400.RData')

obj400_subtree_wHLA_wSex <- max.subtree(rfsrc_obj400_hla_sex_random, max.order =2, sub.order=T)
print(round(obj400_subtree_wHLA_wSex$order, 3))
save(obj400_subtree_wHLA_wSex, file='AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rfSRC_maxSubtree_wHLA_wSexmm_top400.RData')


#### ranodm 400
source('utils/data.utils.R')
source('utils/rfsrc.utils.R')
all_index  <- 2:dim(G)[2]
rep_times <- 1000
random_performance <- data.frame(brierS = double(rep_times),
                                 aucS = double(rep_times),
                                 overall.err.rate = double(rep_times))
for(id in 1:rep_times){
  G400_random <- G[, c(1, sample(all_index, 400))]
  rfsrc_obj_random400 <- rfsrc(agvhi24 ~ ., data = G400_random, ntree=300, na.action='na.impute', nimpute = 50,
                               importance=TRUE, membership=TRUE)
  x <- rfsrc_obj_random400
  random_performance$brierS[id] <- round(brier(x$yvar, if(!is.null(x$predicted.oob) && !all(is.na(x$predicted.oob))) x$predicted.oob else x$predicted) * 100, 2)
  random_performance$aucS[id] <- round(auc(x$yvar, if(!is.null(x$predicted.oob) && !all(is.na(x$predicted.oob))) x$predicted.oob else x$predicted) * 100, 2)
  err.rate <- cbind(x$err.rate)
  random_performance$overall.err.rate[id] <- round(100 * err.rate[nrow(err.rate), 1], 2) # paste(round(100 * err.rate[nrow(err.rate), 1], 2), "%", sep = "")

}

colMeans(random_performance)

save(random_performance, file = 'AML_aGVHD24_VLSRelief_52225_MaxInd5_pct7_random400feat_performance.RData') # random_performance


#
source('~/tools/rfsrc_tools/data.utils.R')
source('~/tools/rfsrc_tools/rfsrc.utils.R')
G400 <- G[, 1:401]
rep_times <- 1000
G400_performance <- data.frame(brierS = double(rep_times),
                               aucS = double(rep_times),
                               overall.err.rate = double(rep_times))

for(id in 1:rep_times){
  rfsrc_obj400_perf <- rfsrc(agvhi24 ~ ., data = G400, ntree=300, na.action='na.impute', nimpute = 50,
                             importance=TRUE, membership=TRUE)
  x <- rfsrc_obj400_perf
  G400_performance$brierS[id] <- round(brier(x$yvar, if(!is.null(x$predicted.oob) && !all(is.na(x$predicted.oob))) x$predicted.oob else x$predicted) * 100, 2)
  G400_performance$aucS[id] <- round(auc(x$yvar, if(!is.null(x$predicted.oob) && !all(is.na(x$predicted.oob))) x$predicted.oob else x$predicted) * 100, 2)
  err.rate <- cbind(x$err.rate)
  G400_performance$overall.err.rate[id] <- round(100 * err.rate[nrow(err.rate), 1], 2) # paste(round(100 * err.rate[nrow(err.rate), 1], 2), "%", sep = "")

}

save(G400_performance, file = 'AML_aGVHD24_VLSRelief_52225_MaxInd5_pct7_top400feat_performance.RData') # random_performance

G = read.csv('AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75.csv', row.names = 1, stringsAsFactors = F)

Labels <- c(rep("Non-GVHD", dim(G)[1]))
Labels[G$agvhi24==0] <- 'Non-GVHD'
Labels[G$agvhi24==1] <- 'aGVHD'
Labels <- as.factor(Labels)
G$agvhi24 <- Labels

G[is.na(G)]<- 0
G400 <- G[, 1:401]
ranger_feat400 <- ranger(agvhi24 ~ ., data = G400, importance = "impurity_corrected",
                         num.threads=detectCores(),
                         num.trees = 300,
                         #save.memory=TRUE,
                         seed=123,
                         verbose=TRUE,
                         write.forest=TRUE)
print(ranger_feat400)
# Type:                             Classification
# Number of trees:                  300
# Sample size:                      331
# Number of independent variables:  400
# Mtry:                             20
# Target node size:                 1
# Variable importance mode:         impurity_corrected
# Splitrule:                        gini
# OOB prediction error:             24.47 %

save(ranger_feat400,file='AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rangerModel_top400.RData')

feat_imp <- importance(ranger_feat400)
sorted_feat_imp <- sort(feat_imp, decreasing = T)

pval_feat_janitza <- importance_pvalues(ranger_feat400, method = 'janitza')
sorted_feat_imp_janitza <- pval_feat_janitza[order(pval_feat_janitza[,1],decreasing = T),]
# Warning message:
# In importance_pvalues(ranger_feat2, method = "janitza") :
#   Only few negative importance values found, inaccurate p-values. Consider the 'altmann' approach.
# > sorted_feat_imp_janitza <- pval_feat_janitza[order(pval_feat_janitza[,1],decreasing = T),]
write.csv(sorted_feat_imp_janitza, file='AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_janitza_top400feat.csv')

pval_feat_altmann <- importance_pvalues(ranger_feat400, method = "altmann", formula = agvhi24 ~ ., data = G400)
sorted_feat_imp_altmann <- pval_feat_altmann[order(pval_feat_altmann[,1],decreasing = T),]
write.csv(sorted_feat_imp_altmann, file='AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_altmann_top400feat.csv')

# HLA
BMT_cases <-read_feather('AML_patient_donor_autosome_SNP_mat_TXoutcome.feather')
BMT_cases$BMT_case <- as.integer(BMT_cases$BMT_case)
HLA_table <- read.csv('../GWAS_cleaned_BMT_HLA_typing.csv')
meta_info <- read.csv('../GWAS_cleaned_BMT_cases_info.csv')

convert_ind <- sapply(1:dim(BMT_cases)[1], function(x) which(HLA_table$bmt_case == BMT_cases$BMT_case[x]))

G400_hla <- cbind(G400, HLA_table[convert_ind, -c(1,2,3)])

ranger_feat400_hla <- ranger(agvhi24 ~ ., data = G400_hla, importance = "impurity_corrected",
                             num.threads=detectCores(),
                             num.trees = 300,
                             #save.memory=TRUE,
                             seed=123,
                             verbose=TRUE,
                             write.forest=TRUE)
print(ranger_feat400_hla)
# Type:                             Classification
# Number of trees:                  300
# Sample size:                      331
# Number of independent variables:  410
# Mtry:                             20
# Target node size:                 1
# Variable importance mode:         impurity_corrected
# Splitrule:                        gini
# OOB prediction error:             23.26 %

save(ranger_feat400_hla,file='AML_aGVHD24_VLSRelief_52225_MaxInd5_pct75_rangerModel_top400_wHla.RData')2