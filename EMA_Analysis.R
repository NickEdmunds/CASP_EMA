######################################################################################################################################################
# 01-Nov-23. New program to load CASP15 results for EMA predictions.
# # Program is located in: /home/nick/Post_confirmation_projects/CASP15/QA/casp15_ema-main-custom_analysis/casp15_ema-main/custom_analysis
# It loads the my_local_df.csv file which is the result of the full dataset (local_df.csv which was downloaded from GitHub) which was huge - loaded 
# into R and then stripped of all models that my EMA program did NOT predict for. Then my results are added under group X999 and reloaded below.
setwd('/home/nick/Post_confirmation_projects/CASP15/QA/casp15_ema-main-custom_analysis/casp15_ema-main/custom_analysis/')
my_local_df=read.csv('my_local_df.csv', h=TRUE)
#######################################################################################################################################################
# The way this has been done (according to the ipynb downloaded from GitHub) is to remove any rows that don't have an entry for either the group or the
# observed value. Then program a ROC AUC figure (define the 0.75 quartile as the a threshold for the binary variable). If the ROC AUC is less than 0.5 it
# returns a value of 0.5 as below this is worse than random. To do this pROC package is needed.
# Needs to be run for: X999, X365, X002, X283, X494, X398, X426, X041, X121, X101, X120, X282, X248, X158, X083, X266, X298, X126, X248_2, X468, X168, X169, X245, X089,
# X275, X086, X234, X098.

calculate_roc_auc <- function(my_local_df, pred_col, target_col) {
  thresh <- quantile(my_local_df[[target_col]], 0.75, na.rm = TRUE)                             # Calculate the 75th percentile of target_col while removing NAs
  sub_df <- my_local_df[!is.na(my_local_df[[target_col]]) & !is.na(my_local_df[[pred_col]]), ]  # Filter rows with non-null target_col and pred_col
  target_classes <- as.integer(sub_df[[target_col]] > thresh)                                   # Binarize target_col based on the threshold
  return(max(0.5, pROC::roc(target_classes, sub_df[[pred_col]])$auc))                           # Calculate ROC AUC with a minimum of 0.5
}
# Here I've just cycled through lddt, cad, patch_qs and patch_dockq as target_col
calculate_roc_auc(my_local_df, pred_col = "X999", target_col = "lddt")  # my EMA program
calculate_roc_auc(my_local_df, pred_col = "X266", target_col = "patch_dockq")  # MFDR (ranked 1)
calculate_roc_auc(my_local_df, pred_col = "X089", target_col = "patch_dockq")  # RocketX (ranked 2)
calculate_roc_auc(my_local_df, pred_col = "X121", target_col = "patch_dockq")  # VoroIF (ranked 3)
calculate_roc_auc(my_local_df, pred_col = "X083", target_col = "patch_dockq")  # MFDS (ranked 5)
calculate_roc_auc(my_local_df, pred_col = "X041", target_col = "patch_dockq")  # MFD (ranked 6)
calculate_roc_auc(my_local_df, pred_col = "X365", target_col = "patch_dockq")
calculate_roc_auc(my_local_df, pred_col = "X494", target_col = "patch_dockq")
calculate_roc_auc(my_local_df, pred_col = "X426", target_col = "patch_dockq")
calculate_roc_auc(my_local_df, pred_col = "X101", target_col = "patch_dockq")
calculate_roc_auc(my_local_df, pred_col = "X248", target_col = "patch_dockq")
calculate_roc_auc(my_local_df, pred_col = "X158", target_col = "patch_dockq")
calculate_roc_auc(my_local_df, pred_col = "X248_2", target_col = "patch_dockq")
calculate_roc_auc(my_local_df, pred_col = "X468", target_col = "patch_dockq")
calculate_roc_auc(my_local_df, pred_col = "X168", target_col = "patch_dockq")
calculate_roc_auc(my_local_df, pred_col = "X245", target_col = "patch_dockq")

calculate_pearson_r <- function(my_local_df, pred_col, target_col) {
  sub_df <- my_local_df[!is.na(my_local_df[[target_col]]) & !is.na(my_local_df[[pred_col]]), ]  # Filter rows with non-null target_col and pred_col
  pearson_r <- cor(sub_df[[pred_col]], sub_df[[target_col]])                                    # Calculate Pearson R
  return(pearson_r)
}
calculate_pearson_r(my_local_df, pred_col = "X999", target_col = "patch_dockq")

calculate_spearman <- function(my_local_df, pred_col, target_col) {
  sub_df <- my_local_df[!is.na(my_local_df[[target_col]]) & !is.na(my_local_df[[pred_col]]), ]  # Filter rows with non-null target_col and pred_col
  spearman <- cor(sub_df[[pred_col]], sub_df[[target_col]], method = "spearman")                # Calculate Spearman correlation
  return(spearman)
}

calculate_spearman(my_local_df, pred_col = "X999", target_col = "patch_dockq")
calculate_spearman(my_local_df, pred_col = "X266", target_col = "lddt")  # MFDR (ranked 1)
calculate_spearman(my_local_df, pred_col = "X089", target_col = "lddt")  # RocketX (ranked 2)
calculate_spearman(my_local_df, pred_col = "X121", target_col = "lddt")  # VoroIF (ranked 3)
calculate_spearman(my_local_df, pred_col = "X083", target_col = "lddt")  # MFDS (ranked 5)
calculate_spearman(my_local_df, pred_col = "X041", target_col = "lddt")  # MFD (ranked 6)
calculate_spearman(my_local_df, pred_col = "X365", target_col = "lddt")
calculate_spearman(my_local_df, pred_col = "X494", target_col = "lddt")
calculate_spearman(my_local_df, pred_col = "X426", target_col = "lddt")
calculate_spearman(my_local_df, pred_col = "X101", target_col = "lddt")
calculate_spearman(my_local_df, pred_col = "X248", target_col = "lddt")
calculate_spearman(my_local_df, pred_col = "X158", target_col = "lddt")
calculate_spearman(my_local_df, pred_col = "X248_2", target_col = "lddt")
calculate_spearman(my_local_df, pred_col = "X468", target_col = "lddt")
calculate_spearman(my_local_df, pred_col = "X168", target_col = "lddt")
calculate_spearman(my_local_df, pred_col = "X245", target_col = "lddt")

# The results have been manually added to the following dataset
results_df=read.csv('results_df.csv', h=TRUE)

results_df$lddt_roc_auc_z  <- scale(results_df$ROC_AUC_lddt)
results_df$lddt_pearson_z  <- scale(results_df$pearson_lddt)
results_df$lddt_spearman_z <- scale(results_df$spearman_lddt)

results_df$cad_roc_auc_z  <- scale(results_df$ROC_AUC_cad)
results_df$cad_pearson_z  <- scale(results_df$pearson_cad)
results_df$cad_spearman_z <- scale(results_df$spearman_cad)

results_df$patchqs_roc_auc_z  <- scale(results_df$ROC_AUC_patch_qs)
results_df$patchqs_pearson_z  <- scale(results_df$pearson_patch_qs)
results_df$patchqs_spearman_z <- scale(results_df$spearman_patch_qs)

results_df$patchdockq_roc_auc_z  <- scale(results_df$ROC_AUC_patch_dockq)
results_df$patchdockq_pearson_z  <- scale(results_df$pearson_patch_dockq)
results_df$patchdockq_spearman_z <- scale(results_df$spearman_patch_dockq)

results_df$lddt_z       <- ((0.5 * results_df$lddt_pearson_z) + (0.5 * results_df$lddt_spearman_z) + results_df$lddt_roc_auc_z)
results_df$cad_z        <- ((0.5 * results_df$cad_pearson_z) + (0.5 * results_df$cad_spearman_z) + results_df$cad_roc_auc_z)
results_df$patchqs_z    <- ((0.5 * results_df$patchqs_pearson_z) + (0.5 * results_df$patchqs_spearman_z) + results_df$patchqs_roc_auc_z)
results_df$patchdockq_z <- ((0.5 * results_df$patchdockq_pearson_z) + (0.5 * results_df$patchdockq_spearman_z) + results_df$patchdockq_roc_auc_z)

results_df$final_score  <- (results_df$lddt_z + results_df$cad_z + results_df$patchqs_z + results_df$patchdockq_z)

png("EMA_Beside_Z-scores_subset.png", width = 16, height = 15, units = 'cm', res = 600)
# Create a data frame with the contributions
results_df <- results_df[order(results_df$final_score), ]
contributions <- results_df[, c("lddt_z", "cad_z", "patchqs_z", "patchdockq_z")]
# Create the stacked horizontal bar graph
par(mar = c(5, 8, 4, 4))  # Increase the left margin
barplot(t(contributions), beside = T, col = c("black", "orange", "lightblue", "green"), xlim=c(-5, 5),
        names.arg = results_df$Gname, horiz = TRUE, las = 1 )
# Add a legend
legend("bottomright", legend = colnames(contributions), fill = c("black", "orange", "lightblue", "green"))
# Add labels and a title
mtext("Group", side = 2, line = 5)
title(main = "Bar Graph of Z-scores", xlab = "Cumulative Z-Score Value", ylab = "")
dev.off()

# The results were slightly different than the CASP official reults and I think they did the thing of setting -ve values to 0. So this does that.
results_df$lddt_z <-       ifelse(results_df$lddt_pearson_z < 0, 0, 0.5 * results_df$lddt_pearson_z) + ifelse(results_df$lddt_spearman_z < 0, 0, 0.5 * results_df$lddt_spearman_z) + ifelse(results_df$lddt_roc_auc_z < 0, 0, results_df$lddt_roc_auc_z)
results_df$cad_z <-        ifelse(results_df$cad_pearson_z < 0, 0, 0.5 *  results_df$cad_pearson_z) +  ifelse(results_df$cad_spearman_z < 0, 0, 0.5 *  results_df$cad_spearman_z) +  ifelse(results_df$cad_roc_auc_z < 0, 0, results_df$cad_roc_auc_z)
results_df$patchqs_z <-    ifelse(results_df$patchqs_pearson_z < 0, 0, 0.5 * results_df$patchqs_pearson_z) +       ifelse(results_df$patchqs_spearman_z < 0, 0, 0.5 * results_df$patchqs_spearman_z) + ifelse(results_df$patchqs_roc_auc_z < 0, 0, results_df$patchqs_roc_auc_z)
results_df$patchdockq_z <- ifelse(results_df$patchdockq_pearson_z < 0, 0, 0.5 * results_df$patchdockq_pearson_z) + ifelse(results_df$patchdockq_spearman_z < 0, 0, 0.5 * results_df$patchdockq_spearman_z) + ifelse(results_df$patchdockq_roc_auc_z < 0, 0, results_df$patchdockq_roc_auc_z)

results_df$final_score  <- (results_df$lddt_z + results_df$cad_z + results_df$patchqs_z + results_df$patchdockq_z)

png("EMA_Stacked_Z-scores_subset_NonNeg.png", width = 16, height = 15, units = 'cm', res = 600)
# Create a data frame with the contributions
results_df <- results_df[order(results_df$final_score), ]
contributions <- results_df[, c("lddt_z", "cad_z", "patchqs_z", "patchdockq_z")]
# Create the stacked horizontal bar graph
par(mar = c(5, 8, 4, 4))  # Increase the left margin
barplot(t(contributions), beside = F, col = c("black", "orange", "lightblue", "green"), xlim=c(0, 15),
        names.arg = results_df$Gname, horiz = TRUE, las = 1 )
# Add a legend
legend("bottomright", legend = colnames(contributions), fill = c("black", "orange", "lightblue", "green"))
# Add labels and a title
mtext("Group", side = 2, line = 5)
title(main = "Stacked Bar Graph of Z-scores", xlab = "Cumulative Positive Z-Score Value", ylab = "")
dev.off()
################################################################################################################################
# This bit analyses the results of the whole set of scores which was downloaded from the CASP website (and also from the GitHub site).
# This is called QA_AVE_predictioncenter_oligo-local.csv which was copied from QA_predictioncenter_oligo-local.csv which is in
# /home/nick/Post_confirmation_projects/CASP15/QA/casp15_ema-main-custom_analysis/ but has only the data for average results (individual target data)
# deleted. Also col names have been changed to match my program above.
# The point of this bit is basically to see if I come up with the same graph as the official CASP results.
full_results_df=read.csv('QA_AVE_predictioncenter_oligo-local.csv', h=TRUE)

full_results_df$lddt_roc_auc_z  <- scale(full_results_df$ROC_AUC_lddt)
full_results_df$lddt_pearson_z  <- scale(full_results_df$pearson_lddt)
full_results_df$lddt_spearman_z <- scale(full_results_df$spearman_lddt)

full_results_df$cad_roc_auc_z  <- scale(full_results_df$ROC_AUC_cad)
full_results_df$cad_pearson_z  <- scale(full_results_df$pearson_cad)
full_results_df$cad_spearman_z <- scale(full_results_df$spearman_cad)

full_results_df$patchqs_roc_auc_z  <- scale(full_results_df$ROC_AUC_patch_qs)
full_results_df$patchqs_pearson_z  <- scale(full_results_df$pearson_patch_qs)
full_results_df$patchqs_spearman_z <- scale(full_results_df$spearman_patch_qs)

full_results_df$patchdockq_roc_auc_z  <- scale(full_results_df$ROC_AUC_patch_dockq)
full_results_df$patchdockq_pearson_z  <- scale(full_results_df$pearson_patch_dockq)
full_results_df$patchdockq_spearman_z <- scale(full_results_df$spearman_patch_dockq)

full_results_df$lddt_z       <- ((0.5 * full_results_df$lddt_pearson_z) + (0.5 * full_results_df$lddt_spearman_z) + full_results_df$lddt_roc_auc_z)
full_results_df$cad_z        <- ((0.5 * full_results_df$cad_pearson_z)  + (0.5 * full_results_df$cad_spearman_z)  + full_results_df$cad_roc_auc_z)
full_results_df$patchqs_z    <- ((0.5 * full_results_df$patchqs_pearson_z) + (0.5 * full_results_df$patchqs_spearman_z) + full_results_df$patchqs_roc_auc_z)
full_results_df$patchdockq_z <- ((0.5 * full_results_df$patchdockq_pearson_z) + (0.5 * full_results_df$patchdockq_spearman_z) + full_results_df$patchdockq_roc_auc_z)

full_results_df$final_score  <- (full_results_df$lddt_z + full_results_df$cad_z + full_results_df$patchqs_z + full_results_df$patchdockq_z)

png("EMA_Full_beside_Z-scores.png", width = 16, height = 15, units = 'cm', res = 600)
# Create a data frame with the contributions
full_results_df <- full_results_df[order(full_results_df$final_score), ]
contributions <- full_results_df[, c("lddt_z", "cad_z", "patchqs_z", "patchdockq_z")]
# Create the stacked horizontal bar graph
par(mar = c(5, 8, 4, 4))  # Increase the left margin
barplot(t(contributions), beside = T, col = c("black", "orange", "lightblue", "green"), xlim=c(-5, 5),
        names.arg = full_results_df$Gname, horiz = TRUE, las = 1 )
# Add a legend
legend("bottomright", legend = colnames(contributions), fill = c("black", "orange", "lightblue", "green"))
# Add labels and a title
mtext("Group", side = 2, line = 5)
title(main = "Bar Graph of full Z-scores", xlab = "Cumulative Z-Score Value", ylab = "")
dev.off()

# The results were slightly different than the CASP official reults and I think they did the thing of setting -ve values to 0. So this does that.
full_results_df$lddt_z <- ifelse(full_results_df$lddt_pearson_z < 0, 0, 0.5 * full_results_df$lddt_pearson_z) + ifelse(full_results_df$lddt_spearman_z < 0, 0, 0.5 * full_results_df$lddt_spearman_z) +
  ifelse(full_results_df$lddt_roc_auc_z < 0, 0, full_results_df$lddt_roc_auc_z)
full_results_df$cad_z <- ifelse(full_results_df$cad_pearson_z < 0, 0, 0.5 * full_results_df$cad_pearson_z) + ifelse(full_results_df$cad_spearman_z < 0, 0, 0.5 * full_results_df$cad_spearman_z) +
  ifelse(full_results_df$cad_roc_auc_z < 0, 0, full_results_df$cad_roc_auc_z)
full_results_df$patchqs_z <- ifelse(full_results_df$patchqs_pearson_z < 0, 0, 0.5 * full_results_df$patchqs_pearson_z) + ifelse(full_results_df$patchqs_spearman_z < 0, 0, 0.5 * full_results_df$patchqs_spearman_z) +
  ifelse(full_results_df$patchqs_roc_auc_z < 0, 0, full_results_df$patchqs_roc_auc_z)
full_results_df$patchdockq_z <- ifelse(full_results_df$patchdockq_pearson_z < 0, 0, 0.5 * full_results_df$patchdockq_pearson_z) + ifelse(full_results_df$patchdockq_spearman_z < 0, 0, 0.5 * full_results_df$patchdockq_spearman_z) +
  ifelse(full_results_df$patchdockq_roc_auc_z < 0, 0, full_results_df$patchdockq_roc_auc_z)

full_results_df$final_score  <- (full_results_df$lddt_z + full_results_df$cad_z + full_results_df$patchqs_z + full_results_df$patchdockq_z)

png("EMA_Full_Stacked_Z-scores_NonNeg.png", width = 16, height = 15, units = 'cm', res = 600)
# Create a data frame with the contributions
full_results_df <- full_results_df[order(full_results_df$final_score), ]
contributions <- full_results_df[, c("lddt_z", "cad_z", "patchqs_z", "patchdockq_z")]
# Create the stacked horizontal bar graph
par(mar = c(5, 8, 4, 4))  # Increase the left margin
barplot(t(contributions), beside = F, col = c("black", "orange", "lightblue", "green"), xlim=c(0, 15),
        names.arg = full_results_df$Gname, horiz = TRUE, las = 1 )
# Add a legend
legend("bottomright", legend = colnames(contributions), fill = c("black", "orange", "lightblue", "green"))
# Add labels and a title
mtext("Group", side = 2, line = 5)
title(main = "Stacked Bar Graph of full Z-scores", xlab = "Cumulative Positive Z-Score Value", ylab = "")
dev.off()

# 18-Nov-23. Local interface Rankings - trying to reproduce the ranking graphs from the official analysis.
png("EMA_local_interface_Ranking.png", width = 16, height = 15, units = 'cm', res = 600)
# Create a data frame with the contributions replaced by Av. per target ROC AUC
full_results_df$Av_ROC <- ((full_results_df$ROC_AUC_patch_dockq+full_results_df$ROC_AUC_patch_qs+full_results_df$ROC_AUC_cad+full_results_df$ROC_AUC_lddt)/4)
full_results_df <- full_results_df[order(full_results_df$Av_ROC), ]
# Create the stacked horizontal bar graph
par(mar = c(5, 8, 4, 4))  # Increase the left margin
barplot(t(full_results_df$Av_ROC), beside = F, xlim=c(0.0, 0.8),
        names.arg = full_results_df$Gname, horiz = TRUE, las = 1 )
# Add labels and a title
mtext("Group", side = 2, line = 5)
title(main = "Local interface residue indentification ranking", xlab = "Ave. per-target ROC AUC", ylab = "")
dev.off()

######### 23-Nov-23 ########
# Adding a little bit of modelling analysis for chapter 4.
setwd('/home/nick/Post_confirmation_projects/CASP15/QA/casp15_ema-main-custom_analysis/casp15_ema-main/custom_analysis/')
CASP15_oligo=read.csv('Composite_CASPoligo_results.csv', h=TRUE)

png("CASP15_local_cor.png", width = 18, height = 15, units = 'cm', res = 600)
plot(main='CASP15 predicted versus CASP local scores (assembly models)',xlab='CASP15 calculated local score', ylab='ModFOLDdockR local score', CASP15_oligo$local, CASP15_oligo$localRscore, col=c("blue"),pch=16,abline(lm(CASP15_oligo$localRscore ~ CASP15_oligo$local)))
legend(x='bottomright', legend=paste('Pearson =',round(cor(CASP15_oligo$local, CASP15_oligo$localRscore),2)))
dev.off()

png("CASP15_Global_cor.png", width = 18, height = 15, units = 'cm', res = 600)
plot(main='CASP15 predicted versus CASP Global scores (assembly models)',xlab='CASP15 calculated Global score', ylab='ModFOLDdockR Global score', CASP15_oligo$Global, CASP15_oligo$GlobalRscore, col=c("blue"),pch=16,abline(lm(CASP15_oligo$GlobalRscore ~ CASP15_oligo$Global)))
legend(x='bottomright', legend=paste('Pearson =',round(cor(CASP15_oligo$Global, CASP15_oligo$GlobalRscore),2)))
dev.off()

CASP15_oligo$Total       <- (CASP15_oligo$local + CASP15_oligo$Global)/2
CASP15_oligo$TotalRscore <- (CASP15_oligo$localRscore + CASP15_oligo$GlobalRscore)/2
png("CASP15_Total_cor.png", width = 18, height = 15, units = 'cm', res = 600)
plot(main='CASP15 predicted versus CASP Total scores (assembly models)',xlab='CASP15 calculated Total score', ylab='ModFOLDdockR Total score', CASP15_oligo$Total, CASP15_oligo$TotalRscore, col=c("blue"),pch=16,abline(lm(CASP15_oligo$TotalRscore ~ CASP15_oligo$Total)))
legend(x='bottomright', legend=paste('Pearson =',round(cor(CASP15_oligo$Total, CASP15_oligo$TotalRscore),2)))
dev.off()

png("CASP15_lDDT_cor.png", width = 18, height = 15, units = 'cm', res = 600)
plot(main='Predicted Global versus CASP15 lDDT-oligo score (assembly models)',xlab='CASP15 lDDT-oligo score', ylab='ModFOLDdockR Global score', CASP15_oligo$lDDToligo, CASP15_oligo$GlobalRscore, col=c("blue"),pch=16,abline(lm(CASP15_oligo$GlobalRscore ~ CASP15_oligo$lDDToligo)))
legend(x='bottomright', legend=paste('Pearson =',round(cor(CASP15_oligo$lDDToligo, CASP15_oligo$GlobalRscore),2)))
dev.off()

png("CASP15_TM_cor.png", width = 18, height = 15, units = 'cm', res = 600)
plot(main='Predicted Global versus CASP15 TM-score (assembly models)',xlab='CASP15 TM-score', ylab='ModFOLDdockR Global score', CASP15_oligo$TMscore, CASP15_oligo$GlobalRscore, col=c("blue"),pch=16,abline(lm(CASP15_oligo$GlobalRscore ~ CASP15_oligo$TMscore)))
legend(x='bottomright', legend=paste('Pearson =',round(cor(CASP15_oligo$TMscore, CASP15_oligo$GlobalRscore),2)))
dev.off()

png("CASP15_QS_cor.png", width = 18, height = 15, units = 'cm', res = 600)
plot(main='Predicted local versus CASP15 QS-score (assembly models)',xlab='CASP15 QS-score', ylab='ModFOLDdockR local score', CASP15_oligo$QSglob, CASP15_oligo$localRscore, col=c("blue"),pch=16,abline(lm(CASP15_oligo$localRscore ~ CASP15_oligo$QSglob)))
legend(x='bottomright', legend=paste('Pearson =',round(cor(CASP15_oligo$QSglob, CASP15_oligo$localRscore),2)))
dev.off()

# CAMEO data  - min-max normalisation (24/3/24)
lDDT_data <- data.frame(
  Server = c('Server 1', 'Server 76', 'Server 4', 'Server 2'),
  Total_lDDT_Score = c(81.6, 61.6, 43.3, 16.9),
  Structures_Modeled = c(127, 127, 80, 40))

# Normalize the scores for Servers 2 and 4 using actual scores for Server 1
lDDT_data$Normalized_Score[lDDT_data$Server == 'Server 1'] <- 1
lDDT_data$Normalized_Score[lDDT_data$Server == 'Server 2'] <-  16.9/23.7
lDDT_data$Normalized_Score[lDDT_data$Server == 'Server 4'] <-  43.3/49.2
lDDT_data$Normalized_Score[lDDT_data$Server == 'Server 76'] <- 61.6/81.6

# Plot the bar chart
setwd('/home/nick/Documents/Thesis_writup')
png("CAMEO_lDDT_scores.png", width = 18, height = 15, units = 'cm', res = 600)
barplot(lDDT_data$Normalized_Score, names.arg = lDDT_data$Server, col = "skyblue",
        main = "Normalised cumulative lDDT-oligo score for CAMEO servers",
        xlab = "Server", ylab = "Normalised lDDT-oligo Score", ylim = c(0, 1), 
        border = "black", space = 0.5)
dev.off()

QS_data <- data.frame(
  Server = c('Server 1', 'Server 76', 'Server 4', 'Server 2'),
  Total_QS_Score = c(44.7, 33.7, 23.3, 8.0),
  Structures_Modeled = c(127, 127, 80, 40))

# Normalize the scores for Servers 2 and 4 using actual scores for Server 1
QS_data$Normalized_Score[QS_data$Server == 'Server 1'] <- 44.7/44.7
QS_data$Normalized_Score[QS_data$Server == 'Server 2'] <- 8.0/11.4
QS_data$Normalized_Score[QS_data$Server == 'Server 4'] <- 23.3/27.3
QS_data$Normalized_Score[QS_data$Server == 'Server 76'] <- 33.7/44.7

# Plot the bar chart
png("CAMEO_QS_scores.png", width = 18, height = 15, units = 'cm', res = 600)
barplot(QS_data$Normalized_Score, names.arg = QS_data$Server, col = "darkgreen",
        main = "Normalised cumulative QS-score for CAMEO servers",
        xlab = "Server", ylab = "Normalised QS-score", ylim = c(0, 1), 
        border = "black", space = 0.5)
dev.off()
