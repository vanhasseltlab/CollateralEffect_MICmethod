#T-test for collateral sensitivity on data from data preparation script

#Libraries and functions
library(tidyverse)
source("scripts/functions_CR.R")

####Import data####
#Pick species
species <- "Escherichia-coli"

load(paste0("data/clean/MIC_clean_", species,".Rdata"))

####Test data####
#Pick direction
CResponse <- "both" #possibilities: "CR", "CS", or "both"

criteria_quant <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
results <- vector(mode = "list", length = length(criteria_quant))
names(results) <- criteria_quant
direction <- c("two.sided", "less", "greater")[c("both", "CS", "CR") == CResponse]

antibiotics <- colnames(MIC_clean)
n <- length(antibiotics)

t_test <- as.data.frame(matrix(0, nrow = n*(n - 1), ncol = 8))
names(t_test) <- c("A", "B", "t", "p", "n", "d", "effect_size", "n_cond")
t_test$p <- 1
t_test$mean_Y <- 999
t_test$mean_X <- 999

###added to test conditional distribution
t_test$t_cond <- t_test$p_cond <- t_test$effect_size_cond <- 0
###

for (crit in 1:length(results)) {
  counter <- 0
  
  for (dep in 1:n) {
    for (indep in 1:n) {
      if (dep == indep) {
        next
      }
      counter <- counter + 1
      
      dat <- log2(MIC_clean[, c(dep, indep)])
      dat <- dat[!is.na(dat[, 1]), ]
     
      t_test[counter, c("A", "B")] <- antibiotics[c(dep, indep)]
      t_test[counter, "n"] <- nrow(dat)


      criterium <- quantile(dat[, 2], criteria_quant[crit], na.rm = TRUE)
      
      vals <- sort(unique(dat[, 2]))
      B <- dat[!is.na(dat[, 2]), 2]
      if (sum(B > criterium) < sum(B < criterium)) {
        criterium <- mean(c(vals[which(vals == criterium) - 1], criterium))
      } else if (sum(B > criterium) > sum(B < criterium)) {
        criterium <- mean(c(vals[which(vals == criterium) + 1], criterium))
      }
      
      t_test[counter, "d"] <- criterium
      
      X_w <- dat[dat[, 2] < criterium | is.na(dat[, 2]), 1]
      X_w_cond <-  dat[dat[, 2] < criterium & !is.na(dat[, 2]), 1]
      Y <- dat[dat[, 2] >= criterium & !is.na(dat[, 2]), 1]
      
      t_test[counter, "mean_Y"] <- mean(2^c(X_w, Y))
      t_test[counter, "mean_X"] <- mean(2^Y)
      t_test[counter, "n_cond"] <- length(Y)
     
      if (any(length(X_w) < 2, length(Y) < 2)) next()
      
      t_test_i <- t.test(Y, X_w, var.equal = TRUE, alternative = direction)
      t_test[counter, 3:4] <- c(t_test_i$statistic, t_test_i$p.value)
      t_test[counter, "effect_size"] <-  mean(Y) - mean(c(X_w, Y))
      
      t_test_i_cond <- t.test(Y, X_w_cond, var.equal = TRUE, alternative = direction)
      t_test[counter, "t_cond"] <- t_test_i_cond$statistic
      t_test[counter, "p_cond"] <- t_test_i_cond$p.value
      t_test[counter, "effect_size_cond"] <-  mean(Y) - mean(X_w_cond)

    }
  }
  result_crit <- data.frame(t_test, p_BY = p.adjust(t_test$p, method = "BY"), p_BY_cond = p.adjust(t_test$p_cond, method = "BY"))
  
  results[[crit]] <- result_crit
}

save(results, file = paste0('results/results', species, "_cond.Rdata"))
write.table(results[["0.5"]], paste0('results/results', species, ".txt"))

#Plot results
#Overview plot

pdf(file = paste0("results/figures/t_allcomparisons_", species, ".pdf"), height = 7*0.8, width = 9*0.8)
print(PlotSignificantEffect(results, 0.4, FDR_crit = 0.05, species = species))
print(PlotSignificantEffect(results, 0.5, FDR_crit = 0.05, species = species))
print(PlotSignificantEffect(results, 0.6, FDR_crit = 0.05, species = species))
print(PlotSignificantEffect(results, 0.7, FDR_crit = 0.05, species = species))
print(PlotSignificantEffect(results, 0.8, FDR_crit = 0.05, species = species))
print(PlotSignificantEffect(results, 0.9, FDR_crit = 0.05, species = species))
dev.off()

pdf(file = paste0("results/figures/t_allcomparisons_", species, "_noFDR.pdf"), height = 7*0.8, width = 9*0.8)
print(PlotSignificantEffect(results, 0.4, FDR_crit = 1, species = species))
print(PlotSignificantEffect(results, 0.5, FDR_crit = 1, species = species))
print(PlotSignificantEffect(results, 0.6, FDR_crit = 1, species = species))
print(PlotSignificantEffect(results, 0.7, FDR_crit = 1, species = species))
print(PlotSignificantEffect(results, 0.8, FDR_crit = 1, species = species))
print(PlotSignificantEffect(results, 0.9, FDR_crit = 1, species = species))
dev.off()


#Significant findings plot
library(gridExtra)

pdf(file = paste0("results/figures/distribution_significant_CS_", species, ".pdf"), height = 7, width = 9)
for (ra in 1:sum(results[['0.5']]$p_BY < 0.05 & results[['0.5']]$t < 0)) {
  PlotCRDistributions(MIC_clean, results, 0.5, t_rank = ra, one_direction = TRUE, CResponse = "CS")
}
dev.off()

pdf(file = paste0("results/figures/distribution_significant_CR_", species, ".pdf"), height = 7, width = 9)
for (ra in 1:sum(results[['0.5']]$p_BY < 0.05 & results[['0.5']]$t > 0)) {
  PlotCRDistributions(MIC_clean, results, 0.5, t_rank = ra, one_direction = TRUE, CResponse = "CR")
}
dev.off()


pdf(file = paste0("results/figures/distribution_twoway_significant_CS_", species, ".pdf"), height = 7, width = 9)
for (ra in 1:sum(results[['0.5']]$p_BY < 0.05 & results[['0.5']]$t < 0))  {
  PlotCRDistributions(MIC_clean, results, 0.5, t_rank = ra, one_direction = FALSE, CResponse = "CS")
}
dev.off()

pdf(file = paste0("results/figures/distribution_stacked_", species, ".pdf"), height = 5, width = 9)
  PlotStackDistribution(MIC_clean, results, 0.5, t_rank = 1, one_direction = FALSE, CResponse = "CS")
  PlotStackDistribution(MIC_clean, results, 0.5, t_rank = 1, one_direction = FALSE, CResponse = "CR")
dev.off()
