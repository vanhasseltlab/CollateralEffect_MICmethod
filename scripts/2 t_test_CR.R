#T-test for collateral sensitivity on data from data preparation script

#Libraries
library(tidyverse)

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
t_test <- as.data.frame(matrix(0, nrow = n*(n - 1), ncol = 6))
names(t_test) <- c("A", "B", "t", "p", "n", "d")
t_test$p <- 1

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
     
      t_test[counter, 1:2] <- antibiotics[c(dep, indep)]
      t_test[counter, 5] <- nrow(dat)
      
      criterium <- quantile(dat[, 2], criteria_quant[crit], na.rm = TRUE)
      t_test[counter, 6] <- criterium
      
      X_w <- dat[dat[, 2] < criterium | is.na(dat[, 2]), 1]
      Y <- dat[dat[, 2] >= criterium & !is.na(dat[, 2]), 1]
      if (any(length(X_w) < 2, length(Y) < 2)) next()
      t_test_i <- t.test(Y, X_w, var.equal = TRUE, alternative = direction)
      t_test[counter, 3:4] <- c(t_test_i$statistic, t_test_i$p.value)
    }
  }
  result_crit <- data.frame(t_test, p_BY = p.adjust(t_test$p, method = "BY"))
  
  results[[crit]] <- result_crit
}


#Plot results
#Overview plot

pdf(file = paste0("results/figures/t_allcomparisons_", species, ".pdf"), height = 7, width = 9)
print(PlotSignificantTvalue(results, 0.4), FDR_crit = 0.05)
print(PlotSignificantTvalue(results, 0.5), FDR_crit = 0.05)
print(PlotSignificantTvalue(results, 0.6), FDR_crit = 0.05)
print(PlotSignificantTvalue(results, 0.7), FDR_crit = 0.05)
print(PlotSignificantTvalue(results, 0.8), FDR_crit = 0.05)
print(PlotSignificantTvalue(results, 0.9), FDR_crit = 0.05)
dev.off()


#Significant findings plot
library(gridExtra)

pdf(file = paste0("results/figures/distribution_significant_CS_", species, ".pdf"), height = 7, width = 9)
for (ra in 1:sum(results[['0.5']]$p_BY < 0.15)) {
  PlotCRDistributions(MIC_clean, results, 0.5, t_rank = ra, one_direction = TRUE, CResponse = "CS")
}
dev.off()

pdf(file = paste0("results/figures/distribution_significant_CR_", species, ".pdf"), height = 7, width = 9)
for (ra in 1:sum(results[['0.5']]$p_BY < 0.15)) {
  PlotCRDistributions(MIC_clean, results, 0.5, t_rank = ra, one_direction = TRUE, CResponse = "CR")
}
dev.off()


pdf(file = paste0("results/figures/distribution_twoway_significant_CS_", species, ".pdf"), height = 7, width = 9)
for (ra in 1:sum(results[['0.5']]$p_BY < 0.15)) {
  PlotCRDistributions(MIC_clean, results, 0.5, t_rank = ra, one_direction = FALSE, CResponse = "CS")
}
dev.off()


View(results[['0.5']])

###experimenting

quantile(log2(MIC_clean$AMK), c(0.4, 0.5, 0.6, 0.8, 0.9), na.rm = T)
quantile(log2(MIC_clean$TZP), c(0.4, 0.5, 0.6, 0.8, 0.9), na.rm = T)

plot(log2(MIC_clean$CFZ), log2(MIC_clean$TZP), cex = seq(0.1, 3, length.out = length(MIC_clean$CFZ)))

hist(log2(MIC_clean$CFZ))
hist(log2(MIC_clean$CFZ[log2(MIC_clean$TZP) > 6]), xlim = c(0, 6))


