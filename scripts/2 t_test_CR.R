#T-test for collateral sensitivity on data from data preparation script

#Libraries and functions
library(tidyverse)
source("scripts/functions_CR.R")

####Import data####
#Pick species
species <- "Escherichia-coli"

load(paste0("data/clean/MIC_clean_", species,".Rdata"))

####collateral effect T-test####
#Pick direction
CResponse <- "both" #possibilities: "CR", "CS", or "both"
#Pick dichotomisation criterium
dich_crit <- "median" #can be an quantile value or "median" where the data is split in (close to) equal halfs 

antibiotics <- colnames(MIC_clean)
m <- length(antibiotics)

t_test_results <- as.data.frame(matrix(NA, nrow = m*(m - 1), ncol = 9))
names(t_test_results) <- c("A", "B", "n", "n_resi", "tau", "mean_resi", "mean_sens", "t", "p")

counter <- 0
for (dep in 1:m) {
  for (indep in 1:m) {
    if (dep == indep) {
      next
    }
    counter <- counter + 1
    antibiotics <- names(MIC_clean)[c(dep, indep)]
    t_test <- CETTest(A = MIC_clean[, dep], B = MIC_clean[, indep], crit_type = "median", CE_type = "both")
    t_test_results[counter, ] <- data.frame(A = antibiotics[1], B = antibiotics[2], 
                                    n = sum(lengths(t_test$data)), n_resi = length(t_test$data$`A|B = r`),
                                    tau = t_test$tau, mean_resi = t_test$estimate[1], 
                                    mean_sens = t_test$estimate[2], t = t_test$statistic,
                                    p = t_test$p.value, stringsAsFactors = F)
    }
}
t_test_results <- t_test_results %>% 
  mutate(effect_size = mean_resi - mean_sens, 
         p_BY = p.adjust(p, method = "BY"),
         effect_type = ifelse(t > 0, "CR", "CS"))

save(t_test_results, file = paste0('results/results_', species, ".Rdata"))

PlotSignificantEffect(t_test_results, FDR_crit = 0.05, species = species)
PlotCondDistribution(MIC_clean, t_test_results, t_rank = 1, one_direction = FALSE, CResponse = "CS", FDR_crit = 0.05, whichAB = NULL,
                     colours = c("#BFC6B8", "#4A5242"))
PlotCondDistribution(MIC_clean, t_test_results, t_rank = 1, one_direction = FALSE, CResponse = "CR", FDR_crit = 0.05, whichAB = NULL,
                     colours = c("#BFC6B8", "#4A5242"))


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
