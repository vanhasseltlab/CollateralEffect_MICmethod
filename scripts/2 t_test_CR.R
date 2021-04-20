#T-test for collateral sensitivity on data from data preparation script

#Libraries and functions
library(tidyverse)
source("scripts/functions_CR.R")

####Import data####
#Pick species
species <- "Escherichia coli"
data_source <- "NIH"

#load(paste0("data/clean/MIC_clean_", species, "_", data_source, ".Rdata"))

load("data/clean/MIC_clean_Ecoli_ARES_NIH_PATRIC.Rdata")

MIC_clean <- NIH_df

#E coli from paper:
load("~/@Work/Shared/zweplb/LZ3_collateral_sensitivity/data/clean/MIC_clean_Escherichia-coli.Rdata")

####collateral effect T-test####
#Pick direction
CResponse <- "both" #possibilities: "CR", "CS", or "both"
#Pick dichotomisation criterium
dich_crit <- NULL #can be a quantile or MIC value or NULL for option median
crit_type <- "median" #can be "quant" or "MIC" or  "median"

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
    
    if (length(unique(MIC_clean[, dep])) < 2 | length(unique(MIC_clean[, indep])) < 2) {
      next
    }
    
    #antibiotics <- c(dep, indep)
    t_test <- CETTest(A = MIC_clean[, dep], B = MIC_clean[, indep], crit_type = crit_type, CE_type = "both", 
                      criterium = dich_crit)
    if (is.null(t_test$statistic)) {
      t_test_results[counter,] <- data.frame(A = antibiotics[1], B = antibiotics[2], 
                                                 n = sum(lengths(t_test$data)), n_resi = length(t_test$data$`A|B = r`),
                                                 tau = t_test$tau, mean_resi = NA, 
                                                 mean_sens = NA, t = NA,
                                                 p = NA, stringsAsFactors = F)
      next
    }
    
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
         effect_type = ifelse(t > 0, "CR", "CS")) %>% 
  filter(!is.na(A))

save(t_test_results, file = paste0('results/results_', species, "_", data_source, "_", dich_crit, ".Rdata"))


#Plot results
PlotSignificantEffect(t_test_results, FDR_crit = 0.05, species = "PATRIC: E. coli")
PlotCondDistribution(MIC_clean, t_test_results, t_rank = 1, one_direction = FALSE, CResponse = "CS", FDR_crit = 0.05,
                     colours = c("#BFC6B8", "#4A5242"))
PlotCondDistribution(MIC_clean, t_test_results, t_rank = 1, one_direction = FALSE, CResponse = "CR", FDR_crit = 0.05,
                     colours = c("#BFC6B8", "#4A5242"))





cairo_pdf("results/figures/ARES_results.pdf", width = 9, height =  7, onefile = T)
for (species in c("Escherichia coli", "Pseudomonas aeruginosa", "Klebsiella pneumoniae", "Stenotrophomonas maltophilia")) {
  load(paste0('results/results_', species, "_ARES_", dich_crit, ".Rdata"))
  load(paste0("data/clean/MIC_clean_", species, "_ARES.Rdata"))
  print(PlotSignificantEffect(t_test_results, FDR_crit = 0.05, species = species))
  Figure4A <- PlotCondDistribution(MIC_clean, t_test_results, t_rank = 1, one_direction = FALSE, CResponse = "CS", FDR_crit = 0.05,
                                   colours = c("#BFC6B8", "#4A5242"), separate_plots = TRUE, ran = c(-7, 9)) 
  
  Figure4B <- PlotCondDistribution(MIC_clean, t_test_results, t_rank = 1, one_direction = FALSE, CResponse = "CR", FDR_crit = 0.05,
                                   colours = c("#BFC6B8", "#4A5242"), separate_plots = TRUE, ran = c(-7, 9))
  
  grid.arrange(Figure4A[[1]] + labs(tag = "CS"), Figure4A[[2]] + labs(tag = " "), 
               Figure4B[[1]] + labs(tag = "CR"), Figure4B[[2]] + labs(tag = " "), ncol = 2, nrow = 2)
  
}
dev.off()




