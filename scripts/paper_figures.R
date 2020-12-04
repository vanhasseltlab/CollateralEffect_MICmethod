#Create plots for manuscript

#Libraries and functions
library(tidyverse)
source("scripts/functions_CR.R")



load("data/clean/MIC_clean_Escherichia-coli.Rdata")

#load results median E. coli
load("results/results_Escherichia-coli_median.Rdata")
CE_results_median <- t_test_results



#Main figures
pdf("paper/figures/Figure1_overview_all_tests.pdf", width = 9, height =  7)
PlotSignificantEffect(t_test_results, FDR_crit = 0.05, species = species)
dev.off()

pdf("paper/figures/Figure2_CS_t-test.pdf", width = 7, height = 4)
PlotCondDistribution(MIC_clean, t_test_results, t_rank = 1, one_direction = FALSE, CResponse = "CS", FDR_crit = 0.05,
                     colours = c("#BFC6B8", "#4A5242"))
dev.off()

pdf("paper/figures/Figure3_CR_t-test.pdf", width = 7, height = 4)
PlotCondDistribution(MIC_clean, t_test_results, t_rank = 1, one_direction = FALSE, CResponse = "CR", FDR_crit = 0.05,
                     colours = c("#BFC6B8", "#4A5242"))
dev.off()






#load and combine results over different quantiles
#ran all dichotomization criteria
quants <- c(4:9*0.1)
columns <- c("A", "B", "t", "tau")

visualize_df <- data.frame()
for (quant in quants) {
  load(paste0("results/results_Escherichia-coli_", quant, ".Rdata"))
  visualize_df <- rbind(visualize_df, data.frame(t_test_results[, columns], quant = quant))
}

#Figure S1
visualize_df_med <- visualize_df %>% 
  right_join(CE_results_median %>% select(A, B, tau) %>% mutate(tau = tau + 0.5))

over_quants_med <- visualize_df %>% 
  ggplot(aes(x = quant, y = t, colour = B)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_point(data = visualize_df_med) +
  geom_line(alpha = 0.7) +
  labs(x = "Dichotomization quantile", y = "T-value") +
  facet_wrap(~A, nrow = 5) +
  theme_bw() + theme(panel.grid.minor.x = element_blank())

pdf("results/figures/effect_of_dichotomization_new.pdf", width = 8, height = 10)
print(over_quants_med)
dev.off()