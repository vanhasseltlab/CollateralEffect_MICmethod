#Create plots for manuscript

#Libraries and functions
library(tidyverse)
source("scripts/functions_CR.R")



load("data/clean/MIC_clean_Escherichia-coli.Rdata")

#load results median E. coli
load("results/results_Escherichia-coli_median.Rdata")
CE_results_median <- t_test_results



#Main figures
meta_antibio_abbr <- read.csv("data/antibiotics_lookup_table.csv", header = T, stringsAsFactors = F, sep = ";")

t_test_results_class <- t_test_results %>% 
  left_join(meta_antibio_abbr %>% 
              mutate(A = abbreviation, class_A = class) %>% 
              select(c(A, class_A))) %>% 
  left_join(meta_antibio_abbr %>% 
              mutate(B = abbreviation, class_B = class) %>% 
              select(c(B, class_B))) %>% 
  mutate(A = factor(A, levels = unique(A[order(class_A)]))) %>% 
  mutate(B = factor(B, levels = unique(B[order(class_B)])))


supltable1 <- data.frame(n_strains = colSums(!is.na(MIC_clean)), abbreviation = names(MIC_clean)) %>% 
  left_join(meta_antibio_abbr %>% select(abbreviation, antibiotic)) %>% 
  mutate(antibiotic = str_to_sentence(antibiotic)) %>% 
  select(c(antibiotic, abbreviation,n_strains)) %>% 
  rename(Abbreviation = abbreviation, Antibiotic = antibiotic, `Number of strains` = n_strains) 
  
rownames(supltable1) <- supltable1$Abbreviation


supltable1 <- supltable1[levels(t_test_results_class$A), ]
write.csv(supltable1, file = "paper/figures/SupplTable1.csv", row.names = F)


# x$V1 <- factor(x$V1, levels=(x$V1)[order(x$V3)])
t_test_results_class %>% mutate(A = factor(A, levels = unique(A[order(class_A)])))
cairo_pdf("paper/figures/Figure3_overview_all_tests.pdf", width = 9, height =  7)
PlotSignificantEffect(t_test_results_class, FDR_crit = 0.05, species = "E. coli")
dev.off()

pdf("paper/figures/Figure4A_CS_t-test.pdf", width = 7, height = 4)
PlotCondDistribution(MIC_clean, t_test_results, t_rank = 1, one_direction = FALSE, CResponse = "CS", FDR_crit = 0.05,
                     colours = c("#BFC6B8", "#4A5242"))
dev.off()

pdf("paper/figures/Figure4B_CR_t-test.pdf", width = 7, height = 4)
PlotCondDistribution(MIC_clean, t_test_results, t_rank = 1, one_direction = FALSE, CResponse = "CR", FDR_crit = 0.05,
                     colours = c("#BFC6B8", "#4A5242"))
dev.off()


#panelplot

Figure4A <- PlotCondDistribution(MIC_clean, t_test_results, t_rank = 1, one_direction = FALSE, CResponse = "CS", FDR_crit = 0.05,
                     colours = c("#BFC6B8", "#4A5242"), separate_plots = TRUE, ran = c(-6, 7))

Figure4B <- PlotCondDistribution(MIC_clean, t_test_results, t_rank = 1, one_direction = FALSE, CResponse = "CR", FDR_crit = 0.05,
                                 colours = c("#BFC6B8", "#4A5242"), separate_plots = TRUE, ran = c(-6, 7))

pdf("paper/figures/Figure4_CSCR_t-test.pdf", width = 7, height = 6.5)
grid.arrange(Figure4A[[1]] + labs(tag = "A"), Figure4A[[2]] + labs(tag = "B"), 
             Figure4B[[1]] + labs(tag = "C"), Figure4B[[2]] + labs(tag = "D"), ncol = 2, nrow = 2)
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

visualize_df_med <- visualize_df %>% 
  mutate(tau_10 = round(tau*10)) %>% 
  select(-tau) %>% 
  right_join(CE_results_median %>% select(A, B, tau) %>% mutate(tau_10 = round((tau + 0.5)*10)))

visualize_df_med <- visualize_df %>% 
  mutate(t_round_100 = round(t*10000)) %>% 
  select(-tau) %>% 
  right_join(CE_results_median %>% 
               mutate(t_round_100 = round(t*10000), tau_between = tau) %>% 
               select(A, B, t_round_100, tau_between), na_matches = "never")



over_quants_med <- visualize_df %>% 
  ggplot(aes(x = quant, y = t, colour = B)) +
  geom_hline(yintercept = 0, colour = "black") +
  
  geom_line(alpha = 0.7) +
  geom_point(data = visualize_df_med, alpha = 0.9, stroke = 0, size = 1.8) +
  labs(x = "Dichotomization quantile", y = "T-value", colour = "Splitting antibiotic (B)") +
  facet_wrap(~A, nrow = 5) +
  theme_bw() + theme(panel.grid.minor.x = element_blank())

pdf("results/figures/effect_of_dichotomization_new.pdf", width = 8.5, height = 10)
print(over_quants_med)
dev.off()


over_quants_med_flipped <- visualize_df %>% 
  ggplot(aes(x = quant, y = t, colour = A)) +
  geom_hline(yintercept = 0, colour = "black") +
  
  geom_line(alpha = 0.7) +
  geom_point(data = visualize_df_med, alpha = 0.9, stroke = 0, size = 1.8) +
  labs(x = "Dichotomization quantile", y = "T-value") +
  facet_wrap(~B, nrow = 5) +
  theme_bw() + theme(panel.grid.minor.x = element_blank())
