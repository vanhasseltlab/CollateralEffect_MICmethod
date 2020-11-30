library(tidyverse)

#ran all dichotomization criteria
quants <- c(4:9*0.1)
columns <- c("A", "B", "t", "tau")

visualize_df <- data.frame()
for (quant in quants) {
  load(paste0("results/results_Escherichia-coli_", quant, ".Rdata"))
  visualize_df <- rbind(visualize_df, data.frame(t_test_results[, columns], quant = quant))
}


visualize_df %>% 
  ggplot(aes(x = quant, y = t, colour = B)) +
  geom_line() +
  facet_wrap(~A, nrow = 5) +
  theme_bw()


load(paste0("results/results_Escherichia-coli_median.Rdata"))


visualize_df_med <- visualize_df %>% 
  right_join(t_test_results %>% select(A, B, tau) %>% mutate(tau = tau + 0.5))

View(visualize_df_med[is.na(visualize_df_med$t), ])


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


over_tau_med <- rbind(visualize_df, t_test_results[, columns] %>%
                        mutate(quant = 0.45, tau = tau + 0.5)) %>% 
  
  ggplot(aes(x = tau, y = t, colour = B)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_line(alpha = 0.7) +
  geom_point(data = (t_test_results %>% mutate(tau = tau + 0.5))) +
  scale_x_continuous(breaks = -2:8) +
  labs(x = expression("Dichotomization criterium ("*log[2]*"(MIC))"), y = "T-value") +
  facet_wrap(~A, nrow = 5) +
  theme_bw() + theme(panel.grid.minor.x = element_blank())

pdf("results/figures/effect_of_dichotomization_tau.pdf", width = 8, height = 10)
print(over_tau_med)
dev.off()


over_tau_med_flipped <- rbind(visualize_df, t_test_results[, columns] %>%
                        mutate(quant = 0.45, tau = tau + 0.5)) %>% 
  
  ggplot(aes(x = tau, y = t, colour = A)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_line(alpha = 0.7) +
  geom_point(data = (t_test_results %>% mutate(tau = tau + 0.5))) +
  scale_x_continuous(breaks = -2:8) +
  labs(x = expression("Dichotomization criterium ("*tau*")"), y = "T-value") +
  facet_wrap(~B, nrow = 5, scales = "free_x") +
  theme_bw() + theme(panel.grid.minor.x = element_blank())

pdf("results/figures/effect_of_dichotomization_tau_flipped.pdf", width = 8, height = 10)
print(over_tau_med_flipped)
dev.off()



over_quants_med_flipped <- visualize_df %>% 
  ggplot(aes(x = quant, y = t, colour = A)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_point(data = visualize_df_med) +
  geom_line(alpha = 0.7) +
  labs(x = "Dichotomization quantile", y = "T-value") +
  facet_wrap(~B, nrow = 5) +
  theme_bw() + theme(panel.grid.minor.x = element_blank())

pdf("results/figures/effect_of_dichotomization_flipped.pdf", width = 8, height = 10)
print(over_quants_med_flipped)
dev.off()
