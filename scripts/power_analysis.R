#power analysis CS test
library(pwr)
library(tidyverse)
library(RColorBrewer)
library(egg)
#columns n_total, n1, n2, frac_disbalance, power, effect size

total_n <- c(20, 30, 50, 70, 100, 200, 500, 1000, 2000)
fracs <- seq(0.1, 0.5, by = 0.05)
sds <- seq(0.6, 2.4, by = 0.2) #based on s in MIC clean E coli PATRIC
effect_sizes <- seq(0, 2, by = 0.05) #based on t results E coli PATRIC

power_results <- expand.grid(total_n, fracs, effect_sizes, sds)

colnames(power_results) <- c("n", "balance", "effect_size", "sd")
power_results[, "power"] <- 0
power_results[, "d"] <- power_results$effect_size/power_results$sd
power_results[, "n1"] <- ceiling(power_results$n * power_results$balance - 0.001)
power_results[, "n2"] <- floor(power_results$n * (1 - power_results$balance) + 0.001)

for (i in 1:nrow(power_results)) {
  row_i <- power_results[i, ]
  power_disb <- pwr.t2n.test(n1 = row_i$n1, n2 = row_i$n2, d = row_i$effect_size/row_i$sd, sig.level = 0.05)
  power_results[i, "power"] <- power_disb$power
}

View(power_results %>% filter(balance == 0.5 & n == 200 & sd > 1.7 & sd < 1.9))

#facets on n
power_plot_n <- power_results %>% 
  mutate(disbalance = paste(balance*100, "/", 100*(1 - balance))) %>% 
  filter(abs(sd - 1.8) < 0.01) %>%
  ggplot(aes(y = power, x = effect_size, group = balance, colour = as.factor(disbalance))) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey65")  +
  geom_line() +
  #scale_x_continuous(breaks = seq(0, 2, by = 0.25)) +
  scale_color_manual(values = c(brewer.pal(n = 9, name = "BuPu")[-c(1)], "black"), name = "disbalance (%)") +
  labs(title = "Power analysis", x = "Effect size (log(FC))", y = "Power") +
  facet_wrap(~n) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank())

#facets on disbalance
power_plot_balance <- power_results %>% 
  filter(abs(sd - 1.8) < 0.01) %>%
  mutate(disbalance = paste(balance*100, "/", 100*(1 - balance))) %>% 
  filter(balance %in% seq(0.1, 0.5, by = 0.1)) %>% 
  filter(!n %in% c(30)) %>% 
  ggplot(aes(y = power, x = effect_size, group = n, colour = as.factor(n))) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey65")  +
  geom_line() +
 # scale_x_continuous(breaks = seq(0, 2, by = 0.25)) +
  scale_color_manual(values = c(brewer.pal(n = 9, name = "BuPu")[-c(1, 2)], "black"), name = "n") +
 
  labs(title = "Power analysis", x = "Effect size (log(FC))", y = "Power") +
  facet_wrap(~disbalance) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank())

pdf("results/figures/power_analysis.pdf", width = 9)
print(power_plot_n)
print(power_plot_balance)
dev.off()

#figures final
bl <- colorRampPalette(c("#BFD3E6", "#8C96C6", "#88419D", "black"))(9)
plot_n <- power_results %>% 
  filter(abs(sd - 1.8) < 0.01) %>% 
  filter(balance == 0.5) %>% 
  mutate(disbalance = paste(balance*100, "/", 100*(1 - balance))) %>% 
  ggplot(aes(y = power, x = effect_size, group = n, colour = as.factor(n))) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey65")  +
  #scale_color_manual(values = c(brewer.pal(n = 9, name = "BuPu")[-c(1)], "black"), name = "n") +
  scale_color_manual(values = bl, name = "n") +
  labs(x = "", y = "", title = "sd = 1.8, disbalance = 50/50") +
  geom_line() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.title=element_blank())

re <- colorRampPalette(c("#FEC6C6", "#D60404", "black"))(9)
plot_disbalance <- power_results %>% 
  filter(abs(sd - 1.8) < 0.01) %>% 
  filter(n == 200) %>% 
  mutate(disbalance = paste(balance*100, "/", 100*(1 - balance))) %>% 
  ggplot(aes(y = power, x = effect_size, group = disbalance, colour = disbalance)) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey65")  +
  #scale_color_manual(values = c(brewer.pal(n = 9, name = "BuPu")[-c(1)], "black"), name = "n") +
  scale_color_manual(values = re, name = "disbalance (%)") +
  labs(x = "", y = "", title = "n = 200, sd = 1.8") +
  geom_line() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.title=element_blank())

gr <- colorRampPalette(c("#BCE4AE", "#5FBB3F", "#3A8420", "black"))(10)
plot_sd <- power_results %>% 
  filter(balance == 0.5) %>% 
  filter(n == 200) %>% 
  mutate(disbalance = paste(balance*100, "/", 100*(1 - balance))) %>% 
  ggplot(aes(y = power, x = effect_size, group = sd, colour = as.factor(sd))) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey65")  +
  #scale_color_manual(values = c(brewer.pal(n = 9, name = "BuPu")[-c(1)], "black"), name = "n") +
  scale_color_manual(values = gr, name = "sd") +
  labs(x = "", y = "", title = "n = 200, disbalance = 50/50") +
  geom_line() + 
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.title=element_blank())


#ggarrange(plot_n, plot_disbalance, plot_sd, nrow = 1, bottom = "Effect size", left = "Power")

pdf("results/figures/power_analysis_summary.pdf", width = 12, height = 5, onefile = FALSE)
# print(plot_n)
# print(plot_disbalance)
# print(plot_sd)
ggarrange(plot_n, plot_disbalance, plot_sd, nrow = 1, bottom = "Effect size (log2(FC))", left = "Power", 
          labels = LETTERS[1:3], label.args = list(gp=grid::gpar(font=1), hjust=1, vjust = 1.5))
dev.off()

power_results %>% 
  ggplot(aes(y = power, x = 2^effect_size, group = balance, colour = as.factor(balance))) +
  geom_line() +
  #scale_x_continuous(breaks = seq(0, 2, by = 0.25)) +
  scale_colour_brewer(palette="BuPu") +
  geom_hline(yintercept = 0.8) +
  labs(title = "Power analysis", x = "Effect size (FC)", y = "Power") +
  facet_wrap(~n) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank())


power_results %>% 
  mutate(disbalance = paste(balance*100, "/", 100*(1 - balance))) %>% 
  filter(balance == 0.5) %>% 
  filter(!n %in% c(30)) %>% 
  ggplot(aes(y = power, x = effect_size, group = n, colour = as.factor(n))) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 2, by = 0.25)) +
  scale_color_manual(values = c(brewer.pal(n = 9, name = "BuPu")[-c(1, 2)], "black")) +
  geom_hline(yintercept = 0.8) +
  labs(title = "Power analysis", x = "Effect size (log(FC))", y = "Power") +
  facet_wrap(~sd) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank())

