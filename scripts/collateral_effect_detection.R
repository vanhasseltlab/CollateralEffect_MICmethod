
################################################################################
# Data analysis for replication of "Identification of collateral effects of    #
# antibiotic resistance in population surveillance data"                       #
# Author: Laura B. Zwep                                                        #
################################################################################

######  Libraries #####
library(tidyverse)
library(collatRal)
library(pwr)
library(egg)
library(grid)

#####  Import data  #####
#Load information files antibiotics
meta_antibio_abbr <- read.csv("data/meta_antibio_abbr.csv", header = T,
                              stringsAsFactors = F)
suppl_table_1 <- read.csv("data/suppl_table_1.csv", header = T, sep = ",",
                          stringsAsFactors = F)
#Choose minimal number of MIC measurements to include drug
na_remove <- 200

#Load the raw data (long format) 
raw_MIC <- read.table("data/raw/Escherichia-coli.txt", header = T, sep = "\t", 
                      dec = ",", stringsAsFactors = F, na.strings = "")
colnames(raw_MIC) <- sapply(strsplit(colnames(raw_MIC), "\\."), function(x) x[2])

#####  Data cleaning  #####
#Create clean MIC values data frame
MIC_df <- raw_MIC %>%
  #remove empty MIC measurements
  filter(!is.na(measurement_value)) %>% 
  #clean strings
  mutate(antibiotic = str_replace_all(antibiotic, "-", "/"),
         genome_id = str_replace_all(genome_id, "\\.", "_")) %>% 
  #add antibiotic abbreviations
  left_join(meta_antibio_abbr, by = "antibiotic") %>% 
  #take MIC value from antibiotic (not beta lactamase blocker)
  mutate(MIC = sapply(strsplit(measurement_value, "/"), 
                      function(x) as.numeric(x[1]))) %>% 
  select(genome_id, abbreviation, MIC)

#Create wide format of MIC values
MIC_table <- MIC_df %>% 
  pivot_wider(names_from = abbreviation, values_from = MIC) %>% 
  column_to_rownames(var = "genome_id")

#Remove antibiotics with to many na's
n_antibiotic <- apply(MIC_table, 2, function(x) sum(!is.na(x)))
MIC_clean <- MIC_table %>% 
  select(which(n_antibiotic > na_remove)) %>% 
  filter(rowSums(!is.na(.)) > 1)

#Save cleaned data for further use
save(MIC_clean, file = paste0("data/clean/MIC_clean_Escherichia-coli.Rdata"))

#####  Data analysis #####
t_test_results <- collateral_mult_test(MIC_clean)


#####  Create figures  #####
# Figure 3
#create_ordering according to antibiotic class
suppl_table_1 <- read.csv("data/suppl_table_1.csv", header = T, stringsAsFactors = F, sep = ",")
suppl_table_1 <- suppl_table_1 %>% arrange(target, class, ab)
order_abs <- suppl_table_1$ab
t_test_results_class <- t_test_results %>% 
  mutate(A = factor(A, levels = order_abs)) %>% 
  mutate(B = factor(B, levels = order_abs))

var_i <- suppl_table_1$class
bold_lines <- cumsum(as.vector(table(var_i)[unique(var_i)]))[-length(unique(var_i))] + 0.5
figure3 <- plot_heatmap_CE(t_test_results_class, sign_criterium = 0.05) +
  geom_hline(yintercept = bold_lines, size = 1) +
  geom_vline(xintercept = bold_lines, size = 1)

cairo_pdf("paper/figures/Figure3_overview_all_tests.pdf", width = 9, height =  7)
print(figure3)
dev.off()

# Figure 4
figure4a <- plot_histogram_CE(MIC_clean$CFZ, MIC_clean$ETP, MIC_range = c(-6, 7))
figure4b <- plot_histogram_CE(MIC_clean$ETP, MIC_clean$CFZ, MIC_range = c(-6, 7))
figure4c <- plot_histogram_CE(MIC_clean$ETP, MIC_clean$MEM, MIC_range = c(-6, 7))
figure4d <- plot_histogram_CE(MIC_clean$MEM, MIC_clean$ETP, MIC_range = c(-6, 7))

pdf("paper/figures/Figure4_CSCR_t-test.pdf", width = 7, height = 6.5)
gridExtra::grid.arrange(figure4a + labs(tag = "A"), figure4b + labs(tag = "B"), 
             figure4c + labs(tag = "C"), figure4d + labs(tag = "D"), ncol = 2, 
             nrow = 2)
dev.off()




################################################################################
# Power analysis                                                               #
################################################################################

#####  Power analysis  #####
#power analysis CS test

library(tidyverse)
#columns n_total, n1, n2, frac_disbalance, power, effect size

total_n <- c(20, 30, 50, 70, 100, 200, 500, 1000, 2000)
fracs <- seq(0.1, 0.5, by = 0.05)
sds <- seq(0.6, 2.4, by = 0.2) #based on s in MIC clean E coli PATRIC
effect_sizes <- seq(0, 2, by = 0.05) #based on t results E coli PATRIC

power_results <- expand.grid(total_n, fracs, effect_sizes, sds) %>% 
  rename(n = Var1, balance = Var2, effect_size = Var3, sd = Var4) %>% 
  mutate(power = 0, d = effect_size/sd, 
         n1 = ceiling(n * balance - 0.001),
         n2 = floor(n * (1 - balance) + 0.001))

for (i in 1:nrow(power_results)) {
  row_i <- power_results[i, ]
  power_disb <- pwr.t2n.test(n1 = row_i$n1, n2 = row_i$n2, 
                             d = row_i$effect_size/row_i$sd, sig.level = 0.05)
  power_results[i, "power"] <- power_disb$power
}



#####  Create figures  #####
bl <- colorRampPalette(c("#BFD3E6", "#8C96C6", "#88419D", "black"))(9)
plot_n <- power_results %>% 
  filter(abs(sd - 1.8) < 0.01 & balance == 0.5) %>% 
  ggplot(aes(y = power, x = effect_size, group = n, colour = as.factor(n))) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey65")  +
  scale_color_manual(values = bl, name = "n") +
  labs(x = "", y = "", title = "sd = 1.8, disbalance = 50/50") +
  geom_line() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.title = element_blank())

re <- colorRampPalette(c("#FEC6C6", "#D60404", "black"))(9)
plot_disbalance <- power_results %>% 
  filter(abs(sd - 1.8) < 0.01 & n == 200) %>% 
  mutate(disbalance = paste(balance*100, "/", 100*(1 - balance))) %>% 
  ggplot(aes(y = power, x = effect_size, group = disbalance, colour = disbalance)) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey65")  +
  scale_color_manual(values = re, name = "disbalance (%)") +
  labs(x = "", y = "", title = "n = 200, sd = 1.8") +
  geom_line() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.title = element_blank())

gr <- colorRampPalette(c("#BCE4AE", "#5FBB3F", "#3A8420", "black"))(10)
plot_sd <- power_results %>% 
  filter(balance == 0.5 & n == 200) %>% 
  ggplot(aes(y = power, x = effect_size, group = sd, colour = as.factor(sd))) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey65")  +
  scale_color_manual(values = gr, name = "sd") +
  labs(x = "", y = "", title = "n = 200, disbalance = 50/50") +
  geom_line() + 
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.title = element_blank())


pdf("paper/figures/Figure5_power_analysis.pdf", width = 12, height = 5, onefile = FALSE)
egg::ggarrange(plot_n, plot_disbalance, plot_sd, nrow = 1, bottom = "Effect size (log2(FC))", 
          left = "Power", labels = LETTERS[1:3], 
          label.args = list(gp = grid::gpar(font = 1), hjust = 1, vjust = 1.5))
dev.off()

################################################################################
# Test different values of tau                                                 #
################################################################################

t_results_quants <- data.frame()
criteria <- c(seq(0.2, 0.9, 0.1), 0.95)
for (i in 1:length(criteria)) {
  t_test_crit <- collateral_mult_test(MIC_clean, crit_type = "quant", criterium = criteria[i])
  t_test_crit$quantile <- criteria[i]
  t_test_crit$is_median <- t_test_crit$tau == t_test_results$tau
  t_results_quants <- rbind(t_results_quants, t_test_crit)
}

figure6 <- t_results_quants %>% 
  ggplot(aes(x = quantile, y = t, colour = B)) +
    geom_hline(yintercept = 0, colour = "black") +
    geom_line(alpha = 0.7) +
    geom_point(data = t_results_quants[t_results_quants$is_median, ], 
               alpha = 0.9, stroke = 0, size = 1.8) +
    labs(x = "Dichotomization quantile", y = "T-value", 
         colour = "Splitting antibiotic (B)") +
    facet_wrap(~ A, nrow = 5) +
    theme_bw() + 
   theme(panel.grid.minor.x = element_blank())
pdf("paper/figures/Figure6_effect_of_dichotomization.pdf", width = 8.5, height = 10)
print(figure6)
dev.off()

### try another way! #####

t_results_taus <- data.frame()
m <- length(names(MIC_clean)) -1
for (b in names(MIC_clean)) {
  tau_values <- log2(sort(unique(MIC_clean[, b]))[-1])
  t_results_b <- as.data.frame(matrix(NA, nrow = m*length(tau_values), ncol = 12))
  for (i in 1:length(tau_values)) {
    t_test_crit <- collateral_mult_test(MIC_clean, crit_type = "log2_MIC", criterium = tau_values[i])
    t_results_b[(((i - 1)*m) + 1):(i*m), ] <- t_test_crit[t_test_crit$B == b, ]
  }
  colnames(t_results_b) <- colnames(t_test_crit)
  
  t_results_b$is_median <- FALSE
  t_results_taus <- rbind(t_results_taus, t_results_b)
}

t_results_taus <- rbind(t_results_taus, t_test_results %>% mutate(is_median = TRUE))

figure6 <- t_results_taus %>% 
  ggplot(aes(x = tau, y = t, colour = A)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_line(alpha = 0.7) +
  geom_point(data = t_results_taus[t_results_taus$is_median, ], 
             alpha = 0.9, stroke = 0, size = 1.8) +
  labs(x = "Dichotomization criterium", y = "T-value", 
       colour = "Testing antibiotic (A)") +
  facet_wrap(~ B, nrow = 5, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank())
pdf("paper/figures/Figure6_effect_of_dichotomization.pdf", width = 8.5, height = 10)
print(figure6)
dev.off()