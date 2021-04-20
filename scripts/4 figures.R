###analyze results
library(tidyverse)
meta_antibio_abbr <- read.csv("data/antibiotics_lookup_table.csv", header = T, stringsAsFactors = F, sep = ";")
load("data/clean/MIC_clean_Escherichia-coli.Rdata")
load("results/results_Escherichia-coli_median.Rdata")
source("scripts/functions_CR.R")

r5 <- results$`0.5`

#CS
sum(r5$p_BY < 0.05 &  r5$t < 0)

#CR
sum(r5$p_BY < 0.05 &  r5$t > 0)
sum(r5$p_BY < 0.05)

top_results <- rbind(r5[order(r5$t, decreasing = T),][1:5, ], r5[order(r5$t),][1:5, ])


#figure 1
dat <- r5
dat$effect_size[dat$p_BY > 0.05] <- 0
dat$t[dat$p_BY > 0.05] <- 0

bl <- colorRampPalette(c("#283c82", "white"))(30) [c(1:10, seq(11, 30, by = 3))] 
re <- colorRampPalette(c("#f54c00", "white"))(30)[c(1:10, seq(11, 30, by = 3))]
limits <- c(-1, 1)*max(dat$effect_size)

bb <- limits
ll <- c("Collateral Sensitivity", "Collateral Resistance") # labels.

dat$Direction <- ifelse(dat$effect_size < 0, ll[1], ll[2])

plot_ <- ggplot(dat, aes(x = B, y = A, color = effect_size, shape = Direction)) +
  geom_point(shape = 15, size = 11) +
  #geom_point(shape = 15, size = 10) +
  scale_colour_gradientn(colours = c(bl, "white", rev(re)), limits = limits) +
  scale_x_discrete(expand = expansion(mult = 0, add = rep(0.5, 2))) +
  scale_y_discrete(expand = expansion(mult = 0, add = rep(0.5, 2))) +
  geom_point(size = 4, colour = "white") +
  scale_shape_manual(values = c("+", "-")) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Antibiotic A", x = paste0("Antibiotic B (conditioned on)"), 
       size = "Difference in means", colour = "Difference between\nmeans", shape = "Direction") +
  geom_vline(xintercept = seq(1.5, length(unique(dat$A)) - 0.5, 1), colour = "grey60") +
  geom_hline(yintercept = seq(1.5, length(unique(dat$B)) - 0.5, 1), colour = "grey60") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12), legend.key = element_rect(fill = "grey60"))

plot_ <- plot_ + geom_abline(slope = 1, intercept = 0, colour = "grey60")

pdf("results/figures/Figure1_significant_Escherichia-coli.pdf", width = 9, height =  7)
print(plot_)
dev.off()



#Figure 2
pdf("results/figures/Figure2_CS_t-test.pdf", width = 7, height = 4)
PlotCondDistribution(MIC_clean, results, 0.5, t_rank = 1, one_direction = FALSE, CResponse = "CS", whichAB = c("CFZ", "ETP"),
                     colours = c("#BFC6B8", "#4A5242"))
dev.off()

#Figure 3
pdf("results/figures/Figure3_CR_t-test.pdf", width = 7, height = 4)
PlotCondDistribution(MIC_clean, results, 0.5, t_rank = 1, one_direction = FALSE, CResponse = "CR", whichAB = c("ETP", "MEM"),
                     colours = c("#BFC6B8", "#4A5242"))
dev.off()


#plots presentation
library(tidyverse)

meta_antibio_abbr <- read.csv("data/meta_antibio_abbr.csv", header = T, stringsAsFactors = F)
load("data/clean/MIC_clean_Escherichia-coli.Rdata")
load("results/resultsEscherichia-coli_cond.Rdata")


#example of marginal distributions
MIC_long <- MIC_clean %>% 
  gather(key = "abbreviation", value = "MIC") %>% 
  left_join(meta_antibio_abbr) %>%
  filter(!is.na(MIC)) %>% 
  filter(abbreviation %in% c("CFZ", "AMC", "FOX", "FEP")) %>% 
  mutate(Antibiotic = paste0(tools::toTitleCase(antibiotic), " (", abbreviation, ")"))

pdf(file = "results/figures/example_MIC_distributions.pdf", height = 4, width = 6)
ggplot(MIC_long, aes(x = log2(MIC), y = ..prop..)) +
  geom_bar(stat = "count", width = 0.6, position = "stack", fill = "#283c82") +
  facet_wrap(~ Antibiotic, ncol = 2) +
  labs(x = expression(log[2]*"(MIC)"), y = "Proportion") +
  scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
  scale_x_continuous(breaks = min(log2(MIC_long$MIC)):max(log2(MIC_long$MIC))) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), strip.background = element_blank(),
        strip.text.x = element_text(angle = 0, hjust = 0))

dev.off()


###plots presentation
#example of joint distribution
test_df <- MIC_clean %>% select(c(TZP, GEN)) %>% filter(!is.na(TZP) & !is.na(GEN))
plot_1 <- ggplot(test_df, aes(x =log2(TZP), y = log2(GEN))) +
  geom_count(aes(size = after_stat(prop)), shape = 15) +
  scale_size_area(max_size = 10) +
  theme_bw()

pdf("results/figures/correlation_plot_AMK_FOX.pdf", width = 5, height = 3)
print(plot_1)
print(plot_1 + geom_smooth(method = "lm", se = F))
dev.off()

mic_dat <- MIC_clean %>% select(selectedAB) %>% 
  filter_all(any_vars(!is.na(.)))


#overview results marginal vs conditional

selectedAB <- c("TZP", "GEN", "LVX", "FOX", "FEP", "CRO", "CIP", "CFZ", "CAZ", "AMC")

dat <- results[["0.5"]]
dat$effect_size[dat$p_BY > 0.05] <- 0
dat$t[dat$p_BY > 0.05] <- 0
dat <- dat %>% 
  filter(A %in% selectedAB & B %in% selectedAB)

bl <- colorRampPalette(c("#283c82", "white"))(30) [c(1:10, seq(11, 30, by = 3))] 
re <- colorRampPalette(c("#f54c00", "white"))(30)[c(1:10, seq(11, 30, by = 3))]
limits <- c(-1, 1)*max(dat$effect_size)

bb <- limits
ll <- c("Collateral Sensitivity", "Collateral Resistance") # labels.

dat$Direction <- ifelse(dat$effect_size < 0, ll[1], ll[2])

plot_ <- ggplot(dat, aes(x = B, y = A, color = effect_size, shape = Direction)) +
  geom_point(shape = 15, size = 14) +
  #geom_point(shape = 15, size = 10) +
  scale_colour_gradientn(colours = c(bl, "white", rev(re)), limits = limits) +
  scale_x_discrete(expand = expansion(mult = 0, add = rep(0.5, 2))) +
  scale_y_discrete(expand = expansion(mult = 0, add = rep(0.5, 2))) +
  geom_point(size = 4, colour = "white") +
  scale_shape_manual(values = c("+", "-")) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Antibiotic A", x = paste0("Antibiotic B (conditioned on)"), 
       size = "Difference in means", colour = "Difference between\nmeans", shape = "Direction") +
  geom_vline(xintercept = seq(1.5, length(unique(dat$A)) - 0.5, 1), colour = "grey60") +
  geom_hline(yintercept = seq(1.5, length(unique(dat$B)) - 0.5, 1), colour = "grey60") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12), legend.key = element_rect(fill = "grey60"))

plot_ <- plot_ + geom_abline(slope = 1, intercept = 0, colour = "grey60")

pdf("results/figures/ppt_significant_selected_Escherichia-coli.pdf", width = 7, height = 5)
print(plot_)
dev.off()


#plot separate results
source("scripts/functions_CR.R")
PlotStackDistribution(MIC_clean, results, 0.5, t_rank = 1, one_direction = FALSE, CResponse = "CS")
PlotStackDistribution(MIC_clean, results, 0.5, t_rank = 1, one_direction = FALSE, CResponse = "CR", whichAB = c("LVX", "CIP"))
PlotStackDistribution(MIC_clean, results, 0.5, t_rank = 1, one_direction = FALSE, CResponse = "CR", whichAB = c("TZP", "CFZ"))

pdf("results/figures/ppt_t-test-distributions.pdf", width = 7, height = 3.5)
PlotCondDistribution(MIC_clean, results, 0.5, t_rank = 1, one_direction = FALSE, CResponse = "CR", whichAB = c("CFZ", "CAZ"), 
                     colours = c("#BFC6B8", "#4A5242"))
dev.off()

pdf("results/figures/ppt_t-test-distributions.pdf", width = 7, height = 3.5)
PlotCondDistribution(MIC_clean, results, 0.5, t_rank = 1, one_direction = FALSE, CResponse = "CR", whichAB = c("TZP", "CRO"), 
                     colours = c("#BFC6B8", "#4A5242"))
dev.off()

View(results$`0.5`)

CRTTest(A = MIC_clean$TZP, B = MIC_clean$CFZ, criterium = 0.5, CResponse = "both", crit_type = "quant")
CRTTest(A = MIC_clean$CFZ, B = MIC_clean$TZP, criterium = 0.5, CResponse = "both", crit_type = "quant")



result <- expand.grid(A = names(MIC_clean), B = names(MIC_clean), stringsAsFactors = F)
result <- result[result$A != result$B, ]

for (i in 1:nrow(result)) {
  t_test <- CRTTest(A = MIC_clean[, result$A[i]], B = MIC_clean[, result$B[i]], criterium = 0.5, CResponse = "both", crit_type = "quant")
  result$t <- t_test$statistic
  result$p <- t_test$p.value
  
  
}




r5$n_sens <- r5$n - r5$n_resi

hist(r5$n_resi/r5$n, breaks = 30)

r5 <- results$`0.5`
for(i in 1:nrow(r5)){
  r5$pair[i] <- paste0(sort(c(r5$A[i], r5$B[i])), collapse = "")
}

for (pair in unique(r5$pair)) {
  es <- r5[r5$pair == pair, "effect_size"]
  r5[r5$pair == pair, "efeect_size_diff"] <- es[1] - es[2]
}
r5$pair

r5 %>% group_by(pair)
  
  
MIC_clean %>% rownames_to_column(var = "strains") %>% 
  pivot_longer(-strains, names_to = "Antibiotic", values_to = "MIC") %>% 
  ggplot(aes(x = Antibiotic, y = log2(MIC))) +
  geom_violin(adjust = 2, fill = "#283C82", colour = "#283C82", trim = T) +
  theme_bw()+
  theme(panel.grid = element_blank())
  
  
MIC_clean$ID <- paste0("I", 1:nrow(MIC_clean))
MIC_clean %>% pivot_longer(-ID, names_to = "antibiotic", values_to = "MIC") %>% 
  filter(!is.na(MIC)) %>% mutate(log2MIC = round(log2(MIC))) %>% 
  ggplot(aes(x = log2MIC)) +
  geom_histogram() +
  facet_wrap(~antibiotic, ncol = 5) +
  theme_bw()