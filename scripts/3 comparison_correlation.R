#Exploration of findings: Would correlation be enough?
library(tidyverse)
#results T-test
load("results/test_resultsPseudomonas-aeruginosa.Rdata")
#Clean MIC data
load("data/clean/MIC_clean_Pseudomonas-aeruginosa.Rdata")
load("data/clean/MIC_clean_Escherichia-coli.Rdata")


r5 <- results$`0.5`
r5$pair <- paste0(r5$A, r5$B)

combination <- c("CST", "TZP")
cs_dat <- MIC_clean[, combination]
names(cs_dat) <- c("A", "B")

result <- r5[r5$pair == paste0(combination, collapse = ""), ]

cs_dat <- cs_dat[!is.na(cs_dat$A) & !is.na(cs_dat$B), ]
d_min <- max(cs_dat$B[cs_dat$B < result$d], na.rm = T)

n <- nrow(cs_dat)

cs_dat$p_sizes <- sample(seq(0.1, 3, length.out = n))
cs_dat$AgivenBr <-  as.factor(ifelse(log2(cs_dat$B) >= result$d & !is.na(cs_dat$B), 
                                     paste0(combination[1],"|", combination[2]," > ", d_min), 
                                     paste0(combination[1])))



#plot MIC vs MIC with every circles representing one data point.
ggplot(cs_dat, aes(x = log2(A), y = log2(B))) +
  labs(x = paste0("log2(", combination[1], ")"), y = paste0("log2(", combination[2], ")")) +
  geom_point(shape = 1, aes(size = p_sizes), show.legend = F) +
  scale_size_continuous(range = c(0.1, 20)) +
#  geom_smooth(method = lm, se = F) +
  geom_hline(yintercept = result$d - 0.5) +
  theme_bw()

ggplot(cs_dat, aes(x = log2(A), y = log2(B))) +
  labs(x = paste0("log2(", combination[1], ")"), y = paste0("log2(", combination[2], ")")) +
  geom_point(shape = 1, aes(size = p_sizes), show.legend = F) +
  scale_size_continuous(range = c(0.1, 20)) +
  #  geom_smooth(method = lm, se = F) +
  geom_hline(yintercept = result$d - 0.5) +
  theme_bw()


ggplot(cs_dat, aes(x = log2(A), fill = AgivenBr, group = AgivenBr)) +
  geom_bar(width = 0.6, position = "stack") + 
  scale_fill_manual(values = c("#283c82", "#F46B2D")) +
  labs(x = expression(log[2]*"(MIC)"), y = "Counts", title = combination[1]) +
  theme_bw()

ggplot(cs_dat, aes(x = log2(A), fill = AgivenBr, group = AgivenBr)) +
  geom_bar(width = 0.6, position = "stack") + 
  scale_fill_manual(values = c("#283c82", "#F46B2D")) +
  labs(x = expression(log[2]*"(MIC)"), y = "Counts", title = combination[1]) +
  theme_bw()

ggplot(cs_dat, aes(x = log2(B), y = log2(A))) +
  labs(x = paste0("log2(", combination[2], ")"), y = paste0("log2(", combination[1], ")")) +
  geom_point(shape = 2, aes(size = p_sizes), show.legend = F) +
  geom_smooth(method = lm, se = F) +
  theme_bw()




cor(log2(cs_dat[, 1:2]), method = "kendall") #ranked correlation (log2 does not matter)
cor(log2(cs_dat[, 1:2]), method = "pearson") #standard correlation (log2 does matter)

############test correlations#########
source("scripts/functions_CR.R")
MIC_t <- log2(MIC_clean)

df_ps <- as.data.frame(expand.grid(names(MIC_t), names(MIC_t), stringsAsFactors = F))
df_ps <- df_ps[df_ps$Var1 != df_ps$Var2, ]

for (i in 1:nrow(df_ps)) {
  test_ <- cor.test(MIC_t[, df_ps$Var1[i]], MIC_t[, df_ps$Var2[i]], method = "kendall")
  df_ps$p_cor[i] <- test_$p.value
  df_ps$cor[i] <- test_$estimate
  test_t <- CRTTest(A = MIC_t[, df_ps$Var1[i]], B = MIC_t[, df_ps$Var2[i]], quant = 0.5)
  df_ps$p_t[i] <- test_t$p.value
  df_ps$t[i] <- test_t$statistic
  df_ps$df[i] <- test_t$parameter
}

#calculate t from correlation
df_ps$cor_t <- qt(df_ps$p_cor, df_ps$df, lower.tail = FALSE)
df_ps$cor_t <- abs(df_ps$cor_t) *sign(df_ps$cor) 


df_ps$sign_t <- df_ps$p_t < 0.05
df_ps$sign_cor <- df_ps$p_cor < 0.05

sum(df_ps$sign_t & df_ps$t < 0)
sum(df_ps$sign_cor & df_ps$cor < 0)
sum(df_ps$sign_t & df_ps$t > 0)
sum(df_ps$sign_cor & df_ps$cor > 0)

View(df_ps[df_ps$sign_cor != df_ps$sign_t, ])



pdf(file = "results/figures/corr/t_values_cor_vs_depttest.pdf", height = 5, width = 5)
plot1 <- ggplot(df_ps, aes(x = cor_t, y = t)) +
  geom_hline(yintercept = qnorm(c(0.975, 0.5, 0.025)), linetype = c(2,1,2), colour = "grey60") +
  geom_vline(xintercept = qnorm(c(0.975, 0.5, 0.025)), linetype = c(2,1,2), colour = "grey60") +
  geom_point(colour = "grey60", shape = 16) +
  geom_point(data = df_ps[df_ps$sign_cor, ], shape = 16, colour = "white") +
  geom_point(data = df_ps[df_ps$sign_t, ], colour = "#283c82", shape = 16) +
  geom_point(data = df_ps[df_ps$sign_cor, ], colour = "#F46B2D", shape = 1) +
  #scale_colour_manual(values = c("#283c82", "#F46B2D")) +

#  geom_smooth(method = "lm", colour = "dark red", linetype = 2) +
  labs(title = "T-statistics from correlation and dependent t-test", x = "t(correlation)", y = "t(t-test)") +
  theme_bw()
print(plot1)
dev.off()

mod_t <- lm(t ~ poly(cor_t, 1, raw = T), data = df_ps)
mod_t <- lm(t ~ cor_t, data = df_ps)

mod_seg <- segmented(mod_t, seg.Z = ~cor_t)
summary(mod_seg)

x <- data.frame(cor_t = seq(-5, 25, by = 0.01))
df_new <- data.frame(x, pred = predict(mod_seg, x))
plot1 +
  geom_line(data = df_new, aes(y = pred, x = cor_t)) +
  lims(y = range(df_new$cor_t))


table(df_ps[, c("sign_cor", "sign_t")])



df_ps$direction_cor <- ifelse(sign(df_ps$cor) < 0, "CS", "CR")
df_ps$direction_t <- ifelse(sign(df_ps$t) < 0, "CS", "CR")

table(df_ps[, c("direction_cor", "direction_t")])

#Number of significant CS in both correlation and T-test
sum(df_ps$sign_t & df_ps$t < 0 & df_ps$sign_cor & df_ps$cor < 0)

#Number of significant CR in both correlation and T-test
sum(df_ps$sign_t & df_ps$t > 0 & df_ps$sign_cor & df_ps$cor > 0)

#Number of significant CS in correlation but not significant in T-test
sum(df_ps$sign_cor & df_ps$cor < 0 & !df_ps$sign_t)

#Number of significant CR in correlation but not significant in T-test
sum(df_ps$sign_cor & df_ps$cor > 0 & !df_ps$sign_t)

#Number of significant CS in T-test but not significant in correlation
sum(df_ps$sign_t & df_ps$t < 0 & !df_ps$sign_cor)

#Number of significant CR in T-test but not significant in correlation
sum(df_ps$sign_t & df_ps$t > 0 & !df_ps$sign_cor)





df_ps$signed_logp_t <- sign(df_ps$t)*log10(df_ps$p_t)
df_ps$signed_logp_cor <- sign(df_ps$cor)*log10(df_ps$p_cor)
plot(df_ps$signed_logp_cor, df_ps$signed_logp_t)
abline(a = 0, b = 1, col = "red")
#abline(h = 0, v = 0, col = "blue")



mod_ps <- lm(df_ps$signed_logp_cor ~ df_ps$signed_logp_t)
abline(mod_ps, col = "blue")
summary(mod_ps)



df_ps$q <- p.adjust(df_ps$p, method = "BY")

#View(df_ps[df_ps$size < 0 & df_ps$q < 0.15, ])
#View(r5[r5$effect_size < 0 & r5$p_BY < 0.15, ])

df_ps$pair <- paste(df_ps[, 1], df_ps[, 2], sep = "_")
df_ps <- df_ps[order(df_ps$q, decreasing = F), ]

r5$pair <-  paste(r5$A, r5$B, sep = "_")

plot_list <- list()

for (pair in df_ps$pair[df_ps$q < 0.15 & df_ps$size < 0]) {
  combination <- strsplit(pair, split = "_")[[1]]
  cs_dat <- MIC_clean[, combination]
  
  names(cs_dat)[1:2] <- c("A", "B")
  
  result <- r5[r5$pair == pair, ]
  
  cs_dat <- cs_dat[!is.na(cs_dat[, 1]) & !is.na(cs_dat[, 2]), ]
  n <- nrow(cs_dat)
  
  cs_dat$p_sizes <- sample(seq(0.1, 3, length.out = n))
  
  #logtransform
  cs_dat[, 1:2] <- log2(cs_dat[, 1:2])

  plot_list[[pair]] <- ggplot(cs_dat, aes(x = A, y = B)) +
    geom_point(shape = 1, aes(size = p_sizes), show.legend = F) +
    scale_size_continuous(range = c(0.1, 20)) +
    geom_smooth(method = lm, se = F) +
    labs(title = paste0(pair, ", \nq_t: ", round(result$p_BY, 4), ", \t\tq_cor: ", round(df_ps[df_ps$pair == pair, "q"], 4),
                        "\np_t: ", round(result$p, 4), ", \t\tp_cor: ", round(df_ps[df_ps$pair == pair, "p"], 4))) +
   # geom_hline(yintercept = result$d - 0.5) +
    theme_bw()

}

pdf(file = "results/figures/corr/significant_kendall.pdf")
  for (i in 1:length(plot_list)) {
    print(plot_list[[i]])
  }
dev.off()

cor.test(cs_dat$A, cs_dat$B, method = "kendall")
cor.test(cs_dat$A, cs_dat$B, method = "pearson")

sum(cs_dat$NIT < 2)

names(results$`0.5`)
cs_dat

PlotStackDistribution(MIC_clean, results, 0.5, one_direction = FALSE, CResponse = "CS", whichAB = c("CFZ", "CIP"))
#plot all significant results correlations
names(df_ps) <- c("A", "B", "p", "effect_size", "p_BY", "pair")
results_cor <- list(`0.5` = df_ps)
PlotSignificantEffect(results_cor, 0.5, FDR_crit = 0.15, species = species)
PlotSignificantEffect(MIC_clean, results_cor, 0.5, one_direction = FALSE, CResponse = "CS", whichAB = c("CFZ", "CIP"))

df_ps$Direction <- ifelse(df_ps$effect_size > 0, "Collateral Sensitivity", "Collateral Resistance")
df_ps$effect_size[df_ps$p_BY > 0.15] <- 0

bl <- colorRampPalette(c("#283c82", "white"))(30) [c(1:10, seq(11, 30, by = 2))] 
#colorRampPalette(c("red","#d53397"))(30)
re <- colorRampPalette(c("#db0087", "white"))(30)[c(1:10, seq(11, 30, by = 2))]

plot_ <- ggplot(df_ps, aes(x = A, y = B, color = effect_size, shape = Direction)) +
  #geom_point(shape = 15, size = 13) +
  geom_point(shape = 15, size = 10) +
  scale_colour_gradientn(colours = c(re, "white", rev(bl)), limits = c(-1, 1)*max(df_ps$effect_size)) +
  scale_x_discrete(expand = expand_scale(mult = 0, add = rep(0.5, 2))) +
  scale_y_discrete(expand = expand_scale(mult = 0, add = rep(0.5, 2))) +
  geom_point(size = 4, colour = "white") +
  scale_shape_manual(values = c("-", "+"), labels = c("Collateral Sensitivity", "Collateral Resistance")) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Anitbiotic A", y = paste0("Antibiotic B"), colour = "Correlation (Kendall)", shape = "Direction",
       title = paste("Significant collateral responses", species)) +
  geom_vline(xintercept = seq(1.5, length(unique(df_ps$A)) - 0.5, 1), colour = "grey60") +
  geom_hline(yintercept = seq(1.5, length(unique(df_ps$B)) - 0.5, 1), colour = "grey60") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12), legend.key = element_rect(fill = "grey60"))

plot_ <- plot_ + geom_abline(slope = 1, intercept = 0, colour = "grey60")

plot_

r5 <- r5[order(r5$p_BY, decreasing = F), ]

#plot all significant from our analysis
plot_list2 <- plot_list <- list()
for (pair in r5$pair[r5$p_BY < 0.15 & r5$effect_size  < 0]) {
  combination <- strsplit(pair, split = "_")[[1]]
  cs_dat <- MIC_clean[, combination]
  
  names(cs_dat)[1:2] <- c("A", "B")
  
  result <- r5[r5$pair == pair, ]
  
  cs_dat <- cs_dat[!is.na(cs_dat[, 1]) & !is.na(cs_dat[, 2]), ]
  n <- nrow(cs_dat)
  
  cs_dat$p_sizes <- sample(seq(0.1, 3, length.out = n))
  
  #logtransform
  cs_dat[, 1:2] <- log2(cs_dat[, 1:2])
  
  plot_list[[pair]] <- ggplot(cs_dat, aes(x = A, y = B)) +
    geom_point(shape = 1, aes(size = p_sizes), show.legend = F) +
    scale_size_continuous(range = c(0.1, 20)) +
    geom_smooth(method = lm, se = F) +
    labs(title = paste0(pair, ", \nq_t: ", round(result$p_BY, 4), ", \t\tq_cor: ", round(df_ps[df_ps$pair == pair, "p_BY"], 4),
                        "\np_t: ", round(result$p, 4), ", \t\tp_cor: ", round(df_ps[df_ps$pair == pair, "p"], 4))) +
    # geom_hline(yintercept = result$d - 0.5) +
    theme_bw()
  dfTab <- as.data.frame(table(cs_dat[, 1:2]))
  dfTab[, 1:2] <- apply(dfTab[, 1:2], 2, function(x) (as.numeric(x)))
  dfTab$lab <- as.character(dfTab$Freq)
  dfTab <- dfTab[dfTab$Freq > 0, ]
  plot_list2[[pair]] <- ggplot(data = dfTab, aes(x = A, y = B, label = lab, colour = Freq)) +
    geom_point(shape = 15, size = 8) +
    labs(title = paste0(pair, ", \nq_t: ", round(result$p_BY, 4), ", \t\tq_cor: ", round(df_ps[df_ps$pair == pair, "p_BY"], 4),
                        "\np_t: ", round(result$p, 4), ", \t\tp_cor: ", round(df_ps[df_ps$pair == pair, "p"], 4))) +
    scale_colour_gradient(low = "white", high = "#d53397") +
    geom_text(size = 4, colour = "black") +
    theme_bw()
}

pdf(file = "results/figures/corr/significant_our_method.pdf")
for (i in 1:length(plot_list)) {
  print(plot_list[[i]])
}
dev.off()

pdf(file = "results/figures/corr/significant_our_method_counts.pdf", height = 6, width = 7)
for (i in 1:length(plot_list)) {
  print(plot_list2[[i]])
}
dev.off()

PlotStackDistribution(MIC_clean, results, 0.5, one_direction = FALSE, CResponse = "CS", whichAB = c("CFZ", "AMC"))



cor.test(cs_dat$A[cs_dat$A < 6 ], cs_dat$B[cs_dat$A < 6])

dfTab <- as.data.frame(table(cs_dat[, 1:2]))
colnames(dfTab)[1] <- "x"
dfTab$lab <- as.character(dfTab$Freq)
dfTab <- dfTab[dfTab$Freq > 0, ]

ggplot(cs_dat, aes(x = A, y = B)) +
  geom_point(shape = 0, size = 13, show.legend = F) +
  #geom_bin2d(binwidth = 1, color = "white") +
  #stat_bin2d(binwidth = 1, show.legend = F, geom = "text", aes(label = ..count..))+
  #scale_fill_gradient(low = "white", high = "#d53397") +
  #scale_size_continuous(range = c(0.1, 20)) +
  geom_smooth(method = lm, se = F, colour = "#283c82") +
  labs(title = paste0(pair, ", \nq_t: ", round(result$p_BY, 4), ", \t\tq_cor: ", round(df_ps[df_ps$pair == pair, "p_BY"], 4),
                      "\np_t: ", round(result$p, 4), ", \t\tp_cor: ", round(df_ps[df_ps$pair == pair, "p"], 4))) +
  # geom_hline(yintercept = result$d - 0.5) +
  geom_text(data = dfTab, aes(x = A, y = B, label = lab), size = 3) +
  theme_bw()


ggplot(data = dfTab, aes(x = A, y = B, label = lab, colour = Freq)) +
  geom_point(shape = 15, size = 6) +
  scale_colour_gradient(low = "white", high = "#d53397") +
  geom_text(size = 3, colour = "black") +
  theme_bw()



###TRY SIMULATION###
set.seed(123)
n_sim <- 250
states <- 0:3
A <- rnorm(n_sim, mean = 3, sd = 2)
A_obs <- ceiling(A)
hist(A_obs)

#BR_ind <- sample(1:length(A_obs), ceiling(n_sim/2)) #if no collateral response
BR_ind <- sample(1:length(A_obs), ceiling(n_sim/2), prob = seq(1, 0, length.out = n_sim)) # if collateral sensitivity


effect_size <- -1
B <- numeric(n_sim)
B[BR_ind] <- rnorm(length(BR_ind), mean = 3 + effect_size, sd = 2)
B[-BR_ind] <- rnorm(n_sim - length(BR_ind), mean = 3, sd = 2)
B_obs <- ceiling(B)

plot(A_obs, B_obs, pch = 1, cex = seq(0.1, 3, length.out = n_sim))
abline(lm(B_obs ~ A_obs))

X_w <- A_obs[B_obs < median(B_obs)]
Y <- A_obs[B_obs >= median(B_obs)]
T_value <- t.test(Y, X_w, var.equal = TRUE, alternative = direction)

print(cor.test(A_obs, B_obs, method = "kendall"))
print(T_value)

par(mfrow = c(2, 1))
hist(A_obs, breaks = 10, xlim = range(A_obs))
hist(A_obs[B_obs >= median(B_obs)], breaks = 10, xlim = range(A_obs))
par(mfrow = c(1, 1))


View(cbind(A_obs, B_obs))


no_col <- sample(states, 100, replace = T)






#Compare with correlations
load("results/test_resultsPseudomonas-aeruginosa.Rdata")

t_results <- results$`0.5`

#calculate correlations
load("data/clean/MIC_clean_Pseudomonas-aeruginosa.Rdata")
cor_matrix <- Hmisc::rcorr(log2(as.matrix(MIC_clean)))


t_results$cor <- t_results$p_cor <- 0
for (i in 1:nrow(t_results)) {
  t_results$cor[i] <- cor_matrix$r[t_results$A[i], t_results$B[i]]
  t_results$p_cor[i] <- cor_matrix$P[t_results$A[i], t_results$B[i]]
}


plot1 <- ggplot(t_results, aes(x = cor, y = t)) +
  geom_hline(yintercept = qnorm(c(0.975, 0.5, 0.025)), linetype = c(2,1,2), colour = "grey60") +
  geom_vline(xintercept = qnorm(c(0.975, 0.5, 0.025)), linetype = c(2,1,2), colour = "grey60") +
  geom_point(colour = "grey60", shape = 16) +
  # geom_point(data = df_ps[df_ps$sign_cor, ], shape = 16, colour = "white") +
  # geom_point(data = df_ps[df_ps$sign_t, ], colour = "#283c82", shape = 16) +
  # geom_point(data = df_ps[df_ps$sign_cor, ], colour = "#F46B2D", shape = 1) +
  #scale_colour_manual(values = c("#283c82", "#F46B2D")) +
  
  #  geom_smooth(method = "lm", colour = "dark red", linetype = 2) +
  labs(title = "T-statistics from correlation and dependent t-test", x = "t(correlation)", y = "t(t-test)") +
  theme_bw()
ggplot(t_results, aes(x = cor, y = t)) +
  geom_point()

plot_ <- ggplot(dat, aes(x = A, y = B, size = abs(t), color = as.factor(sign(t))))+ 
  geom_point(shape = 15) +
  scale_size(range = c(3, 10)) +
  scale_color_manual(values = c("#F46B2D", "#283c82"), labels = c("Collateral Sensitivity", "Collateral Resitance")) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Anitbiotic A", y = "Antibiotic B", size = "T-value", colour = "Direction",
       title = paste("Dichotomised on quantile", dicho_value)) +
  geom_vline(xintercept = seq(1.5, length(unique(dat$A)) - 0.5, 1), colour = "grey60") +
  geom_hline(yintercept = seq(1.5, length(unique(dat$B)) - 0.5, 1), colour = "grey60") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4)), 
         size =  guide_legend(title.theme = element_text(size = 16)))

plot_ <- plot_ + geom_point(data = dat[dat$p_BY < FDR_crit, ], pch = "*", aes(x = A, y = B), colour = "white",
                            size = 7)






