library(tidyverse)
load("test_new_measure_data/MIC_clean.Rdata")

#all_combs <- expand.grid(colnames(MIC_clean), colnames(MIC_clean), stringsAsFactors = F)
#all_combs <- as.matrix(all_combs[all_combs[, 1] != all_combs[, 2], ])

criteria_quant <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)


results <- vector(mode = "list", length = length(criteria_quant))
names(results) <- criteria_quant

antibiotics <- colnames(MIC_clean)
n <- length(antibiotics)
t_test <- as.data.frame(matrix(0, nrow = n*(n - 1), ncol = 5))
names(t_test) <- c("A", "B", "t", "p", "n")


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
      X_w <- dat[dat[, 2] < criterium, 1]
      Y <- dat[dat[, 2] >= criterium , 1]
      t_test_i <- t.test(Y, X_w, var.equal = TRUE, alternative = "less")
      t_test[counter, 3:4] <- c(t_test_i$statistic, t_test_i$p.value)
    }
    
    
    # for (i in 1:n_tests) {
    #   dat <- log2(MIC_clean[, all_combs[i, ]])
    #   dat <- dat[is.na(dat[, 1]) == 0, ]
    #   criteria <- quantile(dat[, 2], c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9), na.rm = TRUE)
    #   n_unique_tests[i] <- length(unique(criteria))
    #   n[i] <- nrow(dat)
    #   X_w <- dat[dat[, 2] < criteria[crit], 1]
    #   Y <- dat[dat[, 2] >= criteria[crit] , 1]
    #   t_test[i, crit] <- t.test(Y, X_w, var.equal = TRUE)$statistic
    # }
  }
  result_crit <- data.frame(t_test, p_BY = p.adjust(t_test$p, method = "BY"))
  
  results[[crit]] <- result_crit
}


quantile(log2(MIC_clean$CFZ), c(0.4, 0.5, 0.6, 0.8, 0.9), na.rm = T)
quantile(log2(MIC_clean$TZP), c(0.4, 0.5, 0.6, 0.8, 0.9), na.rm = T)

plot(log2(MIC_clean$CFZ), log2(MIC_clean$TZP), cex = seq(0.1, 3, length.out = length(MIC_clean$CFZ)))

hist(log2(MIC_clean$CFZ))
hist(log2(MIC_clean$CFZ[log2(MIC_clean$TZP) > 6]), xlim = c(0, 6))

results_t <- data.frame(all_combs, t_test, n_unique_tests, n)

#One sided t-test for collateral sensitivity 
p_value_test <- pt(as.matrix(results_t[, paste0("X", criteria_quant)]), df = n - 2, lower.tail = T)

#Multiple testing correction Benjamini-Yekutieli
p_value_adju <- lapply(results, function(x) p.adjust(x$p, method = "BY"))


#Two sided t-test for plot 
p_value_test_plot <- pt(abs(as.matrix(results_t[, paste0("X", criteria_quant)])), df = n - 2, lower.tail = F)

#Multiple testing correction Benjamini-Yekutieli
p_value_adju_plot <- apply(p_value_test_plot, 2, function(x) p.adjust(x, method = "BY"))

#Summarize in one data frame
results_t_test <- data.frame(results_t, p_value_test, p_value_adju)

#Make clean long dataset with only significant observations



#Plot results
#Overview plot
PlotFunctionTvalue <- function(results_t, dicho_value, p_values = NULL) {
  find_column <- grep(as.character(dicho_value), names(results_t))
  if (length(find_column) == 0) {
    warning("Dichotomisation value has not been found in colnames of given data. First column after antibiotic names is chosen!")
    find_column <- 3
    dicho_value <- names(results_t)[3]
  }
  names(results_t)[c(1, 2, find_column)] <- c("A", "B", "X")
  
  
  plot_ <- ggplot(results_t, aes(x = A, y = B, size = abs(X), color = as.factor(sign(X))))+ 
    geom_point(shape = 15) +
    scale_size(range = c(3, 10)) +
    scale_color_manual(values = c("#d53397", "#283c82"), labels = c("Collateral Sensitivity", "Collateral Resitance")) +
    coord_fixed() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "Anitbiotic A", y = "Antibiotic B", size = "T-value", colour = "Direction",
         title = paste("Dichotomised on quantile", dicho_value)) +
    geom_vline(xintercept = seq(1.5, length(unique(results_t$A)) - 0.5, 1), colour = "grey60") +
    geom_hline(yintercept = seq(1.5, length(unique(results_t$B)) - 0.5, 1), colour = "grey60") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.title = element_text(size = 12)) +
    guides(color = guide_legend(override.aes = list(size = 4)), 
           size =  guide_legend(title.theme = element_text(size = 16)))
  
  if (!is.null(p_values)) {
    significant_dat <- results_t[p_values[, grep(as.character(dicho_value), colnames(p_values))] <= 0.15, ]
    plot_ <- plot_ + geom_point(data = significant_dat, pch = "*", aes(x = A, y = B), colour = "black",
                                size = 5)
  }
  
  return(plot_)
}


pdf(file = paste0("output/t_teststatistic_allcomparisons.pdf"), height = 7, width = 9)
print(PlotFunctionTvalue(results_t, 0.4, p_values = p_value_adju))
print(PlotFunctionTvalue(results_t, 0.5, p_values = p_value_adju))
print(PlotFunctionTvalue(results_t, 0.6, p_values = p_value_adju))
print(PlotFunctionTvalue(results_t, 0.7, p_values = p_value_adju))
print(PlotFunctionTvalue(results_t, 0.8, p_values = p_value_adju))
print(PlotFunctionTvalue(results_t, 0.9, p_values = p_value_adju))
dev.off()


#Significant findings plot
library(gridExtra)
p_value_adju[1:10,]
for (){}


p_A <- MIC_clean %>% 
  filter(!is.na(CFZ)) %>% 
  ggplot((aes(x = log2(CFZ)))) +
  geom_bar(stat = "count", width = 0.5) +
  ggtitle("counts CFZ") +
  scale_y_continuous(limits = c(0, 1600)) +
  theme_bw()

p_AB <- MIC_clean %>% 
  filter(!is.na(CFZ) & AMP > 16) %>% 
  ggplot((aes(x = log2(CFZ)))) +
  geom_bar(stat = "count", width = 0.5) +
  ggtitle("counts CFZ|AMP > 16") +
  scale_y_continuous(limits = c(0, 1600)) +
  theme_bw()

pdf(file = "output/compare_distributions_given_B.pdf", width = 6, height = 6)
grid.arrange(p_A, p_AB, ncol = 1)
dev.off()

