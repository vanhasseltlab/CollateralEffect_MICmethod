###Functions for collateral sensitivity project

##t_test_CR.R
PlotAllTvalues <- function(results, dicho_value, FDR_crit = 0.15) {
  dat <- results[[as.character(dicho_value)]]
  
  plot_ <- ggplot(dat, aes(x = A, y = B, size = abs(t), color = as.factor(sign(t))))+ 
    geom_point(shape = 15) +
    scale_size(range = c(3, 10)) +
    scale_color_manual(values = c("#d53397", "#283c82"), labels = c("Collateral Sensitivity", "Collateral Resitance")) +
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
  
  return(plot_)
}

PlotSignificantTvalue <- function(results, dicho_value, FDR_crit = 0.15) {
  dat <- results[[as.character(dicho_value)]] %>% 
    filter(p_BY < FDR_crit)
  
  
  plot_ <- ggplot(dat, aes(x = A, y = B, size = abs(t), color = as.factor(sign(t))))+ 
    geom_point(shape = 15) +
    scale_size(range = c(3, 10)) +
    scale_color_manual(values = c("#d53397", "#283c82"), labels = c("Collateral Sensitivity", "Collateral Resitance")) +
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
  
  plot_ <- plot_ + geom_abline(slope = 1, intercept = 0, colour = "grey60")
  
  return(plot_)
}

PlotCRDistributions <- function(MIC_clean, results, dicho_value, t_rank = 1, one_direction = TRUE, CResponse = "CS", FDR_crit = 0.15) {
  sign_dat <- results[[as.character(dicho_value)]] %>% 
    filter(p_BY < FDR_crit) %>% 
    arrange(if(CResponse == "CS") {.$t} else {desc(.$t)})
  
  if (nrow(sign_dat) < 1) {
    message("No significant collateral sensitivity, plotting non significant finding with largest T")
    sign_dat <- results[[as.character(dicho_value)]] %>% 
      arrange(t) %>% 
      slice(1)
    t_rank <- 1
  }
  if (nrow(sign_dat) < t_rank) {
    message(paste0("No", t_rank, "th significant collateral sensitivity, plotting significant finding with largest T"))
    t_rank <- 1
  }
  
  
  #for all combinations?
  comb <- unlist(sign_dat[t_rank, 1:2])
  d <- sign_dat[t_rank, "d"]
  dat <- log2(MIC_clean[, comb])
  dat <- dat[!is.na(dat[, 1]), ]
  names(dat) <- c("A", "B")
  
  ran <- range(dat$A) + (0.5 * c(-1, + 1))
  ticks <- floor(seq(ran[1], ran[2] + 2))
  p_A <- dat %>% 
    ggplot(aes(x = A, y = ..prop..)) +
    geom_bar(stat = "count", width = 0.5, position = "dodge") +
    labs(x = "log2(MIC)", y = "Probability mass", title = comb[1]) +
    scale_x_continuous(breaks = ticks, limits = ran) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05)), labels = function(x) sprintf("%.2f", x)) +
    geom_vline(xintercept = mean(dat$A, na.rm = T), colour = "#d53397") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank())
  
  p_AB <- dat %>% 
    filter(B >= d & !is.na(B)) %>% 
    ggplot((aes(x = A, y = ..prop..))) +
    geom_bar(stat = "count", width = 0.5, position = "dodge") +
    labs(x = "log2(MIC)", y = "Probability mass", 
         title = paste0(comb[1], "|", comb[2], " > ", d)) +
    scale_x_continuous(breaks = ticks, limits = ran) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05)), labels = function(x) sprintf("%.2f", x)) +
    geom_vline(xintercept = mean(dat$A[dat$B >= d & !is.na(dat$B)], na.rm = T), colour = "#d53397") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank())
  
  if (one_direction) {
    return(grid.arrange(p_A, p_AB, ncol = 1))
  }
  #Show directionality
  new_dat <- results[[as.character(dicho_value)]]
  rownames(new_dat) <- paste(new_dat$A, new_dat$B, sep = "_")
  comb <- unlist(sign_dat[t_rank, 2:1])
  d <- new_dat[paste(comb, collapse = "_"), "d"]
  dat <- log2(MIC_clean[, comb])
  dat <- dat[!is.na(dat[, 1]), ]
  names(dat) <- c("A", "B")
  
  ran <- range(dat$A) + (0.5 * c(-1, + 1))
  ticks <- floor(seq(ran[1], ran[2] + 2))
  p_B <- dat %>% 
    ggplot(aes(x = A, y = ..prop..)) +
    geom_bar(stat = "count", width = 0.5, position = "dodge") +
    labs(x = "log2(MIC)", y = "Probability mass", title = comb[1]) +
    scale_x_continuous(breaks = ticks, limits = ran) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05)), labels = function(x) sprintf("%.2f", x)) +
    geom_vline(xintercept = mean(dat$A, na.rm = T), colour = "#d53397") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank())
  
  p_BA <- dat %>% 
    filter(B >= d & !is.na(B)) %>% 
    ggplot((aes(x = A, y = ..prop..))) +
    geom_bar(stat = "count", width = 0.5, position = "dodge") +
    labs(x = "log2(MIC)", y = "Probability mass", 
         title = paste0(comb[1], "|", comb[2], " > ", d)) +
    scale_x_continuous(breaks = ticks, limits = ran) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05)), labels = function(x) sprintf("%.2f", x)) +
    geom_vline(xintercept = mean(dat$A[dat$B >= d & !is.na(dat$B)], na.rm = T), colour = "#d53397") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank())
  
  return(grid.arrange(p_A, p_B, p_AB, p_BA, ncol = 2, nrow = 2))
  
  
}
