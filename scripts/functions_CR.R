###Functions for collateral sensitivity project

#Libraries needed
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)

##our test statistic (based on t.test())
CRTTest <- function(A, B, quant = 0.5) {
  
  criterium <- quantile(B, quant, na.rm = TRUE)
  
  Z <- A[B < criterium | is.na(B)]
  Y <- A[B >= criterium & !is.na(B)]
  w <- length(Y)/length(A)
  #TODO: include some if's for possible errors
  #if (any(length(Z) < 2, length(Y) < 2)) warning("Not enough observations for test!")
  
  dep_t_test <- t.test(Y, Z, var.equal = TRUE, alternative = direction)
  
  #adjust estimates
  dep_t_test$estimate <- c(mean(Y), mean(A))
  names(dep_t_test$estimate) <- c("mean of A|B = R", "mean of A")
  
  #adjust standard error
  dep_t_test$stderr <- dep_t_test$stderr*(1-w)
  dep_t_test$conf.int <- dep_t_test$conf.int*(1-w)
  #diff(dep_t_test$estimate[2:1]) + qt(c(0.025, 0.975), dep_t_test$parameter)*dep_t_test$stderr
  
  #adjust method
  dep_t_test$method <- "Overlapping Sample t-test"
  dep_t_test$data.name <- "A|B = R and A"
  
  return(dep_t_test)
}

##t_test_CR.R
PlotAllTvalues <- function(results, dicho_value, FDR_crit = 0.15) {
  dat <- results[[as.character(dicho_value)]]
  
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
  
  return(plot_)
}

PlotSignificantEffect <- function(results, dicho_value, FDR_crit = 0.15, species, selectedAB = NULL) {
  dat <- results[[as.character(dicho_value)]]

  if (!is.null(selectedAB)) {
    dat$effect_size[dat$p_BY > FDR_crit] <- 0
    dat <- dat %>% 
      filter(A %in% selectedAB & B %in% selectedAB)
  } else {
    dat <- dat[dat$p_BY < FDR_crit, ]
  }

  bl <- colorRampPalette(c("#283c82", "white"))(30) [c(1:10, seq(11, 30, by = 2))] 
  #colorRampPalette(c("red","#d53397"))(30)
  re <- colorRampPalette(c("#f54c00", "white"))(30)[c(1:10, seq(11, 30, by = 2))]
  
  limits <- c(-1, 1)*max(dat$effect_size)
  
  bb <- limits
  ll <- c("Collateral Sensitivity", "Collateral Resistance") # labels.
  
  dat$Direction <- ifelse(dat$effect_size > 0, ll[1], ll[2] )
  
  plot_ <- ggplot(dat, aes(x = A, y = B, color = effect_size, shape = Direction)) +
    #geom_point(shape = 15, size = 13) +
    geom_point(shape = 15, size = 10) +
    scale_colour_gradientn(colours = c(re, "white", rev(bl)), limits = limits) +
    scale_x_discrete(expand = expand_scale(mult = 0, add = rep(0.5, 2))) +
    scale_y_discrete(expand = expand_scale(mult = 0, add = rep(0.5, 2))) +
    geom_point(size = 4, colour = "white") +
    scale_shape_manual(values = c("-", "+"), labels = ll) +
    coord_fixed() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "Antibiotic A", y = paste0("Antibiotic B (dichotomized on quantile ", dicho_value, ")"), 
         size = "Difference in means", colour = "Difference between means", shape = "Type",
         title = paste("Significant collateral responses", species)) +
    geom_vline(xintercept = seq(1.5, length(unique(dat$A)) - 0.5, 1), colour = "grey60") +
    geom_hline(yintercept = seq(1.5, length(unique(dat$B)) - 0.5, 1), colour = "grey60") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.title = element_text(size = 12), legend.key = element_rect(fill = "grey60"))
  
  plot_ <- plot_ + geom_abline(slope = 1, intercept = 0, colour = "grey60")
  plot_
  
  return(plot_)
}


PlotCRDistributions <- function(MIC_clean, results, dicho_value, t_rank = 1, one_direction = TRUE, CResponse = "CS", FDR_crit = 0.15) {
  sign_dat <- results[[as.character(dicho_value)]] %>% 
    filter(p_BY < FDR_crit & (sign(t) == c(-1, 1)[CResponse == c("CS", "CR")])) %>% 
    arrange(desc(abs(effect_size)))
  
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
    ggplot(aes(x = A, y = ..count..)) +
    geom_bar(stat = "count", width = 0.5, position = "dodge") +
    labs(x = "log2(MIC)", y = "Probability mass", title = comb[1]) +
    scale_x_continuous(breaks = ticks, limits = ran) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05)), labels = function(x) sprintf("%.2f", x)) +
    geom_vline(xintercept = mean(dat$A), colour = "#d53397") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank())
  
  p_AB <- dat %>% 
    filter(B >= d & !is.na(B)) %>% 
    ggplot((aes(x = A, y = ..count..))) +
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
  
  #remove all labels from the plots
  p_BA <- p_BA + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  p_B <- p_B + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  p_AB <- p_AB + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  p_A <- p_A + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

  return(grid.arrange(p_A, p_B, p_AB, p_BA, ncol = 2, nrow = 2, bottom = textGrob(label = expression(log[2]*"(MIC)")), 
                      left = "Probability mass"))
  
}


PlotStackDistribution <- function(MIC_clean, results, dicho_value, t_rank = 1, one_direction = TRUE, CResponse = "CS", FDR_crit = 0.15, whichAB = NULL) {

  if(!is.null(whichAB)) {
    comb <- whichAB
    sign_dat <- results[[as.character(dicho_value)]] %>% 
      filter(A == comb[1] & B == comb[2])
    row.names(sign_dat) <- paste0(sign_dat$A, sign_dat$B)
    if (whichAB[1] %in% sign_dat$A & whichAB[2] %in% sign_dat$B) {
      d <- sign_dat[paste0(comb, collapse = ""), "d"]
    } else {
      message("Specified antibiotics are not in data, using highest effect size")
      whichAB <- NULL
    }
  }
  
  
  if(is.null(whichAB)) {
    sign_dat <- results[[as.character(dicho_value)]] %>% 
      filter(p_BY < FDR_crit & (sign(t) == c(-1, 1)[CResponse == c("CS", "CR")])) %>% 
      arrange(desc(abs(effect_size)))
    row.names(sign_dat) <- paste0(sign_dat$A, sign_dat$B)
    
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
    comb <- unlist(sign_dat[t_rank, 1:2])
    d <- sign_dat[t_rank, "d"]
  }
  

  
  #for all combinations?
  dat <- log2(MIC_clean[, comb])
  dat <- dat[!is.na(dat[, 1]), ]
  names(dat) <- c("A", "B")
  d_min <- max(dat$B[dat$B < d], na.rm = T)
  
  dat$Condition <- as.factor(ifelse(dat$B >= d & !is.na(dat$B), 
                                    paste0(comb[1],"|", comb[2]," > ", d_min), 
                                    paste0(comb[1])))
  means <- data.frame(mean = c(mean(dat$A), mean(dat$A[dat$Condition == paste0(comb[1],"|", comb[2]," > ", d_min)])), 
                      Means = sort((unique(dat$Condition))))
  ran <- range(dat$A) + (0.5 * c(-1, + 1))
  ticks <- floor(seq(ran[1], ran[2] + 2))
  
  plotPanel <- function(dat) {
    
    plotje <- dat %>% 
      ggplot(aes(x = A, y = ..count.., group = Condition, fill = Condition)) +
      geom_bar(stat = "count", width = 0.6, position = "stack") +
      labs(x = expression(log[2]*"(MIC)"), y = "Counts", title = comb[1]) +
      #scale_y_continuous(expand = expand_scale(mult = c(0, .05)), labels = function(x) sprintf("%.2f", x)) +
      scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
      scale_x_continuous(breaks = ticks, limits = ran) +
      geom_vline(data = means, aes(xintercept = mean), colour = "white") +
      geom_vline(data = means, aes(xintercept = mean, colour = Means), show.legend  = TRUE, linetype = 2) +
      scale_fill_manual(values = c("#283c82", "#F46B2D")) +
      scale_colour_manual(labels = c(bquote(hat(mu)[.(as.character(means$Means[1]))]),
      #                               bquote(paste(hat(mu)[paste(comb[1], "|", comb[2],>=, d)]))),
                                     bquote(hat(mu)[.(as.character(means$Means[2]))])), 
                          values = c("#283c82", "#F46B2D")) +
      theme_bw() +
      theme(panel.grid.minor.x = element_blank(), legend.position = "bottom",
       #     legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "black"),
            legend.box = "horizontal", legend.direction = "vertical") + 
      guides(fill = guide_legend(override.aes = list(linetype = 0), order = 1, title = NULL), 
             colour = guide_legend(order = 2, title = NULL, label.theme = element_text(size = 12)))
      
  
    return(plotje)
  }
  

  p_A <- plotPanel(dat)
  
  if (one_direction) {
    return(p_A)
  }
  
  new_dat <- results[[as.character(dicho_value)]]
  rownames(new_dat) <- paste(new_dat$A, new_dat$B, sep = "_")
  comb <- comb[2:1]
  d <- new_dat[paste(comb, collapse = "_"), "d"]
  dat <- log2(MIC_clean[, comb])
  dat <- dat[!is.na(dat[, 1]), ]
  names(dat) <- c("A", "B")
  d_min <- max(dat$B[dat$B < d], na.rm = T)
  
  dat$Condition <- as.factor(ifelse(dat$B >= d & !is.na(dat$B), 
                                    paste0(comb[1],"|", comb[2]," > ", d_min), 
                                    paste0(comb[1])))
  means <- data.frame(mean = c(mean(dat$A), mean(dat$A[dat$Condition == paste0(comb[1],"|", comb[2]," > ", d_min)])), 
                      Means = sort((unique(dat$Condition))))
  ran <- range(dat$A) + (0.5 * c(-1, + 1))
  ticks <- floor(seq(ran[1], ran[2] + 2))
  p_B <- plotPanel(dat) 
  
  p_A <- p_A + theme(axis.title.y = element_blank())
  p_B <- p_B + theme(axis.title.y = element_blank())
  
  return(grid.arrange(p_A, p_B, ncol = 2, left = textGrob("Counts", hjust = -0.45, rot = 90)))
}
