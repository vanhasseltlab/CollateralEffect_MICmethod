###Functions for collateral sensitivity project

#Libraries required
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)

#data preparation
RemoveDuplicateMICs <- function(MIC_df, key = NULL) {
  
  if (!is.null(key)) {
    MIC_df$key <- apply(MIC_df[, key], 1, paste, collapse = "_")
  }
  
  # two unique measurements for the same key wil be summarized: max if ">" and mean if "<=" 
  if (length(unique(MIC_df$key)) == nrow(MIC_df)) {
    return(MIC_df)
  }
  find_unique <- ave(MIC_df$key, MIC_df$key, FUN = length) == 1
  new_MIC <- MIC_df[find_unique, ]
  search_keys <- unique(MIC_df$key[!find_unique])
  for (i in search_keys) {
    cat("\r", "At (%): ", round(which(search_keys == i)/length(search_keys)*100, 2), ", key = ", i)
    dat <- MIC_df[MIC_df$key == i, ]
    if (nrow(dat) == 1) {
      new_MIC <- rbind(new_MIC, dat)
    }
    if (nrow(dat) > 1) {
      if (dat$measurement_sign[1] == dat$measurement_sign[2]) {
        if(dat$measurement_sign[1] == "<=") {
          MIC_i <- mean(dat$MIC)
        }
        if(dat$measurement_sign[1] == ">") {
          MIC_i <- max(dat$MIC)
        }
      } else {
        MIC_i <- mean(dat$MIC[dat$measurement_sign == "<="])
        dat$measurement_sign <- "<="
      }
      dat <- dat[1, ]
      dat$MIC <- MIC_i
      new_MIC <- rbind(new_MIC, dat)
    }
  }
  return(new_MIC)
}

##CE test statistic (based on t.test())
CETTest <- function(A, B, CE_type = "both", crit_type = "median", criterium = NULL) {
  # Performs T-test on vector A, which is split by a dichotomization of vector B  
  # 
  # Args:
  #   A: MIC values antibiotic A, numeric vector
  #   B: MIC values antibiotic B, numeric vecter, same length as A
  #   CE_type: type of collateral effect to be evaluated, string: "both", "CR" or "CS"
  #   crit_type: type of dichotomization criterion, string: "quant", "median" or "log2(MIC)"
  #   criterium: value of dichotomization criterium, single numeric, overwritten if crit_type = "median"
  # 
  # Returns: 
  #   List with class "htest" including statistic, parameter, p.value, etc.

  #Error handling
  if (length(A) != length(B)) {
    stop("Arguments A and B have different lengths: ", length(A), " and ", length(B), ".")
  }
  
  #Remove NA's
  ind_na <- is.na(A) | is.na(B)
  A <- log2(A[!ind_na])
  B <- log2(B[!ind_na])
  
  #Translate CE_type to t-test direction
  direction <- c("two.sided", "less", "greater")[c("both", "CS", "CR") == CE_type]
  
  #Calculate dichotomization criterium tau based on criterium type
  if (crit_type == "quant") {
    tau <- quantile(B, criterium)
    
  } else if (crit_type == "median") {

    tau <- quantile(B, 0.5)
    B_values <- sort(unique(B))
    #Adjust tau to create most equal split
    if (sum(B > tau) < sum(B < tau)) {
      tau <- mean(c(B_values[which(B_values == tau) - 1], tau))
    } else if (sum(B > tau) > sum(B < tau)) {
      tau <- mean(c(B_values[which(B_values == tau) + 1], tau))
    }
    
  } else {
    tau <- criterium
    if (!crit_type %in% c("quant", "median", "log2(MIC)")) {
      warning("crit_type is not specified (correctly), defaults to \"log2(MIC)\"")
    }
  }
  
  A_Bsens <- A[B < tau]
  A_Bres <- A[B >= tau]
  
  #Error handling
  if (length(A_Bsens) < 2 | length(A_Bres) < 2) {
    warning("Not enough observations for test! Returning list with limited information")
    return(list(A = A, B = B, tau = tau, data = list(`A|B = r` = A_Bres, `A|B != r` = A_Bsens)))
  }
  
  #Perform t.test
  A_t_test <- t.test(A_Bres, A_Bsens, var.equal = TRUE, alternative = direction)
  
  #Adjust names of groups
  names(A_t_test$estimate) <- c("mean of A|B = r", "mean of A|B != r")
  A_t_test$data.name <- "A|B = r and A|B != r"
  
  #Add data to output
  A_t_test$data <- list(`A|B = r` = A_Bres, `A|B != r` = A_Bsens)
  
  #Add tau to object
  A_t_test$tau <- tau
  
  return(A_t_test)
}

#plotting the results
PlotSignificantEffect <- function(t_result, FDR_crit = 0.15, species, selectedAB = NULL, t_or_effect = "effect") {
  
  t_result[is.na(t_result$t), c("effect_size", "t")] <- 0
  
  if (t_or_effect == "t") {
    t_result$effect_size <- t_result$t
    effect_label <- "T-value"
  } else { 
    effect_label <- expression(log[2]*"(FC)")
  }
  

  t_result$effect_size[t_result$q > FDR_crit] <- 0
  
  
  # find reciprocal AB combinations
  t_result$reciprocal <- "one-directional"
  combinations <- t(combn(unique(t_result$A), 2))
  
  for (i in 1:nrow(combinations)) {
    comb <- combinations[i, ]
    ind_comb <- t_result$A %in% comb & t_result$B %in% comb
    df_comb <- t_result[ind_comb, ]
    if (all(df_comb$effect_size > 0) | all(df_comb$effect_size < 0)){
      t_result[ind_comb, "reciprocal"] <- "reciprocal"
    }
  }
  levels(t_result$reciprocal) <- c("one-directional", "reciprocal")
  
  if (!is.null(selectedAB)) {
    t_result <- t_result %>% 
      filter(A %in% selectedAB & B %in% selectedAB)
  }
  bl <- colorRampPalette(c("#283c82", "white"))(30) [c(1:10, seq(11, 30, by = 3))] 
  re <- colorRampPalette(c("#f54c00", "white"))(30)[c(1:10, seq(11, 30, by = 3))]
  
  limits <- c(-1, 1)*max(t_result$effect_size)
  
  ll <- c("Collateral Sensitivity", "Collateral Resistance") # labels.
  
  t_result$Direction <- ifelse(t_result$effect_size < 0, ll[1], ll[2] )
  
  
  plot_ <- ggplot(t_result, aes(x = B, y = A)) +
    geom_tile(aes(fill = effect_size)) +
    scale_fill_gradientn(colours = c(bl, "white", rev(re)), limits = limits) +
    scale_x_discrete(expand = expansion(mult = 0, add = rep(0.5, 2))) +
    scale_y_discrete(expand = expansion(mult = 0, add = rep(0.5, 2))) +
    geom_point(aes(shape = reciprocal), size = 5, colour = "white") +
    scale_shape_manual(values = c("", "\u2194"), drop = FALSE) +
    coord_fixed() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "Splitting antibiotic (B)", y = paste0("Testing antibiotic (A)"), 
         fill = effect_label, shape = "", 
         title = bquote("Significant collateral responses"~italic(.(species)))) +
    geom_vline(xintercept = seq(1.5, length(unique(t_result$A)) - 0.5, 1), colour = "grey60") +
    geom_hline(yintercept = seq(1.5, length(unique(t_result$B)) - 0.5, 1), colour = "grey60") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.title = element_text(size = 12), legend.key = element_rect(fill = "grey60"),
          axis.title = element_text(size = 14)) + 
    geom_abline(slope = 1, intercept = 0, colour = "grey60")
  
  return(plot_)
}

PlotCondDistribution <- function(MIC_clean, results, t_rank = 1, one_direction = TRUE, CResponse = "CS", FDR_crit = 0.15, whichAB = NULL,
                                  colours = c("#283c82", "#F46B2D"), separate_plots = NULL, ran = NULL, t_or_effect = "t") {
  
  if (!is.null(whichAB)) {
    comb <- whichAB
    sign_dat <- results %>% 
      filter(A == comb[1] & B == comb[2])
    row.names(sign_dat) <- paste0(sign_dat$A, sign_dat$B)
    if (whichAB[1] %in% sign_dat$A & whichAB[2] %in% sign_dat$B) {
      d <- sign_dat[paste0(comb, collapse = ""), "tau"]
    } else {
      message("Specified antibiotics are not in data, using highest effect size")
      whichAB <- NULL
    }
  }
  
  
  if (is.null(whichAB)) {
    if (t_or_effect == "effect"){
      results$t <- results$effect_size
    }
    
    sign_dat <- results %>% 
      filter(q < FDR_crit & (sign(t) == c(-1, 1)[CResponse == c("CS", "CR")])) %>% 
      arrange(desc(abs(t)))
    row.names(sign_dat) <- paste0(sign_dat$A, sign_dat$B)
    
    if (nrow(sign_dat) < 1) {
      message("No significant ", CResponse, " effect, plotting non significant finding with largest T")
      if (CResponse == "CS") {
        sign_dat <- results %>% 
          arrange(t) %>% 
          slice(1)
      } else {
        sign_dat <- results %>% 
          arrange(desc(t)) %>% 
          slice(1)
      }
      t_rank <- 1
    }
    
    if (nrow(sign_dat) < t_rank) {
      message(paste0("No", t_rank, "th significant ", CResponse, " effect, plotting significant finding with largest T"))
      t_rank <- 1
    }
    comb <- unlist(sign_dat[t_rank, 1:2])
    d <- sign_dat[t_rank, "tau"]
  }
  
  
  
  #for all combinations?
  dat <- log2(MIC_clean[, comb])
  dat <- dat[!is.na(dat[, 1]) & !is.na(dat[, 2]), ]
  names(dat) <- c("A", "B")
  
  dat$Condition <- as.factor(ifelse(dat$B >= d, 
                                    paste0(comb[1],"|", comb[2]," > ", round(d, 2)), 
                                    paste0(comb[1],"|", comb[2]," < ", round(d, 2))))
  means <- data.frame(mean = c(mean(dat$A[dat$B < d]), mean(dat$A[dat$B >= d])), 
                      Means = sort((unique(dat$Condition))))
  change_range <- FALSE
  if (is.null(ran)){
    ran <- range(dat$A) + (0.5 * c(-1, + 1))
    change_range <- TRUE
  }
  ticks <- floor(seq(ran[1], ran[2] + 2))
  
  plotPanel <- function(dat) {
    
    plotje <- dat %>% 
      ggplot(aes(x = A, y = ..count.., group = Condition, fill = Condition)) +
      geom_bar(stat = "count", width = 0.6, position = "stack") +
      labs(x = bquote(log[2]*"(MIC) " ~ .(comb[1])), y = "Counts") +
      #labs(x = expression(log[2]*"(MIC) "*comb[1]), y = "Counts") +
      #scale_y_continuous(expand = expansion(mult = c(0, .05)), labels = function(x) sprintf("%.2f", x)) +
      scale_y_continuous(expand = expansion(mult = c(0, .05))) +
      scale_x_continuous(breaks = ticks, limits = ran) +
      geom_vline(data = means, aes(xintercept = mean), colour = "white") +
      geom_vline(data = means, aes(xintercept = mean, colour = Means), show.legend  = TRUE, linetype = 2) +
      scale_fill_manual(values = colours) +
      scale_colour_manual(labels = c(bquote(hat(mu)[.(as.character(means$Means[1]))]),
                                     #                               bquote(paste(hat(mu)[paste(comb[1], "|", comb[2],>=, d)]))),
                                     bquote(hat(mu)[.(as.character(means$Means[2]))])), 
                          values = colours) +
      theme_bw() +
      theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom",
            #     legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "black"),
            legend.box = "horizontal", legend.direction = "vertical",
            legend.margin = margin(0,0,0,0),
            legend.box.margin = margin(-15,0,0,0),
            legend.background = element_blank()) + 
      guides(fill = guide_legend(override.aes = list(linetype = 0), order = 1, title = NULL), 
             colour = guide_legend(order = 2, title = NULL, label.theme = element_text(size = 12)))
    
    
    return(plotje)
  }
  
  
  p_A <- plotPanel(dat)
  
  if (one_direction) {
    return(p_A)
  }
  
  new_dat <- results
  rownames(new_dat) <- paste(new_dat$A, new_dat$B, sep = "_")
  comb <- comb[2:1]
  d <- new_dat[paste(comb, collapse = "_"), "tau"]
  dat <- log2(MIC_clean[, comb])
  dat <- dat[!is.na(dat[, 1]) & !is.na(dat[, 2]), ]
  names(dat) <- c("A", "B")
  d_min <- max(dat$B[dat$B < d], na.rm = T)
  
  dat$Condition <- as.factor(ifelse(dat$B >= d, 
                                    paste0(comb[1],"|", comb[2]," > ", round(d, 2)), 
                                    paste0(comb[1],"|", comb[2]," < ", round(d, 2))))
  means <- data.frame(mean = c(mean(dat$A[dat$B < d]), mean(dat$A[dat$B >= d])), 
                      Means = sort((unique(dat$Condition))))
  
  if (change_range == TRUE){
    ran <- range(dat$A) + (0.5 * c(-1, + 1))
  }
  ticks <- floor(seq(ran[1], ran[2] + 2))
  p_B <- plotPanel(dat) 
  

  
  if (is.null(separate_plots)) {
    p_A <- p_A + theme(axis.title.y = element_blank())
    p_B <- p_B + theme(axis.title.y = element_blank())
    return(grid.arrange(p_A, p_B, ncol = 2, left = textGrob("Counts", hjust = -0.45, rot = 90)))
  } else {
    return(list(p_A, p_B + ylab("  ")))
  }
}
  

