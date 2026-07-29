plot_kobe <- function(stk, refpts_obj, use_ssb = TRUE, title = "Kobe Stock Status Phase Plot") {  
  # 1. Extract Reference Points
  b_msy <- as.numeric(refpts_obj['msy', if (use_ssb) 'ssb' else 'biomass'])[1]
  f_msy <- as.numeric(refpts_obj['msy', 'harvest'])[1]
  
  b_flq <- if (use_ssb) ssb(stk) else stock(stk)
  years <- as.numeric(dimnames(stk)$year)
  terminal_year <- max(years)
  
  # 2. Extract Medians for the Historical Trajectory
  b_vec <- as.numeric(iterMedians(b_flq))
  f_vec <- as.numeric(iterMedians(fbar(stk)))
  
  kobe_df <- data.frame(
    year = years,
    b_r  = b_vec / b_msy,
    f_r  = f_vec / f_msy
  )
  
  # 3. Extract Iterations for the Terminal Year (The Cloud)
  b_term_iters <- as.numeric(b_flq[, ac(terminal_year)]) / b_msy
  f_term_iters <- as.numeric(fbar(stk)[, ac(terminal_year)]) / f_msy
  
  iters_df <- data.frame(
    b_r = b_term_iters,
    f_r = f_term_iters
  )
  
  # 4. Calculate Probabilities per Quadrant
  p_green <- mean(iters_df$b_r >= 1 & iters_df$f_r <= 1) * 100 # Not overfished, no overfishing
  p_red   <- mean(iters_df$b_r < 1  & iters_df$f_r > 1) * 100  # Overfished, overfishing
  p_obl   <- mean(iters_df$b_r < 1  & iters_df$f_r <= 1) * 100 # Overfished, no overfishing
  p_otr   <- mean(iters_df$b_r >= 1 & iters_df$f_r > 1) * 100  # Not overfished, overfishing
  
  # 5. Dynamically Scale Plot Limits to Encompass All Iterations
  max_x <- max(c(2.2, kobe_df$b_r * 1.1, iters_df$b_r * 1.05), na.rm = TRUE)
  max_y <- max(c(2.2, kobe_df$f_r * 1.1, iters_df$f_r * 1.05), na.rm = TRUE)
  x_label <- if (use_ssb) expression(SSB / SSB[MSY]) else expression(B / B[MSY])
  
  df_first <- kobe_df[1, ]
  df_last  <- kobe_df[nrow(kobe_df), ]
  
  # 6. Build the Plot
  ggplot(kobe_df, aes(x = b_r, y = f_r)) +
    ## Quadrants
    annotate("rect", xmin = 1, xmax = max_x, ymin = 0, ymax = 1, fill = "#a6d96a", alpha = 0.8) + 
    annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = max_y, fill = "#d7191c", alpha = 0.8) + 
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "#fdae61", alpha = 0.8) +    
    annotate("rect", xmin = 1, xmax = max_x, ymin = 1, ymax = max_y, fill = "#fdae61", alpha = 0.8) +
    
    ## Reference lines
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray30") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
    
    ## Cloud of Points: Terminal year iterations
    geom_point(data = iters_df, aes(x = b_r, y = f_r), color = "gray30", alpha = 0.35, size = 1.2) +

    ## Trajectory path & neutral intermediate points
    geom_path(color = "gray30", linewidth = 0.6, arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_point(color = "gray30", size = 1.8) +
    
    ## First year (Blue square + year label)
    geom_point(data = df_first, color = "#2b83ba", size = 3.5, shape = 15) +
    geom_text(data = df_first, aes(label = year), vjust = -1.2, fontface = "bold", color = "#2b83ba", size = 3.2) +
    
    ## Last year (Yellow circle with black outline + year label)
    geom_point(data = df_last, fill = "yellow", color = "black", size = 4, shape = 21, stroke = 1.2) +
    geom_text(data = df_last, aes(label = year), vjust = -1.3, fontface = "bold", color = "black", size = 3.2) +
    
    ## Probability Labels (Anchored to corners)
    annotate("text", x = 0.05, y = max_y - 0.1, label = sprintf("%.1f%%", p_red), fontface = "bold", hjust = 0, color = "darkred") +
    annotate("text", x = max_x - 0.05, y = 0.1, label = sprintf("%.1f%%", p_green), fontface = "bold", hjust = 1, color = "darkgreen") +
    annotate("text", x = 0.05, y = 0.1, label = sprintf("%.1f%%", p_obl), fontface = "bold", hjust = 0, color = "darkorange") +
    annotate("text", x = max_x - 0.05, y = max_y - 0.1, label = sprintf("%.1f%%", p_otr), fontface = "bold", hjust = 1, color = "darkorange") +

    coord_cartesian(xlim = c(0, max_x), ylim = c(0, max_y), expand = FALSE) +
    labs(x = x_label, y = expression(F / F[MSY]), title = title) +
    theme_bw(base_size = 10)
}