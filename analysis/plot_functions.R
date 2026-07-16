plot_status <- function(stk, fit, a4a.stk){
  aa.stk <- stk + fit
  landings.n(aa.stk)<-catch.n(aa.stk)
  landings.wt(aa.stk)<-catch.wt(aa.stk)
  landings(aa.stk)<-catch(aa.stk)
  discards(aa.stk)<-0
  discards.n(aa.stk)<-0
  discards.wt(aa.stk)<-0

  model = "bevholt"
  srr_upd <- fmle(as.FLSR(aa.stk, model = model),control=list(trace=0),fixed=list(b=100))
  brp_upd<- brp(FLBRP(aa.stk, sr = srr_upd))
  ref_points_upd <- refpts(brp_upd)

  fmsy_upd <- as.numeric(ref_points_upd['msy','harvest'])
  ssbmsy_upd <- as.numeric(ref_points_upd['msy','ssb'])
  bbmsy_upd <- as.numeric(ref_points_upd['msy','biomass'])
  p1 <- plot(a4a.stk, metrics=list(SSB=ssb, F=fbar, B = stock)) +
    geom_flpar(data=FLPars(SSB=FLPar(SSBmsy = ssbmsy_upd),
    F=FLPar(FMSY=fmsy_upd),B=FLPar(BMSY=bbmsy_upd)), 
    x=c(1990, 2015,1997))
  p1
}

plot_catch_err <- function(flq_sim, title_text) {
      sim_data <- sum_sim(flq_sim)
      ggplot() +
        geom_ribbon(data = sim_data, aes(x = year, ymin = lwr, ymax = upr), 
                    fill = "#2b83ba", alpha = 0.25) +
        geom_line(data = sim_data, aes(x = year, y = med), 
                  color = "#2b83ba", linetype = "dashed", linewidth = 0.7) +
        geom_line(data = est_df, aes(x = year, y = data), 
                  color = "black", linewidth = 0.8) +
        geom_point(data = obs_df, aes(x = year, y = data), 
                   color = "#d7191c", size = 1.8) +
        labs(title = title_text, x = "Year", y = "Catch (t)") +
        theme_bw(base_size = 9)
  }

sum_sim <- function(flq) {
      df <- as.data.frame(flq)
      agg <- aggregate(data ~ year, data = df, FUN = function(x) c(
        med = median(x, na.rm = TRUE),
        lwr = unname(quantile(x, 0.05, na.rm = TRUE)),
        upr = unname(quantile(x, 0.95, na.rm = TRUE))
      ))
      mat <- as.data.frame(agg$data)
      names(mat) <- c("med", "lwr", "upr")
      data.frame(year = agg$year, mat)
  }

    plot_kobe <- function(stk, refpts_obj, use_ssb = TRUE, title = "Kobe Stock Status Phase Plot") {  
      b_msy <- as.numeric(refpts_obj['msy', if (use_ssb) 'ssb' else 'biomass'])[1]
      f_msy <- as.numeric(refpts_obj['msy', 'harvest'])[1]
      
      b_flq <- if (use_ssb) ssb(stk) else stock(stk)
      b_vec <- as.numeric(iterMedians(b_flq))
      f_vec <- as.numeric(iterMedians(fbar(stk)))
      years <- as.numeric(dimnames(stk)$year)
      
      kobe_df <- data.frame(
        year = years,
        b_r  = b_vec / b_msy,
        f_r  = f_vec / f_msy
      )
      
      max_x <- max(c(2.2, kobe_df$b_r * 1.1), na.rm = TRUE)
      max_y <- max(c(2.2, kobe_df$f_r * 1.1), na.rm = TRUE)
      x_label <- if (use_ssb) expression(SSB / SSB[MSY]) else expression(B / B[MSY])
      
      df_first <- kobe_df[1, ]
      df_last  <- kobe_df[nrow(kobe_df), ]
      
      ggplot(kobe_df, aes(x = b_r, y = f_r)) +
        ## Quadrants
        annotate("rect", xmin = 1, xmax = max_x, ymin = 0, ymax = 1, fill = "#a6d96a", alpha = 0.6) + 
        annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = max_y, fill = "#d7191c", alpha = 0.4) + 
        annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "#fdae61", alpha = 0.4) +    
        annotate("rect", xmin = 1, xmax = max_x, ymin = 1, ymax = max_y, fill = "#fdae61", alpha = 0.4) +
        ## Reference lines
        geom_vline(xintercept = 1, linetype = "dashed", color = "gray30") +
        geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
        ## Trajectory path & neutral intermediate points
        geom_path(color = "gray30", linewidth = 0.6, arrow = arrow(length = unit(0.2, "cm"), type =   
  "closed")) +
        geom_point(color = "gray30", size = 1.8) +
        ## First year (Blue square + year label)
        geom_point(data = df_first, color = "#2b83ba", size = 3.5, shape = 15) +
        geom_text(data = df_first, aes(label = year), vjust = -1.2, fontface = "bold", color =        
  "#2b83ba", size = 3.2) +
        ## Last year (Yellow circle with black outline + year label)
        geom_point(data = df_last, fill = "yellow", color = "black", size = 4, shape = 21, stroke = 1.2) +
        geom_text(data = df_last, aes(label = year), vjust = -1.3, fontface = "bold", color = "black",
  size = 3.2) +
        coord_cartesian(xlim = c(0, max_x), ylim = c(0, max_y), expand = FALSE) +
        labs(x = x_label, y = expression(F / F[MSY]), title = title) +
        theme_bw(base_size = 10)
    }
    
                                                                  
                                                                                                      
plotSels <- function(stk){
  sel <- harvest(stk)%/%apply(harvest(stk), c(2:6),max)
  ss <- as.data.frame(sel)
  ss$unique = 'Selectivity'
  p_sel <- ggplot(data = ss, aes(x = age, y = data, group = year,color = year)) +
    geom_line() + scale_colour_viridis_c(option = 'C') + facet_wrap(~unique) +
    theme(legend.position = 'bottom', axis.title.y=element_blank())+
    scale_x_continuous(expand = c(0,0), 
                       breaks = seq(as.numeric(range(stk)['min']), as.numeric(range(stk)['max']),1))+
    labs(colour = "")+ggtitle(name(stk))
  return(p_sel)
}                                                                                
                                                                                                                     
    

