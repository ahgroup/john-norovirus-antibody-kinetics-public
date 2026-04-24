###
# 06-figures.R
###

# manuscript ====
func_figure_setup <- function() {
  
  dose_levels <- c("4.8","48","4800",
                   "5","15","50","150")
  
  my_palette <- ggokabeito::palette_okabe_ito()
  
  color_scale <- c(
    "4.8" = my_palette[1],
    "48" = my_palette[7],
    "4800" = my_palette[6],
    "5" = my_palette[2],
    "15" = my_palette[5],
    "50" = my_palette[3],
    "150" = my_palette[8]
  )
  
  shape_scale <- c(
    "4.8" = 2, "48" = 6, "4800" = 11,
    "5" = 0, "15" = 5,
    "50" = 7, "150" = 9
  )
  
  list(
    dose_levels = dose_levels,
    color_scale = color_scale,
    shape_scale = shape_scale
  )
}

func_plot_group_model_traj <- function(predictions, data, filepath, height, width) {
  
  library(ggplot2)
  
  if (unique(data$study)[1] == "NI") {
    max_day <- max(data$day, na.rm = TRUE) + 8
  } else {
    max_day <- max(data$day, na.rm = TRUE)
  }
  
  p <- ggplot(predictions) +
    geom_line(aes(x=day, y=y,
                  col=model,
                  linetype=model),
              size=1) +
    geom_ribbon(aes(x=day, ymin=lo, ymax=up, fill=model),
                alpha=0.1, show.legend = FALSE) +
    scale_y_log10(breaks=c(1,10,100,1000),
                  expand=c(0,0)) +
    coord_cartesian(xlim=c(0,max_day), ylim=c(0.5,5000)) +
    geom_point(data=data, 
               mapping=aes(x=day, y=ifelse(cens==0, y, y2),
                           group=id, shape=factor(ifelse(cens==0, "circle", "tri"))),
               size=1) +
    scale_shape_manual(
      values = c("circle" = 1, "tri" = 2),
      guide = "none"
    ) +
    facet_grid(dose~antibody_clean, switch="y") +
    ggokabeito::scale_color_okabe_ito(
      labels=c("exp"="Exponential decay", "pow"="Power-Law decay")
    ) +
    ggokabeito::scale_fill_okabe_ito(
      labels=c("exp"="Exponential decay", "pow"="Power-Law decay"), 
      guide="legend") +  
    scale_linetype_discrete(
      labels=c("exp"="Exponential decay", "pow"="Power-Law decay"), 
      guide="legend") +
    labs(x="Days post-inoculation",
         y="Titer",
         fill="Model",
         linetype="Model",
         col="Model") +
    theme_bw() +
    theme(legend.position="top",
          axis.text = element_text(size=11),
          axis.title = element_text(size=16),
          strip.text = element_text(size=12),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14))
  
  ggsave(
    filename = filepath,
    plot = p,
    width = width,
    height = height
  )
  
  return(p)
  
}

func_plot_measures <- function(dat_group, model_name, study_name,
                               filepath, height, width) {
  
  library(ggplot2)
  library(patchwork)
  
  measures = c(
    "growth_rate",
    "peak_response",
    "peak_time",
    "decay_rate_D90",
    "response_1year"
  )
  
  dat_group_process <- dat_group %>%
    filter(
      model_func==model_name,
      study==study_name,
      measure %in% measures
    ) %>%
    mutate(
      plot_type = factor(case_when(
        measure == "peak_response" ~ "Peak response",
        measure == "peak_time" ~ "Peak time",
        measure == "growth_rate" ~ "Growth rate",
        measure == "decay_rate_D90" & model_name == "pow" ~ "Decay rate D90",
        measure == "decay_rate_D90" & model_name == "exp" ~ "Decay rate",
        measure == "response_1year" ~ "1-year response",
      ), levels = c("Growth rate",
                    "Peak response", 
                    "Peak time",
                    "Decay rate",
                    "Decay rate D90",
                    "1-year response"
      )),
      antibody = factor(antibody),
      dose = factor(dose, 
                    levels = c("4.8","48","4800",
                               "5","15","50","150"))
    )
  
  # get color scale
  full_palette <- ggokabeito::palette_okabe_ito()
  if (study_name == "NI") {
    color_subset <- full_palette[c(1,7,6)]
    growth_y_max <- 0.8
    peak_time_y_max <- 35
    x_axis_size <- 10
  } else {
    color_subset <- full_palette[c(2,5,3,8)]
    growth_y_max <- 3.5 
    peak_time_y_max <- 18
    x_axis_size <- 7
  }
  
  p <- ggplot() +
    geom_pointrange(
      data = dat_group_process,
      aes(x=dose,
          y=est,
          ymin = lo, 
          ymax = up, 
          color = dose,
          shape = dose,
          fill=dose),
      fatten = 3
    ) +
    scale_color_manual(values = color_subset) +
    scale_shape_manual(values = c(21, 24, 22, 25)) +
    scale_fill_manual(values = color_subset) +
    ggh4x::facet_grid2(vars(plot_type), vars(antibody), switch="y",
                       scales="free"
    ) +
    ggh4x::facetted_pos_scales(
      y = list(
        plot_type == "Peak time" ~ scale_y_continuous(limits = c(0, peak_time_y_max)),
        plot_type == "1-year response" ~ scale_y_continuous(),
        plot_type == "Growth rate" ~ scale_y_continuous(limits = c(0, growth_y_max)),
        plot_type == "Peak response" ~ scale_y_continuous(),
        plot_type %in% c("Decay rate D90", "Decay rate") ~ scale_y_continuous()
      )
    ) +
    labs(x = "Dose",
         color="Dose",
         shape="Dose",
         fill="Dose") +
    theme_bw() + 
    theme(axis.title.y = element_blank(),
          legend.position = "bottom",
          axis.text.x = element_text(size=x_axis_size))
  
  ggsave(
    filename = filepath,
    plot = p,
    width = width,
    height = height
  )
  
  return(p)
}

func_plot_measures_vs_1year_draws <- function(setup, 
                                              dat_measures_ni, dat_measures_nv,
                                              ni_draws, nv_draws,
                                              decay_day, main=NULL, 
                                              baseline=FALSE, model_name, 
                                              file_path, height, width) {
  browser()
  
  library(ggplot2)
  library(patchwork)
  
  # process decay day
  decay_name  <- paste0("decay_rate_D", decay_day)
  decay_label <- paste0("Decay rate D", decay_day)
  
  
  antibodies <- if (!is.null(main)) {
    c("GII.4 IgA","GII.4 IgG")
  } else {
    c("GI.1 HBGA", "GI.1 IgA", "GI.1 IgG")
  }
  
  x_measure <- if (baseline) "peak_response" else "peak"
  y_measure <- if (baseline) "response_1year" else "titer_1year"
  y_label <- if (baseline) "Predicted response at 1 year" else "Predicted titer at 1 year"
  
  measure_labels <- c(
    peak = "Peak titer",
    peak_response = "Peak response",
    growth_rate = "Growth rate",
    decay_rate_D90 = "Decay rate D90"
  )
  
  # combine data
  dat_draws <- bind_rows(ni_draws, nv_draws)
  dat_measures <- bind_rows(dat_measures_ni, dat_measures_nv)
  
  # process data
  dat_draws_process <- dat_draws %>%
    filter(
      model_func == model_name,
      antibody %in% antibodies,
    ) %>%
    tidyr::pivot_longer(
      cols = c(x_measure,growth_rate,decay_rate_D90),
      names_to = "measure_name",
      values_to = "measure"
    ) %>%
    mutate(
      measure_name_clean = factor(
        measure_labels[measure_name],
        levels = c("Growth rate", measure_labels[x_measure], decay_label)
      ),
      antibody = factor(antibody),
      dose = factor(dose, levels = setup$dose_levels)
    )
  
  dat_measures_process <- dat_measures %>%
    filter(
      model_func == model_name,
      antibody %in% antibodies,
      measure %in% c(x_measure, "growth_rate", decay_name, y_measure)
    ) %>%
    select(antibody, dose, measure, est) %>%
    tidyr::pivot_wider(
      names_from = measure,
      values_from = est
    ) %>%
    tidyr::pivot_longer(
      cols = c(x_measure,growth_rate,!!decay_name),
      names_to = "measure_name",
      values_to = "measure"
    ) %>%
    mutate(
      measure_name_clean = factor(
        case_when(
          measure_name == decay_name ~ decay_label,
          TRUE ~ measure_labels[measure_name]
        ),
        levels = c("Growth rate", measure_labels[x_measure], decay_label)
      ),
      antibody = factor(antibody),
      dose = factor(dose, levels = setup$dose_levels)
    )
  
  plots <- list()
  
  for (i in seq_along(levels(dat_measures_process$measure_name_clean))) {
    
    m <- levels(dat_measures_process$measure_name_clean)[i]
    dat_measures_m <- dat_measures_process %>% filter(measure_name_clean == m)
    dat_draws_m <- dat_draws_process %>% filter(measure_name_clean == m)
    
    # growth rate scale truncate
    if (m == "Growth rate") {
      dat_draws_m <- dat_draws_m %>%
        group_by(antibody) %>%
        filter(measure <= stats::quantile(measure, 0.95, na.rm = TRUE)) %>%
        ungroup()
    }
    
    # truncate for baseline plot
    if (baseline) {
      y_scale <- scale_y_continuous()
      
      y_col <- rlang::sym(y_measure)
      
      dat_draws_m <- dat_draws_m %>%
        group_by(antibody) %>%
        filter(
          !!y_col >= quantile(!!y_col, 0.10, na.rm = TRUE)
        ) %>%
        ungroup()
    } else {
      y_scale <- scale_y_log10()
    }
    
    x_scale <- if (m == "Peak titer") scale_x_log10() else scale_x_continuous(limits = c(0, NA))
    facet_scales <- if (m == measure_labels[x_measure]) "free_y" else "free"
    
    p <- ggplot(
      dat_draws_m,
      aes(
        x = measure,
        y = !!rlang::sym(y_measure),
        color = dose,
        shape = dose,
        fill = dose
      )
    ) +
      geom_point(alpha = 0.1, size = 0.5) +
      geom_point(
        data = dat_measures_m,
        size = 4,
        stroke = 1
      ) +
      ggh4x::facet_grid2(
        .~antibody, 
        scales = facet_scales, 
        switch = "y",
        independent = "y"
      ) +
      scale_color_manual(values = setup$color_scale) +
      scale_shape_manual(values = setup$shape_scale) +
      scale_fill_manual(values = setup$color_scale) +
      guides(
        color = guide_legend(nrow = 1),
        shape = guide_legend(nrow = 1),
        fill  = guide_legend(nrow = 1)
      ) +
      x_scale +
      y_scale +
      labs(
        x = m,
        y = y_label,
        color="Dose",
        shape="Dose",
        fill="Dose"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        strip.text = element_text(size = 12)
      )
    
    if (i > 1) {
      p <- p + theme(strip.text.x = element_blank())
    } else {
      p <- p + theme(strip.text = element_text(size = 12))
    }
    
    plots[[m]] <- p
    
  }
  
  p_measure <- patchwork::wrap_plots(plots, ncol = 1, guides = "collect") +
    plot_layout(axes = "collect_y") &
    plot_annotation(tag_levels = "A") &
    theme(legend.position="bottom")
  
  ggsave(
    filename = file_path,
    plot = p_measure,
    width = width,
    height = height
  )
  
  return(p_measure)
}

func_plot_group_peak_decay_draws <- function(setup, 
                                             ni_draws, nv_draws, 
                                             ni_measures, nv_measures, 
                                             decay_day, model_type,
                                             main=NULL, baseline=FALSE) {
  library(ggplot2)
  library(patchwork)
  
  if (!is.null(main)) {
    antibodies <- c("GII.4 IgA","GII.4 IgG")
    ncols <- 2
    if (model_type == "pow") {
      limits_peak <- c(NA, 300)
      breaks_peak <- c(30, 100, 300)
    }
    if (model_type == "exp") {
      limits_peak <- c(NA, NA)
      breaks_peak <- c(30, 100, 300)
    }
  } else {
    antibodies <- c("GI.1 HBGA", "GI.1 IgA", "GI.1 IgG")
    ncols <- 3
    if (model_type == "pow") {
      limits_peak <- c(100, 10000)
      breaks_peak <- c(100, 300, 1000, 3000, 10000)
    }
    if (model_type == "exp") {
      limits_peak <- c(100, NA)
      breaks_peak <- c(100, 300, 1000, 3000)
    }
  }
  
  decay_col  <- paste0("decay_rate_D", decay_day)
  decay_sym  <- rlang::sym(decay_col)
  decay_lab  <- paste0("Decay rate D", decay_day)
  
  y_measure <- if (baseline) "peak_response" else "peak"
  y_sym <- rlang::sym(y_measure)
  y_label <- if (baseline) "Peak response" else "Peak titer"
  
  # process data
  dat_draws <- bind_rows(ni_draws, nv_draws)
  dat_measures <- bind_rows(ni_measures, nv_measures)
  
  dat_draws_process <- dat_draws %>%
    filter(
      model_func == model_type,
      antibody %in% antibodies
    ) %>%
    mutate(
      antibody = factor(antibody),
      dose = factor(dose, levels = setup$dose_levels)
    )
  
  dat_measures_process <- dat_measures %>%
    filter(
      model_func == model_type,
      antibody %in% antibodies
    ) %>%
    distinct(study, model_func, antibody, dose, measure, est) %>%
    tidyr::pivot_wider(
      names_from = measure,
      values_from = est
    ) %>%
    mutate(
      antibody = factor(antibody),
      dose = factor(dose, levels = setup$dose_levels)
    )
  
  plots <- list()
  
  for (a in unique(dat_draws_process$antibody)) {
    
    dat_draws_a <- filter(dat_draws_process, antibody == a)
    dat_measures_a <- filter(dat_measures_process, antibody == a)
    
    p <- ggplot(
      dat_draws_a, 
      aes(
        x = !!decay_sym,
        y = !!y_sym,
        color = dose,
        fill  = dose,
        shape = dose
      )
    ) +
      geom_point(alpha = 0.1, size = 0.5) +
      geom_point(
        data = dat_measures_a,
        size = 4,
        stroke = 1
      ) +
      {if (!baseline) scale_y_log10(limits = limits_peak,
                                    breaks = breaks_peak)} +
      {if (baseline) scale_y_continuous(limits = c(NA, 11))} +
      scale_x_continuous(limits = c(0, NA)) +
      scale_color_manual(values = setup$color_scale) +
      scale_shape_manual(values = setup$shape_scale) +
      scale_fill_manual(values = setup$color_scale) +
      labs(
        title=a,
        x = decay_lab,
        y = y_label,
        color = "Dose",
        fill = "Dose",
        shape = "Dose"
      ) +
      theme_bw() +
      guides(
        color = guide_legend(nrow = 1),
        fill  = guide_legend(nrow = 1),
        shape = guide_legend(nrow = 1)
      ) +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=4)
      )
    
    plots[[a]] <- p
    
  }
  
  p_all <- patchwork::wrap_plots(plots, ncol = ncols) + 
    plot_layout(
      axis_titles = "collect",
      guides = "collect"
    ) & 
    theme(
      legend.position = "bottom",
      axis.title = element_text(size = 14),
      legend.text = element_text(size=11),
      strip.text = element_text(size = 12)
    )
  
  return(p_all)
  
}

func_plot_group_peak_decay_1year_draws <- function(setup,
                                                   ni_draws, nv_draws, 
                                                   ni_measures, nv_measures,
                                                   decay_day, model_type,
                                                   main=NULL, baseline=FALSE) {
  library(ggplot2)
  library(patchwork)
  
  if (!is.null(main)) {
    antibodies <- c("GII.4 IgA","GII.4 IgG")
    ncols <- 2
    if (model_type == "pow") {
      limits_peak <- c(NA, 300)
      breaks_peak <- c(30, 100, 300)
    }
    if (model_type == "exp") {
      limits_peak <- c(NA, NA)
      breaks_peak <- c(30, 100, 300)
    }
  } else {
    antibodies <- c("GI.1 HBGA", "GI.1 IgA", "GI.1 IgG")
    ncols <- 3
    if (model_type == "pow") {
      limits_peak <- c(100, 10000)
      breaks_peak <- c(100, 300, 1000, 3000, 10000)
    }
    if (model_type == "exp") {
      limits_peak <- c(100, NA)
      breaks_peak <- c(100, 300, 1000, 3000)
    }
  }
  
  decay_col  <- paste0("decay_rate_D", decay_day)
  decay_sym  <- rlang::sym(decay_col)
  decay_lab  <- paste0("Decay rate D", decay_day)
  
  y_measure <- if (baseline) "peak_response" else "peak"
  y_sym <- rlang::sym(y_measure)
  y_label <- if (baseline) "Peak response" else "Peak titer"
  
  color_var <- if (baseline) "response_1year" else "titer_1year"
  color_sym <- rlang::sym(color_var)
  color_label <- if (baseline) "1 year response" else "1 year titer"
  
  # process data
  dat_draws <- bind_rows(ni_draws, nv_draws)
  dat_measures <- bind_rows(ni_measures, nv_measures)
  
  dat_draws_process <- dat_draws %>%
    filter(
      model_func == model_type,
      antibody %in% antibodies
    ) %>%
    mutate(
      antibody = factor(antibody),
      dose = factor(dose, levels = setup$dose_levels)
    )
  
  dat_measures_process <- dat_measures %>%
    filter(
      model_func == model_type,
      antibody %in% antibodies
    ) %>%
    distinct(model_func, antibody, dose, measure, est) %>%
    tidyr::pivot_wider(
      names_from = measure,
      values_from = est
    ) %>%
    mutate(
      antibody = factor(antibody),
      dose = factor(dose, levels = setup$dose_levels)
    )
  
  plots <- list()
  
  for (a in unique(dat_draws_process$antibody)) {
    
    dat_draws_a <- dat_draws_process %>% filter(antibody == a)
    dat_measures_a <- dat_measures_process %>% filter(antibody == a)
    
    # conditional color/fill scale
    color_scale <- if (!baseline) {
      list(
        scale_color_viridis_c(
          begin = 0, end = 0.9, trans = "log10",
          guide = guide_colorbar(
            title.position = "top", direction = "horizontal", title.hjust = 0.5,
            barwidth = unit(2.5, "cm"), barheight = unit(0.3, "cm"), title = color_label
          )
        ),
        scale_fill_viridis_c(
          begin = 0, end = 0.9, trans = "log10",
          guide = guide_colorbar(
            title.position = "top", direction = "horizontal", title.hjust = 0.5,
            barwidth = unit(2.5, "cm"), barheight = unit(0.3, "cm"), title = color_label
          )
        )
      )
    } else {
      list(
        scale_color_viridis_c(
          begin = 0, end = 0.9,
          guide = guide_colorbar(
            title.position = "top", direction = "horizontal", title.hjust = 0.5,
            barwidth = unit(2.5, "cm"), barheight = unit(0.3, "cm"), title = color_label
          )
        ),
        scale_fill_viridis_c(
          begin = 0, end = 0.9,
          guide = guide_colorbar(
            title.position = "top", direction = "horizontal", title.hjust = 0.5,
            barwidth = unit(2.5, "cm"), barheight = unit(0.3, "cm"), title = color_label
          )
        )
      )
    }
    
    p <- ggplot(
      dat_draws_a, 
      aes(
        x = !!decay_sym,
        y = !!y_sym,
        color = !!color_sym,
        fill  = !!color_sym,
        shape = dose
      )
    ) +
      geom_point(alpha = 0.15, size = 0.5) +
      geom_point(
        data = dat_measures_a,
        size = 4,
        stroke = 1
      ) +
      color_scale +
      {if (!baseline) scale_y_log10(limits = limits_peak,
                                    breaks = breaks_peak)} +
      {if (baseline) scale_y_continuous(limits = c(NA, 11))} +
      scale_x_continuous(limits = c(0, NA)) +
      scale_shape_manual(values = setup$shape_scale) +
      guides(shape = "none") +
      labs(
        title = if (decay_day == 60) a,
        y = y_label,
        x = decay_lab,
        color = color_label,
        fill = color_label,
        shape = "Dose"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        legend.justification = "center",
        legend.box.margin = margin(0,0,0,0),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 13)
      )
    
    if (decay_day != 90) {
      p <- p + theme(legend.position = "none")
    }
    
    plots[[a]] <- p
    
  }
  
  p_all <- wrap_plots(plots, ncol = ncols) +
    plot_layout(axis_titles = "collect")
  
  return(p_all)
  
}

func_plot_comb_peak_decay_vs <- function(p1, p2, filepath, height, width) {
  
  library(patchwork)
  library(ggplot2)
  
  p <- patchwork::wrap_plots(
    list(
      wrap_elements(p1), wrap_elements(p2)
    ),
    ncol = 1
  ) +
    plot_annotation(tag_levels = "A")
  
  ggsave(
    filename = filepath,
    plot = p,
    width = width,
    height = height
  )
  
}

# supplement ====
func_plot_schematic <- function(dat_pred, file_path, width, height) {
  
  library(ggplot2)
  
  dat_schematic <- dat_pred %>% filter(dose == "4800", 
                                       antibody == "GI.1_IgG", 
                                       model == "pow")
  
  peak <- dat_schematic %>% filter(y==max(y))
  peak_day <- peak %>% pull(day)
  half_peak_day <- peak_day * 0.5 # technically growth rate is using 70% peak titer but thats complicated
  baseline <- dat_schematic %>% filter(day == 0)
  response_1year <- dat_schematic %>% filter(day == 365)
  decay_85 <- dat_schematic %>% filter(day %in% c(85))
  decay_95 <- dat_schematic %>% filter(day %in% c(95))
  growth1 <- dat_schematic %>% slice(which.min(abs(day - half_peak_day)))
  growth2 <- dat_schematic %>% slice(which.min(abs(day - (half_peak_day + 1))))
  
  p <- ggplot() +
    geom_line(data = dat_schematic, aes(x = day, y = exp(logy))) +
    # peak 
    geom_point(data = peak, aes(x = day, y = exp(logy)), color = "black", size = 3) +
    geom_text(data = peak, aes(x = day, y = exp(logy), label = "Peak titer"), vjust = -1, size = 3) +
    # peak increase
    geom_text(data = peak, aes(x = day, y = exp(logy), label = "Peak response"), vjust=15, hjust=-0.1, size = 3) +
    geom_segment(
      aes(x = peak$day, xend = peak$day,
          y = exp(baseline$logy), yend = exp(peak$logy)),
      arrow = arrow(ends = "both", length = unit(0.3, "cm"), type = "closed"),
      color = "#56B4E9", size = 0.5
    ) +
    # time of peak
    geom_segment(
      data = peak,
      aes(x=day, xend = day, y = 0, yend=exp(baseline$logy)-1.5),
      color="black", linetype = "solid"
    ) +
    geom_segment(
      x=peak$day, xend = peak$day,
      y = 0, yend=0.65,
      color="#E69F00"
    ) +
    geom_text(data = peak, aes(x = day, y = 3.3, label = "Time of peak"), vjust=0, hjust=-0.05, size = 3) +
    
    # baseline 
    geom_point(data = baseline, aes(x = day, y = exp(logy)), color = "black", size = 3) +
    geom_text(data = baseline, aes(x = day, y = exp(logy), label = "Baseline"), vjust = 1.9, hjust=0.6, size = 3) +
    geom_hline(yintercept = exp(baseline$logy), linetype="dashed") +
    
    # titer 1-year
    geom_text(data = response_1year, aes(x = day, y = exp(logy), label = "1-year titer"), vjust=-2, hjust=1.1, size = 3) +
    
    # response_1year
    geom_point(data = response_1year, aes(x = day, y = exp(logy)), color = "black", size = 3) +
    geom_text(data = response_1year, aes(x = day, y = exp(logy), label = "1-year response"), vjust=4.5, hjust=1.1, size = 3) +
    geom_segment(
      aes(x = 365, xend = 365,
          y = exp(baseline$logy), yend = exp(response_1year$logy)),
      arrow = arrow(ends = "both",length = unit(0.3, "cm"),type = "closed"), 
      color = "#56B4E9", size = 0.5
    ) +
    
    # decay
    geom_text(data = decay_85, aes(x = day, y = exp(logy), label = "Decay rate"), vjust=0, hjust=-0.2, size = 3) +
    geom_segment(
      aes(x = 85, xend = 95,
          y = exp(decay_85$logy), yend = exp(decay_95$logy)),
      color = "#E69F00", size = 1
    ) +
    
    # growth
    geom_text(data = growth1, aes(x = day, y = exp(logy), label = "Growth rate"), vjust=-0.5, hjust=0.5, angle = 85, size = 3) +
    geom_segment(
      aes(x = growth1$day, xend = growth2$day,
          y = exp(growth1$logy), yend = exp(growth2$logy)),
      color = "#E69F00", size = 1
    ) +
    
    scale_y_log10(breaks = c(1, 10, 100, 1000), expand = c(0, 0)) +
    scale_x_continuous(breaks = c(0, 100, 200, 300, 365)) +
    coord_cartesian(xlim = c(-10, 380), ylim = c(3, 1000)) +
    labs(x = "Days post-inoculation", y = "Titer") +
    theme_bw()
  
  ggsave(filename = file_path, plot = p, width = width, height = height)
  
  return(p)
}

func_plot_prior_posterior_draws <- function(draws, file_path, width, height) {
  browser()
  library(ggplot2)
  
  # process draws 
  draws_processed <- draws %>%
    mutate(
      dose = factor(
        case_when(
          stringr::str_detect(parameter,"dose4\\.8") ~ "4.8",
          stringr::str_detect(parameter,"dose4800") ~ "4800",
          stringr::str_detect(parameter,"dose48") ~ "48",
          stringr::str_detect(parameter,"dose50") ~ "50",
          stringr::str_detect(parameter,"dose150") ~ "150",
          stringr::str_detect(parameter,"dose5") ~ "5",
          stringr::str_detect(parameter,"dose15") ~ "15",
          TRUE ~ NA_character_
        ),
        levels = c("4.8","48","4800","5","15","50","150")
      ),
      parameter = case_when(
        stringr::str_detect(parameter,"b_d") ~ "d",
        stringr::str_detect(parameter,"b_k") ~ "k",
        stringr::str_detect(parameter,"b_g") ~ "g",
        stringr::str_detect(parameter,"b_p") ~ "p",

        stringr::str_detect(parameter,"sigma") ~ "sigma",
        stringr::str_detect(parameter,"nu") ~ "nu",
        stringr::str_detect(parameter,"p_Intercept") ~ "p_Intercept",
        stringr::str_detect(parameter,"g_Intercept") ~ "g_Intercept",
        stringr::str_detect(parameter,"k_Intercept") ~ "k_Intercept",
        stringr::str_detect(parameter,"d_Intercept") ~ "d_Intercept",
        TRUE ~ NA_character_
      )
    )
  
  # store plots
  plots <- list()
  
  for (m in unique(draws_processed$model)) {
    
    name <- paste(
      stringr::str_sub(m, 1, 2),
      ifelse(stringr::str_detect(m, "exp"), "Exponential",
             ifelse(stringr::str_detect(m, "pow"), "Power-Law", "")),
      stringr::str_split_fixed(m, "_", 3)[,3] %>% 
        stringr::str_replace_all("_", " ")
    )
    
    draws_sub <- draws_processed %>% filter(model == m)
    
    plots[[m]] <- ggplot(draws_sub, aes(x = value, fill = type)) +
      geom_density(alpha = 0.5) +
      ggh4x::facet_grid2(
        rows = vars(dose), 
        cols = vars(parameter), 
        scales = "free", 
        independent = "y", 
        switch = "y"
      ) +
      theme_bw() +
      labs(title = name, fill = "") + 
      theme(
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm")
      )
  }
  
  for (m in names(plots)) {
    ggsave(
      filename = paste0(file_path, "/", m, ".png"),
      plot = plots[[m]],
      width = width,
      height = height
    )
  }
  
  return(plots)
  
}

func_plot_fitted_vs_observed <- function(residuals, file_path, width, height) {
  
  library(ggplot2)
  library(patchwork)
  
  plots <- list()
  
  for (m in unique(residuals$model)) {
    
    name <- paste(
      stringr::str_sub(m, 1, 2),
      ifelse(stringr::str_detect(m, "exp"), "Exponential",
             ifelse(stringr::str_detect(m, "pow"), "Power-Law", "")),
      stringr::str_split_fixed(m, "_", 3)[,3] %>% 
        stringr::str_replace_all("_", " ")
    )
    
    residuals_sub <- residuals %>% filter(model == m)
    
    # custom colors
    if (unique(residuals_sub$study) == "NI") {
      dose_colors <- ggokabeito::palette_okabe_ito(c(1,7,6))
    } else {
      dose_colors <- ggokabeito::palette_okabe_ito(c(2,5,3,8))
    }
    
    names(dose_colors) <- sort(unique(residuals_sub$dose))
    
    max <- max(max(residuals_sub$y_obs), max(residuals_sub$y_pred)) + 5
    
    plots[[m]] <- ggplot(residuals_sub, aes(x=y_obs, y=y_pred, color=dose)) +
      geom_point(alpha=0.7) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
      coord_cartesian(xlim = c(0,max), ylim = c(0,max)) +
      scale_color_manual(values = dose_colors) +
      theme_bw() +
      labs(title=name,
           y="Predicted value",
           x="Observed value") +
      theme(
        plot.title = element_text(size=10)
      )
  }
  
  combined <- wrap_plots(plots, ncol = 3, guides = "collect") & 
    theme(
      legend.position = "bottom"
    )
  
  ggsave(
    filename = file_path,
    plot = combined,
    width = width,
    height = height
  )
  
  return(plots)
  
}

func_plot_id_model <- function(predictions, data, ncols, zoom=FALSE,
                               filepath, height, width) {
  
  library(ggplot2)
  library(patchwork)
  
  data_sort <- data %>% mutate(
    dose = factor(
      dose, 
      levels = c("4.8","48","4800","5","15","50","150")
    )
  )
  
  for (ab in unique(data_sort$antibody_clean)) {
    
    predictions_ab <- predictions %>% filter(antibody_clean == ab)
    data_ab <- data_sort %>% filter(antibody_clean == ab)
    
    plots_id <- list()
    
    for (d in levels(data_sort$dose)) {
      
      predictions_d <- predictions_ab %>% filter(dose == d)
      data_d <- data_ab %>% filter(dose == d)
      
      for (idx in unique(data_d$id)) {
        
        predictions_id <- predictions_d %>% filter(id == idx)
        data_id <- data_d %>% filter(id == idx)
        
        p_id <- ggplot(predictions_id) +
          geom_line(aes(x=day, y=y,
                        col=model,
                        linetype=model),
                    size=1) +
          geom_ribbon(aes(x=day, ymin=lo, ymax=up, fill=model),
                      alpha=0.1, show.legend = FALSE) +
          scale_y_log10(breaks=c(1,10,100,1000),
                        expand=c(0,0)) +
          coord_cartesian(xlim=c(0,max(data_ab$day)+5), ylim=c(0.5,max(predictions_ab$up)+5)) +
          geom_point(data=data_id, 
                     mapping=aes(x=day, y=ifelse(cens==0, y, y2),
                                 group=idx, shape=factor(ifelse(cens==0, "circle", "tri"))),
                     size=1) +
          scale_shape_manual(
            values = c("circle" = 1, "tri" = 2),
            guide = "none"
          ) +
          ggokabeito::scale_color_okabe_ito(
            labels=c("exp"="Exponential decay", "pow"="Power-Law decay")
          ) +
          ggokabeito::scale_fill_okabe_ito(
            labels=c("exp"="Exponential decay", "pow"="Power-Law decay"), 
            guide="legend") +  
          scale_linetype_discrete(
            labels=c("exp"="Exponential decay", "pow"="Power-Law decay"), 
            guide="legend") +
          labs(fill="Model",
               linetype="Model",
               col="Model",
               title=paste0(idx," ", paste(unique(data_id$dose)))) +
          theme_bw() +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                legend.text = element_text(size=12),
                legend.title = element_text(size=14))
        
        if (zoom) {
          
          p_id <- p_id +
            coord_cartesian(xlim=c(0,60), 
                            ylim=c(0.5,max(predictions_ab$up)+5)) +
            geom_vline(xintercept=28, linetype="dashed") +
            geom_point(data=data_id %>% filter(y==max(y)),
                       aes(x=day, y=y),
                       shape=16, size=2)
          
        }
        
        plots_id[[paste0(d, "_", idx)]] <- p_id
        
      }
      
    }
    
    p_ab <- patchwork::wrap_plots(plots_id, ncol = ncols, guides = "collect") + 
      plot_annotation(title = paste(ab),
                      theme = theme(
                        plot.title = element_text(size = 20)
                      )) & 
      theme(legend.position="bottom")
    
    ggsave(
      filename = file.path(filepath, paste0(ab, ".png")),
      plot = p_ab,
      width = width,
      height = height
    )
    
  }
}

func_plot_decay_rates <- function(setup, dat_group, model_name, 
                                  study_name, filepath, height, width) {
  library(ggplot2)
  library(patchwork)
  
  measures = c(
    "decay_rate_D60",
    "decay_rate_D90",
    "decay_rate_D180",
    "decay_rate_D365"
  )
  
  measure_labels <- c(
    decay_rate_D60  = "Decay rate D60",
    decay_rate_D90  = "Decay rate D90",
    decay_rate_D180 = "Decay rate D180",
    decay_rate_D365 = "Decay rate D365"
  )
  
  dat_group_process <- dat_group %>%
    filter(
      model_func==model_name,
      study==study_name,
      measure %in% measures
    ) %>%
    mutate(
      plot_type = factor(
        measure_labels[measure],
        levels = measure_labels
      ),
      antibody = factor(antibody),
      dose = factor(dose, levels = setup$dose_levels),
    )
  
  p <- ggplot() +
    geom_pointrange(
      data = dat_group_process,
      aes(x=dose,
          y=est,
          ymin = lo, 
          ymax = up, 
          color = dose,
          shape = dose,
          fill=dose),
      fatten = 3
    ) +
    scale_color_manual(values = setup$color_scale) +
    scale_fill_manual(values = setup$color_scale) +
    scale_shape_manual(values = c(21, 24, 22, 25)) +
    ggh4x::facet_grid2(vars(plot_type), vars(antibody), 
                       switch="y", scales="free") +
    labs(
      x = "Dose",
      color="Dose",
      shape="Dose",
      fill="Dose"
    ) +
    theme_bw() + 
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.justification = "center",
      legend.box.margin = margin(0, 0, 0, 0),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10)
    )
  
  ggsave(
    filename = filepath,
    plot = p,
    width = width,
    height = height
  )
  
  return(p)
}

func_plot_absolute_measures <- function(setup, dat_group, 
                                        model_name, study_name, 
                                        filepath, height, width) {
  library(ggplot2)
  library(patchwork)
  
  measures = c(
    "baseline",
    "peak",
    "peak_response",
    "titer_1year",
    "response_1year"
  )
  
  dat_group_process <- dat_group %>%
    filter(
      model_func==model_name,
      study==study_name,
      measure %in% measures
    ) %>%
    mutate(
      plot_type = factor(case_when(
        measure == "baseline" ~ "Baseline",
        measure == "peak" ~ "Peak titer",
        measure == "peak_response" ~ "Peak response",
        measure == "titer_1year" ~ "1-year titer",
        measure == "response_1year" ~ "1-year response",
      ), levels = c("Baseline",
                    "Peak titer",
                    "Peak response",
                    "1-year titer",
                    "1-year response")),
      antibody = factor(antibody),
      dose = factor(dose, 
                    levels = c("4.8","48","4800",
                               "5","15","50","150"))
    )
  
  # get color scale
  full_palette <- ggokabeito::palette_okabe_ito()
  if (study_name == "NI") {
    color_subset <- full_palette[c(1,7,6)]
    growth_y_max <- .05
    peak_time_y_max <- 30
    x_axis_size <- 10
  } else {
    color_subset <- full_palette[c(2,5,3,8)]
    growth_y_max <- .5
    peak_time_y_max <- 18
    x_axis_size <- 7
  }
  
  p <- ggplot() +
    geom_pointrange(
      data = dat_group_process,
      aes(x=dose,
          y=est,
          ymin = lo, 
          ymax = up, 
          color = dose,
          shape = dose,
          fill=dose),
      fatten = 3
    ) +
    scale_color_manual(values = color_subset) +
    scale_fill_manual(values = color_subset) +
    scale_shape_manual(values = c(21, 24, 22, 25)) +
    ggh4x::facet_grid2(vars(plot_type), vars(antibody), switch="y",
                       scales="free"
    ) +
    ggh4x::facetted_pos_scales(
      y = list(
        plot_type == "Baseline" ~ scale_y_log10(breaks=c(0.1, 1, 10, 100)),
        plot_type %in% c("Baseline","Peak titer","1-year titer") ~ scale_y_log10(),
        plot_type %in% c("Peak response","1-year response") ~ scale_y_continuous()
      )
    ) +
    labs(x = "Dose",
         color="Dose",
         shape="Dose",
         fill="Dose") +
    theme_bw() + 
    theme(axis.title.y = element_blank(),
          legend.position = "bottom",
          axis.text = element_text(size=7))
  
  ggsave(
    filename = filepath,
    plot = p,
    width = width,
    height = height
  )
  
  return(p)
}

func_plot_comb_group_peak_decay_1year_draws <- function(p1,p2,p3,
                                                        filepath, height, width) {
  
  library(patchwork)
  library(ggplot2)
  
  plots <- list(p1, p2, p3)
  plots <- lapply(plots, function(p) p + theme(legend.position = "none"))
  
  p <- wrap_plots(plots, ncol = 1, guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(
    filename = filepath,
    plot = p,
    width = width,
    height = height
  )
  
}

func_plot_pairwise_compare <- function(probabilities, model, 
                                       filepath, height, width) {
  browser()
  library(ggplot2)
  
  dat_prob <- probabilities %>% 
    filter(model_func == model) %>%
    mutate(
      measure = factor(measure, 
                       levels = rev(c("growth_rate","peak_response","peak_time",
                                      "decay_rate_D90","response_1year",
                                      "peak","titer_1year")),
                       labels = rev(c("Growth rate","Peak response","Peak time",
                                      "Decay rate D90","1-year response",
                                      "Peak titer","1-year titer"))),
      comparison = factor(comparison,
                          levels = c("P(4.8 > 48)","P(4.8 > 4800)","P(48 > 4800)",
                                     "P(5 > 15)","P(5 > 50)","P(5 > 150)",
                                     "P(15 > 50)","P(15 > 150)",
                                     "P(50 > 150)")),
      prob_label = case_when(
        round(prob, 2) == 1 ~ ">0.99",
        round(prob, 2) == 0 ~ "<0.01",
        TRUE ~ as.character(round(prob, 2))
      )
    )
  
  p <- dat_prob %>%
    ggplot(aes(x=comparison, y=measure, fill=prob)) +
    geom_tile(color = "white") +
    geom_text(aes(label = prob_label), size = 3) +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white",
      midpoint = 0.5, limit = c(0,1), name = "Probability"
    ) +
    facet_wrap(~ antibody, ncol=1) +
    labs(
      x = NULL,
      y = "Measure"
    ) +
    theme_minimal()
  
  ggsave(
    filename = filepath,
    plot = p,
    height = height,
    width = width
  )
  
}

# END OF SCRIPT ====