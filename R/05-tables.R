###
# 05-tables.R
###

func_tab_model_diagnostics <- function(models, filepath) {
  
  for (m in names(models)) {
    
    title <-  m |>
      sub("_exp_", "_Exponential_", x = _) |>
      sub("_pow_", "_Power-Law_", x = _) |>
      gsub("_", " ", x = _)
    
    fixed_summary <- as.data.frame(summary(models[[m]])$fixed) %>%
      tibble::rownames_to_column(var = "Parameter") %>%
      mutate(
        across(
          where(is.numeric) & !matches("ESS|Rhat"),
          ~ round(.x, 2)
        ),
        across(
          matches("ESS|n_eff"),
          ~ round(.x, 0)
        ),
        Rhat = sprintf("%.2f", Rhat)
      )
    
    tab <- fixed_summary %>%
      flextable::flextable() %>%
      flextable::add_header_lines(
        values = paste(title)
      ) %>%
      flextable::autofit()
    
    saveRDS(
      tab,
      file = paste0(filepath, m, ".rds")
    )
    
  }
  
}


func_tab_model_metrics <- function(metrics_ni, metrics_nv, filepath) {
  
  metrics <- bind_rows(metrics_ni, metrics_nv) %>%
    mutate(
      antibody = stringr::str_replace(antibody,"_"," "),
      across(c(delta_elpd, se_diff, p_adj_pow_better),~ round(.x, 3))
    )
  
  tab_metrics <- metrics %>%
    flextable::flextable(
      col_keys = c(
        "study", "antibody",
        "delta_elpd", "se_diff", "p_adj_pow_better"
      )
    ) %>%
    flextable::set_header_labels(
      study = "Study",
      antibody = "Antibody",
      
      delta_elpd = "ΔELPD",
      se_diff = "ΔSE",
      p_adj_pow_better = "P(ΔELPD>0)"
    ) %>%
    flextable::align(align = "center", part = "all")
  
  saveRDS(tab_metrics, filepath)
  
}

func_tab_random_effect_sd <- function(models, filepath) {
  browser()
  results <- list()
  
  for (m in names(models)) {
    
    fit <- models[[m]]
    
    draws <- brms::as_draws_df(fit)
    
    # random effect SD parameters
    sd_pars <- c(
      p = "sd_id__p_Intercept",
      g = "sd_id__g_Intercept",
      k = "sd_id__k_Intercept",
      d = "sd_id__d_Intercept"
    )
    
    df <- purrr::map_dfr(names(sd_pars), function(par) {
      
      x <- draws[[sd_pars[[par]]]]
      
      tibble(
        parameter = par,
        estimate = median(x),
        lower = quantile(x, 0.025),
        upper = quantile(x, 0.975)
      )
      
    })
    
    # add model information
    df <- df %>%
      mutate(
        study = stringr::str_extract(m, "^[^_]+"),
        model = case_when(
          stringr::str_detect(m, "_exp_") ~ "Exponential",
          stringr::str_detect(m, "_pow_") ~ "Power-law",
          TRUE ~ NA_character_
        ),
        antibody = stringr::str_remove(
          m,
          "^[^_]+_(exp|pow)_"
        ) %>%
          stringr::str_replace_all("_", " "),
        model_name = m
      ) %>%
      select(
        study,
        model,
        antibody,
        parameter,
        estimate,
        lower,
        upper
      )
    
    results[[m]] <- df
    
  }
  
  bind <- bind_rows(results)
  
  tab <- bind %>%
    mutate(
      SD = paste0(
        round(estimate, 2),
        " (",
        round(lower, 2),
        "-",
        round(upper, 2),
        ")"
      )
    ) %>%
    select(
      study,
      model,
      antibody,
      parameter,
      SD
    )
  
  tab_split <- split(tab, tab$model)
  
  for (t in names(tab_split)) {
    
    dat <- tab_split[[t]]
    
    title <- paste(unique(dat$study), t)
    
    tab <- dat %>%
      select(antibody, parameter, SD) %>%
      flextable::flextable() %>%
      flextable::add_header_lines(
        values = title
      ) %>%
      flextable::autofit()
    
    saveRDS(
      tab,
      file.path(filepath, paste0("tab-", unique(dat$study),"-", t, "-random-effect-SD.rds"))
    )
  }
  
}

func_tab_n_cens <- function(data, filepath) {
  
  library(gtsummary)
  
  # get clean time points for ni study
  dat_timepoint <- data %>%
    mutate(
      timepoint = ifelse(
        study == "NI",
        case_when(
          day == 0 ~ 0,
          day == 2 ~ 2,
          day == 7 ~ 7,
          day %in% c(13,14,15) ~ 14,
          day %in% c(27,28,29,30,32) ~ 28,
          day %in% c(163,174,175,179,180,181,182,183,187) ~ 180,
          TRUE ~ day
        ),
        day
      ),
      cens = factor(
        cens, 
        levels = c(0,1,2),
        labels = c("Not censored","Above limit","Below limit")
      )
    ) %>%
    # get ordered group variable
    mutate(
      antibody_clean = factor(antibody_clean, levels = sort(unique(antibody_clean))),
      dose = factor(dose, levels = sort(unique(dose)))
    )
  
  tab_list <- list()
  
  # by antibody-dose combination
  for (ab in levels(dat_timepoint$antibody_clean)) {
    for (d in levels(dat_timepoint$dose)) {
      
      dat_sub <- dat_timepoint %>%
        filter(antibody_clean == ab, dose == d)
      
      tab <- dat_sub %>%
        tbl_summary(
          by = cens,
          include = timepoint,
          percent = "row"
        ) %>%
        add_n() %>%
        as_flex_table() %>%
        flextable::add_header_lines(
          values = paste(unique(data$study), ab, d)
        ) %>%
        flextable::autofit()
      
      saveRDS(tab, paste0(filepath, ab, "-", d, ".rds"))
      
      tab_list[[paste0(ab,"_",d)]] <- tab
    }
  }
  
  return(tab_list)
  
}

func_tab_correlations <- function(correlations, model, by, baseline, filepath) {
  
  dat_filter <- correlations %>%
    filter(
      model_func == model
    ) %>%
    mutate(
      across(
        starts_with("cor"),
        ~ round(.x, 2)
      )
    )
  
  measure_vars <- if (baseline) {
    c(
      "Growth rate and 1-year response" = "cor_resp_growth",
      "Peak response and 1-year response" = "cor_resp_peakresp",
      "Decay rate D90 and 1-year response" = "cor_resp_D90",
      "Peak response and Growth rate" = "cor_peakresp_growth",
      "Peak response and Decay rate D90" = "cor_peakresp_D90"
    )
  } else {
    c(
      "Growth rate and 1-year titer" = "cor_titer_growth",
      "Peak and 1-year titer" = "cor_titer_peak",
      "Decay rate D90 and 1-year titer" = "cor_titer_D90",
      "Peak and Growth rate" = "cor_peak_growth",
      "Peak and Decay rate D90" = "cor_peak_D90"
    )
  }
  
  if (by == "dose") {
    
    studies   <- unique(dat_filter$study)
    antibodies <- unique(dat_filter$antibody)
    
    for (s in studies) {
      for (ab in antibodies) {
        
        dat_sub <- dat_filter %>%
          filter(study == s, antibody == ab)
        
        if (nrow(dat_sub) == 0) next
        
        dat_wide <- dat_sub %>%
          select(dose, all_of(unname(measure_vars))) %>%
          rename(!!!setNames(unname(measure_vars), names(measure_vars))) %>%
          tidyr::pivot_longer(
            cols = all_of(names(measure_vars)),
            names_to = "Measure",
            values_to = "value"
          ) %>%
          tidyr::pivot_wider(
            names_from  = dose,
            values_from = value
          )
        
        tab <- flextable::flextable(dat_wide) %>%
          flextable::set_header_labels(Measure = "Measure") %>%
          flextable::add_header_lines(values = paste0(s," ", ab)) %>%
          flextable::autofit()
        
        safe_s  <- gsub("[^A-Za-z0-9]", "_", s)
        safe_ab <- gsub("[^A-Za-z0-9]", "_", ab)
        fpath   <- paste0(filepath, "_", safe_s, "_", safe_ab)
        
        saveRDS(tab, paste0(fpath, ".rds"))
        flextable::save_as_image(tab, paste0(fpath, ".png"))
      }
    }
  } else if (by == "overall") {
    
    dat_wide <- dat_filter %>%
      select(antibody, all_of(unname(measure_vars))) %>%
      rename(!!!setNames(unname(measure_vars), names(measure_vars))) %>%
      tidyr::pivot_longer(
        cols      = all_of(names(measure_vars)),
        names_to  = "Measure",
        values_to = "value"
      ) %>%
      tidyr::pivot_wider(
        names_from  = antibody,
        values_from = value
      )
    
    tab <- flextable::flextable(dat_wide) %>%
      flextable::set_header_labels(Measure = "Measure") %>%
      flextable::autofit()
    
    saveRDS(tab, paste0(filepath, ".rds"))
    flextable::save_as_image(tab, paste0(filepath, ".png"))
    
  }
}

func_tab_percent_loss <- function(predictions, filepath) {
  
  # get percent ab loss
  antibody_loss <- predictions %>%
    group_by(study, model, antibody, antibody_clean, dose) %>%
    summarise(
      peak_day = day[which.max(y)],
      peak_ab = max(y, na.rm = TRUE),
      ab_365 = y[day == 365],
      pct_loss = 100 * (peak_ab - ab_365) / peak_ab,
      .groups = "drop"
    ) %>%
    select(study, antibody_clean, dose, model, pct_loss) %>%
    tidyr::pivot_wider(
      names_from = model,
      values_from = pct_loss,
      names_glue = "{model}_pct_loss"
    ) %>%
    mutate(
      exp_pct_loss = round(exp_pct_loss, 1),
      pow_pct_loss = round(pow_pct_loss, 1)
    ) %>%
    arrange(study, antibody_clean, dose)
  
  tab <- antibody_loss %>%
    rename(
      Study = study,
      Antibody = antibody_clean,
      Dose = dose,
      `Exponential decay\n(% loss)` = exp_pct_loss,
      `Power-law decay\n(% loss)` = pow_pct_loss
    ) %>%
    flextable::flextable() %>%
    flextable::autofit()
  
  saveRDS(tab, filepath)
  
}

func_tab_dose_median_iqr <- function(data, filepath) {
  browser()
  dat_timepoint <- data %>%
    mutate(
      timepoint = ifelse(
        study == "NI",
        case_when(
          day == 0 ~ 0,
          day == 2 ~ 2,
          day == 7 ~ 7,
          day %in% c(13,14,15) ~ 14,
          day %in% c(27,28,29,30,32) ~ 28,
          day %in% c(163,174,175,179,180,181,182,183,187) ~ 180,
          TRUE ~ day
        ),
        day
      )
    )
  
  tab_group <- dat_timepoint %>%
    group_by(antibody_clean, dose, timepoint) %>%
    summarise(
      value = paste0(
        round(median(y2, na.rm = TRUE), 1),
        " (",
        round(quantile(y2, 0.25, na.rm = TRUE), 1),
        "–",
        round(quantile(y2, 0.75, na.rm = TRUE), 1),
        ")"
      ),
      .groups = "drop"
    )
  
  # Create one flextable per antibody
  tabs <- tab_group %>%
    split(.$antibody_clean) %>%
    lapply(function(x) {
      
      antibody <- unique(x$antibody_clean)
      study <- unique(data$study)
      
      x %>%
        select(timepoint, dose, value) %>%
        tidyr::pivot_wider(
          names_from = dose,
          values_from = value
        ) %>%
        flextable::flextable() %>%
        flextable::autofit() %>%
        flextable::set_header_labels(
          timepoint = "Study Day"
        ) %>%
        flextable::add_header_lines(
          values = paste(study, antibody)
        ) %>%
        flextable::align(
          i = 1,
          part = "header",
          align = "left"
        )
      
    })
  
  for (antibody in names(tabs)) {
    
    clean_ab <- stringr::str_replace(antibody, " ", "_")
    
    saveRDS(
      tabs[[antibody]],
      file = paste0(filepath, clean_ab, ".rds")
    )
  }
}

func_tab_iga_igg <- function(draws, filepath) {
  
  dat_ratio <- draws %>%
    filter(antibody %in% c("GI.1 IgA","GI.1 IgG","GII.4 IgA","GII.4 IgG")) %>%
    select(study, model_func, dose, draw, antibody, peak_response) %>%
    tidyr::pivot_wider(
      names_from = antibody,
      values_from = peak_response
    ) %>%
    mutate(
      g1_ratio = `GI.1 IgA` / `GI.1 IgG`
    )
  
  summarize_ratio <- function(x) {
    s <- brms::posterior_summary(
      x,
      probs = c(0.025, 0.975),
      robust = TRUE
    )
    
    paste0(
      round(s[1], 2),
      " (",
      round(s[3], 2),
      "-",
      round(s[4], 2),
      ")"
    )
  }
  
  g1_ratio_summary <- dat_ratio %>%
    group_by(model_func, dose) %>%
    summarise(
      `GI.1 IgA:IgG` = summarize_ratio(g1_ratio),
      .groups = "drop"
    )
  
  if (unique(dat_ratio$study) == "NV") {
    dat_ratio <- dat_ratio %>%
      mutate(
        g2_ratio = `GII.4 IgA` / `GII.4 IgG`
      )
    
    g2_ratio_summary <- dat_ratio %>%
      group_by(model_func, dose) %>%
      summarise(
        `GII.4 IgA:IgG` = summarize_ratio(g2_ratio),
        .groups = "drop"
      )
    
    tab <- g1_ratio_summary %>%
      left_join(
        g2_ratio_summary,
        by = c("model_func", "dose")
      )
    
  } else {
    
    tab <- g1_ratio_summary
    
  }
  
  tab_format <- tab %>% 
    mutate(
      model_func = recode(
        model_func,
        pow = "Power-law",
        exp = "Exponential"
      )
    ) %>%
    select(
      Model = model_func,
      Dose = dose,
      everything()
    ) %>%
    flextable::flextable() %>%
    flextable::add_header_lines(
      values = paste(unique(dat_ratio$study))
    ) %>%
    flextable::autofit()
  
  saveRDS(tab_format, filepath)
  
}

func_tab_prior_sens <- function(measures_ni_main, measures_nv_main, 
                                measures_ni_sens, measures_nv_sens, filepath) {
  browser()
  # combine all data
  dat <- bind_rows(
    measures_ni_main %>% mutate(prior_spec = "Main"),
    measures_nv_main %>% mutate(prior_spec = "Main"),
    measures_ni_sens %>% mutate(prior_spec = "Sensitivity"),
    measures_nv_sens %>% mutate(prior_spec = "Sensitivity")
  ) %>%
    filter(
      measure %in% c("peak_response", "peak_time", "growth_rate",
                     "decay_rate_D90", "response_1year")
    ) %>%
    mutate(
      digits   = if_else(measure == "decay_rate_D90", 4, 2),
      estimate = paste0(
        round(est, digits), " (",
        round(lo,  digits), ", ",
        round(up,  digits), ")"
      ),
      model_func = case_when(
        model_func == "exp" ~ "Exponential",
        model_func == "pow" ~ "Power-law"
      )
    ) %>%
    select(-digits) %>%
    arrange(dose)
  
  studies    <- unique(dat$study)
  antibodies <- unique(dat$antibody)
  models     <- unique(dat$model_func)
  
  for (s in studies) {
    for (ab in antibodies) {
      for (m in models) {
        
        dat_sub <- dat %>%
          filter(study == s, antibody == ab, model_func == m)
        
        if (nrow(dat_sub) == 0) next
        
        dat_wide <- dat_sub %>%
          select(measure, dose, prior_spec, estimate) %>%
          tidyr::pivot_wider(
            names_from  = prior_spec,
            values_from = estimate
          ) %>%
          rename(
            Measure = measure,
            Dose    = dose
          )
        
        tab <- flextable::flextable(dat_wide) %>%
          flextable::merge_v(j = "Measure") %>%
          flextable::valign(j = "Measure", valign = "top") %>%
          flextable::add_header_row(
            values    = c("", "Prior specification"),
            colwidths = c(2, 2)
          ) %>%
          flextable::add_header_lines(
            values = paste(s, ab, m)
          ) %>%
          flextable::autofit()
        
        safe_ab <- gsub(" ", "_", ab)
        fpath   <- paste0(filepath, s, "-", safe_ab, "-", m)
        
        saveRDS(tab, paste0(fpath, ".rds"))
      }
    }
  }
  
}

func_tab_nv_second_dose <- function(data, filepath1, filepath2) {
  
  dat_increase <- data %>%
    filter(day %in% c(28, 35)) %>%
    select(study, dose, antibody_clean, id, day, run, y) %>%
    tidyr::pivot_wider(
      id_cols = c(dose, antibody_clean, id, run),
      names_from = day,
      values_from = y,
      names_prefix = "day_"
    ) %>%
    # drop NA either days
    filter(!is.na(day_28), 
           !is.na(day_35)) %>%
    #for multiple runs TRUE if at least 1 greater
    group_by(id, dose, antibody_clean) %>%
    summarise(
      increase = any(day_35 > day_28, na.rm = TRUE),
      .groups = "drop"
    )
  
  # checked for duplicate days with same max and none found
  dat_peak <- data %>%
    group_by(id, dose, antibody_clean) %>%
    slice_max(y, n = 1, with_ties = FALSE) %>%
    select(id, dose, antibody_clean, day)
  
  tab_increase <- dat_increase %>%
    group_by(antibody_clean, dose) %>%
    rename(Antibody=antibody_clean) %>%
    summarise(
      cell = paste0(
        sum(increase), "/", n(),
        " (", round(100 * mean(increase)), "%)"
      ),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = dose,
      values_from = cell,
      names_prefix = "Dose "
    ) %>%
    flextable::flextable() %>%
    flextable::add_header_lines("Participants with increased antibody response after second vaccination") %>%
    flextable::autofit()
  
  tab_peak <- list()
  
  for (ab in unique(dat_peak$antibody_clean)) {
    
    tab_peak[[ab]] <- dat_peak %>%
      filter(antibody_clean == ab) %>%
      rename(Dose=dose) %>%
      gtsummary::tbl_summary(
        include = Dose,
        by = day,
        percent = "row",
        missing = "no"
      ) %>%
      gtsummary::modify_header(label ~ "**Day**") %>%
      gtsummary::as_flex_table() %>%
      flextable::add_header_lines(values = paste0(ab, "; Timing of maximum observed antibody response"))
  }
  
  saveRDS(tab_increase, filepath1)
  
  for (ab in names(tab_peak)) {
    safe_ab <- gsub(" ", "_", ab)
    saveRDS(tab_peak[[ab]], paste0(filepath2, "-", safe_ab, ".rds"))
  }
  
}

# END OF SCRIPT ====