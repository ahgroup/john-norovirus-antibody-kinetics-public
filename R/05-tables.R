###
# 05-tables.R
###

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
        as_flex_table()
      
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
          flextable::bold(part = "header") %>%
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
      flextable::bold(part = "header") %>%
      flextable::autofit()
    
    saveRDS(tab, paste0(filepath, ".rds"))
    flextable::save_as_image(tab, paste0(filepath, ".png"))
    
  }
}

# END OF SCRIPT ====