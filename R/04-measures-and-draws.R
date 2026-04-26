###
# 04-measures-and-draws.R
###

func_get_measures_group <- function(data, models, preds) {
  
  library(brms)
  
  measures <- list()
  ci_probs <- c(0.025, 0.975)
  
  for (model_name in names(models)) {
    
    measures[[model_name]] <- list()
    
    # get current model
    model <- models[[model_name]]
    model_preds <- preds[[model_name]]
    
    # get doses in current model
    doses <- unique(model$data$dose)
    
    # dose loop
    for (dose_current in doses) {
      
      epred_log <- model_preds[[dose_current]]$epred
      # exponentiate to get non-scaled
      epred <- exp(epred_log)
      timestamps <- model_preds[[dose_current]]$predictions$day
      
      # time stamp in cols, draws in rows
      rownames(epred) <- make.unique(paste0("ndraws_", seq_len(nrow(epred))))
      colnames(epred) <- make.unique(paste0("timestamp_", seq_len(ncol(epred))))
      
      # peak
      peak_draws <- apply(epred, 1, max)
      peak_summary <- posterior_summary(
        peak_draws,
        probs = ci_probs,
        robust = TRUE
      )
      
      # time of peak
      ## get index
      peak_index <- apply(epred, 1, which.max)
      ## get timestamp
      peak_time_draws <- timestamps[peak_index]
      peak_time_summary <- posterior_summary(
        peak_time_draws,
        probs = ci_probs,
        robust = TRUE
      )
      
      # growth rate
      ## using log epred
      ## at "50% of peak titer"
      ## get 50% peak values
      peak_frac <- 0.5
      frac_peak_log <- log(peak_draws * peak_frac)
      # initalize output vector
      n_draws <- nrow(epred_log)
      growth_draws <- numeric(n_draws)
      growth_draws[] <- NA
      
      # loop over each draw
      for (i in seq_len(n_draws)) {
        y <- epred_log[i, ]
        t <- timestamps
        pk_idx <- peak_index[i]
        # restict to "pre-peak" values
        y_pre <- y[1:(pk_idx-1)]
        t_pre <- t[1:(pk_idx-1)]
        if (length(y_pre) < 2) next # if fewer than 2 pts, skip
        # find interval that contains 70% of peak
        below <- which(y_pre < frac_peak_log[i])
        above <- which(y_pre >= frac_peak_log[i])
        
        if (length(below) == 0 | length(above) == 0) next
        i1 <- max(below)
        i2 <- min(above)
        
        # linear interpolate to exact 70% peak time
        t_frac <- t_pre[i1] + (frac_peak_log[i] - y_pre[i1]) / (y_pre[i2] - y_pre[i1]) * (t_pre[i2] - t_pre[i1])
        
        # convert to log values
        log_y1 <- y_pre[i1]
        log_y2 <- y_pre[i2]
        
        # compute growth rate (slope btwn two pts)
        growth_draws[i] <- (log_y2 - log_y1) / (t_pre[i2] - t_pre[i1])
        
      }
      
      growth_summary <- posterior_summary(
        growth_draws,
        probs = ci_probs,
        robust = TRUE
      )
      
      # baseline
      baseline <- epred[cbind(seq_len(nrow(epred)), 1)]
      baseline_summary <- posterior_summary(
        baseline,
        probs = ci_probs,
        robust = TRUE
      )
      
      # peak response
      peak_response_draws <- log(peak_draws) - log(baseline) 
      peak_response_summary <- posterior_summary(
        peak_response_draws,
        probs = ci_probs,
        robust = TRUE
      )
      
      # decay rate
      ## check multiple time points
      decay_check = c(60, 90, 180, 365)
      decay_list <- list()
      decay_draws_list <- list()
      
      # decay check loop
      for (t in decay_check) {
        
        ## get indexes for decay_check val
        index1 <- which(timestamps == t)
        index2 <- which(timestamps == (t+1))
        
        ## get estimates
        decay_val1 <- epred_log[,index1]
        decay_val2 <- epred_log[,index2]
        
        ## get decay rate
        decay_rate <- 
          (decay_val1 - decay_val2) / 
          abs((timestamps[index1] - timestamps[index2]))
        
        decay_summary <- posterior_summary(
          decay_rate,
          probs = ci_probs,
          robust = TRUE
        )
        
        decay_draws_list[[paste0("rate_D", t)]] <- decay_rate
        decay_list[[paste0("rate_D", t)]] <- decay_summary
        
      }
      
      # response 1 year
      ## get estimate at 1 year
      index_1year <- which(timestamps == 365)
      titer_1year_draws <- epred[,index_1year]
      titer_1year_summary <- posterior_summary(
        titer_1year_draws,
        probs = ci_probs,
        robust = TRUE
      )
      
      response_1year_draws <- log(titer_1year_draws) - log(baseline) 
      response_1year_summary <- posterior_summary(
        response_1year_draws,
        probs = ci_probs,
        robust = TRUE
      )
      
      # store all measures
      measures[[model_name]][[paste0(dose_current)]] <- list(
        min_peak_index = min(peak_index),
        min_growth = min(growth_draws),
        peak = peak_summary,
        peak_time = peak_time_summary,
        growth_rate = growth_summary,
        peak_response = peak_response_summary,
        baseline = baseline_summary,
        decay_rate = decay_list,
        titer_1year = titer_1year_summary,
        response_1year = response_1year_summary,
        
        peak_draws = peak_draws,
        growth_draws = growth_draws,
        peak_time_draws = peak_time_draws,
        peak_response_draws = peak_response_draws,
        decay_rate_draws = decay_draws_list,
        titer_1year_draws = titer_1year_draws,
        response_1year_draws = response_1year_draws
      )
      
    }
    
  }
  
  return(measures)
  
}

func_process_measures <- function(measures, use_id=NULL, 
                                  dat_decay_constants=NULL) {
  
  process_measures <- list()
  
  # helper function
  .process_summary <- function(x) {
    tibble::tibble(
      est = x[1],
      lo  = x[3],
      up  = x[4]
    )
  }
  
  # model loop
  for (model_name in names(measures)) {
    
    # get study, model function, antibody
    study <- stringr::str_split(model_name, "_", simplify = TRUE)[, 1]
    model_func <- stringr::str_split(model_name, "_", simplify = TRUE)[, 2]
    antibody <- paste(
      stringr::str_split(model_name, "_", simplify = TRUE)[, 3:4],collapse = "_"
    )
    antibody_clean <- stringr::str_replace(antibody,"_"," ")
    dose_list <- measures[[model_name]]
    
    # dose loop
    for (dose_val in names(dose_list)) {
      
      dose_data <- dose_list[[dose_val]]
      
      # group-level
      if (use_id == FALSE) {
        
        dat_new <- tibble(
          study = study,
          model_func = model_func,
          antibody = antibody_clean,
          dose = dose_val
        )
        
        scalar_measures <- c("peak","peak_time","growth_rate","peak_response",
                             "baseline","titer_1year","response_1year")
        
        for (m in scalar_measures) {
          
          summary <- .process_summary(dose_data[[m]])
          
          process_measures[[length(process_measures)+1]] <- bind_cols(
            dat_new, tibble(measure=m), summary
          )
        }
        
        # decay processing (multiple days)
        for (decay_days in names(dose_data$decay_rate)) {
          
          summary <- .process_summary(dose_data$decay_rate[[decay_days]])
          
          process_measures[[length(process_measures)+1]] <- bind_cols(
            dat_new, tibble(measure = paste0("decay_", decay_days)),
            summary
          )
          
        }
        
        # id-level
      } else{
        
        for (id_i in names(dose_data)) {
          
          id_data <- dose_data[[id_i]]
          
          dat_new <- tibble(
            study = study,
            model_func = model_func,
            antibody = antibody_clean,
            dose = dose_val,
            id = id_i
          )
          
          scalar_measures <- c("peak","peak_time","growth_rate","peak_response",
                               "baseline","titer_1year","response_1year")
          
          for (m in scalar_measures) {
            
            summary <- .process_summary(id_data[[m]])
            
            process_measures[[length(process_measures)+1]] <- bind_cols(
              dat_new, tibble(measure=m), summary
            )
          }
          
          # decay processing (multiple days)
          for (decay_days in names(id_data$decay_rate)) {
            
            summary <- .process_summary(id_data$decay_rate[[decay_days]])
            
            process_measures[[length(process_measures)+1]] <- bind_cols(
              dat_new, tibble(measure = paste0("decay_", decay_days)),
              summary
            )
          }
        }
      }
    }
  }
  
  dat_all <- bind_rows(process_measures) %>%
    mutate(
      dose = factor(dose, levels=c("4.8","48","4800",
                                   "5","15","50","150"))
    )
  
  return(dat_all)
  
}

func_process_draws <- function(measures, use_id) {
  
  dose_levels <- c("4.8","48","4800","5","15","50","150")
  base_measure_cols <- c("peak", "peak_response", "growth_rate", "peak_time", 
                         "titer_1year", "response_1year")
  
  process_draws <- list()
  
  for (model_name in names(measures)) {
    
    parts <- stringr::str_split(model_name, "_", simplify = TRUE)
    study <- parts[,1]
    model_func <- parts[,2]
    antibody <- paste(parts[,3:4], collapse = "_")
    antibody_clean <- stringr::str_replace(antibody, "_", " ")
    
    dose_list <- measures[[model_name]]
    
    for (dose_val in names(dose_list)) {
      
      dose_data <- dose_list[[dose_val]]
      
      make_draws <- function(draw_obj, id_val = NULL) {
        n_draws <- length(draw_obj$peak_response_draws)
        
        dat_new <- tibble::tibble(
          study = study,
          model_func = model_func,
          antibody = antibody_clean,
          dose = dose_val,
          id = id_val,
          draw = seq_len(n_draws),
          peak = draw_obj$peak_draws,
          growth_rate = draw_obj$growth_draws,
          peak_time = draw_obj$peak_time_draws,
          peak_response = draw_obj$peak_response_draws,
          titer_1year = draw_obj$titer_1year_draws,
          response_1year = draw_obj$response_1year_draws
        )
        
        for (decay_day in names(draw_obj$decay_rate_draws)) {
          col_name <- paste0("decay_", decay_day)
          decay_vec <- draw_obj$decay_rate_draws[[decay_day]]
          dat_new[[col_name]] <- decay_vec
        }
        
        return(dat_new)
        
      }
      
      if (!use_id) {
        # group-level
        process_draws[[length(process_draws) + 1]] <- make_draws(dose_data)
      } else {
        # id-level
        for (id_i in names(dose_data)) {
          process_draws[[length(process_draws) + 1]] <- make_draws(dose_data[[id_i]], id_i)
        }
      }
    }
  }
  
  dat_all <- dplyr::bind_rows(process_draws) %>%
    dplyr::mutate(dose = factor(dose, levels = dose_levels))
  
  # filter for 500 draws per dose, remove the lowest/highest 5% to control for outliers
  set.seed(333)
  
  measure_cols <- intersect(
    names(dat_all),
    c(base_measure_cols, grep("^decay_rate_D", names(dat_all), value = TRUE))
  )
  
  dat_random <- dat_all %>%
    group_by(model_func, antibody, dose) %>%
    filter(
      if_all(all_of(measure_cols), ~ . >= quantile(., 0.05, na.rm = TRUE) &
               . <= quantile(., 0.95, na.rm = TRUE))
    ) %>%
    filter(draw %in% sample(unique(draw), size = min(500, n_distinct(draw)))) %>%
    ungroup()
  
  list(
    draws_random=dat_random,
    draws_all=dat_all
  )
}

func_pairwise_compare <- function(draws) {
  
  measures <- c("peak", "growth_rate", "peak_response", "peak_time", 
                "titer_1year", "response_1year", "decay_rate_D90")
  
  results <- draws %>%
    tidyr::pivot_longer(
      cols = all_of(measures), 
      names_to = "measure", 
      values_to = "value"
    ) %>%
    group_by(study, model_func, antibody, measure) %>%
    group_modify(~ {
      doses <- unique(.x$dose)
      pairs <- tidyr::expand_grid(dose_a = doses, dose_b = doses) %>%
        filter(dose_a != dose_b)
      
      purrr::map_dfr(seq_len(nrow(pairs)), function(i) {
        paired <- inner_join(
          .x %>% filter(dose == pairs$dose_a[i], !is.na(value)) %>% select(draw, value),
          .x %>% filter(dose == pairs$dose_b[i], !is.na(value)) %>% select(draw, value),
          by = "draw"
        )
        
        tibble(
          dose_a = pairs$dose_a[i],
          dose_b = pairs$dose_b[i],
          comparison = paste0("P(", pairs$dose_a[i], " > ", pairs$dose_b[i], ")"),
          prob = mean(paired$value.x > paired$value.y)
        )
      })
    }) %>%
    ungroup()
  
  results_filter <- results %>%
    filter(
      (study == "NI" & comparison %in% c("P(4.8 > 48)", "P(4.8 > 4800)", 
                                         "P(48 > 4800)")) |
        (study != "NI" & comparison %in% c("P(5 > 15)","P(5 > 50)","P(5 > 150)",
                                           "P(15 > 50)","P(15 > 150)",
                                           "P(50 > 150)"))
    )
  
  return(results_filter)
  
}

func_get_correlations <- function(draws_ni, draws_nv, by) {
  
  dat_bind <- bind_rows(draws_ni, draws_nv)
  
  # overall
  dat_overall <- dat_bind %>%
    group_by(antibody, model_func) %>%
    summarise(
      ## relative
      # compute Spearman correlations with response_1year
      cor_resp_growth = cor(response_1year, growth_rate, method = "spearman", use = "pairwise.complete.obs"),
      cor_resp_peakresp = cor(response_1year, peak_response, method = "spearman"),
      cor_resp_D90 = cor(response_1year, decay_rate_D90, method = "spearman"),
      
      # peak_response and decay
      cor_peakresp_D90 = cor(peak_response, decay_rate_D90, method = "spearman"),
      # peak_response and growth
      cor_peakresp_growth = cor(peak_response, growth_rate, method = "spearman", use = "pairwise.complete.obs"),
      
      ## absolute
      # compute Spearman correlations with titer_1year
      cor_titer_growth = cor(titer_1year, growth_rate, method = "spearman", use = "pairwise.complete.obs"),
      cor_titer_peak = cor(titer_1year, peak, method = "spearman"),
      cor_titer_D90 = cor(titer_1year, decay_rate_D90, method = "spearman"),
      
      # peak and decay
      cor_peak_D90 = cor(peak, decay_rate_D90, method = "spearman"),
      # peak and growth
      cor_peak_growth = cor(peak, growth_rate, method = "spearman", use = "pairwise.complete.obs"),
      .groups = "drop"
    )
  
  # by dose
  dat_dose <- dat_bind %>%
    group_by(antibody, model_func, study, dose) %>%
    summarise(
      ## relative
      # compute Spearman correlations with response_1year
      cor_resp_growth = cor(response_1year, growth_rate, method = "spearman", use = "pairwise.complete.obs"),
      cor_resp_peakresp = cor(response_1year, peak_response, method = "spearman"),
      cor_resp_D90 = cor(response_1year, decay_rate_D90, method = "spearman"),
      
      # peak_response and decay
      cor_peakresp_D90 = cor(peak_response, decay_rate_D90, method = "spearman"),
      # peak_response and growth
      cor_peakresp_growth = cor(peak_response, growth_rate, method = "spearman", use = "pairwise.complete.obs"),
      
      ## absolute
      # compute Spearman correlations with titer_1year
      cor_titer_growth = cor(titer_1year, growth_rate, method = "spearman", use = "pairwise.complete.obs"),
      cor_titer_peak = cor(titer_1year, peak, method = "spearman"),
      cor_titer_D90 = cor(titer_1year, decay_rate_D90, method = "spearman"),
      
      # peak and decay
      cor_peak_D90 = cor(peak, decay_rate_D90, method = "spearman"),
      # peak and growth
      cor_peak_growth = cor(peak, growth_rate, method = "spearman", use = "pairwise.complete.obs"),
      .groups = "drop"
    )
  
  if (by == "overall") {
    
    return(dat_overall)
  
    } else if (by == "dose") {
    
    return(dat_dose)
    
  }
}

# END OF SCRIPT ====