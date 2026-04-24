###
# 02-model-predictions.R
###

func_get_preds <- function(data, models, use_id=NULL) {
  
  library(brms)
  
  predictions <- list()
  
  # define sequence of time points for predictions
  pred_seq <- unique(c(seq(0,20,1/6), seq(20,400,1), 1826)) # 5 years
  
  # get antibodies
  antibodies <- unique(data$antibody)
  
  # antibody loop
  for (antibody_current in antibodies) {
    
    # get current antibody
    dat_pred <- data %>%
      filter(
        antibody == antibody_current
      )
    
    # get doses in current model
    doses <- unique(dat_pred$dose)
    
    # get both models for current antibody
    model_indices <- grep(antibody_current, names(models))
    
    # model (exp/pow) loop
    for (model_index in model_indices) {
      
      model <- models[[model_index]]
      model_name <- names(models)[model_index]
      
      # dose loop
      for (dose_current in doses) {
        
        # group-level predictions
        if (use_id == FALSE) {
          
          dat_new <- tidyr::crossing(
            study = unique(data$study),
            model = stringr::str_split(model_name, "_")[[1]][2],
            antibody = antibody_current,
            antibody_clean = stringr::str_replace(antibody_current,"_"," "),
            day = pred_seq,
            dose = dose_current
          ) %>%
            mutate(dose = factor(dose, levels = levels(dat_pred$dose)))
          
          # get posterior point predictions
          dat_epred <- posterior_epred(
            object = model,
            newdata = dat_new,
            re_formula = NA # marginalization step
          ) 
          
          # get CIs
          dat_post <- dat_epred %>%
            posterior_summary(
              probs = c(0.025, 0.975),
              robust = TRUE
            ) # robust = TRUE assigns median and mad
          
          # process preds
          dat_pred_process <- tibble::as_tibble(dat_post)
          colnames(dat_pred_process) <- c("logy","se","log_lo","log_up")
          
          # merge with data
          dat_all <- bind_cols(dat_new, dat_pred_process)
          
          predictions[[model_name]][[as.character(dose_current)]] <- list(
            type = "group",
            epred = dat_epred,
            predictions = dat_all
          )
          
        } else {
          
          # id-level predictions
          ids <- dat_pred %>%
            filter(dose == dose_current) %>%
            pull(id) %>%
            unique()
          
          for (id_i in ids) {
            
            dat_new <- tidyr::crossing(
              study = unique(data$study),
              model = stringr::str_split(model_name, "_")[[1]][2],
              antibody = antibody_current,
              antibody_clean = stringr::str_replace(antibody_current,"_"," "),
              day = pred_seq,
            ) %>%
              mutate(
                dose = factor(dose_current, levels = levels(dat_pred$dose)),
                id = id_i
              )
            
            dat_epred <- posterior_epred(
              object = model,
              newdata = dat_new,
              re_formula = NULL # includes id random effects
            )
            
            # get CIs
            dat_post <- posterior_summary(
              dat_epred,
              probs = c(0.025, 0.975),
              robust = TRUE
            )
            
            dat_pred_process <- tibble::as_tibble(dat_post)
            colnames(dat_pred_process) <- c("logy","se","log_lo","log_up")
            
            dat_all <- bind_cols(dat_new, dat_pred_process)
            
            predictions[[model_name]][[as.character(dose_current)]][[id_i]] <- list(
              type = "id",
              epred = dat_epred,
              predictions = dat_all
            )
          }
        }
      }
    }
  }
  
  return(predictions)
  
}

func_process_preds <- function(predictions, use_id=NULL) {
  
  process_data <- list()
  
  # model loop
  for (model_name in names(predictions)) {
    
    dose_list <- predictions[[model_name]]
    
    # dose loop
    for (dose_val in names(dose_list)) {
      
      dose_data <- dose_list[[dose_val]]
      
      # group-level
      if (use_id == FALSE) {
        
        dat_process <- dose_data$predictions
        
        process_data[[paste(model_name, dose_val, sep = "_")]] <- dat_process %>%
          mutate(
            y = exp(logy),
            lo = exp(log_lo),
            up = exp(log_up)
          )
        
      } else{
        
        for (id in names(dose_data)) {
          
          id_data <- dose_data[[id]]
          
          dat_process_id <- id_data$predictions
          
          process_data[[paste(model_name, dose_val, id, sep = "_")]] <- dat_process_id %>%
            mutate(
              y = exp(logy),
              lo = exp(log_lo),
              up = exp(log_up)
            )
          
        }
      }
    }
  }
  
  dat_all <- bind_rows(process_data)
  
  return(dat_all)
  
}

func_get_residuals <- function(dat_pred, dat_obs) {
  
  # filter predictions to reduce computational cost
  days <- unique(dat_obs$day)
  dat_pred_red <- dat_pred %>%
    filter(day %in% days) %>%
    select(id, model_func=model, antibody, day, dose, y_pred=y)
  
  # remove censored obs
  dat_obs_red <- dat_obs %>%
    filter(cens == 0) %>%
    select(id, study, dose, day, antibody, antibody_clean, y_obs=y)
  
  residuals <- left_join(
    dat_obs_red, dat_pred_red, by=c("id","day","antibody","dose")
  ) %>%
    mutate(
      residual = y_obs - y_pred,
      model = paste(
        stringr::str_to_upper(study), model_func, antibody, 
        sep="_"
      )
    )
  
  return(residuals)
  
}

# END OF SCRIPT ====