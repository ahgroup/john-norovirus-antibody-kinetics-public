###
# 01-model-fits.R
###

func_setup_model <- function() {
  
  library(brms)
  
  total_cores <- parallelly::availableCores(omit = 2)
  
  # brms settings
  set_brms <- list(
    myseed=333,
    mywarmup=1000,
    myiter=3000,
    mychains=4,
    mycores=max(1, floor(total_cores / 4)),
    my_adapt_delta=0.99,
    my_max_tree=15
  )
  
  # priors
  priors_exp <- c(
    brms::prior(student_t(3, 0, 2), class = "sigma", lb = 0),
    brms::prior(constant(4), class = "nu"),
    
    brms::prior(normal(7, 2), class = "b", nlpar = "p"),  # fixed antibody peak
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "p", lb = 0),  # sd
    
    brms::prior(normal(0, 1.5), class = "b", nlpar = "g"),  # fixed growth rate
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "g", lb = 0),  # sd
    
    brms::prior(normal(log(10), 1), class = "b", nlpar = "k"),  # time of half peak
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "k", lb = 0),  # sd
    
    brms::prior(normal(-4, 1.5), class = "b", nlpar = "d"),  # fixed waning parameter
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "d", lb = 0)  # sd
  )
  
  priors_pow <- c(
    brms::prior(student_t(3, 0, 2), class = "sigma", lb = 0),
    brms::prior(constant(4), class = "nu"),
    
    brms::prior(normal(7, 2), class = "b", nlpar = "p"),  # fixed antibody peak
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "p", lb = 0),  # sd
    
    brms::prior(normal(0, 1.5), class = "b", nlpar = "g"),  # fixed growth rate
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "g", lb = 0),  # sd
    
    brms::prior(normal(log(10), 1), class = "b", nlpar = "k"),  # half-peak time
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "k", lb = 0),  # sd
    
    brms::prior(normal(-1, 1.5), class = "b", nlpar = "d"),  # fixed waning parameter
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "d", lb = 0)  # sd
  )
  
  priors <- list(
    exp = priors_exp,
    pow = priors_pow
  )
  
  # equations
  ## exponential decay model
  eqn_exp <- brms::bf(
    logy | cens(cens, logy2) ~ log((exp(p)/(1+exp(-exp(g)*(day-exp(k))))) * exp(-exp(d)*day)),
    p ~ 0 + dose + (1|id),
    g ~ 0 + dose + (1|id),
    k ~ 0 + dose + (1|id),
    d ~ 0 + dose + (1|id),
    nl = TRUE)
  
  ## power-law decay model
  eqn_pow <- brms::bf(
    logy | cens(cens, logy2) ~ log((exp(p)/(1+exp(-exp(g)*((day+1)-exp(k))))) * (day+1)^(-exp(d))),
    p ~ 0 + dose + (1|id),
    g ~ 0 + dose + (1|id),
    k ~ 0 + dose + (1|id),
    d ~ 0 + dose + (1|id),
    nl = TRUE)
  
  equations <- list("exp" = eqn_exp, "pow" = eqn_pow)
  
  list(
    set_brms=set_brms,
    priors=priors,
    equations=equations
  )
}

func_fit_models <- function(data, setup) {
  
  library(brms)
  
  models <- list()
  
  antibodies <- unique(data$antibody)
  study_id   <- unique(data$study)
  
  # fit models loop
  for (k in seq_along(setup$equations)) {
    
    equation <- names(setup$equations)[k]
    
    priors_use <- setup$priors[[equation]]
    
    for (i in seq_along(antibodies)) {
      
      dat_model <- data %>% 
        filter(
          antibody == antibodies[i]
        )
      
      model_name <- paste0(
        study_id, "_",
        equation, "_",
        antibodies[i]
      )
      
      models[[model_name]] <- brms::brm(
        formula = setup$equations[[k]],
        data = dat_model,
        family = student(),
        prior = priors_use,
        init = 0,
        backend = "cmdstanr",
        silent = 2,
        refresh = setup$set_brms$myiter/5,
        warmup = setup$set_brms$mywarmup,
        iter = setup$set_brms$myiter,
        chains = setup$set_brms$mychains,
        cores = setup$set_brms$mycores,
        seed = setup$set_brms$myseed,
        save_pars = save_pars(all = TRUE),
        control = list(
          adapt_delta = setup$set_brms$my_adapt_delta,
          max_treedepth = setup$set_brms$my_max_tree
        )
      )
      
    }
  }
  
  return(models)
  
}

func_get_model_metrics <- function(models) {
  
  library(brms)
  
  # add loo 
  model_metrics <- lapply(models, brms::add_criterion, "loo")
  
  # get study, model function, and antibody
  info <- tibble::tibble(
    model_name = names(model_metrics)
  ) %>%
    mutate(
      study = stringr::str_extract(model_name, "^[^_]+"),
      model_func  = stringr::str_extract(model_name, "(exp|pow)"),
      antibody = stringr::str_remove(model_name, "^[^_]+_(exp|pow)_")
    )
  
  # helper to extract pointwise elpd
  get_pointwise <- function(mod) {
    mod$criteria$loo$pointwise[, "elpd_loo"]
  }
  
  metrics_long <- purrr::map2_dfr(model_metrics, names(model_metrics), function(mod, name) {
    tibble::tibble(
      model_name = name,
      elpd_pointwise = list(get_pointwise(mod)),
      elpd_loo = mod$criteria$loo$estimates["elpd_loo", "Estimate"],
      se_loo   = mod$criteria$loo$estimates["elpd_loo", "SE"]
    )
  }) %>%
    left_join(info, by = "model_name")
  
  metrics_wide <- metrics_long %>%
    select(study, antibody, model_func, elpd_pointwise, elpd_loo) %>%
    tidyr::pivot_wider(
      names_from = model_func,
      values_from = c(elpd_pointwise, elpd_loo)
    )
  
  metrics_compare <- metrics_wide %>%
    mutate(
      # difference in total ELPD
      delta_elpd = elpd_loo_pow - elpd_loo_exp,
      
      # pointwise differences
      diff_i = purrr::map2(elpd_pointwise_pow, elpd_pointwise_exp, ~ .x - .y),
      
      # SE of the difference
      se_diff = purrr::map_dbl(diff_i, ~ sqrt(length(.x) * var(.x, na.rm = TRUE))),
      
      # probability that power-law is better than exponential
      p_pow_better = 1 - pnorm(0, mean = delta_elpd, sd = se_diff),
      p_adj_pow_better = 1 - pnorm(0, mean = delta_elpd, sd = 2*se_diff)
    ) %>%
    ungroup() %>%
    select(study, antibody, delta_elpd, se_diff, p_pow_better, p_adj_pow_better)
  
  return(metrics_compare)
  
}

func_get_prior_posterior_draws <- function(models) {
  
  library(brms)
  
  all_draws <- list()
  
  for (m in names(models)) {
    
    fit <- models[[m]]
    
    # posterior draws
    draws_posterior <- as_draws_df(fit) %>%
      select(starts_with(c("b_"))) %>%
      tidyr::pivot_longer(
        everything(),
        names_to = "parameter",
        values_to = "value"
      ) %>%
      mutate(
        type = "Posterior",
        model = m
      )
    
    # prior draws
    fit_prior <- update(
      fit,
      sample_prior = "only",
      refresh = 0
    )
    
    draws_prior <- as_draws_df(fit_prior) %>%
      select(starts_with(c("b_"))) %>%
      tidyr::pivot_longer(
        everything(),
        names_to = "parameter",
        values_to = "value"
      ) %>%
      mutate(
        type = "Prior",
        model =  m
      )
    
    all_draws[[m]] <- bind_rows(draws_prior, draws_posterior)
    
  }
  
  draws <- bind_rows(all_draws)
  return(draws)
  
}

# sensitivity analysis ====
# model setup function for prior sensitivity analysis in supplement
# doubles the standard deviations
func_setup_model_prior_sens <- function() {
  
  library(brms)
  
  total_cores <- parallelly::availableCores(omit = 2)
  
  # brms settings
  set_brms <- list(
    myseed=333,
    mywarmup=1000,
    myiter=3000,
    mychains=4,
    mycores=max(1, floor(total_cores / 4)),
    my_adapt_delta=0.99,
    my_max_tree=15
  )
  
  # priors
  priors_exp <- c(
    brms::prior(student_t(3, 0, 2), class = "sigma", lb = 0),
    brms::prior(constant(4), class = "nu"),
    
    brms::prior(normal(7, 4), class = "b", nlpar = "p"),  # fixed antibody peak
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "p", lb = 0),  # sd
    
    brms::prior(normal(0, 3), class = "b", nlpar = "g"),  # fixed growth rate
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "g", lb = 0),  # sd
    
    brms::prior(normal(log(30), 2), class = "b", nlpar = "k"),  # time of half peak
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "k", lb = 0),  # sd
    
    brms::prior(normal(-4, 3), class = "b", nlpar = "d"),  # fixed waning parameter
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "d", lb = 0)  # sd
  )
  
  priors_pow <- c(
    brms::prior(student_t(3, 0, 2), class = "sigma", lb = 0),
    brms::prior(constant(4), class = "nu"),
    
    brms::prior(normal(7, 4), class = "b", nlpar = "p"),  # fixed antibody peak
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "p", lb = 0),  # sd
    
    brms::prior(normal(0, 3), class = "b", nlpar = "g"),  # fixed growth rate
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "g", lb = 0),  # sd
    
    brms::prior(normal(log(30), 2), class = "b", nlpar = "k"),  # half-peak time
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "k", lb = 0),  # sd
    
    brms::prior(normal(-1, 3), class = "b", nlpar = "d"),  # fixed waning parameter
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "d", lb = 0)  # sd
  )
  
  priors <- list(
    exp = priors_exp,
    pow = priors_pow
  )
  
  # equations
  ## exponential decay model
  eqn_exp <- brms::bf(
    logy | cens(cens, logy2) ~ log((exp(p)/(1+exp(-exp(g)*(day-exp(k))))) * exp(-exp(d)*day)),
    p ~ 0 + dose + (1|id),
    g ~ 0 + dose + (1|id),
    k ~ 0 + dose + (1|id),
    d ~ 0 + dose + (1|id),
    nl = TRUE)
  
  ## power-law decay model
  eqn_pow <- brms::bf(
    logy | cens(cens, logy2) ~ log((exp(p)/(1+exp(-exp(g)*((day+1)-exp(k))))) * (day+1)^(-exp(d))),
    p ~ 0 + dose + (1|id),
    g ~ 0 + dose + (1|id),
    k ~ 0 + dose + (1|id),
    d ~ 0 + dose + (1|id),
    nl = TRUE)
  
  equations <- list("exp" = eqn_exp, "pow" = eqn_pow)
  
  list(
    set_brms=set_brms,
    priors=priors,
    equations=equations
  )
}

# flexible power-law model
# allows slower long-term waning after the initial post-peak decline
func_setup_model_flexible <- function() {
  
  library(brms)
  
  total_cores <- parallelly::availableCores(omit = 2)
  
  # brms settings
  set_brms <- list(
    myseed=333,
    mywarmup=1000,
    myiter=3000,
    mychains=4,
    mycores=max(1, floor(total_cores / 4)),
    my_adapt_delta=0.99,
    my_max_tree=15
  )
  
  priors_pow <- c(
    brms::prior(student_t(3, 0, 2), class = "sigma", lb = 0),
    brms::prior(constant(4), class = "nu"),
    
    brms::prior(normal(7, 2), class = "b", nlpar = "p"),  # fixed antibody peak
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "p", lb = 0),  # sd
    
    brms::prior(normal(0, 1.5), class = "b", nlpar = "g"),  # fixed growth rate
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "g", lb = 0),  # sd
    
    brms::prior(normal(log(10), 1), class = "b", nlpar = "k"),  # time of half peak
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "k", lb = 0),  # sd
    
    brms::prior(normal(-1, 1.5), class = "b", nlpar = "d"),  # fixed waning parameter
    brms::prior(student_t(3, 0, 2), class = "sd", nlpar = "d", lb = 0),  # sd
    
    brms::prior(normal(0, 0.5), class = "b", nlpar = "lambda")  # decay slowing parameter
    
  )
  
  priors <- list(
    pow = priors_pow
  )
  
  # equations
  
  ## flexible power-law decay model
  eqn_pow <- brms::bf(
    logy | cens(cens, logy2) ~ 
      log(
        (exp(p) / (1+exp(-exp(g)*((day+1)-exp(k))))) * 
          (1+exp(lambda)*day)^(-exp(d) / exp(lambda))
      ),
    p ~ 0 + dose + (1|id),
    g ~ 0 + dose + (1|id),
    k ~ 0 + dose + (1|id),
    d ~ 0 + dose + (1|id),
    lambda ~ 0 + dose,
    nl = TRUE)
  
  equations <- list("pow" = eqn_pow)
  
  list(
    set_brms=set_brms,
    priors=priors,
    equations=equations
  )
}

func_get_model_metrics_flexible <- function(models_original, models_flexible) {
  browser()
  library(brms)
  
  models_original <- models_original[
    !grepl("_exp_", names(models_original))
  ]
  
  info_original <- tibble::tibble(
    model_name = names(models_original),
    model_type = "original"
  )
  
  info_flexible <- tibble::tibble(
    model_name = names(models_flexible),
    model_type = "flexible"
  )
  
  info <- bind_rows(
    info_original,
    info_flexible
  ) %>%
    mutate(
      study = stringr::str_extract(model_name, "^[^_]+"),
      antibody = stringr::str_remove(model_name, "^[^_]+_(pow|flex)_")
    )
  
  models <- c(models_original, models_flexible)
  
  model_id <- paste0(
    rep(c("original", "flexible"),
        c(length(models_original), length(models_flexible))),
    "_",
    names(models)
  )
  
  names(models) <- model_id
  
  # add loo 
  model_metrics <- lapply(models, brms::add_criterion, "loo")
  
  # helper to extract pointwise elpd
  get_pointwise <- function(mod) {
    mod$criteria$loo$pointwise[, "elpd_loo"]
  }
  
  metrics_long <- purrr::map2_dfr(
    model_metrics, 
    names(model_metrics), 
    function(mod, name) {
      tibble::tibble(
        model_name = name,
        elpd_pointwise = list(get_pointwise(mod)),
        elpd_loo = mod$criteria$loo$estimates["elpd_loo", "Estimate"],
        se_loo   = mod$criteria$loo$estimates["elpd_loo", "SE"]
      )
    }) %>%
    left_join(
      info %>%
        mutate(
          model_id = paste0(model_type, "_", model_name)
        ), 
      by = c("model_name" = "model_id")
    )
  
  metrics_wide <- metrics_long %>%
    select(study, antibody, model_type, elpd_pointwise, elpd_loo) %>%
    tidyr::pivot_wider(
      names_from = model_type,
      values_from = c(elpd_pointwise, elpd_loo)
    )
  
  metrics_compare <- metrics_wide %>%
    mutate(
      # difference in total ELPD
      delta_elpd = elpd_loo_flexible - elpd_loo_original,
      
      # pointwise differences
      diff_i = purrr::map2(
        elpd_pointwise_flexible, 
        elpd_pointwise_original, 
        ~ .x - .y
        ),
      
      # SE of the difference
      se_diff = purrr::map_dbl(diff_i, ~ sqrt(length(.x) * var(.x, na.rm = TRUE))),
      
      # probability that power-law is better than exponential
      p_flexible_better = 1 - pnorm(0, mean = delta_elpd, sd = se_diff),
      p_adj_flexible_better = 1 - pnorm(0, mean = delta_elpd, sd = 2*se_diff)
    ) %>%
    ungroup() %>%
    select(study, antibody, delta_elpd, se_diff, p_flexible_better, p_adj_flexible_better)
  
  return(metrics_compare)
  
}

# END OF SCRIPT ====