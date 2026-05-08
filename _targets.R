###
# _targets.R
#
# Reproducible pipeline for:
# - Data processing (NI and NV studies)
# - Bayesian model fitting (power-law and exponential)
# - Posterior prediction (group and individual levels)
# - Derivation of peak, decay, and durability measures
# - Figure and table generation (main + supplement)
# - Sensitivity analyses (alternative censoring bounds)
#
# To reproduce entire analysis, run targets::tar_make()
# To reproduce the manuscript pipeline, run targets::tar_make("reproduce_manuscript")
###

###
# Setup ====
###

# load packages
library(targets)
library(brms)

# number of cores to use for parallel processing
n_local_cores <- parallelly::availableCores(
  omit=2, 
  constraints="connections"
)

# define local parallel controller, store a log for each worker in /logs/
controller_local <- crew::crew_controller_local(
  name="local",
  workers=n_local_cores,
  options_local=crew::crew_options_local(log_directory="logs")
)

# set global target options
tar_option_set(
  format="rds",
  controller=controller_local,
  resources=targets::tar_resources(
    crew=tar_resources_crew(controller="local")
  ),
  packages=c("dplyr"),
  seed=333
)

# source all functions in /R/
tar_source(
  files=list.files(
    "R", 
    pattern="\\.R$", 
    full.names=TRUE, 
    recursive=TRUE
  )
)

###
# Targets pipeline ====
###

list(
  
  ###
  ## Run manuscript results only ====
  # reproduces only the results in the main text
  # in console, run targets::tar_make("reproduce_manuscript")
  ###
  
  tar_target(
    reproduce_manuscript,
    list(
      tab_model_metrics,
      plot_ni_group_model,
      plot_nv_group_model,
      plot_ni_pow_measures,
      plot_nv_pow_measures,
      plot_measures_vs_1year_draws,
      plot_pow_comb_peak_decay_vs
    )
  ),
  
  ###
  ## Run all results ====
  # reproduces all results in the main text and supplementary materials
  # in console, run targets::tar_make()
  ###
  
  ###
  ### 01 read processed data ====
  ###
  
  tar_target(
    process_ni_data,
    func_process_ni_data(
      data_path = here::here("data/dat_ni.rds")
    )
  ),
  
  tar_target(
    process_nv_data,
    func_process_nv_data(
      data_path = here::here("data/dat_nv.rds")
    )
  ),
  
  ###
  ### 02 model fitting ====
  # fit Bayesian models to NI and NV data
  # compare exponential vs power-law model performance metrics
  # get prior and posterior draws
  ###
  
  tar_target(
    setup_model,
    func_setup_model()
  ),
  
  ###
  # fit models
  ###
  
  tar_target(
    fit_ni_models,
    func_fit_models(
      data=process_ni_data,
      setup=setup_model
    )
  ),
  
  tar_target(
    fit_nv_models,
    func_fit_models(
      data=process_nv_data,
      setup=setup_model
    )
  ),
  
  ###
  # get model metrics
  ###
  
  tar_target(
    get_ni_model_metrics,
    func_get_model_metrics(
      models=fit_ni_models
    )
  ),
  
  tar_target(
    get_nv_model_metrics,
    func_get_model_metrics(
      models=fit_nv_models
    )
  ),
  
  ###
  # get prior and posterior draws
  ###
  
  tar_target(
    get_ni_prior_posterior_draws,
    func_get_prior_posterior_draws(
      models=fit_ni_models
    )
  ),
  
  tar_target(
    get_nv_prior_posterior_draws,
    func_get_prior_posterior_draws(
      models=fit_nv_models
    )
  ),
  
  ###
  ### 03 posterior predictions ====
  # generate predicted trajectories at population- and individual-level
  # process into plotting-ready format
  # get residuals
  ###
  
  ###
  # population-level predictions
  ###
  
  tar_target(
    get_ni_group_preds,
    func_get_preds(
      data=process_ni_data,
      models=fit_ni_models,
      use_id=FALSE
    )
  ),
  
  tar_target(
    get_nv_group_preds,
    func_get_preds(
      data=process_nv_data,
      models=fit_nv_models,
      use_id=FALSE
    )
  ),
  
  ###
  # individual-level predictions
  ###
  
  tar_target(
    get_ni_id_preds,
    func_get_preds(
      data=process_ni_data,
      models=fit_ni_models,
      use_id=TRUE
    )
  ),
  
  tar_target(
    get_nv_id_preds,
    func_get_preds(
      data=process_nv_data,
      models=fit_nv_models,
      use_id=TRUE
    )
  ),
  
  ###
  # process predictions
  ###
  
  tar_target(
    process_ni_group_preds,
    func_process_preds(
      predictions=get_ni_group_preds,
      use_id=FALSE
    )
  ),
  
  tar_target(
    process_nv_group_preds,
    func_process_preds(
      predictions=get_nv_group_preds,
      use_id=FALSE
    )
  ),
  
  tar_target(
    process_ni_id_preds,
    func_process_preds(
      predictions=get_ni_id_preds,
      use_id=TRUE
    )
  ),
  
  tar_target(
    process_nv_id_preds,
    func_process_preds(
      predictions=get_nv_id_preds,
      use_id=TRUE
    )
  ),
  
  ###
  # get individual-level residuals
  ###
  
  tar_target(
    get_ni_residuals,
    func_get_residuals(
      dat_pred=process_ni_id_preds,
      dat_obs=process_ni_data
    )
  ),
  
  tar_target(
    get_nv_residuals,
    func_get_residuals(
      dat_pred=process_nv_id_preds,
      dat_obs=process_nv_data
    )
  ),
  
  ###
  ### 04 measures and draws ====
  # from posterior trajectories, get group- and individual-level measures and draws
  # compute pairwise dose probabilities
  # compute measures correlations
  ###
  
  ###
  # population-level measures
  ###
  
  tar_target(
    get_ni_group_measures,
    func_get_measures_group(
      data=process_ni_data,
      models=fit_ni_models,
      preds=get_ni_group_preds
    )
  ),
  
  tar_target(
    get_nv_group_measures,
    func_get_measures_group(
      data=process_nv_data,
      models=fit_nv_models,
      preds=get_nv_group_preds
    )
  ),
  
  ###
  # process measures
  ###
  
  tar_target(
    process_ni_group_measures,
    func_process_measures(
      measures=get_ni_group_measures,
      use_id=FALSE
    )
  ),
  
  tar_target(
    process_nv_group_measures,
    func_process_measures(
      measures=get_nv_group_measures,
      use_id=FALSE
    )
  ),
  
  ###
  # process draws
  ###
  
  tar_target(
    process_nv_group_draws,
    func_process_draws(
      measures=get_nv_group_measures,
      use_id=FALSE
    )
  ),
  
  tar_target(
    process_ni_group_draws,
    func_process_draws(
      measures=get_ni_group_measures,
      use_id=FALSE
    )
  ),
  
  ###
  # get pairwise dose probabilities
  ###
  
  tar_target(
    pairwise_compare_ni,
    func_pairwise_compare(
      draws = process_ni_group_draws$draws_all
    )
  ),
  
  tar_target(
    pairwise_compare_nv,
    func_pairwise_compare(
      draws = process_nv_group_draws$draws_all
    )
  ),
  
  ###
  # get correlations from draws
  ###
  
  tar_target(
    get_correlations_overall,
    func_get_correlations(
      draws_ni=process_ni_group_draws$draws_all,
      draws_nv=process_nv_group_draws$draws_all,
      by = "overall"
    )
  ),
  
  tar_target(
    get_correlations_dose,
    func_get_correlations(
      draws_ni=process_ni_group_draws$draws_all,
      draws_nv=process_nv_group_draws$draws_all,
      by = "dose"
    )
  ),
  
  ###
  ### 05 tables ====
  ###
  
  ###
  # model metrics
  ###
  
  tar_target(
    tab_model_metrics,
    func_tab_model_metrics(
      metrics_ni=get_ni_model_metrics,
      metrics_nv=get_nv_model_metrics,
      filepath=here::here("results/tables/model-metrics/tab-model-metrics.rds")
    )
  ),
  
  ###
  # number censored observations
  ###
  
  tar_target(
    tab_n_ni_cens_obs,
    func_tab_n_cens(
      data=process_ni_data,
      filepath=here::here("results/tables/n-cens/ni/tab-ni-cens-")
    )
  ),
  
  tar_target(
    tab_n_nv_cens_obs,
    func_tab_n_cens(
      data=process_nv_data,
      filepath=here::here("results/tables/n-cens/nv/tab-nv-cens-")
    )
  ),
  
  ###
  # measures correlations
  ###
  
  tar_target(
    tab_correlations_dose,
    func_tab_correlations(
      correlations=get_correlations_dose, 
      model="pow", 
      by="dose", 
      baseline=TRUE,
      filepath=here::here("results/tables/correlations/dose/tab-cor-dose")
    )
  ),
  
  tar_target(
    tab_correlations_overall,
    func_tab_correlations(
      correlations=get_correlations_overall, 
      model="pow", 
      by="overall", 
      baseline=TRUE,
      filepath=here::here("results/tables/correlations/overall/tab-cor-overall")
    )
  ),
  
  ###
  ### 06 manuscript figures ====
  ###
  
  tar_target(
    fig_setup, 
    func_figure_setup()
  ),
  
  ###
  # model trajectories
  ###
  
  tar_target(
    plot_ni_group_model,
    func_plot_group_model_traj(
      predictions=process_ni_group_preds,
      data=process_ni_data,
      main=TRUE,
      filepath=here::here("results/figures/manuscript/plot-ni-model-trajectories.png"),
      height=5,
      width=7
    )
  ),
  
  tar_target(
    plot_nv_group_model,
    func_plot_group_model_traj(
      predictions=process_nv_group_preds,
      data=process_nv_data,
      main=TRUE,
      filepath=here::here("results/figures/manuscript/plot-nv-model-trajectories.png"),
      height=7,
      width=10
    )
  ),
  
  tar_target(
    plot_comb_group_model,
    func_plot_comb_group_model_traj(
      ni_plot=plot_ni_group_model, 
      nv_plot=plot_nv_group_model, 
      filepath=here::here("results/figures/manuscript/plot-comb-model-trajectories.png"),
      width=11, 
      height=14
    )
  ),
  
  ###
  # measures
  ###
  
  tar_target(
    plot_ni_pow_measures,
    func_plot_measures(
      dat_group=process_ni_group_measures, 
      model_name="pow", 
      study_name="NI",
      filepath=here::here("results/figures/manuscript/plot-ni-measures.png"), 
      height=7, width=5
    )
  ),
  
  tar_target(
    plot_nv_pow_measures,
    func_plot_measures(
      dat_group=process_nv_group_measures, 
      model_name="pow", 
      study_name="NV",
      filepath=here::here("results/figures/manuscript/plot-nv-measures.png"), 
      height=7, width=8
    )
  ),
  
  ###
  # measures vs response
  ###
  
  tar_target(
    plot_measures_vs_1year_draws,
    func_plot_measures_vs_1year_draws(
      setup=fig_setup,
      dat_measures_ni=process_ni_group_measures, 
      dat_measures_nv=process_nv_group_measures,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random,
      decay_day=90, 
      model_name="pow",
      baseline=TRUE,
      file_path=here::here("results/figures/manuscript/plot-measures-vs.png"), 
      height=8, 
      width=7
    )
  ),
  
  ###
  # peak vs decay, color response
  ###
  
  tar_target(
    plot_pow_group_peak_decay_draws,
    func_plot_group_peak_decay_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random, 
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      decay_day=90,
      baseline=TRUE,
      model_type="pow"
    )
  ),
  
  tar_target(
    plot_pow_group_peak_decay_1year_draws,
    func_plot_group_peak_decay_1year_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random, 
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      decay_day=90,
      baseline=TRUE,
      model_type="pow"
    )
  ),
  
  tar_target(
    plot_pow_comb_peak_decay_vs,
    func_plot_comb_peak_decay_vs(
      p1=plot_pow_group_peak_decay_draws,
      p2=plot_pow_group_peak_decay_1year_draws,
      filepath=here::here("results/figures/manuscript/plot-peak-decay-vs-comb.png"),
      width=8.5,
      height=8
    )
  ),
  
  ###
  ### 07 supplement figures ====
  ###
  
  ###
  # measures schematic
  ###
  
  tar_target(
    plot_schematic,
    func_plot_schematic(
      dat_pred=process_ni_group_preds, 
      file_path=here::here("results/figures/supplement/measures-schematic/measures-schematic.png"), 
      width=5, 
      height=4
    )
  ),
  
  ###
  # prior posterior check
  ###
  
  tar_target(
    plot_ni_prior_posterior_draws,
    func_plot_prior_posterior_draws(
      draws=get_ni_prior_posterior_draws,
      file_path=here::here("results/figures/supplement/prior-posterior-plots/ni"),
      width=6,
      height=4
    )
  ),
  
  tar_target(
    plot_nv_prior_posterior_draws,
    func_plot_prior_posterior_draws(
      draws=get_nv_prior_posterior_draws,
      file_path=here::here("results/figures/supplement/prior-posterior-plots/nv"),
      width=6,
      height=5
    )
  ),
  
  ###
  # fitted vs observed
  ###
  
  tar_target(
    plot_ni_fitted_vs_observed,
    func_plot_fitted_vs_observed(
      residuals=get_ni_residuals,
      file_path=here::here("results/figures/supplement/fitted-observed-plots/plot-ni-fitted-observed.png"),
      width=7,
      height=5
    )
  ),
  
  tar_target(
    plot_nv_fitted_vs_observed,
    func_plot_fitted_vs_observed(
      residuals=get_nv_residuals,
      file_path=here::here("results/figures/supplement/fitted-observed-plots/plot-nv-fitted-observed.png"),
      width=9,
      height=11
    )
  ),
  
  ###
  # individual-level model trajectories
  ###
  
  tar_target(
    plot_ni_id_model,
    func_plot_id_model(
      predictions=process_ni_id_preds, 
      data=process_ni_data, 
      filepath=here::here("results/figures/supplement/id-model-trajectories/ni/"),
      ncols=4,
      height=9, 
      width=7
    )
  ),
  
  tar_target(
    plot_nv_id_model,
    func_plot_id_model(
      predictions=process_nv_id_preds, 
      data=process_nv_data, 
      filepath=here::here("results/figures/supplement/id-model-trajectories/nv/"),
      ncols=5,
      height=14, 
      width=10
    )
  ),
  
  tar_target(
    plot_nv_id_model_zoom,
    func_plot_id_model(
      predictions=process_nv_id_preds, 
      data=process_nv_data, 
      zoom=TRUE,
      filepath=here::here("results/figures/supplement/id-model-trajectories/nv/zoom/"),
      ncols=5,
      height=14, 
      width=10
    )
  ),
  
  ###
  # absolute measures
  ###
  
  tar_target(
    plot_ni_pow_absolute_measures,
    func_plot_absolute_measures(
      setup=fig_setup,
      dat_group=process_ni_group_measures,
      model_name="pow",
      study_name="NI",
      filepath=here::here("results/figures/supplement/plot-absolute-measures/plot-ni-absolute-measures.png"),
      height=7, width=5
    )
  ),
  
  tar_target(
    plot_nv_pow_absolute_measures,
    func_plot_absolute_measures(
      setup=fig_setup,
      dat_group=process_nv_group_measures,
      model_name="pow",
      study_name="NV",
      filepath=here::here("results/figures/supplement/plot-absolute-measures/plot-nv-absolute-measures.png"),
      height=7, width=8
    )
  ),
  
  tar_target(
    plot_pow_absolute_measures_vs_1year_draws,
    func_plot_measures_vs_1year_draws(
      setup=fig_setup,
      dat_measures_ni=process_ni_group_measures, 
      dat_measures_nv=process_nv_group_measures,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random,
      decay_day=90, 
      baseline=FALSE,
      model_name="pow", 
      file_path=here::here("results/figures/supplement/plot-absolute-measures-vs/plot-absolute-measures-vs.png"), 
      height=8, 
      width=7.5
    )
  ),
  
  tar_target(
    plot_pow_absolute_group_peak_decay_draws,
    func_plot_group_peak_decay_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random, 
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      decay_day=90,
      baseline=FALSE,
      model_type="pow"
    )
  ),
  
  tar_target(
    plot_pow_absolute_group_peak_decay_1year_draws,
    func_plot_group_peak_decay_1year_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random, 
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      decay_day=90,
      model_type="pow",
      baseline=FALSE
    )
  ),
  
  tar_target(
    plot_pow_absolute_comb_peak_decay_vs,
    func_plot_comb_peak_decay_vs(
      p1=plot_pow_absolute_group_peak_decay_draws,
      p2=plot_pow_absolute_group_peak_decay_1year_draws,
      filepath=here::here("results/figures/supplement/plot-absolute-measures-vs/plot-absolute-peak-decay-vs-comb.png"),
      width=8.5,
      height=8
    )
  ),
  
  ###
  # pairwise dose probabilities
  ###
  
  tar_target(
    plot_pairwise_compare_ni,
    func_plot_pairwise_compare(
      probabilities=pairwise_compare_ni,
      model="pow",
      filepath = here::here("results/figures/supplement/plot-pairwise-measure-probabilities/plot-ni-pairwise-measure-probabilities.png"),
      height=6,
      width=5
    )
  ),
  
  tar_target(
    plot_pairwise_compare_nv,
    func_plot_pairwise_compare(
      probabilities=pairwise_compare_nv,
      model="pow",
      filepath = here::here("results/figures/supplement/plot-pairwise-measure-probabilities/plot-nv-pairwise-measure-probabilities.png"),
      height=11,
      width=9
    )
  ),
  
  ###
  # alt decay rate days
  ###
  
  tar_target(
    plot_ni_decay_rates,
    func_plot_decay_rates(
      setup=fig_setup,
      dat_group=process_ni_group_measures,
      model_name="pow",
      study_name="NI",
      filepath=here::here("results/figures/supplement/plot-decay-rates/plot-ni-decay-rates.png"),
      width=5,
      height=6
    )
  ),
  
  tar_target(
    plot_nv_decay_rates,
    func_plot_decay_rates(
      setup=fig_setup,
      dat_group=process_nv_group_measures,
      model_name="pow",
      study_name="NV",
      filepath=here::here("results/figures/supplement/plot-decay-rates/plot-nv-decay-rates.png"),
      width=8,
      height=6
    )
  ),
  
  tar_target(
    plot_pow_group_peak_decay60_1year_draws,
    func_plot_group_peak_decay_1year_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random,
      nv_draws=process_nv_group_draws$draws_random,
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      baseline=TRUE,
      decay_day=60,
      model_type="pow"
    )
  ),
  
  tar_target(
    plot_pow_group_peak_decay180_1year_draws,
    func_plot_group_peak_decay_1year_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random,
      nv_draws=process_nv_group_draws$draws_random,
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      baseline=TRUE,
      decay_day=180,
      model_type="pow"
    )
  ),
  
  tar_target(
    plot_pow_group_peak_decay365_1year_draws,
    func_plot_group_peak_decay_1year_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random,
      nv_draws=process_nv_group_draws$draws_random,
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      baseline=TRUE,
      decay_day=365,
      model_type="pow"
    )
  ),
  
  tar_target(
    plot_comb_pow_group_peak_decay_1year_draws,
    func_plot_comb_group_peak_decay_1year_draws(
      p1=plot_pow_group_peak_decay60_1year_draws,
      p2=plot_pow_group_peak_decay180_1year_draws,
      p3=plot_pow_group_peak_decay365_1year_draws,
      filepath=here::here("results/figures/supplement/plot-decay-rates/plot-decay-days.png"),
      height=9,
      width=8.5
    )
  ),
  
  ###
  # GII.4 plots
  ###
  
  tar_target(
    plot_gii4_measures_vs_1year_draws,
    func_plot_measures_vs_1year_draws(
      setup=fig_setup,
      dat_measures_ni=process_ni_group_measures, 
      dat_measures_nv=process_nv_group_measures,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random,
      decay_day=90, 
      main=FALSE,
      baseline=TRUE,
      model_name="pow", 
      file_path=here::here("results/figures/supplement/plot-gii4-measures-draws-vs/plot-gii4-measures-vs.png"), 
      height=8, 
      width=5
    )
  ),
  
  tar_target(
    plot_gii4_pow_group_peak_decay_draws,
    func_plot_group_peak_decay_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random, 
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      decay_day=90,
      main=FALSE,
      baseline=TRUE,
      model_type="pow"
    )
  ),
  
  tar_target(
    plot_gii4_pow_group_peak_decay_1year_draws,
    func_plot_group_peak_decay_1year_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random, 
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      decay_day=90,
      main=FALSE,
      baseline=TRUE,
      model_type="pow"
    )
  ),
  
  tar_target(
    plot_gii4_pow_comb_peak_decay_vs,
    func_plot_comb_peak_decay_vs(
      p1=plot_gii4_pow_group_peak_decay_draws,
      p2=plot_gii4_pow_group_peak_decay_1year_draws,
      filepath=here::here("results/figures/supplement/plot-gii4-measures-draws-vs/plot-gii4-peak-decay-vs-comb.png"),
      width=6,
      height=8
    )
  ),
  
  tar_target(
    plot_gii4_pow_absolute_measures_vs_1year_draws,
    func_plot_measures_vs_1year_draws(
      setup=fig_setup,
      dat_measures_ni=process_ni_group_measures, 
      dat_measures_nv=process_nv_group_measures,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random,
      decay_day=90, 
      baseline=FALSE,
      main=FALSE,
      model_name="pow", 
      file_path=here::here("results/figures/supplement/plot-absolute-measures-vs/gii4/plot-gii4-absolute-measures-vs.png"), 
      height=8, 
      width=5
    )
  ),
  
  tar_target(
    plot_gii4_pow_absolute_group_peak_decay_draws,
    func_plot_group_peak_decay_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random, 
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      decay_day=90,
      baseline=FALSE,
      main=FALSE,
      model_type="pow"
    )
  ),
  
  tar_target(
    plot_gii4_pow_absolute_group_peak_decay_1year_draws,
    func_plot_group_peak_decay_1year_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random, 
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      decay_day=90,
      model_type="pow",
      baseline=FALSE,
      main=FALSE
    )
  ),
  
  tar_target(
    plot_gii4_pow_absolute_comb_peak_decay_vs,
    func_plot_comb_peak_decay_vs(
      p1=plot_gii4_pow_absolute_group_peak_decay_draws,
      p2=plot_gii4_pow_absolute_group_peak_decay_1year_draws,
      filepath=here::here("results/figures/supplement/plot-absolute-measures-vs/gii4/plot-gii4-absolute-peak-decay-vs-comb.png"),
      width=6,
      height=8
    )
  ),
  
  ###
  # exponential decay model plots
  ###
  
  tar_target(
    plot_ni_exp_measures,
    func_plot_measures(
      dat_group=process_ni_group_measures, 
      model_name="exp", 
      study_name="NI",
      filepath=here::here("results/figures/supplement/exponential-results/plot-model-measures/plot-ni-measures.png"), 
      height=7, width=5
    )
  ),
  
  tar_target(
    plot_nv_exp_measures,
    func_plot_measures(
      dat_group=process_nv_group_measures, 
      model_name="exp", 
      study_name="NV",
      filepath=here::here("results/figures/supplement/exponential-results/plot-model-measures/plot-nv-measures.png"), 
      height=7, width=8
    )
  ),
  
  tar_target(
    plot_exp_measures_vs_1year_draws,
    func_plot_measures_vs_1year_draws(
      setup=fig_setup,
      dat_measures_ni=process_ni_group_measures, 
      dat_measures_nv=process_nv_group_measures,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random,
      decay_day=90, 
      model_name="exp", 
      baseline=TRUE,
      file_path=here::here("results/figures/supplement/exponential-results/plot-measures-draws-vs/plot-measures-vs.png"), 
      height=8, 
      width=7.5
    )
  ),
  
  tar_target(
    plot_exp_group_peak_decay_draws,
    func_plot_group_peak_decay_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random, 
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      decay_day=90,
      model_type="exp",
      baseline=TRUE
    )
  ),
  
  tar_target(
    plot_exp_group_peak_decay_1year_draws,
    func_plot_group_peak_decay_1year_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random, 
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      decay_day=90,
      model_type="exp",
      baseline=TRUE
    )
  ),
  
  tar_target(
    plot_exp_comb_peak_decay_vs,
    func_plot_comb_peak_decay_vs(
      p1=plot_exp_group_peak_decay_draws,
      p2=plot_exp_group_peak_decay_1year_draws,
      filepath=here::here("results/figures/supplement/exponential-results/plot-measures-draws-vs/plot-peak-decay-vs-comb.png"),
      width=8.5,
      height=8
    )
  ),
  
  tar_target(
    plot_gii4_exp_measures_vs_1year_draws,
    func_plot_measures_vs_1year_draws(
      setup=fig_setup,
      dat_measures_ni=process_ni_group_measures, 
      dat_measures_nv=process_nv_group_measures,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random,
      decay_day=90, 
      main=FALSE,
      baseline=TRUE,
      model_name="exp", 
      file_path=here::here("results/figures/supplement/exponential-results/plot-measures-draws-vs/gii4/plot-gii4-measures-vs.png"), 
      height=8, 
      width=5
    )
  ),
  
  tar_target(
    plot_gii4_exp_group_peak_decay_draws,
    func_plot_group_peak_decay_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random, 
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      decay_day=90,
      main=FALSE,
      baseline=TRUE,
      model_type="exp"
    )
  ),
  
  tar_target(
    plot_gii4_exp_group_peak_decay_1year_draws,
    func_plot_group_peak_decay_1year_draws(
      setup=fig_setup,
      ni_draws=process_ni_group_draws$draws_random, 
      nv_draws=process_nv_group_draws$draws_random, 
      ni_measures=process_ni_group_measures,
      nv_measures=process_nv_group_measures,
      decay_day=90,
      main=FALSE,
      model_type="exp",
      baseline=TRUE
    )
  ),
  
  tar_target(
    plot_gii4_exp_comb_peak_decay_vs,
    func_plot_comb_peak_decay_vs(
      p1=plot_gii4_exp_group_peak_decay_draws,
      p2=plot_gii4_exp_group_peak_decay_1year_draws,
      filepath=here::here("results/figures/supplement/exponential-results/plot-measures-draws-vs/gii4/plot-gii4-peak-decay-vs-comb.png"),
      width=6,
      height=8
    )
  ),
  
  ###
  ### 08 supplement analyses ====
  ###
  
  ###
  # alternative censoring bounds
  ###
  
  ###
  ### 01 data processing
  ###
  
  tar_target(
    alt_process_ni_data,
    func_alt_process_ni_data(
      data=process_ni_data
    )
  ),
  
  tar_target(
    alt_process_nv_data,
    func_alt_process_nv_data(
      data=process_nv_data
    )
  ),
  
  ###
  # increase bound (10^-3)
  ###
  
  ###
  ### 02 model fitting
  ###
  
  tar_target(
    alt1_fit_ni_models,
    func_fit_models(
      data=alt_process_ni_data$dat_ni_increase,
      setup=setup_model
    )
  ),
  
  tar_target(
    alt1_fit_nv_models,
    func_fit_models(
      data=alt_process_nv_data$dat_nv_increase,
      setup=setup_model
    )
  ),
  
  tar_target(
    alt1_get_ni_model_metrics,
    func_get_model_metrics(
      models=alt1_fit_ni_models
    )
  ),
  
  tar_target(
    alt1_get_nv_model_metrics,
    func_get_model_metrics(
      models=alt1_fit_nv_models
    )
  ),
  
  ###
  ### 03 posterior predictions
  ###
  
  tar_target(
    alt1_get_ni_group_preds,
    func_get_preds(
      data=alt_process_ni_data$dat_ni_increase,
      models=alt1_fit_ni_models,
      use_id=FALSE
    )
  ),
  
  tar_target(
    alt1_get_nv_group_preds,
    func_get_preds(
      data=alt_process_nv_data$dat_nv_increase,
      models=alt1_fit_nv_models,
      use_id=FALSE
    )
  ),
  
  tar_target(
    alt1_process_ni_group_preds,
    func_process_preds(
      predictions=alt1_get_ni_group_preds,
      use_id=FALSE
    )
  ),
  
  tar_target(
    alt1_process_nv_group_preds,
    func_process_preds(
      predictions=alt1_get_nv_group_preds,
      use_id=FALSE
    )
  ),
  
  ###
  ### 04 measures
  ###
  
  tar_target(
    alt1_get_ni_group_measures,
    func_get_measures_group(
      data=alt_process_ni_data$dat_ni_increase,
      models=alt1_fit_ni_models,
      preds=alt1_get_ni_group_preds
    )
  ),
  
  tar_target(
    alt1_get_nv_group_measures,
    func_get_measures_group(
      data=alt_process_nv_data$dat_nv_increase,
      models=alt1_fit_nv_models,
      preds=alt1_get_nv_group_preds
    )
  ),
  
  tar_target(
    alt1_process_ni_group_measures,
    func_process_measures(
      measures=alt1_get_ni_group_measures,
      use_id=FALSE
    )
  ),
  
  tar_target(
    alt1_process_nv_group_measures,
    func_process_measures(
      measures=alt1_get_nv_group_measures,
      use_id=FALSE
    )
  ),
  
  ###
  ### 05 tables 
  ###
  
  tar_target(
    alt1_tab_model_metrics,
    func_tab_model_metrics(
      metrics_ni=alt1_get_ni_model_metrics,
      metrics_nv=alt1_get_nv_model_metrics,
      filepath=here::here("results/tables/cens-bound-sensitivity/increase/tab-model-metrics.rds")
    )
  ),
  
  ###
  ### 06 figures
  ###
  
  tar_target(
    alt1_plot_ni_group_model,
    func_plot_group_model_traj(
      predictions=alt1_process_ni_group_preds,
      data=alt_process_ni_data$dat_ni_increase,
      filepath=here::here("results/figures/supplement/cens-bound-sensitivity/increase/plot-ni-model-trajectories.png"),
      height=5,
      width=7
    )
  ),
  
  tar_target(
    alt1_plot_nv_group_model,
    func_plot_group_model_traj(
      predictions=alt1_process_nv_group_preds,
      data=alt_process_nv_data$dat_nv_increase,
      filepath=here::here("results/figures/supplement/cens-bound-sensitivity/increase/plot-nv-model-trajectories.png"),
      height=7,
      width=10
    )
  ),
  
  tar_target(
    alt1_plot_ni_pow_measures,
    func_plot_measures(
      dat_group=alt1_process_ni_group_measures, 
      model_name="pow", 
      study_name="NI",
      filepath=here::here("results/figures/supplement/cens-bound-sensitivity/increase/plot-ni-measures.png"), 
      height=7, width=5
    )
  ),
  
  tar_target(
    alt1_plot_nv_pow_measures,
    func_plot_measures(
      dat_group=alt1_process_nv_group_measures, 
      model_name="pow", 
      study_name="NV",
      filepath=here::here("results/figures/supplement/cens-bound-sensitivity/increase/plot-nv-measures.png"), 
      height=7, width=8
    )
  ),
  
  ###
  # decrease bound (10^-15)
  ###
  
  ###
  ### 02 model fitting
  ###
  
  tar_target(
    alt2_fit_ni_models,
    func_fit_models(
      data=alt_process_ni_data$dat_ni_decrease,
      setup=setup_model
    )
  ),
  
  tar_target(
    alt2_fit_nv_models,
    func_fit_models(
      data=alt_process_nv_data$dat_nv_decrease,
      setup=setup_model
    )
  ),
  
  tar_target(
    alt2_get_ni_model_metrics,
    func_get_model_metrics(
      models=alt2_fit_ni_models
    )
  ),
  
  tar_target(
    alt2_get_nv_model_metrics,
    func_get_model_metrics(
      models=alt2_fit_nv_models
    )
  ),
  
  ###
  ### 03 posterior predictions
  ###
  
  tar_target(
    alt2_get_ni_group_preds,
    func_get_preds(
      data=alt_process_ni_data$dat_ni_decrease,
      models=alt2_fit_ni_models,
      use_id=FALSE
    )
  ),
  
  tar_target(
    alt2_get_nv_group_preds,
    func_get_preds(
      data=alt_process_nv_data$dat_nv_decrease,
      models=alt2_fit_nv_models,
      use_id=FALSE
    )
  ),
  
  tar_target(
    alt2_process_ni_group_preds,
    func_process_preds(
      predictions=alt2_get_ni_group_preds,
      use_id=FALSE
    )
  ),
  
  tar_target(
    alt2_process_nv_group_preds,
    func_process_preds(
      predictions=alt2_get_nv_group_preds,
      use_id=FALSE
    )
  ),
  
  ###
  ### 04 measures
  ###
  
  tar_target(
    alt2_get_ni_group_measures,
    func_get_measures_group(
      data=alt_process_ni_data$dat_ni_decrease,
      models=alt2_fit_ni_models,
      preds=alt2_get_ni_group_preds
    )
  ),
  
  tar_target(
    alt2_get_nv_group_measures,
    func_get_measures_group(
      data=alt_process_nv_data$dat_nv_decrease,
      models=alt2_fit_nv_models,
      preds=alt2_get_nv_group_preds
    )
  ),
  
  tar_target(
    alt2_process_ni_group_measures,
    func_process_measures(
      measures=alt2_get_ni_group_measures,
      use_id=FALSE
    )
  ),
  
  tar_target(
    alt2_process_nv_group_measures,
    func_process_measures(
      measures=alt2_get_nv_group_measures,
      use_id=FALSE
    )
  ),
  
  ###
  ### 05 tables
  ###
  
  tar_target(
    alt2_tab_model_metrics,
    func_tab_model_metrics(
      metrics_ni=alt2_get_ni_model_metrics,
      metrics_nv=alt2_get_nv_model_metrics,
      filepath=here::here("results/tables/cens-bound-sensitivity/decrease/tab-model-metrics.rds")
    )
  ),
  
  ###
  ### 06 figures
  ###
  
  tar_target(
    alt2_plot_ni_group_model,
    func_plot_group_model_traj(
      predictions=alt2_process_ni_group_preds,
      data=alt_process_ni_data$dat_ni_decrease,
      filepath=here::here("results/figures/supplement/cens-bound-sensitivity/decrease/plot-ni-model-trajectories.png"),
      height=5,
      width=7
    )
  ),
  
  tar_target(
    alt2_plot_nv_group_model,
    func_plot_group_model_traj(
      predictions=alt2_process_nv_group_preds,
      data=alt_process_nv_data$dat_nv_decrease,
      filepath=here::here("results/figures/supplement/cens-bound-sensitivity/decrease/plot-nv-model-trajectories.png"),
      height=7,
      width=10
    )
  ),
  
  tar_target(
    alt2_plot_ni_pow_measures,
    func_plot_measures(
      dat_group=alt2_process_ni_group_measures, 
      model_name="pow", 
      study_name="NI",
      filepath=here::here("results/figures/supplement/cens-bound-sensitivity/decrease/plot-ni-measures.png"), 
      height=7, width=5
    )
  ),
  
  tar_target(
    alt2_plot_nv_pow_measures,
    func_plot_measures(
      dat_group=alt2_process_nv_group_measures, 
      model_name="pow", 
      study_name="NV",
      filepath=here::here("results/figures/supplement/cens-bound-sensitivity/decrease/plot-nv-measures.png"), 
      height=7, width=8
    )
  )
)
