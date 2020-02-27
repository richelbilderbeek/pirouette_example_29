# pirouette example 29:
# Check the effect of MCMC chain length on ESS estimation
# Similar to pirouette example 28, except multiple MCMC settings
# From https://github.com/richelbilderbeek/pirouette_article/issues/55 :
#
# Write script that shows the true and twin error for a fixed tree prior and
# reasonably small trees (say in the magnitudes of hundreds of them)
# with 10:40 taxa
suppressMessages(library(pirouette))
suppressMessages(library(ggplot2))
library(testthat)
expect_true(mcbette::can_run_mcbette())

root_folder <- getwd()
example_no <- 29
n_replicates <- 5
unique_mcmc_chain_lengths <- c(1e6, 1e7, 1e8)
crown_age <- 10
n_taxa <- 6
is_testing <- is_on_travis()

# Number of replicates per number of taxa
if (is_testing) {
  n_replicates <- 2
  root_folder <- tempdir()
  unique_mcmc_chain_lengths <- c(3000, 4000)
}
rng_seeds <- seq(314, 314 - 1 + length(unique_mcmc_chain_lengths) * n_replicates)
mcmc_chain_lengths <- rep(unique_mcmc_chain_lengths, each = n_replicates)
expect_equal(length(rng_seeds), length(mcmc_chain_lengths))
################################################################################
# Create phylogenies
################################################################################
phylogenies <- list()
for (i in seq_along(rng_seeds)) {
  rng_seed <- rng_seeds[i]

  set.seed(rng_seed)

  phylogenies[[i]] <- create_yule_tree(
    n_taxa = n_taxa,
    crown_age = crown_age
  )
}
expect_equal(length(phylogenies), length(rng_seeds))
################################################################################
# Create pirouette parameter sets
################################################################################
pir_paramses <- list()
for (i in seq_along(phylogenies)) {

  alignment_params <- create_alignment_params(
    sim_tral_fun = get_sim_tral_with_std_nsm_fun(
      mutation_rate = 1.0 / crown_age
    ),
    root_sequence = create_blocked_dna(length = 1000)
  )

  # Hand-pick a generating model
  # By default, this is JC69, strict, Yule
  generative_experiment <- create_gen_experiment()
  # Create the set of candidate birth-death experiments
  candidate_experiments <- create_all_bd_experiments(
    exclude_model = generative_experiment$inference_model
  )
  # Combine all experiments
  experiments <- c(list(generative_experiment), candidate_experiments)

  twinning_params <- create_twinning_params(
    sim_twin_tree_fun = get_sim_bd_twin_tree_fun(),
    sim_twal_fun = get_sim_twal_same_n_muts_fun(
      mutation_rate = 1.0 / crown_age,
      max_n_tries = 10000
    ),
    twin_evidence_filename = get_temp_evidence_filename()
  )

  pir_params <- create_pir_params(
    alignment_params = alignment_params,
    experiments = experiments,
    twinning_params = twinning_params,
    evidence_filename = get_temp_evidence_filename()
  )

  pir_paramses[[i]] <- pir_params
}
expect_equal(length(pir_paramses), length(phylogenies))
################################################################################
# Shorter run on Travis
################################################################################
if (is_testing) {
  for (i in seq_along(pir_paramses)) {
    pir_paramses[[i]]$experiments <- shorten_experiments(
      pir_paramses[[i]]$experiments
    )
  }
}

################################################################################
# Set the RNG seeds
################################################################################
pir_paramses <- renum_rng_seeds(
  pir_paramses = pir_paramses,
  rng_seeds = seq(314, 314 - 1 + length(pir_paramses))
)

################################################################################
# Rename filenames
################################################################################
for (i in seq_along(pir_paramses)) {
  rng_seed <- pir_paramses[[i]]$alignment_params$rng_seed
  pir_paramses[[i]] <- pir_rename_to_std(
    pir_params = pir_paramses[[i]],
    folder_name = file.path(root_folder, paste0("example_", example_no, "_", rng_seed))
  )
}

################################################################################
# Save tree to files
################################################################################
for (i in seq_along(pir_paramses)) {
  expect_equal(length(pir_paramses), length(phylogenies))
  rng_seed <- pir_paramses[[i]]$alignment_params$rng_seed
  folder_name <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))

  # Create folder, do not warn if it already exists
  dir.create(folder_name, showWarnings = FALSE, recursive = TRUE)
  ape::write.tree(
    phylogenies[[i]],
    file = file.path(folder_name, "true_tree.newick")
  )
}

################################################################################
# Delete previous files
################################################################################
for (pir_params in pir_paramses) {
  check_pir_params(pir_params)
  rm_pir_param_files(pir_params)
}

################################################################################
# Do the runs
################################################################################
pir_outs <- pir_runs(
  phylogenies = phylogenies,
  pir_paramses = pir_paramses
)

################################################################################
# Save
################################################################################
for (i in seq_along(pir_outs)) {
  expect_equal(length(pir_paramses), length(pir_outs))
  rng_seed <- pir_paramses[[i]]$alignment_params$rng_seed
  folder_name <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))

  utils::write.csv(
    x = pir_outs[[i]],
    file = file.path(folder_name, "errors.csv"),
    row.names = FALSE
  )

  pir_plots(pir_outs[[i]]) +
    ggsave(file.path(folder_name, "errors.png"))
}

