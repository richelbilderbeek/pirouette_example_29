# pirouette example 29:
# Check the effect of MCMC chain length on ESS estimation
# Similar to pirouette example 28, except multiple MCMC settings
# From https://github.com/richelbilderbeek/pirouette_article/issues/55 :
#
# Write script that shows the true and twin error for a fixed tree prior and
# reasonably small trees (say in the magnitudes of hundreds of them)
# with 10:40 taxa
library(pirouette)
library(beautier)
library(beastier)
library(testthat)
library(ggplot2)

# Constants
example_no <- 29
rng_seed <- 314
crown_age <- 10
mcmc_chain_lengths <- c(1e5, 1e6, 1e7, 1e8)
n_phylogenies_per_mcmc_chain_length <- 3
folder_name <- paste0("example_", example_no)
is_testing <- is_on_ci()
if (is_testing) {
  mcmc_chain_lengths <- c(3000, 4000)
  n_phylogenies_per_mcmc_chain_length <- 2
}
n_mcmc_chain_lengths <- length(mcmc_chain_lengths)
n_pir_params <- n_mcmc_chain_lengths * n_phylogenies_per_mcmc_chain_length

# Create phylogenies
phylogenies <- list()
for (i in seq_len(n_phylogenies_per_mcmc_chain_length)) {
  set.seed(rng_seed)
  phylogeny <- create_yule_tree(n_taxa = 6, crown_age = 10)
  phylogenies[[i]] <- phylogeny
}
# 1 2 3 1 2 3
phylogenies <- rep(phylogenies, n_mcmc_chain_lengths)
expect_equal(n_pir_params, length(phylogenies))

# Create pirouette parameter sets
pir_paramses <- create_std_pir_paramses(
  n = n_pir_params,
  folder_name = folder_name
)
expect_equal(length(pir_paramses), length(phylogenies))
if (is_testing) {
  pir_paramses <- shorten_pir_paramses(pir_paramses)
}

# Set the alignment lengths
# 1 1 1 2 2 2
mcmc_chain_lengthses <- rep(
  mcmc_chain_lengths, each = n_phylogenies_per_mcmc_chain_length
)
expect_equal(length(mcmc_chain_lengthses), length(pir_paramses))
for (i in seq_along(mcmc_chain_lengthses)) {
  for (j in seq_along(pir_paramses[[i]]$experiments)) {
    pir_paramses[[i]]$experiments[[j]]$inference_model$mcmc$chain_length <-
      mcmc_chain_lengthses[[i]]
  }
}

# Do the runs per MCMC chain length
pir_outs <- vector("list", n_pir_params)
for (i in seq_along(mcmc_chain_lengths)) {
  n <- mcmc_chain_lengths[i]
  from_index <- ((i - 1) * n_phylogenies_per_mcmc_chain_length) + 1
  to_index <- ((i - 1) * n_phylogenies_per_mcmc_chain_length) + n_phylogenies_per_mcmc_chain_length
  pir_outs[from_index:to_index] <- pir_runs(
    phylogenies = phylogenies[from_index:to_index],
    pir_paramses = pir_paramses[from_index:to_index]
  )
}

# Save per MCMC chain length
for (i in seq_along(mcmc_chain_lengths)) {
  n <- mcmc_chain_lengths[i]
  from_index <- ((i - 1) * n_phylogenies_per_mcmc_chain_length) + 1
  to_index <- from_index + n_phylogenies_per_mcmc_chain_length - 1
  pir_plots(
    pir_outs = pir_outs[from_index:to_index]
  ) + ggtitle(
      paste(
        "MCMC chain length:", n, ", number of replicates: ", n_phylogenies_per_mcmc_chain_length
      )
    ) +
    ggsave(filename = paste0("errors_", i, ".png"), width = 7, height = 7)
}

# Save individual runs
expect_equal(length(pir_paramses), length(pir_outs))
expect_equal(length(pir_paramses), length(phylogenies))
for (i in seq_along(pir_outs)) {
  pir_save(
    phylogeny = phylogenies[[i]],
    pir_params = pir_paramses[[i]],
    pir_out = pir_outs[[i]],
    folder_name = dirname(pir_paramses[[i]]$alignment_params$fasta_filename)
  )
}

