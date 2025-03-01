# %%
library(tidyverse)
library(arrow)
library(slendr)
init_env()
check_env()
set.seed(42)

# Define Populations ----
# %%
afr <- population("AFR", time = 6e6, N = 10000)
eur <- population("EUR", parent = afr, time = 70e3, N = 5000)
nea <- population("NEA", parent = afr, time = 600e3, N = 1000, remove = 40e3)

gf <- gene_flow(from = nea, to = eur, rate = 0.05, start = 55000, end = 50000)

model <- compile_model(
  populations = list(afr, eur, nea),
  gene_flow = gf,
  generation_time = 28,
  simulation_length = 6e6
)

# Sample reference panel and targets ----
# %%
reference_panel_size <- 1000
present_samples <- schedule_sampling(
  model,
  times = 0,
  list(eur, reference_panel_size)
)

broad_samples <- seq(from = 0, to = 50000, length.out = 1000) # 5000
detailed_samples <- seq(from = 49000, to = 50000, length.out = 1000)
sampling_times <- sort(unique(c(broad_samples, detailed_samples)))

test_samples <- schedule_sampling(model, times = sampling_times, list(eur, 1))
samples_count <- nrow(test_samples)
ts <- msprime(
  model,
  sequence_length = 100000000,
  recombination_rate = 1e-8,
  samples = rbind(present_samples, test_samples),
  random_seed = 42,
  verbose = TRUE
)
samples <- ts_samples(ts)

# %%
simulate_site_mask <- function(ts, mask_rate) {
  # TRUE for sites to exclude, FALSE to include
  num_sites <- ts$num_sites
  site_mask <- runif(num_sites) < mask_rate

  np <- reticulate::import("numpy")
  np$array(site_mask, dtype = "bool")
}

ts_vcf_damaged <- function(
  ts,
  path,
  chrom = NULL,
  individuals = NULL,
  site_mask = NULL
) {
  if (!attr(ts, "recapitated") && !ts_coalesced(ts)) {
    stop(
      "Tree sequence was not recapitated and some nodes do not ",
      "have parents over some portion of their genome. This is interpreted as ",
      "missing data, which is not currently supported by tskit. For more context, ",
      "take a look at <https://github.com/tskit-dev/tskit/issues/301#issuecomment-520990038>.",
      call. = FALSE
    )
  }

  if (!attr(ts, "mutated")) {
    stop(
      "Attempting to extract genotypes from a tree sequence which has not been mutated",
      call. = FALSE
    )
  }

  data <- ts_nodes(ts) %>%
    dplyr::filter(!is.na(name)) %>%
    dplyr::as_tibble() %>%
    dplyr::distinct(name, ind_id)

  if (is.null(individuals)) individuals <- data$name

  present <- individuals %in% unique(data$name)
  if (!all(present)) {
    stop(
      "",
      paste(individuals[!present], collapse = ", "),
      " not present in the tree sequence",
      call. = FALSE
    )
  }

  if (!is.null(site_mask)) {
    np <- reticulate::import("numpy")
    if (!inherits(site_mask, "numpy.ndarray")) {
      site_mask <- np$array(as.logical(site_mask), dtype = "bool")
    }
    if (length(site_mask) != ts$num_sites) {
      stop(
        "Site mask length (",
        length(site_mask),
        ") does not match number of sites (",
        ts$num_sites,
        ")",
        call. = FALSE
      )
    }
  }

  gzip <- reticulate::import("gzip")
  with(reticulate::`%as%`(gzip$open(path.expand(path), "wt"), vcf_file), {
    ts$write_vcf(
      vcf_file,
      contig_id = chrom,
      individuals = as.integer(data$ind_id),
      individual_names = data$name,
      site_mask = site_mask
    )
  })
}

# Write samples data to disk ----
# %%
unlink("data/*", recursive = TRUE)
write_parquet(samples, "data/sample_ages.parquet")

test_samples <- ts_simplify(
  ts,
  simplify_to = c(paste0("EUR_", 1:samples_count))
)
reference_panel <- ts_simplify(
  ts,
  simplify_to = c(paste0(
    "EUR_",
    (samples_count + 1):(samples_count + reference_panel_size)
  ))
)

ts_vcf_damaged(test_samples, path = "data/test_samples.vcf.gz", chrom = "1")
ts_vcf_damaged(
  reference_panel,
  path = "data/reference_panel.vcf.gz",
  chrom = "1"
)

for (mask_rate in seq(0.1, 0.9, by = 0.1)) {
  site_mask <- simulate_site_mask(test_samples, mask_rate = mask_rate)
  filename <- paste0("data/damaged_test_samples_", mask_rate, ".vcf.gz")
  ts_vcf_damaged(
    test_samples,
    path = filename,
    chrom = "1",
    site_mask = site_mask
  )
}

ts_tracts(ts, census = 55000, source = "NEA", target = "EUR") |>
  write_parquet("data/neanderthal_tracts.parquet")
