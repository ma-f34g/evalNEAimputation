library(conflicted)
library(tidyverse)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
library(duckplyr)
library(glue)

test_samples <- read_parquet_duckdb("data/out/test_samples.parquet")

process_rate <- function(rate) {
  "data/out/imputed_filtered_{rate}.parquet" |>
    glue() |>
    read_parquet_duckdb() |>
    full_join(test_samples, by = c("POS", "Sample")) |>
    mutate(
      Allele1 = if_else(is.na(Allele1.y), 0, Allele1.y),
      Allele2 = if_else(is.na(Allele2.y), 0, Allele2.y),
      Allele1_imp = if_else(is.na(Allele1.x), 0, Allele1.x),
      Allele2_imp = if_else(is.na(Allele2.x), 0, Allele2.x)
    ) |>
    select(POS, Sample, Allele1, Allele2, Allele1_imp, Allele2_imp) |>
    filter(Allele1 != Allele2) |>
    mutate(
      phase_match = if_else(
        Allele1 == Allele1_imp & Allele2 == Allele2_imp,
        1,
        if_else(Allele1 == Allele2_imp & Allele2 == Allele1_imp, 0, NA)
      )
    ) |>
    arrange(Sample, POS) |>
    mutate(prev_phase_match = lag(phase_match), .by = Sample) |>
    filter(!is.na(prev_phase_match)) |>
    summarize(
      switch_count = sum(phase_match != prev_phase_match, na.rm = TRUE),
      total_switch_opportunities = n(),
      .by = Sample
    ) |>
    mutate(
      rate = rate,
      masking_rate = rate / 10,
      switch_error_rate = switch_count / total_switch_opportunities
    )
}

map_dfr(1:9, process_rate) |>
  compute_parquet("data/out/switch_error_rate.parquet")
