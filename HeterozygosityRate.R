library(conflicted)
library(tidyverse)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
library(duckplyr)
library(glue)

process_rate <- function(rate) {
  "data/out/imputed_filtered_{rate}.parquet" |>
    glue() |>
    read_parquet_duckdb() |>
    filter(Allele1 != 0 | Allele2 != 0) |>
    summarise(
      hetero = sum(Allele1 != Allele2),
      homo_nonref = sum(Allele1 == Allele2),
      .by = Sample
    ) |>
    mutate(
      rate = rate,
      masking_rate = rate / 10,
      heterozygosity_ratio = hetero / homo_nonref
    )
}

map_dfr(1:9, process_rate) |>
  compute_csv("data/out/heterozygosity_ratio_imp.csv")

"data/out/test_samples.parquet" |>
  read_parquet_duckdb() |>
  filter(Allele1 != 0 | Allele2 != 0) |>
  summarise(
    hetero = sum(Allele1 != Allele2),
    homo_nonref = sum(Allele1 == Allele2),
    .by = Sample
  ) |>
  mutate(heterozygosity_ratio = hetero / homo_nonref) |>
  compute_csv("data/out/heterozygosity_ratio.csv")
