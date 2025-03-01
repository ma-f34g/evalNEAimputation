library(tidyverse)
library(duckplyr)

tracts <- read_parquet_duckdb("data/out/tract_intervals.parquet")
sample_ages <- read_parquet_duckdb("data/sample_ages.parquet") |>
	select(Sample = name, age = time)

nea <- read_parquet_duckdb("data/neanderthal_tracts.parquet") |>
	select(Sample = name, length) |>
	summarise(NEA = sum(length / 100000000), .by = "Sample")

ser <- read_parquet_duckdb("data/out/switch_error_rate.parquet")

process_rate <- function(rate) {
	original <- read_parquet_duckdb(glue::glue(
		"data/out/masked_samples_{rate}.parquet"
	))
	imputed <- read_parquet_duckdb(glue::glue(
		"data/out/imputed_filtered_{rate}.parquet"
	))

	# Join and replace NAs with 0 in Allele columns
	df <- original |>
		full_join(imputed, by = c("POS", "Sample")) |>
		mutate(
			Allele1.x = if_else(is.na(Allele1.x), 0, Allele1.x),
			Allele2.x = if_else(is.na(Allele2.x), 0, Allele2.x),
			Allele1.y = if_else(is.na(Allele1.y), 0, Allele1.y),
			Allele2.y = if_else(is.na(Allele2.y), 0, Allele2.y),
		)

	flags <- df |>
		select(POS, Sample) |>
		left_join(tracts, join_by(Sample, POS >= left, POS < right)) |>
		summarise(.by = c("Sample", "POS"), is_NEA = !any(is.na(left)))

	result <- df |>
		inner_join(flags)

	q <- result |>
		mutate(
			matching_allele1 = Allele1.x == Allele1.y,
			matching_allele2 = Allele2.x == Allele2.y,
			matching_alleles = if_else(matching_allele1, 1, 0) +
				if_else(matching_allele2, 1, 0),
			matching_pair = matching_alleles == 2,
			homo_orig = Allele1.x == Allele2.x,
			homo_imp = Allele1.y == Allele2.y,
			matching_homozygous = homo_orig & homo_imp,
			matching_heterozygous = !homo_orig & !homo_imp,
			matching_alleles_NEA = if_else(is_NEA, matching_alleles, 0),
			matching_pair_NEA = if_else(is_NEA, matching_pair, FALSE),
			homo_orig_NEA = if_else(is_NEA, homo_orig, FALSE),
			homo_imp_NEA = if_else(is_NEA, homo_imp, FALSE),
			matching_homozygous_NEA = if_else(
				is_NEA,
				matching_homozygous,
				FALSE
			),
			matching_heterozygous_NEA = if_else(
				is_NEA,
				matching_heterozygous,
				FALSE
			),
		)

	sums <- q |>
		summarize(
			total_variants = n(),
			matching_alleles = sum(matching_alleles),
			matching_pair = sum(matching_pair),
			homo_orig = sum(homo_orig),
			homo_imp = sum(homo_imp),
			matching_homozygous = sum(matching_homozygous),
			matching_heterozygous = sum(matching_heterozygous),
			total_variants_NEA = sum(is_NEA),
			matching_alleles_NEA = sum(matching_alleles_NEA),
			matching_pair_NEA = sum(matching_pair_NEA),
			homo_orig_NEA = sum(homo_orig_NEA),
			homo_imp_NEA = sum(homo_imp_NEA),
			matching_homozygous_NEA = sum(matching_homozygous_NEA),
			matching_heterozygous_NEA = sum(matching_heterozygous_NEA),
			.by = "Sample"
		)

	finale <- sums |>
		mutate(
			rate = !!rate,
			masking_rate = rate / 10,
			allele_accuracy = matching_alleles / (2 * total_variants),
			pair_accuracy = matching_pair / total_variants,
			allele_accuracy_NEA = matching_alleles_NEA /
				(2 * total_variants_NEA),
			pair_accuracy_NEA = matching_pair_NEA / total_variants_NEA,
			total_variants_EUR = total_variants - total_variants_NEA,
			matching_alleles_EUR = matching_alleles - matching_alleles_NEA,
			allele_accuracy_EUR = matching_alleles_EUR /
				(total_variants_EUR * 2)
		)

	finale
}

# Process all rates and combine results
results <- map_dfr(1:9, process_rate) |>
	inner_join(sample_ages) |>
	inner_join(nea) |>
	inner_join(ser)

# Save combined results
results |> compute_csv("out/result.csv")
