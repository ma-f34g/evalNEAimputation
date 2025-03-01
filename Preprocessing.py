import polars as pl
from tqdm import tqdm

SAMPLES = 1980
eur_cols = [f"EUR_{i}" for i in range(1, SAMPLES + 1)]


def files(rate):
    return [
        (
            "data/test_samples.vcf.gz",
            f"data/out/masked_samples_{rate}.parquet",
            5,
            f"data/damaged_test_samples_0.{rate}.vcf.gz",
        ),
        (
            f"data/imputed_0.{rate}_filtered.vcf.gz",
            f"data/out/imputed_filtered_{rate}.parquet",
            12,
            None,
        ),
    ]


def rewrite(input_path, output_path, skip_rows, mask_file=None):
    df = pl.read_csv(
        input_path, separator="\t", skip_rows=skip_rows, columns=["POS"] + eur_cols
    )
    if mask_file:
        mask = pl.read_csv(
            mask_file, separator="\t", skip_rows=skip_rows, columns=["POS"]
        )
        df = df.join(mask, on="POS", how="anti")
    df = df.with_columns(pl.col("POS").cast(pl.UInt32()))

    df = df.unpivot(
        index="POS",
        on=eur_cols,
        variable_name="Sample",
        value_name="Genotype",
    )

    df = df.with_columns(
        pl.col("Genotype").str.slice(0, 1).alias("Allele1").cast(pl.Int8()),
        pl.col("Genotype").str.slice(2, 1).alias("Allele2").cast(pl.Int8()),
    ).drop("Genotype")
    df.write_parquet(output_path)


for input_path, output_path, skip_rows, mask_file in tqdm(
    sum([files(rate) for rate in range(1, 10)], [])
):
    rewrite(input_path, output_path, skip_rows, mask_file)
rewrite(
    "data/test_samples.vcf.gz",
    "data/out/test_samples.parquet",
    5,
)

tracts = pl.read_parquet(
    "data/neanderthal_tracts.parquet", columns=["name", "left", "right"]
)
tract_intervals = tracts.select(
    pl.col("name").alias("Sample"),
    pl.col("left").cast(pl.UInt32),
    pl.col("right").cast(pl.UInt32),
)
tract_intervals.write_parquet("data/out/tract_intervals.parquet")
