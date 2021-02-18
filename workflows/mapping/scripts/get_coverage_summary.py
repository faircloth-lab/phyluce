import os
import pandas as pd

# import pdb

temp_summary_across_all = []
summary_output = snakemake.output[-1]
species_output = snakemake.output[:-1]
for coverage_file, species_output_file in zip(snakemake.input, species_output):
    taxon_name = os.path.basename(coverage_file).split(".")[0]
    cov = pd.read_table(
        coverage_file,
        names=["contig", "pos", "cov"],
        dtype={"pos": "Int64", "cov": "Int64"},
    )
    # summary stats in row form
    temp_sum = cov.describe().transpose().iloc[[1]]
    temp_sum.insert(0, "taxon", taxon_name)
    temp_summary_across_all.append(temp_sum)
    # this gets coverage by contig for each organism
    taxon_summary = cov.groupby("contig").describe()["cov"]
    taxon_summary.rename(columns={"count": "base_pairs"}, inplace=True)
    with open(species_output_file, "w") as outf:
        taxon_summary.to_csv(outf, header=True)
summary_across_all = pd.concat(temp_summary_across_all)
summary_across_all.rename(columns={"count": "total_base_pairs"}, inplace=True)
# now write the summary output file
with open(summary_output, "w") as outf:
    summary_across_all.to_csv(outf, header=True, index=False)
