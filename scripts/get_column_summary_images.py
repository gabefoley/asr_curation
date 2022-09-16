import matplotlib.pyplot as plt

import pandas as pd


# Read in the dataframe
df = pd.read_csv(snakemake.input[0])

# Define some BRENDA columns to skip
skip_brenda_cols = [
    "BRENDA_CL",
    "BRENDA_GI",
    "BRENDA_SN",
    "BRENDA_SY",
    "BRENDA_MW",
    "BRENDA_SP",
    "BRENDA_NSP",
    "BRENDA_PM",
    "BRENDA_LO",
    "BRENDA_SU",
    "BRENDA_PU",
    "BRENDA_ST",
    "BRENDA_CR",
    "BRENDA_CF",
    "BRENDA_RN",
    "BRENDA_RT",
    "BRENDA_ME",
    "BRENDA_REFERENCES",
]


# Only get the columns that start with BRENDA and that we don't want to skip
brenda_cols = [
    x
    for x in df
    if x.startswith("BRENDA")
    and "COMMENT" not in x
    and "REFS" not in x
    and "UNITS" not in x
    and not any(skip in x for skip in skip_brenda_cols)
]

# Drop rows if they don't have at least one entry in one of the BRENDA columns
b_df = df[brenda_cols].dropna(thresh=1)

# Create and save a plot of the BRENDA column counts

fig, ax = plt.subplots(figsize=(20, len(brenda_cols) / 3))
plot = b_df.count().sort_values(ascending=False).plot.barh(ax=ax)
fig_for_save = plot.get_figure()

print(snakemake.output.img)
# fig.savefig(f"{snakemake.output.img}")
fig_for_save.savefig(f"{snakemake.output.img}", bbox_inches="tight")
