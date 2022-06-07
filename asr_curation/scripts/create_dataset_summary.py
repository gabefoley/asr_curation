import pandas as pd

# Read in the dataframe
df = pd.read_csv(snakemake.input[0])

with open(snakemake.output.summary, "w+") as summary:
    summary.write("# Markdown file \n")
    summary.write("This is the summary file")
