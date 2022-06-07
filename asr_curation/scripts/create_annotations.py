import seqcurate as sc

print("Creating the initial annotation file")

# Create the sequence dataframe and save it
annotations = sc.get_sequence_df(snakemake.input[0])
annotations.to_csv(snakemake.output[0], index=False)
