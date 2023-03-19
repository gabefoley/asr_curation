
with open(snakemake.input[0], "r") as in_file:
    buf = in_file.readlines()

with open(snakemake.output[0], "w+") as out_file:

    for line in buf:
        if line.strip() == '"snakemake-job-properties"':
            line += ',\n\t"remove-cell"\n'
        out_file.write(line)
