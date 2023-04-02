## **Snakefile** 
The main workflow is contained within the `snakefile` file in the top folder of the `asr_curation` repository.

A `snakefile` contains a series of rules that define -

- an input
- an output
- how to get to the output from the input, typically either a Python script or a shell command

## **Rules in snakefiles**
For example, this rule defines -  

- the input as an alignment, 
- the output as a tree file (in Newick format) 
- a shell command to run in order to generate this - in this case calling `FastTree`

```py
rule infer_tree:
    input:
        WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.aln"

    output:
        WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.nwk"

    shell:
        "FastTree {input} > {output}"
```

Anything defined within curly brackets `{}` is a wildcard. 

The wildcards in the shell command let us use the whatever we define as `input` and `output` within the shell command.

The wildcards in the input and output rules let as define different `dataset` and `subset` wildcards which makes it possible to run the
entire snakemake pipeline on different sets of data.

## **Triggering a snakemake pipeline**
Once you trigger a snakemake pipeline, it will look for the first command defined in the `snakefile` file and try to create the output.

By convention, this is a rule called `all` which we place at the top of the `snakefile`

In `asr_curation` the following rule is defined as the first rule -

```py
rule all:
        input:
            annotations = [f'{WORKDIR}/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_annotations.txt' for dataset in DATASETS for subset in subsets[dataset]],
            ancestors = [f'{WORKDIR}/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_ancestors.csv' for dataset in DATASETS for subset in subsets[dataset]]
```

Snakemake will look to see if these `input` files already exists for each `dataset` and `subset` in terms of both the `annotations.txt` and the `ancestors.csv`

If these files do not exist, snakemake finds a rule defined in the `snakefile` that it would need to run in order to generate these files

Ignoring the other detail, focus here on the output line - this rule will generate the output of `ancestors.csv` which is one of the files we needed as input to the `all` rule.

If this file didn't already exist, then we know we must run this rule in order to generate it!

```python
rule add_annotations_from_ancestors:
    input:
        csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_alignment.csv",
        aln = WORKDIR + "/{dataset}/subsets/{subset}/grasp_results/GRASP_ancestors.fa",
    output:
        csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_ancestors.csv"
    script:
        CUSTOM_ANCESTOR_DIR + "/add_annotations_from_ancestors_for_use_in_example_workflow.py"
```

## **Snakemake pipelines only rerun what is needed**
Snakemake backtracks in this way, checking to see which rules need to be run in order to generate the required input for the `all` rule.

If no data at all has been generated it will start with the very first step in the workflow, but if some steps have been completed then
it does not need to regenerate them. This is an important concept in workflow management and saves us a lot of time.

For example, in the `asr_curation` pipeline, all of the subset files use the `annotation_files` generated from going to the UniProt and BRENDA databases.
So only one call to the databases needs to be made, and then all downstream analysis simply uses the generated output from the annotation steps.


## **A closer look at snakemake rules**

Here are two `snakemake` rules. We can compare the differences in order to get a better sense of how to read `snakemake` rules in the `snakefile`

```py
rule validate_ids:
   input:
       WORKDIR + "/{dataset}/csv/original/{dataset}_original.csv"
   output:
       WORKDIR + "/{dataset}/csv/validated/{dataset}_validated.csv"
   script:
       "scripts2/validate_ids.py"

```


```py
rule infer_tree:
    input:
        WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.aln"

    output:
        WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.nwk"

    shell:
        "FastTree {input} > {output}"
```

- `validate_ids` only uses the `dataset` wildcard as it is only run once per `dataset`
- `infer_tree` uses both the `dataset` and the `subset` wildcard as it is run once per `subset`
- Both rules use the wildcards to name the output files - for example, in the example workflow, the phylogenetic tree generated for the dataset `kari_example_ec_1_1_1_86` and subset `eukaryotic_kari` is written to `/{dataset}/subsets/{subset}/{dataset}_{subset}.nwk` or `/kari_example_ec_1_1_1_86/subsets/eukaryotic_kari/kari_example_ec_1_1_1_86_eukaryotic_kari.nwk`. This means that files can easily be identified and will not overwrite each other. 
- `validate_ids` defines a `script` command - which tells `snakemake` to run the Python script `scripts/validate_ids.py`. You can look in the scripts folder to see exactly what is being run for each step.
- `infer_tree` defines a `shell` commmand - which tells `snakemake` to run the defined command on the command line.
