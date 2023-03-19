## **Running example data**
The repository contains example data that can be inspected and rerun.

You can see the config file in `config/example_config.yaml`

The first few lines of this config file define where data will be stored and where the required input files are

`config/example_config.yaml`
```
# Main working directory
workdir: "workflows/example_workflow"

# Store all FASTA files here
fastadir: "workflows/example_workflow/fasta"

# Each FASTA file needs an according .subset file in this folder (test1.fasta -> test1.subset)
subdir: "workflows/example_workflow/subset_rules"
```

All of these are directories within the main `asr_curation` repository.

## **FASTA directory**

You can see that the FASTA directory `fastadir` at `workflows/example_workflow/fasta` contains two FASTA files 

- als_example_ec_2_2_1_6.fasta
- kari_example_ec_1_1_1_86.fasta

## **Subset directory** 
And the subset directory `subdir` at `workflows/example_workflow/subset_rules` contains a matching set of subset_rules for each dataset


- als_example_ec_2_2_1_6.subset
- kari_example_ec_1_1_1_86.subset

The ALS subset file contains a single subset -

```
als_interpro_IPR012782 = Non_AA_Character : False $ protein_name : NOT Deleted, Merged $ xref_interpro : IPR012782

```

And the KARI subset file contains three subsets - 

One of these will only contain Eukaryotic KARI sequences, one will only contains Class I KARI sequences, and one will 
contain all of the KARI sequences

```
eukaryotic_kari = lineage_superkingdom : Eukaryota $ protein_name : NOT Deleted, Merged $ Non_AA_Character : FALSE

classI_kari = KARI_Class : Class_1 $ protein_name : NOT Deleted, Merged $ Non_AA_Character : FALSE

all_kari = *
```

Therefore if we look at the generated data, we will see -




## **Running it via the tests**

You can validate that everything is working correctly by running a small example workflow locally. 

There is a test within `test/test_snakemake.py` called `test_snakemake_pipeline`

This file only contains a small number of sequences so that it can be run quickly.


You can run this test by - 

`` conda activate asr_curation``

`` pytest test/test_snakemake.py``

This will run the entire `snakemake` pipeline from end to end, including any database retrieval.

The output folders stay after the test is run so that you can inspect them, but are deleted every time the test is rerun.

See `test_snakemake_pipeline`  in `test/test_snakemake.py` for more details about this test and `test/files/config/test_config.yaml`
for the full details of where output folders, fasta directories, and subset directories are set.


There is also another test within `test/test_snakemake.py` called `test_snakemake_pipeline_with_existing_annotation_file` 
that does the same test without deleting the original annotations file - meaning that no calls to external databases need to be made
and only the subset generation, alignment, tree inference, and ancestor prediction is rerun - making this even quicker to test the
pipeline.

## **Rerunning the entire example_workflow**

You can also delete the entire output of the example_workflow locally and then run it again to check that you can
regenerate the data.

Note that these files contain substantially more data than the tests, so will take longer to run.

Remove the `datasets` directory stored in `workflows/example_workflow`. Make sure to keep the `fastadir` and `subdir` directories
```
rm -rf workflows/example_workflow/datasets/
```

Rerun the snakemake pipeline to re-generate the data

!!! warning
    Two important points here - 
    This will retrigger calls to the UniProt and BRENDA databases, so the annotations may not be identical to those stored in this repository as these are live databases. This should be treated as just a test to check that it is working. 
    Any changes you make to this files are excluded from the Git repository - so a fresh git pull will overwrite any regeneration of data or deleting of folders that you do.