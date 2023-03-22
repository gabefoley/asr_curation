# Defining input files

## **Defining the config files**

Every snakemake pipeline needs to be run with a specific config file - a YAML file with a .yaml extension.

For example this call uses the `--configfile` flag to set the location to the `example_config.yaml` file stored in the asr_curation repository


```
snakemake --cores 1 --configfile ./asr_curation/config/example_config.yaml
```


*A config file must define the following parameters*

- **workdir** - where the `datasets` folder and all of the created data will be stored

- **fastadir** - where the fasta files are stored

- **subdir** - where the `subset` files are stored - they must match the name of the .fasta file they refer to but contain a .subset extension

- **annotation_cols** - this specifies which columns in the annotation file should be written to the FigTree annotation file

- **blocked_datasets** - a list of datasets that will not be run - meaning you can keep the FASTA and .subset file in the fastadir and subdir and just skip over the generation of these files. An empty list `[]` will run every fasta and matching .subset file

*A config file can contain the following optional parameters*

These define the custom annotations that can be added. 

See [Customising the annotations scripts](customising_workflow.md#customising-the-annotations-scripts) for more details

- **custom_dir** - if you want to use a custom `add_custom_annotations.py` file then set the directory location here and place a custom file name `add_custom_annotations.py` into this folder
- **custom_align_dir** - if you want to use a custom `add_annotations_from_alignment.py` file then set the directory location here and place a custom file name `add_annotations_from_alignment.py` into this folder
- **custom_align_dir** - if you want to use a custom `add_annotations_from_ancestors.py` file then set the directory location here and place a custom file name `add_annotations_from_ancestors.py` into this folder

## Example config file

Here is an example config file from the example_workflow in the `asr_curation` repository


```
# Main working directory
workdir: "workflows/example_workflow"

# Store all FASTA files here
fastadir: "workflows/example_workflow/fasta"

# Each FASTA file needs an according .subset file in this folder (test1.fasta -> test1.subset)
subdir: "workflows/example_workflow/subset_rules"

# Which columns should we plot in the summary documents?
annotation_cols: ['lineage_superkingdom', 'xref_supfam', 'xref_panther', 'KARI_Class', 'ec']

# Directories for custom annotations
custom_dir : "scripts/example_custom_annotations"
custom_align_dir: "scripts/example_custom_annotations"


# Block a data set from being run (you can still keep the FASTA and .subset file in the fastadir and subdir)
blocked_datasets : []
```


## **Defining the subset files**

Subset files need to be in the format 

`<subset_name> = <column_name> : <string> $ <other_column_name> : NOT <string>`

Where `subset_name` can be any name you wish to call the subset, `column_name` and `other_column_name` are columns that appear in the final annotation file (see /csv/custom/*.csv) and `string` is a term that you wish to include or exclude.

For each row / sequence in your annotation file, if the given column contains the string (or substring) then it will be included.

If you use the `NOT` modifier, for each row / sequence in your annotation file, if the given column contains the string (or substring) then it will not be included.

You can define multiple subsets within the one file.



```
# THIS IS A COMMENT LINE.

eukaryotic_kari = lineage_superkingdom : Eukaryota $ protein_name : NOT Deleted, Merged $ Non_AA_Character : FALSE

classI_kari = KARI_Class : Class_1 $ protein_name : NOT Deleted, Merged $ Non_AA_Character : FALSE

all_kari = *
```


## Including all sequences in a subset
The default way to include every row / sequence is 

`<subset_name> = *`

For example,

`all = *`  

