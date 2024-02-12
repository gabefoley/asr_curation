# Customising workflows

## **Customising the annotations scripts**

The following rules / Python files stored in scripts are run at three different steps during the `asr_curation` pipeline 

- **add_custom_annotations / scripts/add_custom_annotations.py'** - run after all of the database annotations are retrived and before subsets are made

- **add_annotations_from_alignment / scripts/add_annotations_from_alignment.py'** - run after the sequences for a given subset have been aligned

- **add_annotations_from_ancestors / scripts/add_annotations_from_ancestors.py'** - run after the ancestors for a given subset have been predicted

Which script to use for which custom annotation is a decision as to what information you require to build the rules. For example, if your custom rules make use of aligned positions, this will obviously need to come after alignments have been generated.

Currently only `add_custom_annotations` comes before the creation of subsets, so only annotations added from this custom script (or annotations from the database retrievals) can be used to make subsets. 

The files stored in `scripts` are template files - they will read in the data and write it out correctly, but not add any annotations.

You can override these files by copying them to a new location and writing code to be annotated into your data. 

Within your config file you then specify the location of any custom script. 

See [Defining the config files](defining_files.md#defining-the-config-files) for more information.


## Writing code for new annotation scripts


You must keep the structure for the input and output the same but to help with this all three files define the following line -

```
# IF CUSTOMISING THIS FILE PLACE CUSTOM CODE HERE
```

An example of how to write custom code is given in the `scripts/configs/example_custom_annotations` folder, and the example_workflow uses these files to add custom annotations to the datasets in the example_workflow

There is a lot of code in `scripts/annot_functions` that is designed to do common annotation tasks.


# Customising annotations

## **Provided annotation functions**

The following is a description of some of the provided annotation functions and how to incorporate them.

All changes should be made in your `add_custom_annotations.py` file and you should add the custom location of this file
to you `config` file.


## Create a top column

Creates a new column where if a value within a column is more than a given percentage, it will create a new column,
labelled TOP_<column_name>_<value_name> with `True` or `False` values depending on whether a sequence has
that annotation.

This can be useful for grouping together multiple annotations that differ from the dominant value in a column.

`annot_df = an.create_top_column(annot_df, 80)`

# Separate out the note values from the feature columns

annot_df = an.separate_notes(annot_df)


## Generating embeddings and running DBSCAN

# Add embeddings

Generate embeddings using TM Vec. You need to download the checkpoint and config.json and point to their path in your `add_custom_annotations.py` code


To download the files -
`!wget https://users.flatironinstitute.org/thamamsy/public_www/tm_vec_cath_model.ckpt -q gwpy\n"`
`!wget https://users.flatironinstitute.org/thamamsy/public_www/tm_vec_cath_model_params.json -q gwpy"`

Then set their path in your file - 

```
model_checkpoint_path = <path_to_ckpt>
model_config_path = <path_to_json>
```

And then create the embeddings. This will store the generated embeddings in a local file, so that embeddings can be reused over subsequent runs of asr_curation.

The default path is in the custom_annotations folder, in a pickle file - "embeddings.pkl", but this can be changed with the parameter `embedding_df_path` when using the following
`process_and_store_embeddings` function

annot_df = tm_vec_embed.process_and_store_embeddings(annot_df, 'Prot_T5', model_checkpoint_path, model_config_path)

# Generate DBSCAN images

Generate a DBSCAN coverage image

```
db.generate_dbscan_coverage(annot_df, 'Prot_T5 Embed Encoded', 
	f"{snakemake.input.custom_dir}/{snakemake.wildcards.dataset}_dbscan_coverage" )
```


Here we provide the `dataframe` name, the name of the column to find the embeddings in (by default 'Prot_T5 Embed Encoded' if you are using Prot_T5 in the previous step) and a prefix for the output files


You can also pass columns that you do not wish to include in the DBSCAN analysis using the parameter `skip_cols`
```
db.generate_dbscan_coverage(annot_df, 'Prot_T5 Embed Encoded', 
	f"{snakemake.input.custom_dir}/{snakemake.wildcards.dataset}_dbscan_coverage_no_ft",
	skip_cols = ['ft_var_seq||', 'ft_variant||', 'ft_conflict||', 'ft_chain||', 'ft_crosslnk||', 'ft_carbohyd||', 'ft_init_met||' , 'ft_mod_res||', 'ft_lipid||', 'ft_transit||', 'ft_compbias||', 'ft_domain||', 'ft_motif||',  'ft_region||', 'ft_repeat||', 'ft_zn_fing||', 'ft_binding||', 'ft_topo_dom||', 'ft_act_site||'])

```

# Full minimal example of generating embeddings and running DBSCAN outlier detection

The following can be created as `add_custom_annotations.py` to generate embeddings. Make sure to specify that you are using a custom annotations file in your `config` file 
and to update the paths to the TM-Vec checkpoint / config.

```python
import os
import annot_functions as an
import get_funfams as ff
import map_to_cdd as m2c
import seqcurate as sc
import pandas as pd
import numpy as np
import add_embeddings as embed
import add_tm_vec_embeddings as tm_vec_embed
import create_itol_files as itol
import create_dbscan_coverage as db


annot_df = pd.read_csv(snakemake.input.csv)

# Create TOP column for high scoring values

annot_df = an.create_top_column(annot_df, 80)

# Separate out the note values from the feature columns

annot_df = an.separate_notes(annot_df)

# Add embeddings

model_checkpoint_path = <path_to_checkpoint>
model_config_path = <path_to_config>

annot_df = tm_vec_embed.process_and_store_embeddings(annot_df, 'Prot_T5', model_checkpoint_path, model_config_path)

# Generate DBSCAN images
db.generate_dbscan_coverage(annot_df, 'Prot_T5 Embed Encoded', 
	f"{snakemake.input.custom_dir}/{snakemake.wildcards.dataset}_dbscan_coverage" )


db.generate_dbscan_coverage(annot_df, 'Prot_T5 Embed Encoded', 
	f"{snakemake.input.custom_dir}/{snakemake.wildcards.dataset}_dbscan_coverage_no_ft",
	skip_cols = ['ft_var_seq||', 'ft_variant||', 'ft_conflict||', 'ft_chain||', 'ft_crosslnk||', 'ft_carbohyd||', 'ft_init_met||' , 'ft_mod_res||', 'ft_lipid||', 'ft_transit||', 'ft_compbias||', 'ft_domain||', 'ft_motif||',  'ft_region||', 'ft_repeat||', 'ft_zn_fing||', 'ft_binding||', 'ft_topo_dom||', 'ft_act_site||'])

annot_df.to_csv(snakemake.output[0], index=False)

```