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

Creates a new column 

annot_df = an.create_top_column(annot_df, 80)

# Separate out the note values from the feature columns

annot_df = an.separate_notes(annot_df)

# Add embeddings

Generate embeddings using TM Vec.

annot_df = tm_vec_embed.process_and_store_embeddings(annot_df, 'Prot_T5')

# Generate DBSCAN images

Generate a DBSCAN coverage image

generate_dbscan_coverage(annot_df, 'Prot_T5 Embed Encoded', f"{snakemake.input.custom_dir}/{snakemake.wildcards.dataset}_dbscan_coverage")

Here we provide the `dataframe` name, the name of the column to find the embeddings in (by default 'Prot_T5 Embed Encoded' if you are using Prot_T5 in the previous step) and a prefix for the output files

Change to docs

