# Customising annotations

## **Provided annotation functions**

The following is a description of some of the provided annotation functions and how to incorporate them.

All changes should be made in your `add_custom_annotations.py` file and you should add the custom location of this file
to you `config` file.


## Create a top column

Because it can be useful for generating rules that 

annot_df = an.create_top_column(annot_df, 80)

# Separate out the note values from the feature columns

annot_df = an.separate_notes(annot_df)

# Add embeddings

Generate embeddings using TM Vec.

annot_df = tm_vec_embed.process_and_store_embeddings(annot_df, 'Prot_T5')

# Generate DBSCAN images

Generate a DBSCAN coverage image

generate_dbscan_coverage(df, 'Prot_T5 Embed Encoded', snakemake.input.custom_dir + "/dbscan_coverage.png" )

## Writing code for new annotation scripts


You must keep the structure for the input and output the same but to help with this all three files define the following line -

```
# IF CUSTOMISING THIS FILE PLACE CUSTOM CODE HERE
```

An example of how to write custom code is given in the `scripts/configs/example_custom_annotations` folder, and the example_workflow uses these files to add custom annotations to the datasets in the example_workflow

There is a lot of code in `scripts/annot_functions` that is designed to do common annotation tasks.
