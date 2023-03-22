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
