1. Generate the datasets with the notebook in /notebooks/BRENDA_Analysis/1_Create_input_files.ipynb

This allows you select which columns you wish to set a cutoff for (brenda_cols), and what that cutoff should be (cutoff).
The final cell puts the files where they need to be /workflows/brenda_workflow/fasta and /workflows/brenda_workflow/subset_rules

2. Run the snakemake workflow using the brenda_config config file.

snakemake  --cores 1 --configfile brenda_config.yaml

The brenda config file tells snakemake to just use the folders determined in the brenda_config.yaml file, defined at the top of the file

If you want to change this, change the top 3 variables -
workdir
fastadir
subdir

Note you'll also need to change where 1_Create_input_files.ipynb is writing out to.


If you want to make just brenda annotations (not alignments, trees, and ancestors)

Comment out line 109 'tree_images ='

3. Look at the datasets with the notebook in  /notebooks/BRENDA_Analysis/2_Check_annotations.ipynb

You can change fasta_name = 'uniprot_ec_3_5_2_6' at the top of the notebook and it'll change the file you're looking at.

