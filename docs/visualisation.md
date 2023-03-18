# ASR curation

ASR curation is a snakemake pipeline for annotating sequences to be used in ancestral sequence reconstruction.


## Super quick overview

* Input is a FASTA file of all the sequences you're interested in
* ASR curation goes through all of the basic steps of setting up a reconstruction - sequence alignment, tree inference, ancestral prediction
* First it retives annotations from the UniProt and BRENDA databases
* These annotations can be used to build rules to create subsets of your data
* The annotations can be further used to visualise your data on phylogenetic trees
* This becomes an iterative process whereby you build trees, look at how data is annotated, create rules for excluding sequences based on these annotations
* The focus on subset rule files means that you can 1) create many different subsets at once and 2) have a list of understandable and reproducible rules for curating your data 

## Installing ASR curation

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.

## Cool features