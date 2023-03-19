## **ASR curation overview**
<figure markdown>
  ![Image title](images/logo.png){ width="400" }
</figure>


ASR curation is a snakemake pipeline for annotating sequences to be used in ancestral sequence reconstruction.

It implements a re-runnable workflow of all the basic steps of a reconstruction -

- sequence alignment 
- tree inference
- ancestral prediction

The focus is on building up a series of annotation files for all of the sequences in your dataset.

The columns in these alignment files are used to subset the data by defining a series of rules including or excluding sequences.

The annotations can be further used to visualise your data on phylogenetic trees and within sequence alignments.

The focus on subset rule files means that you can -

- create many different subsets at once and 
- have a list of understandable and reproducible rules for curating your data 

## **Documentation overview**

- **Installation and quickstart** 
     - *Installation -* Shows how to install `asr_curation` locally and what files and commands run the pipeline. 
     - *Running example data -* Explains how to run the example data to verify the pipeline is working correctly 
     - *Explanation of output -* Shows what each of the output folders contain after running the example data
- **Detailed explanation** 
     - *Basic snakemake -* A quick introduction to key snakemake concepts and how to read the steps defined in the `snakefile`  
     - *Explanation of steps -* A detailed rundown of what each step within the `asr_curation` pipeline is doing 
- **Customisation** 
     - *Defining files -* How to define the `config` and `subset` files so that you can run the `asr_curation` pipeline on your data 
     - *Customising workflows -* How to customise the workflow by adding custom annotations to be added to each dataset

## **Cool features**