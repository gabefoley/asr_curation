![Coverage](coverage.svg)
[![Twitter](https://badgen.net/badge/icon/gabefoley?icon=twitter&label)](https://twitter.com/gabefoley)
# ASR Curation workflow

<p align="center">
	<img src="docs/images/logo.png" alt="ASR Curation" width="400" />

</p>
 
Snakemake pipeline for annotating sequences to be used in ancestral sequence reconstruction.



# Documentation

[Read the full documentation here](http://gabefoley.github.io/asr_curation)

# Basic concept

A common task in phylogenetics and ancestral sequence reconstruction to have a large set of data you are interested in with the need to curate this data to include only relevant sequences. 

This pipeline starts by allowing the user to submit a large set of sequences (in a single FASTA file) and then to also create a set of 'rules' by which the data should be split up - which sequences to include in downstream phylogenetic analyses, based on the metadata that is retrieved, by this pipeline, from the UniProt and BRENDA databases.

As phylogenetic analysis is an iterative process that benefits from a deep understanding of the underlying sequences, these annotations can be viewed in interactive notebooks generated by Jupyter-book, and further sets of rules can be created in .subset files to create alternative subsets of the data.

Because this is executed within a snakemake pipeline, it has the added advantage of keeping all the iterations of subsets available and the rules that exclude sequences clearly defined and therefore entirely reproducible. 


# Install instructions 


1. Clone this repository to your desktop

```
git clone https://github.com/gabefoley/asr_curation.git
```


2. Create a conda environment

```
conda create -n asr_curation python=3.9
```

3. Activate the conda environment

```
conda activate asr_curation
```

4. Install the required Python packages

```
pip install -r requirements.txt
```


5. Install the following so that they are callable from the command line
- [mafft](https://mafft.cbrc.jp/alignment/software/) - callable as `mafft`
- [FastTree](http://www.microbesonline.org/fasttree/) - callable as `FastTree`
- [GRASP](https://bodenlab.github.io/GRASP-suite/project/graspcmd/) - callable as `grasp`


Optional (for viewing trees with generated annotation files)

- [FigTree](http://tree.bio.ed.ac.uk/software/figtree/)


