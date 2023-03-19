# **Installation and quick start**


## **Installation guide**


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


    1. [mafft](https://mafft.cbrc.jp/alignment/software/) - callable as `mafft`

    2.  [FastTree](http://www.microbesonline.org/fasttree/) - callable as `FastTree`

    3. [GRASP](https://bodenlab.github.io/GRASP-suite/project/graspcmd/) - callable as `grasp`
   
    4. [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) (Optional - for viewing trees with generated annotation files)

## **Quick start**

To run an `asr_curation` pipeline you need to define -

- a `config.yaml` file that defines where to look for the `fasta` and `subset` files. See [Defining the config files](defining_files.md#defining-the-config-files) for further details.
)
- a `fasta` file containing all the sequences you wish to include in the pipeline
- a `subset` file containing rules for subsetting the data. See [Defining the subset files](defining_files.md#defining-the-subset-files) for further details. 


You can then run snakemake directly by calling it from within the `asr_curation` environment we created and activated earlier.

You need to specify the number of cores to run `snakemake` with and point to the `config.yaml` file.
```
snakemake --cores 1 --configfile ./asr_curation/config/example_config.yaml
```