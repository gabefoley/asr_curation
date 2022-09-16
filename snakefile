import glob
from collections import defaultdict


# Collect config info from config file
WORKDIR = config['workdir'] + "/datasets"
FASTADIR = config['fastadir']
SUBDIR = config['subdir']
VERBOSE = True

# Collect config info from parameters passed from command line
if 'BRENDA_RUN' in config.keys():
    BRENDA_RUN = config['BRENDA_RUN']
else:
    BRENDA_RUN = False



DATASETS = expand(os.path.basename(x).split('.')[0] for x in glob.glob(FASTADIR + "/*.fasta") if os.path.basename(x).split('.')[0] not in config['blocked_datasets'])

if VERBOSE:

    print ("Running ASR curation pipeline with the following datasets:")
    print (DATASETS)

# for blocked in config['blocked_datasets']:
#     if blocked in DATASETS:
#         DATASETS.remove(blocked)


def get_col_val_name(col_val_dict):
    name =  "_".join(str(k) + "_" +  str(v) for k,v in col_val_dict.items())

    return name


def get_subset_names(subset_rules_dir):
    subset_dict = defaultdict(list)

    subsets = expand(os.path.basename(x).split('.')[0] for x in glob.glob(SUBDIR + "/*.subset"))



    for subset in subsets:

        print ('subset is ')
        print (subset)

        # If a subset rule doesn't exist, then just write out the full data set
        if not open(f"{SUBDIR}/{subset}.subset").read().splitlines():

            subset_dict[subset].append('full_set')


        for line in open(f"{SUBDIR}/{subset}.subset").read().splitlines():

            if line and not line.startswith("#"):

                print ('line is ')
                print (line)


                # Check to see if custom name exists
                name = line.split("=")[0].strip()

                # If no custom name make name based on values of dictionary
                if len(name) < 3:
                    name = get_col_val_name(col_val_dict)

                subset_dict[subset].append(name)
    print ('subset dict is')

    print (subset_dict)
    return subset_dict

def get_ancestor_col_dict(ancestor_cols):
    ancestor_col_dict = {}
    for col_val in ancestor_cols:
        print (col_val)
        col = col_val.split(":")[0].strip()
        val = col_val.split(":")[1].strip()
        ancestor_col_dict[col] = val

    print (ancestor_col_dict)
    return ancestor_col_dict


# subsets = get_subset_names(SUBDIR) {'uniprot-ec_1_1_1_86-filtered-reviewed_yes' : ['subset1']}

subsets = get_subset_names(SUBDIR)


print (subsets)

print (DATASETS)


rule all:
        input:
            # files = [f'{WORKDIR}/{dataset}/subsets/{subset}/grasp_results/GRASP_ancestors.fasta' for dataset in DATASETS for subset in subsets[dataset]]
            annotations = [f'{WORKDIR}/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_annotations.txt' for dataset in DATASETS for subset in subsets[dataset]],
            summary_document = [f'{WORKDIR}/{dataset}/dataset_summary/{subset}/_build/html/index.html' for dataset in DATASETS for subset in subsets[dataset]],
            ancestors = [f'{WORKDIR}/{dataset}/subsets/{subset}/grasp_results/GRASP_ancestors.fa' for dataset in DATASETS for subset in subsets[dataset] for column in config['annotation_cols']]

            # expand(WORKDIR + "/{dataset}/grasp_results/GRASP_ancestors.fasta",
            #     dataset=DATASETS)



# Create the initial annotation file from the FASTA file or list of IDs
rule create_annotations:
    input:
        FASTADIR + "/{dataset}.fasta"
    output:
        WORKDIR + "/{dataset}/csv/original/{dataset}_original.csv"
    script:
        "scripts/create_annotations.py"

# Map the IDs to UniProt and NCBI in order to validate the type of IDs they are
rule validate_ids:
   input:
       WORKDIR + "/{dataset}/csv/original/{dataset}_original.csv"
   output:
       WORKDIR + "/{dataset}/csv/validated/{dataset}_validated.csv"
   script:
       "scripts/validate_ids.py"

# Map to UniProt to get all of the known UniProt annotations
rule get_uniprot_annotations:
    input:
        WORKDIR + "/{dataset}/csv/validated/{dataset}_validated.csv"
    output:
        WORKDIR + "/{dataset}/csv/uniprot/{dataset}_uniprot.csv"
    script:
        "scripts/get_uniprot_annotations_pagination.py"

# # # Map the UniProt sequences to IDs we can query IPG with, and then query IPG to add genome information
# rule get_genomes:
#     input:
#         WORKDIR + "/{dataset}/csv/uniprot/{dataset}_uniprot.csv"
#     output:
#         WORKDIR + "/{dataset}/csv/genome/{dataset}_genome.csv"
#         #WORKDIR + "/{dataset}/csv/gtdb_processed/{dataset}_gtdb_processed.csv"
#     script:
#         "scripts/get_genomes.py"

# # Map the Genomes to GTDB
# rule get_gtdb_annotations:
#    input:
#        WORKDIR + "/{dataset}/csv/genome/{dataset}_genome.csv"
#    output:
#        WORKDIR + "/{dataset}/csv/gtdb/{dataset}_gtdb.csv"
#    script:
#        "scripts/get_gtdb_annotations.py"

# #Use GTDB to cross-check CheckM scores and GTDB taxonomy
# rule process_gtdb_annotations:
#    input:
#        WORKDIR + "/{dataset}/csv/gtdb/{dataset}_gtdb.csv"
#    output:
#        WORKDIR + "/{dataset}/csv/gtdb_processed/{dataset}_gtdb_processed.csv"
#    script:
#        "scripts/process_gtdb_annotations.py"


# Map to BRENDA database to get all of the known BRENDA annotations
rule get_brenda_annotations:
    input:
        WORKDIR + "/{dataset}/csv/uniprot/{dataset}_uniprot.csv"
        # WORKDIR + "/{dataset}/csv/gtdb_processed/{dataset}_gtdb_processed.csv"
    output:
        WORKDIR + "/{dataset}/csv/brenda/{dataset}_brenda.csv"
    script:
        "scripts/get_brenda_annotations.py"

# Add any custom annotations
rule add_custom_annotations:
    input:
        WORKDIR + "/{dataset}/csv/brenda/{dataset}_brenda.csv"
    output:
        WORKDIR + "/{dataset}/csv/custom/{dataset}_annotated.csv"
    script:
        "scripts/custom_annotations.py"


# Create the su
rule create_column_summary_images:
    input:
        WORKDIR + "/{dataset}/csv/custom/{dataset}_annotated.csv"
    output:
        img = WORKDIR + "/{dataset}/col_images/{dataset}_brenda.png"
    script:
        "scripts/get_column_summary_images.py"

rule create_subsets:
    input:
        rules=SUBDIR + "/{dataset}.subset",
        csv=WORKDIR + "/{dataset}/csv/custom/{dataset}_annotated.csv"
    output:
        fasta = WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.fasta",
        csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}.csv"
    script:
        "scripts/create_subsets.py"

rule align_seqs:
    input:
        WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.fasta"
    output:
        WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.aln"
    shell:
        "mafft --reorder {input} > {output}"

rule infer_tree:
    input:
        WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.aln"

    output:
        WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.nwk"

    shell:
        "FastTree {input} > {output}"

rule run_grasp:
    input:
        aln= WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.aln",
        tree= WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.nwk"

    output:
        dir= directory(WORKDIR + "/{dataset}/subsets/{subset}/grasp_results/"),
        aln= WORKDIR + "/{dataset}/subsets/{subset}/grasp_results/GRASP_ancestors.fa",
        tree = WORKDIR + "/{dataset}/subsets/{subset}/grasp_results/GRASP_ancestors.nwk"

    shell:
        "grasp -a {input.aln} -n {input.tree} -s LG -o {output.dir} -i BEP -j --save-as FASTA TREE -pre GRASP -t 2"

rule add_annotations_from_alignment:
    input:
        aln= WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.aln",
        csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}.csv"
    output:
        csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_alignment.csv"

    script:
        "scripts/add_annotations_from_alignment.py"

# rule add_annotations_from_ancestors:
#     input:
#         csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_alignment.csv",
#         aln = WORKDIR + "/{dataset}/subsets/{subset}/grasp_results/GRASP_ancestors.fasta",
#     output:
#         csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_ancestors.csv"
#     script:
#         "scripts/add_annotations_from_ancestors.py"


rule create_annotation_file:
    input:
        csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_alignment.csv"

    params:
        annotation_cols = config['annotation_cols']
    output:
        tsv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_annotations.txt"
    script:
        "scripts/create_annotation_file.py"



# Rules for creating the PDF summary files

rule create_dataset_summary:
    input:
       WORKDIR + "/{dataset}/csv/custom/{dataset}_annotated.csv",
    params:
        annotation_cols = config['annotation_cols']
    output:
       summary = WORKDIR + "/{dataset}/dataset_summary/temp/{dataset}_summary.ipynb",
       # dir = directory(WORKDIR + "/{dataset}/dataset_summary")

    log:
       # optional path to the processed notebook
       notebook=WORKDIR + "/{dataset}/dataset_summary/temp/{dataset}_summary.ipynb"
    notebook:
       "notebooks/create_summary_document.py.ipynb"


#Hide the first cell that snakemake adds to the notebook
rule clean_summary_document:
   input:
       WORKDIR + "/{dataset}/dataset_summary/temp/{dataset}_summary.ipynb"
   output:
       WORKDIR + "/{dataset}/dataset_summary/{dataset}_summary_cleaned.ipynb"

   script:
       "scripts/clean_summary_document.py"

rule create_subset_summary:
    input:
       aln= WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.aln",
       csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_alignment.csv"

    params:
        annotation_cols = config['annotation_cols']

    output:
       summary = WORKDIR + "/{dataset}/dataset_summary/{dataset}/subsets/{subset}/temp/{subset}_subset_summary.ipynb",
    log:
        #optional path to the processed notebook
       notebook=WORKDIR + "/{dataset}/dataset_summary/{dataset}/subsets/{subset}/temp/{subset}_subset_summary.ipynb"
    notebook:
       "notebooks/create_subset_document.py.ipynb"


#Hide the first cell that snakemake adds to the notebook
rule clean_subset_summary:
   input:
       summary = WORKDIR + "/{dataset}/dataset_summary/{dataset}/subsets/{subset}/temp/{subset}_subset_summary.ipynb",
   output:
       summary = WORKDIR + "/{dataset}/dataset_summary/{dataset}/subsets/{subset}/{subset}_subset_summary_cleaned.ipynb",

   script:
       "scripts/clean_summary_document.py"

        #expand("indelible_summaries/{{taxon}}/{rep}.csv", rep = [x for x in range(1, config['REPS'] + 1)])
rule create_subset_document:
   input:
       pdf = WORKDIR + "/{dataset}/dataset_summary/{dataset}_summary_cleaned.ipynb",
       subsets = WORKDIR + "/{dataset}/dataset_summary/{dataset}/subsets/{subset}/{subset}_subset_summary_cleaned.ipynb"
       # subsets = expand(WORKDIR + "/{{dataset}}/dataset_summary/{{dataset}}/subsets/{subset}/temp/{subset}_summary.ipynb", subset = [subset for subset in subsets[{{wildcards.dataset}}]])
   output:
       config = WORKDIR + "/{dataset}/dataset_summary/{subset}/_config.yml",
       toc = WORKDIR + "/{dataset}/dataset_summary/{subset}/_toc.yml"

   script:
       "scripts/create_summary_document.py"

rule compile_summary_document:
   input:
       # dir = WORKDIR + "/{dataset}/dataset_summary",
       config = WORKDIR + "/{dataset}/dataset_summary/{subset}/_config.yml"

   params:
       dir = WORKDIR + "/{dataset}/dataset_summary/{subset}",
   output:
      config = WORKDIR + "/{dataset}/dataset_summary/{subset}/_build/html/index.html"
   shell:
       "jb build {params.dir}"



rule concat_ancestor_alignment:
    input:
        extants = WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.aln",
        ancestors = WORKDIR + "/{dataset}/subsets/{subset}/grasp_results/GRASP_ancestors.fa",

    output:
        WORKDIR + "/{dataset}/subsets/{subset}/concatenated_seqs/{dataset}_{subset}_ancestors.aln"
    shell:
        "cat {input.extants} {input.ancestors} > {output}"

# rule create_ancestor_image:
#     input:
#         tree = WORKDIR + "/{dataset}/subsets/{subset}/grasp_results/GRASP_ancestors.nwk",
#         aln = WORKDIR + "/{dataset}/subsets/{subset}/concatenated_seqs/{dataset}_{subset}_ancestors.aln",
#         csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_ancestors.csv"
#     output:
#         img = "{WORKDIR}/{dataset}/subsets/{subset}/ancestor_images/{anc_col}/{anc_val}.png"
#     script:
#         "scripts/create_ancestor_images.py"


# rule create_tree_image:
#     input:
#         tree = WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.nwk",
#         aln = WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.aln",
#         csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_ancestors.csv"

#     output:
#         WORKDIR + "/{dataset}/subsets/{subset}/tree_images/{column}.png"
#     # output:
#     #     expand(WORKDIR + "/{dataset}/tree_images/{dataset}_{column}.png",
#     #         dataset=DATASETS,
#     #         column=config['annotation_cols'])
#     shell:
#         "python scripts/tree_annot.py -t {input.tree} -a {input.aln} -c {input.csv} --col '{wildcards.column}' --match_from 'Entry' -r 4425 -o {output}"