import glob
from collections import defaultdict
import os

# Set SNAKEMAKE environment variable (so that we can run scripts from the command line as well)

os.environ["SNAKEMAKE"] = 'True'

# Add the scripts folder to PYTHONPATH. This is needed if we use custom annotations located elsewhere (so they can access the scripts folder)
if "PYTHONPATH" in os.environ:
    os.environ["PYTHONPATH"] += f':{os.getcwd()}/scripts'

else:
    os.environ["PYTHONPATH"] = f'{os.getcwd()}/scripts'

# Collect config info from config file
WORKDIR = config['workdir'] + "/datasets"
FASTADIR = config['fastadir']
SUBDIR = config['subdir']
VERBOSE = True

# Path for custom annotation scripts. Defaults to scripts folder

try:
    CUSTOM_DIR = config['custom_dir']
except:
    CUSTOM_DIR = "scripts"

try:
    CUSTOM_ALIGN_DIR = config['custom_align_dir']
except:
    CUSTOM_ALIGN_DIR = "scripts"

try:
    CUSTOM_ANCESTOR_DIR = config['custom_ancestor_dir']
except:
    CUSTOM_ANCESTOR_DIR = "scripts"


try:
    ALIGNMENT_TOOL = config['alignment_tool']
except:
    ALIGNMENT_TOOL = "mafft"

try:
    KEY_SEQUENCES = config['key_sequences']
except:
    KEY_SEQUENCES = "scripts/blank_key_seqs.fasta"


# Collect config info from parameters passed from command line
if 'BRENDA_RUN' in config.keys():
    BRENDA_RUN = config['BRENDA_RUN']
else:
    BRENDA_RUN = False



# cluster_threshes = ["1", "0.9", "0.7"]
# cluster_threshes = ["0.65", "0.7"]
cluster_threshes = ["1"]


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

subsets = get_subset_names(SUBDIR)
print (WORKDIR)
print (FASTADIR)


rule all:
        input:
            # brenda =     [f'{WORKDIR}/{dataset}/csv/brenda/{dataset}_brenda.csv' for cluster_thresh in cluster_threshes for dataset in DATASETS for subset in subsets[dataset]]
            annotations = [f'{WORKDIR}/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_alignment_annotations.txt' for cluster_thresh in cluster_threshes for dataset in DATASETS for subset in subsets[dataset]],
            reordered_annotations = [f'{WORKDIR}/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_alignment_reordered.csv' for cluster_thresh in cluster_threshes for dataset in DATASETS for subset in subsets[dataset]],

            trees = [f'{WORKDIR}/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.nwk' for cluster_thresh in cluster_threshes for dataset in DATASETS for subset in subsets[dataset]],
            # ancestors = [f'{WORKDIR}/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_ancestors.csv' for cluster_thresh in cluster_threshes for dataset in DATASETS for subset in subsets[dataset]],
            # extants_and_ancestors = [f'{WORKDIR}/{dataset}/subsets/{subset}/{cluster_thresh}/concatenated_seqs/{dataset}_{subset}_{cluster_thresh}_ancestors.aln' for cluster_thresh in cluster_threshes for dataset in DATASETS for subset in subsets[dataset]]

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
#        "scripts/validate_ids.py"
       "scripts/validate_ids_not_uniprot.py"

# Map to UniProt to get all of the known UniProt annotations
rule get_uniprot_annotations:
    input:
        WORKDIR + "/{dataset}/csv/validated/{dataset}_validated.csv"
    output:
        WORKDIR + "/{dataset}/csv/uniprot/{dataset}_uniprot.csv"
    script:
        "scripts/get_uniprot_annotations.py"

# Map to BRENDA database to get all of the known BRENDA annotations
rule get_brenda_annotations:
    input:
        WORKDIR + "/{dataset}/csv/uniprot/{dataset}_uniprot.csv"
        # WORKDIR + "/{dataset}/csv/gtdb_processed/{dataset}_gtdb_processed.csv"
    output:
        WORKDIR + "/{dataset}/csv/brenda/{dataset}_brenda.csv"
    script:
        "scripts/get_brenda_annotations.py"

# Add any generic annotations
rule add_generic_annotations:
    input:
        WORKDIR + "/{dataset}/csv/brenda/{dataset}_brenda.csv"
    output:
        WORKDIR + "/{dataset}/csv/custom/{dataset}_generic_annotated.csv"
    script:
        "scripts/add_generic_annotations.py"

# Add any custom annotations
rule add_custom_annotations:
    input:
        csv = WORKDIR + "/{dataset}/csv/custom/{dataset}_generic_annotated.csv",
        custom_dir = CUSTOM_DIR
    output:
        WORKDIR + "/{dataset}/csv/custom/{dataset}_annotated.csv"
    script:
        CUSTOM_DIR + "/add_custom_annotations.py"

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
        csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}.csv",
        # subset_log = WORKDIR + "/{dataset}/subsets/{subset}/{subset}.log"
    script:
        "scripts/create_subsets.py"

rule cluster_sequences:
        input:
            WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.fasta"

        output:
            WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.fasta"
        shell:
            "cd-hit -i {input} -o {output} -c {wildcards.cluster_thresh}"

# rule create_clustered_subset
#         input:
#             WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.fasta",
#         output:
#             WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}.csv"
#         script:
#         "scripts/create_clustered_subset.py"
#


rule add_key_sequences:
    input:
        fasta = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.fasta",
        key_sequences = KEY_SEQUENCES,
        full_csv = WORKDIR + "/{dataset}/csv/custom/{dataset}_annotated.csv",
        align_csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}.csv"

    output:
        fasta = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}_key_sequences_added.fasta",
        csv = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_key_sequences_added.csv"

    script:
        "scripts/add_key_sequences.py"




if ALIGNMENT_TOOL == 'mafft':
    rule align_seqs:
        input:
            WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}_key_sequences_added.fasta"
        output:
            WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.aln"
        shell:
            "mafft --reorder {input} > {output}"

elif ALIGNMENT_TOOL == 'mafft-dash':
    rule align_seqs_dash:
        input:
            WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}_key_sequences_added.fasta"
        output:
            WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.aln"
        shell:
            "mafft --dash --reorder --originalseqonly {input} > {output}"

rule infer_tree:
    input:
            WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.aln"

    output:
            WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.nwk"

    shell:
        "FastTree {input} > {output}"

rule run_grasp:
    input:
        aln=  WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.aln",
        tree= WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.nwk"


    output:
        dir= directory(WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/grasp_results/"),
        aln= WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/grasp_results/GRASP_ancestors.fa",
        tree = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/grasp_results/GRASP_ancestors.nwk"

    shell:
        "grasp -a {input.aln} -n {input.tree} -s LG -o {output.dir} -i BEP -j --save-as FASTA TREE -pre GRASP -t 2"

rule add_annotations_from_alignment:
    input:
        aln = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.aln",
        csv = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_key_sequences_added.csv"
    output:
        csv = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_alignment.csv"


    script:
        CUSTOM_ALIGN_DIR + "/add_annotations_from_alignment.py"

rule reorder_alignment_annotations:
    input:
        aln = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.aln",
        csv = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_alignment.csv"
    output:
        csv = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_alignment_reordered.csv"

    script:
        CUSTOM_ALIGN_DIR + "/reorder_alignment_annotations.py"

rule create_alignment_annotation_file:
    input:
        csv = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_alignment.csv"

    params:
        annotation_cols = config['annotation_cols']
    output:
        tsv = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_alignment_annotations.txt",
    script:
        "scripts/create_annotation_file.py"


rule add_annotations_from_ancestors:
    input:
        csv = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_alignment.csv",
        aln= WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/grasp_results/GRASP_ancestors.fa",
        tree = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/grasp_results/GRASP_ancestors.nwk"
    output:
        csv = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_ancestors.csv"
    script:
        CUSTOM_ANCESTOR_DIR + "/add_annotations_from_ancestors.py"


rule create_ancestor_annotation_file:
    input:
        csv = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_ancestors.csv"

    params:
        annotation_cols = config['annotation_cols']
    output:
        tsv = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/csv/{dataset}_{subset}_{cluster_thresh}_ancestor_annotations.txt",
    script:
        "scripts/create_annotation_file.py"

rule concat_ancestor_alignment:
    input:
        extants =  WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/{dataset}_{subset}_{cluster_thresh}.aln",
        ancestors = WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/grasp_results/GRASP_ancestors.fa",

    output:
        WORKDIR + "/{dataset}/subsets/{subset}/{cluster_thresh}/concatenated_seqs/{dataset}_{subset}_{cluster_thresh}_ancestors.aln"
    shell:
        "cat {input.extants} {input.ancestors} > {output}"
