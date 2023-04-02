# output for making the summary documents
# summary_document = [f'{WORKDIR}/{dataset}/dataset_summary/{subset}/_build/html/index.html' for dataset in DATASETS for subset in subsets[dataset]],


## Rules for creating the PDF summary files
#
# rule create_dataset_summary:
#     input:
#        WORKDIR + "/{dataset}/csv/custom/{dataset}_annotated.csv",
#     params:
#         annotation_cols = config['annotation_cols']
#     output:
#        summary = WORKDIR + "/{dataset}/dataset_summary/temp/{dataset}_summary.ipynb",
#        # dir = directory(WORKDIR + "/{dataset}/dataset_summary")
#
#     log:
#        # optional path to the processed notebook
#        notebook=WORKDIR + "/{dataset}/dataset_summary/temp/{dataset}_summary.ipynb"
#     notebook:
#        "notebooks/create_summary_document.py.ipynb"
#
#
# #Hide the first cell that snakemake adds to the notebook
# rule clean_summary_document:
#    input:
#        WORKDIR + "/{dataset}/dataset_summary/temp/{dataset}_summary.ipynb"
#    output:
#        WORKDIR + "/{dataset}/dataset_summary/{dataset}_summary_cleaned.ipynb"
#
#    script:
#        "scripts2/clean_summary_document.py"
#
# rule create_subset_summary:
#     input:
#        aln= WORKDIR + "/{dataset}/subsets/{subset}/{dataset}_{subset}.aln",
#        csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_alignment.csv"
#
#     params:
#         annotation_cols = config['annotation_cols']
#
#     output:
#        summary = WORKDIR + "/{dataset}/dataset_summary/{dataset}/subsets/{subset}/temp/{subset}_subset_summary.ipynb",
#     log:
#         #optional path to the processed notebook
#        notebook=WORKDIR + "/{dataset}/dataset_summary/{dataset}/subsets/{subset}/temp/{subset}_subset_summary.ipynb"
#     notebook:
#        "notebooks/create_subset_document.py.ipynb"
#
#
# #Hide the first cell that snakemake adds to the notebook
# rule clean_subset_summary:
#    input:
#        summary = WORKDIR + "/{dataset}/dataset_summary/{dataset}/subsets/{subset}/temp/{subset}_subset_summary.ipynb",
#    output:
#        summary = WORKDIR + "/{dataset}/dataset_summary/{dataset}/subsets/{subset}/{subset}_subset_summary_cleaned.ipynb",
#
#    script:
#        "scripts2/clean_summary_document.py"
#
#         #expand("indelible_summaries/{{taxon}}/{rep}.csv", rep = [x for x in range(1, config['REPS'] + 1)])
# rule create_subset_document:
#    input:
#        pdf = WORKDIR + "/{dataset}/dataset_summary/{dataset}_summary_cleaned.ipynb",
#        subsets = WORKDIR + "/{dataset}/dataset_summary/{dataset}/subsets/{subset}/{subset}_subset_summary_cleaned.ipynb"
#        # subsets = expand(WORKDIR + "/{{dataset}}/dataset_summary/{{dataset}}/subsets/{subset}/temp/{subset}_summary.ipynb", subset = [subset for subset in subsets[{{wildcards.dataset}}]])
#    output:
#        config = WORKDIR + "/{dataset}/dataset_summary/{subset}/_config.yml",
#        toc = WORKDIR + "/{dataset}/dataset_summary/{subset}/_toc.yml"
#
#    script:
#        "scripts2/create_summary_document.py"
#
# rule compile_summary_document:
#    input:
#        # dir = WORKDIR + "/{dataset}/dataset_summary",
#        config = WORKDIR + "/{dataset}/dataset_summary/{subset}/_config.yml"
#
#    params:
#        dir = WORKDIR + "/{dataset}/dataset_summary/{subset}",
#    output:
#       config = WORKDIR + "/{dataset}/dataset_summary/{subset}/_build/html/index.html"
#    shell:
#        "jb build {params.dir}"
#


# # # Map the UniProt sequences to IDs we can query IPG with, and then query IPG to add genome information
# rule get_genomes:
#     input:
#         WORKDIR + "/{dataset}/csv/uniprot/{dataset}_uniprot.csv"
#     output:
#         WORKDIR + "/{dataset}/csv/genome/{dataset}_genome.csv"
#         #WORKDIR + "/{dataset}/csv/gtdb_processed/{dataset}_gtdb_processed.csv"
#     script:
#         "scripts2/get_genomes.py"

# # Map the Genomes to GTDB
# rule get_gtdb_annotations:
#    input:
#        WORKDIR + "/{dataset}/csv/genome/{dataset}_genome.csv"
#    output:
#        WORKDIR + "/{dataset}/csv/gtdb/{dataset}_gtdb.csv"
#    script:
#        "scripts2/get_gtdb_annotations.py"

# #Use GTDB to cross-check CheckM scores and GTDB taxonomy
# rule process_gtdb_annotations:
#    input:
#        WORKDIR + "/{dataset}/csv/gtdb/{dataset}_gtdb.csv"
#    output:
#        WORKDIR + "/{dataset}/csv/gtdb_processed/{dataset}_gtdb_processed.csv"
#    script:
#        "scripts2/process_gtdb_annotations.py"


# rule create_ancestor_image:
#     input:
#         tree = WORKDIR + "/{dataset}/subsets/{subset}/grasp_results/GRASP_ancestors.nwk",
#         aln = WORKDIR + "/{dataset}/subsets/{subset}/concatenated_seqs/{dataset}_{subset}_ancestors.aln",
#         csv = WORKDIR + "/{dataset}/subsets/{subset}/csv/{dataset}_{subset}_ancestors.csv"
#     output:
#         img = "{WORKDIR}/{dataset}/subsets/{subset}/ancestor_images/{anc_col}/{anc_val}.png"
#     script:
#         "scripts2/create_ancestor_images.py"


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
#         "python scripts2/tree_annot.py -t {input.tree} -a {input.aln} -c {input.csv} --col '{wildcards.column}' --match_from 'Entry' -r 4425 -o {output}"
