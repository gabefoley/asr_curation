import pandas as pd
import click
import os


# annotation_cols = ['Entry', 'Taxonomic_lineage_PHYLUM', 'Taxonomic_lineage_SUPERKINGDOM', 'Loop_Length', 'Binding_positions_extracted', 'Binding_positions_character', 'Acidic_Binding', 'Cross_reference_InterPro_2', 'KARI_Class',
# 'Dimer_alignment_residues', 'Binding_positions_alignment_residues', 'Acidic_Binding_alignment', 'thermo_bacteria_split', 'thermo_bacteria', 'thermo_terms']

# annotation_cols = ['Entry', 'Taxonomic_lineage_PHYLUM', 'Taxonomic_lineage_SUPERKINGDOM', 'Taxonomic_lineage_CLASS', 'Cross_reference_OMA', 'Cross_reference_InterPro', 'Cross_reference_Pfam']

@click.command()
@click.option('--df', help='Dataframe with all annotations')
@click.option('--annot', help='Columns to add to annotation file')
@click.option('--outpath', default='./annotation_cols.txt', help='Outpath for annotation file')
def create_annotations(df, annot, outpath):

    if 'SNAKEMAKE' in os.environ:
        df = pd.read_csv(snakemake.input.csv)
        annotation_cols = snakemake.params.annotation_cols
        outpath = snakemake.output.tsv
    
    else:

        df = pd.read_csv(df)

        with open(annot) as annot_file:
            annotation_cols = [line.strip() for line in annot_file]

    
    # Accession needs to be the first column so if it isn't requested, add it in.
    if annotation_cols[0] != "accession":
        annotation_cols = ["accession"] + annotation_cols


    # Subset the columns
    subset_df = df[[x for x in annotation_cols if x in df.columns]]
    subset_df = subset_df.fillna("None")
    
    # Write out the annotation columns file
    subset_df.to_csv(outpath, sep="\t", index=False)


if __name__ == "__main__":
    create_annotations()

