import pandas as pd

align_df = pd.read_csv(snakemake.input.csv)


# annotation_cols = ['Entry', 'Taxonomic_lineage_PHYLUM', 'Taxonomic_lineage_SUPERKINGDOM', 'Loop_Length', 'Binding_positions_extracted', 'Binding_positions_character', 'Acidic_Binding', 'Cross_reference_InterPro_2', 'KARI_Class',
# 'Dimer_alignment_residues', 'Binding_positions_alignment_residues', 'Acidic_Binding_alignment', 'thermo_bacteria_split', 'thermo_bacteria', 'thermo_terms']

# annotation_cols = ['Entry', 'Taxonomic_lineage_PHYLUM', 'Taxonomic_lineage_SUPERKINGDOM', 'Taxonomic_lineage_CLASS', 'Cross_reference_OMA', 'Cross_reference_InterPro', 'Cross_reference_Pfam']

if snakemake.params.annotation_cols[0] != "accession":

    annotation_cols = ["accession"] + snakemake.params.annotation_cols

subset_df = align_df[[x for x in annotation_cols if x in align_df.columns]]

subset_df = subset_df.fillna("None")

subset_df.to_csv(snakemake.output.tsv, sep="\t", index=False)
