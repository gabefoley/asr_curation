import annot_functions as an
import seqcurate as sc
import pandas as pd

annot_df = pd.read_csv(snakemake.input[0])

annot_df = an.annotate_nonAA(annot_df)
annot_df = an.annotate_AA(annot_df)

annot_df = an.annotate_sp_tr(annot_df)
# annot_df = an.annotate_fragment(annot_df)

if 'uniprot_ec_4_2_1_9' in snakemake.wildcards.dataset:
	annot_df = an.annotate_motif(annot_df, '[AG][AQNP][KTAYG][LVI][PG][RK][HD][NS]H')


	annot_df = sc.add_from_csv(annot_df, './additional_data/dhad/experimental_values.csv')

# Check if the ALS EC number is in the dataset we're currently looking at
if 'uniprot_ec_2_2_1_6' in snakemake.wildcards.dataset:
	print ("Adding ALS specific annotations")


	# Add the catabolic and anabolic ALS motifs
	annot_df = an.annotate_motif(annot_df, 'SPVEY')
	annot_df = an.annotate_motif(annot_df, 'RFDDR')

	# List of truncated sequences
	truncated_seqs = ['T0TXU9','A0A4U9YPB6','A0A139QHE1','K2NT36','C2E671','K0N695','Q3EXB3','A0A0M8WM32','X8FB41','W0BGQ1','A0A377NDB6','A0A7H4P2E2','A0A4P0UPF2','A0A485ALN6','A0A485BQI2','A0A3P8M0G1','A0A6J5DDJ5']
	
	truncated_seqs2 = ['C2E671','K0N695','A0A4U9YPB6','T0TXU9','A0A139QHE1','X0PAD2','Q3EXB3','X8FB41','A0A485BQI2','A0A485ALN6','A0A377NDB6','A0A0S2SI96','A0A377X9X4','W0BGQ1','T2JQD0','M1WSK7','R6CCA3','A0A259KA15','A0A0G1FUD6']

	extended_seqs = ['A0A1Q7NUF9','A0A536GDZ0','A0A3D1RW85','A0A523F6C3']

	truncated_seqs += truncated_seqs2 + extended_seqs

	annot_df['truncated_seqs'] = annot_df.apply(lambda row: True if row['Entry'] in truncated_seqs else False, axis=1)


	# Add an additional Length column (this can be removed once we can call more complex, multiple queries on a single column)
	annot_df['Length_2'] = annot_df['Length']
	annot_df['Cross_reference_InterPro_2'] = annot_df['Cross_reference_InterPro']


	# annot_df = an.summarise_motifs(annot_df, 'SPVEY', 'RFDDR')


# Check if the KARI EC number is in the dataset we're currently looking at
if 'uniprot_ec_1_1_1_86' in snakemake.wildcards.dataset:
	print ("Adding KARI specific annotations")

	
	print ("Adding KARI Class")
	annot_df['KARI_Class'] = annot_df.apply(lambda row : an.classify_KARI(row['Features'] if pd.notnull(row['Features']) else ''), axis = 1)
	
	print ("Adding loop length")
	annot_df['Loop_Length'] = annot_df.apply(lambda row : 
		an.classify_loop_length(
			an.get_binding_pos(row['Binding_site']) 
			if pd.notnull(row['Binding_site']) 
			else 'No_binding_positions'), axis = 1)

	print ("Adding Binding positions extracted")
	annot_df['Binding_positions_extracted'] = annot_df.apply(lambda row : "No_binding_positions" if pd.isnull(row['Binding_site']) else an.get_binding_pos(row['Binding_site']), axis = 1)
	# annot_df['Binding_positions_extracted'] = annot_df.apply(lambda row : "No_binding_positions" if pd.notnull(row['Binding_site']) else an.get_binding_pos(row['Binding_site']), axis = 1)


	# df['Testing']=df.apply(lambda x: 1 if x['Liq_Factor']=='Nan'  else min(x['Use']/x['Tw'],1), axis=1)

	print ("Adding Binding positions character")
	annot_df['Binding_positions_character'] = annot_df.apply(lambda row : "No_binding_positions" if pd.isnull(row['Binding_site']) else  an.get_amino_acids(row['Sequence'], *row['Binding_positions_extracted']), axis = 1)
	
	print ("Adding acidic binding")
	# annot_df['Acidic_Binding'] = annot_df.apply(lambda row : "No_binding_positions" if pd.notnull(row['Binding_site']) else an.check_binding_for_acidic(row['Sequence'], *row['Binding_positions_extracted']), axis = 1)

	print ("Add a tag to NADH preferring sequences")

	import os

	print (os.getcwd())
	psuedomonas_aeruginosa = sc.get_entry_ids_from_fasta("./additional_data/kari/pseudomonas_aeruginosa.fasta")


	annot_df['psuedomonas_aeruginosa'] = annot_df.apply(lambda row: True if row['Entry'] in psuedomonas_aeruginosa else False, axis=1)


	nadh_pref = sc.get_entry_ids_from_fasta("./additional_data/NADH_preferring.fasta")

	print (nadh_pref)

	annot_df['NADH_pref'] = annot_df.apply(lambda row: True if row['Entry'] in nadh_pref else False, axis=1)

	general_approach = sc.get_entry_ids_from_fasta("./additional_data/general_approach_class_I.fasta")
	print ("Add a tag to sequences in the general approach paper that are also Class I")

	annot_df['General_approach_Class_I'] = annot_df.apply(lambda row: True if row['Entry'] in general_approach else False, axis=1)



annot_df.to_csv(snakemake.output[0], index=False)