# 2021 July Note: 'comment(ENZYME REGULATION)', 'context', 'tools'  are not returning a value correctly
# 2022 Sep Note: ' commented fields not supported by uniprot api ' - SANJANA


full_uniprot_cols = [
    "accession",
    "id",
    "gene_names",
    "gene_primary",
    "gene_synonym",
    "gene_oln",
    "gene_orf",
    "organism_name",
    "organism_id",
    "protein_name",
    "sequence",
    "xref_proteomes",
    "lineage",
    "virus_hosts",
    "fragment",
    "organelle",
    "cc_alternative_products",
    "error_gmodel_pred",
    # "comment(ERRONEOUS INITIATION)",
    # "comment(ERRONEOUS TERMINATION)",
    # "comment(ERRONEOUS TRANSLATION)",
    # "comment(FRAMESHIFT)",
    "cc_mass_spectrometry",
    "cc_polymorphism",
    "cc_rna_editing",
    "cc_sequence_caution",
    "length",
    "mass",
    "ft_var_seq",
    "ft_variant",
    "ft_non_cons",
    "ft_non_std",
    "ft_non_ter",
    "ft_conflict",
    "ft_unsure",
    "sequence_version",
    "ec",
    "absorption",
    # "comment(CATALYTIC ACTIVITY)",
    # "ft_ca_bind",
    "cc_catalytic_activity",
    "cc_cofactor",
    # "chebi-id",
    # "comment(COFACTOR)",
    "cc_function",
    "kinetics",
    "cc_pathway",
    "redox_potential",
    "temp_dependence",
    "ph_dependence",
    "ft_act_site",
    "ft_binding",
    "ft_dna_bind",
    # "ft_metal",
    # "ft_np_bind",
    "ft_site",
    "annotation_score",
    "feature_count",
    "cc_caution",
    "cc_miscellaneous",
    "keyword",
    "protein_existence",
    "reviewed",
    "cc_subunit",
    "cc_interaction",
    "cc_developmental_stage",
    "cc_induction",
    "cc_tissue_specificity",
    "go",
    "go_p",
    "go_f",
    "go_c",
    "go_id",
    "cc_allergen",
    "cc_biotechnology",
    "cc_disruption_phenotype",
    "cc_disease",
    "cc_pharmaceutical",
    "cc_toxic_dose",
    "cc_subcellular_location",
    "ft_intramem",
    "ft_topo_dom",
    "ft_transmem",
    "cc_ptm",
    "ft_chain",
    "ft_crosslnk",
    "ft_disulfid",
    "ft_carbohyd",
    "ft_init_met",
    "ft_lipid",
    "ft_mod_res",
    "ft_peptide",
    "ft_propep",
    "ft_signal",
    "ft_transit",
    "structure_3d",
    "ft_strand",
    "ft_helix",
    "ft_turn",
    # "citationmapping",
    "lit_pubmed_id",
    "date_created",
    "date_modified",
    "date_sequence_modified",
    "version",
    "cc_domain",
    # "comment(SIMILARITY)",
    "protein_families",
    "ft_coiled",
    "ft_compbias",
    "ft_domain",
    "ft_motif",
    "ft_region",
    "ft_repeat",
    "ft_zn_fing",
    "lineage",
    # "lineage(SUPERKINGDOM)",
    # "lineage(KINGDOM)",
    # "lineage(SUBKINGDOM)",
    # "lineage(SUPERPHYLUM)",
    # "lineage(PHYLUM)",
    # "lineage(SUBPHYLUM)",
    # "lineage(SUPERCLASS)",
    # "lineage(CLASS)",
    # "lineage(SUBCLASS)",
    # "lineage(INFRACLASS)",
    # "lineage(SUPERORDER)",
    # "lineage(ORDER)",
    # "lineage(SUBORDER)",
    # "lineage(INFRAORDER)",
    # "lineage(PARVORDER)",
    # "lineage(SUPERFAMILY)",
    # "lineage(FAMILY)",
    # "lineage(SUBFAMILY)",
    # "lineage(TRIBE)",
    # "lineage(SUBTRIBE)",
    # "lineage(GENUS)",
    # "lineage(SUBGENUS)",
    # "lineage(SPECIES GROUP)",
    # "lineage(SPECIES SUBGROUP)",
    # "lineage(SPECIES)",
    # "lineage(SUBSPECIES)",
    # "lineage(VARIETAS)",
    # "lineage(FORMA)",
    "xref_abcd",
    "xref_allergome",
    "xref_antibodypedia",
    "xref_arachnoserver",
    "xref_araport",
    "xref_bgee",
    "xref_bindingdb",
    "xref_biocyc",
    "xref_biogrid",
    "xref_biogrid-orcs",
    "xref_biomuta",
    "xref_bmrb",
    "xref_brenda",
    "xref_carbonyldb",
    "xref_cazy",
    "xref_ccds",
    "xref_cdd",
    "xref_cgd",
    "xref_chembl",
    "xref_chitars",
    "xref_clae",
    "xref_collectf",
    "xref_complexportal",
    "xref_compluyeast-2dpage",
    "xref_conoserver",
    "xref_corum",
    "xref_cptac",
    "xref_cptc",
    "xref_ctd",
    "xref_dbsnp",
    # "database(DDBJ)",
    "xref_depod",
    "xref_dictybase",
    "xref_dip",
    "xref_disgenet",
    "xref_disprot",
    "xref_dmdm",
    "xref_dnasu",
    "xref_dosac-cobs-2dpage",
    "xref_drugbank",
    "xref_drugcentral",
    "xref_echobase",
    "xref_eggnog",
    "xref_elm",
    "xref_embl",
    "xref_ensembl",
    "xref_ensemblbacteria",
    "xref_ensemblfungi",
    "xref_ensemblmetazoa",
    "xref_ensemblplants",
    "xref_ensemblprotists",
    # "database(ENZYME)",
    "xref_epd",
    "xref_esther",
    "xref_euhcvdb",
    "xref_evolutionarytrace",
    "xref_expressionatlas",
    "xref_flybase",
    # "database(GenAtlas)",
    # "database(GenBank)",
    "xref_gene3d",
    "xref_genecards",
    # "database(GeneDB)",
    "xref_geneid",
    "xref_genereviews",
    "xref_genetree",
    "xref_genevisible",
    "xref_genewiki",
    "xref_genomernai",
    "xref_glyconnect",
    "xref_glygen",
    # "database(GO)",
    # "database(GPCRDB)",
    "xref_gramene",
    "xref_guidetopharmacology",
    "xref_hamap",
    "xref_hgnc",
    "xref_hogenom",
    "xref_hpa",
    # "database(HUGE)",
    "xref_ideal",
    "xref_imgt_gene-db",
    "xref_inparanoid",
    "xref_intact",
    "xref_interpro",
    "xref_iptmnet",
    "xref_jpost",
    "xref_kegg",
    "xref_legiolist",
    "xref_leproma",
    "xref_maizegdb",
    "xref_malacards",
    "xref_massive",
    "xref_maxqb",
    "xref_merops",
    "xref_metosite",
    "xref_mgi",
    # "database(Micado)",
    "xref_mim",
    "xref_mint",
    # "database(MobiDB)",
    # "database(ModBase)",
    "xref_moondb",
    "xref_moonprot",
    "xref_nextprot",
    "xref_niagads",
    # "database(OGP)",
    "xref_oma",
    "xref_opentargets",
    "xref_orphanet",
    "xref_orthodb",
    "xref_panther",
    "xref_pathwaycommons",
    "xref_patric",
    "xref_paxdb",
    "xref_pcddb",
    "xref_pdb",
    # "database(PDBe-KB)",
    # "database(PDBj)",
    "xref_pdbsum",
    "xref_peptideatlas",
    "xref_peroxibase",
    "xref_pfam",
    "xref_pharmgkb",
    "xref_pharos",
    "xref_phi-base",
    "xref_phosphositeplus",
    "xref_phylomedb",
    "xref_pir",
    "xref_pirsf",
    "xref_plantreactome",
    "xref_pombase",
    "xref_pride",
    "xref_prints",
    "xref_pro",
    "xref_promex",
    "xref_prosite",
    # "database(Proteomes)",
    "xref_proteomicsdb",
    # "database(ProtoNet)",
    "xref_pseudocap",
    # "database(RCSB-PDB)",
    "xref_reactome",
    "xref_rebase",
    "xref_refseq",
    "xref_reproduction-2dpage",
    "xref_wormbase",
    "xref_rnact",
    # "database(Rouge)",
    "xref_sabio-rk",
    "xref_sasbdb",
    "xref_sfld",
    "xref_sgd",
    "xref_signalink",
    "xref_signor",
    "xref_smart",
    "xref_smr",
    # "database(SOURCE)",
    "xref_string",
    "xref_supfam",
    "xref_swiss-2dpage",
    # "database(SWISS-MODEL-Workspace)",
    "xref_swisslipids",
    "xref_swisspalm",
    "xref_tair",
    "xref_tcdb",
    # "xref_tigrfams",
    "xref_topdownproteomics",
    "xref_treefam",
    "xref_tuberculist",
    "xref_ucd-2dpage",
    "xref_ucsc",
    "xref_unilectin",
    "xref_unipathway",
    "xref_veupathdb",
    "xref_vgnc",
    "xref_wbparasite",
    "xref_world-2dpage",
    "xref_wormbase",
    "xref_xenbase",
    "xref_zfin",
    "comment_count",  # extra fields 2022 sep
    "xref_unicarbkb",
    "xref_vectorbase",
    "xref_wbparasitetranscriptprotein",
    "xref_ko",
    "xref_cleanex",
    "xref_prodom",
]
