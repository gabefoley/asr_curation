import scripts.annot_functions as an

def test_get_amino_acids():

    pos = an.get_amino_acids('PPGP', 0)
    assert pos == 'P'


def test_get_amino_acids_multiple():

    pos = an.get_amino_acids('PPGP', 0, 2)
    assert pos == 'PG'

def test_get_binding_pos():
    ft_binding = 'BINDING 123..130; /ligand="NADP(+)"; /ligand_id="ChEBI:CHEBI:58349"; /evidence="ECO:0000250"; BINDING 156..161; /ligand="NADP(+)"; /ligand_id="ChEBI:CHEBI:58349"; /evidence="ECO:0000250"; BINDING 195..199; /ligand="NADP(+)"; /ligand_id="ChEBI:CHEBI:58349"; /evidence="ECO:0000250"; BINDING 309; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="1"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 309; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="2"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 313; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="1"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 486; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="2"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 490; /ligand="Mg(2+)"; /ligand_id="ChEBI:CHEBI:18420"; /ligand_label="2"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198"; BINDING 512; /ligand="substrate"; /evidence="ECO:0000255|PROSITE-ProRule:PRU01198'
