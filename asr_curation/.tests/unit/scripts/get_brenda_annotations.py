from collections import defaultdict

import pandas as pd
from brendapy import BrendaParser

brenda_cols = ['AC', 'AP', 'CF', 'CL', 'CR', 'EN', 'EXP', 'GI', 'GS',
               'ID', 'IN', 'KI', 'KKM', 'KM', 'LO', 'ME', 'MW', 'NSP', 'OS', 'OSS', 'PHO',
               'PHR', 'PHS', 'PI', 'PM', 'PU', 'RE', 'REN', 'RN', 'RT', 'SA', 'SN', 'SP',
               'SS', 'ST', 'SU', 'SY', 'TN', 'TO', 'TR', 'TS', 'IC50']

# brenda_cols = ['AC', 'AP', 'CF', 'CL', 'CR', 'EN', 'EXP', 'GI', 'GS', 'IC50', 
# 'ID', 'IN',  'KKM', 'LO', 'ME', 'MW', 'NSP', 'OS', 'OSS', 'PHO', 
# 'PHR', 'PHS', 'PI', 'PM', 'PU', 'RE', 'REN', 'RN', 'RT', 'SA', 'SN', 'SP', 
# 'SS', 'ST', 'SU', 'SY',  'TO', 'TR', 'TS']

# brenda_cols = ['CL', 'KM']

BRENDA_PARSER = BrendaParser()  # reuse parser


def parse_proteins_for_ec(ec="1.1.1.1"):
    """Parse the protein entries for a given EC number in BRENDA.
    """
    proteins = BRENDA_PARSER.get_proteins(ec)
    return proteins


def get_ec_nums(dataset_name):
    ec_nums = []
    if 'ec' in dataset_name:
        for ec_string in dataset_name.split('ec_')[1:]:
            ec_nums.append(".".join(ec_string.split("_")[0:4]).replace("_", "."))

    return ec_nums


def parse_all_proteins_for_all_ecs():
    ec_dict = {}
    for ec in BRENDA_PARSER.keys():
        proteins = BRENDA_PARSER.get_proteins(ec)
    ec_dict[ec] = proteins

    return ec_dict


def count_uniprot_entries(ec_dict):
    totals = {}

    for ec, proteins in ec_dict.items():

        for k, v in proteins.items():
            if v.data['uniprot']:
                up_count += 1

            totals[ec] = up_count

    return totals


def add_col_from_brenda_dict(df, entry_id, cols_to_add, brenda_dict):
    for name, annots in brenda_dict.items():
        df.loc[df['Entry'].str.contains(entry_id), name] = ";".join(str(x) for x in annots)

    return df


# def add_val(brenda_dict, protein, attrib, attrib_count, val_name, suffix=""):

#     if 'comment' not in attrib or 'mutant' not in attrib['comment']:

#         if 'substrate' in attrib:
#             brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}_{str(attrib['substrate'])}"].append(f"{attrib['value']}_count={attrib_count}")
#             brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}"].append(f"{attrib['value']};{attrib['substrate']}_count={attrib_count}")


#         elif val_name in attrib:
#             brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}{suffix}"].append(f"{attrib[val_name]}_count={attrib_count}")


def add_val(brenda_dict, protein, attrib, attrib_count):
    
    if 'comment' not in attrib or 'mutant' not in attrib['comment']:

        if 'substrate' in attrib:
            
            terms = ['units', 'refs', 'comment']
            
            brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}_{str(attrib['substrate'])}_DATA"].append(f"{attrib['value']}_count={attrib_count}")
            brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}"].append(f"{attrib['value']};{attrib['substrate']}_count={attrib_count}")
                                                                     
            for term in terms:
                if term in attrib:
                    brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}_{str(attrib['substrate'])}_{term.upper()}"].append(f"{attrib[term]}_count={attrib_count}")


        else:
            terms = ['data', 'units', 'refs', 'comment']
                                                    
            for term in terms:
                if term in attrib:
                    brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}_{term.upper()}"].append(f"{attrib[term]}_count={attrib_count}")


# Create a BRENDA dictionary that maps all of the available Uniprot IDs from BRENDA to their annotations
brenda_dict = defaultdict(lambda: defaultdict(list))


# def add_km(brenda_dict, protein, attrib):

#     print (attrib)
#     if 'comment' not in attrib or 'mutant' not in attrib['comment']:
#         print (attrib)
#         if 'substrate' in attrib:
#             brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}_{str(attrib['substrate'])}].append(f"{attrib['value']}_count={attrib_count}")
#             brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}"].append(f"{attrib['value']};{attrib['substrate']}_count={attrib_count}")


# def add_kkm(brenda_dict, protein, attrib):

#     print (attrib)


#     if 'comment' not in attrib or 'mutant' not in attrib['comment']:
#         brenda_dict[protein.uniprot][f"{str(bc)}"].append(f"{attrib['value']};{attrib['units']}")


ec_nums = get_ec_nums(snakemake.wildcards.dataset)

if ec_nums:

    for ec_num in ec_nums:
        ec_dict = parse_proteins_for_ec(ec_num)

        original_df = pd.read_csv(snakemake.input[0])



        print(f"\nEC number is {ec_num}")
        print(f"\nBRENDA columns we are adding are {', '.join(brenda_cols)}")

        # Create a BRENDA dictionary that maps all of the available Uniprot IDs from BRENDA to their annotations
        # Not getting SN (synonyms) or RN (accepted name (IUPAC)) or IC50
        for prot_id, protein in sorted(ec_dict.items()):
            if protein.uniprot:
                print ('here are references')
                print (protein.references)
                
                # For all the references, add them as separate columns so we can search them
                for refno, ref_dict in protein.references.items():
                    print (ref_dict)
                    brenda_dict[protein.uniprot][f"BRENDA_REFERENCES_{refno}"].append(ref_dict['info'])
                    
                    # If there is a pubmed ID add this as well
                    if 'pubmed' in ref_dict.keys():
                        brenda_dict[protein.uniprot][f"BRENDA_REFERENCES_{refno}_PUBMED"].append((ref_dict['pubmed']))
 
                for bc in brenda_cols:
                    attrib_count = 0
                    attribs = getattr(protein, bc)
                    if attribs:
                        for attrib in attribs:
                            attrib_count +=1
                            if bc in [ 'AP', 'AC', 'CF', 'CL', 'CR', 'EXP', 'IC50', 'LO', 'NSP', 'PHO', 'PU', 'PM', 'SP', 'EN',
                                      'IN', 'ME', 'MW', 'SA', 'ST', 'SU', 'SY', 'TO', 'TR', 'TS', 'KKM', 'KM', 'TN', 'KI', 'OS', 'PHR', 'SS']:
                                add_val(brenda_dict, protein, attrib, attrib_count)
                            if bc == 'GI':
                                print(f"WARNING {bc} is not implemented")
                                print(attrib)
                            if bc == 'GS':
                                print(f"WARNING {bc} is not implemented")
                                print(attrib)
                            if bc == 'OSS':
                                print(f"WARNING {bc} is not implemented")
                                print(attrib)
                            if bc == 'PHS':
                                print(f"WARNING {bc} is not implemented")
                                print(attrib)
                            if bc == 'PI':
                                print(f"WARNING {bc} is not implemented")
                                print(attrib)
                            if bc == 'REN':
                                print(f"WARNING {bc} is not implemented")
                                print(attrib)


    

        # Add the annotations from BRENDA dictionary to the annotation file
        for entry_id, bd in brenda_dict.items():
            print ('Getting BRENDA DF')
            brenda_df = add_col_from_brenda_dict(original_df, entry_id, bd.keys(), bd)
        
    print (f"Writing out the BRENDA annotations to {snakemake.output[0]}")
    brenda_df.to_csv(snakemake.output[0], index=False)

    # original_df.to_csv(snakemake.output[0], index=False)


else:
    print("This dataset doesn't have an EC number associated with it")
    original_df.to_csv(snakemake.output[0], index=False)
