from collections import defaultdict
import pandas as pd
from brendapy import BrendaParser
from configs.brenda_cols import brenda_cols

BRENDA_PARSER = BrendaParser()


def parse_proteins_for_ec(ec="1.1.1.1"):
    """Parse the protein entries for a given EC number in BRENDA."""
    proteins = BRENDA_PARSER.get_proteins(ec)
    return proteins


def get_ec_nums(dataset_name):
    ec_nums = []
    if "ec" in dataset_name:
        for ec_string in dataset_name.split("ec_")[1:]:
            ec_nums.append(".".join(ec_string.split("_")[0:4]).replace("_", "."))
    return ec_nums


def count_uniprot_entries(ec_dict):
    totals = {}
    for ec, proteins in ec_dict.items():
        up_count = sum(1 for v in proteins.values() if v.data["uniprot"])
        totals[ec] = up_count
    return totals


def add_col_from_brenda_dict(df, entry_id, cols_to_add, brenda_dict):
    for name, annots in brenda_dict.items():
        df.loc[df["extracted_id"].str.contains(entry_id), name] = ";".join(str(x) for x in annots)
    return df


def add_val(brenda_dict, protein, attrib, attrib_count, bc):
    if "comment" not in attrib or "mutant" not in attrib["comment"]:
        if "substrate" in attrib:
            terms = ["units", "refs", "comment"]
            brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}_{str(attrib['substrate'])}_DATA"].append(
                f"{attrib['value']}_count={attrib_count}")
            brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}"].append(
                f"{attrib['value']};{attrib['substrate']}_count={attrib_count}")
            for term in terms:
                if term in attrib:
                    brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}_{str(attrib['substrate'])}_{term.upper()}"].append(
                        f"{attrib[term]}_count={attrib_count}")
        else:
            terms = ["data", "units", "refs", "comment"]
            for term in terms:
                if term in attrib:
                    brenda_dict[protein.uniprot][f"BRENDA_{str(bc)}_{term.upper()}"].append(
                        f"{attrib[term]}_count={attrib_count}")


def main():
    ec_nums = get_ec_nums(snakemake.wildcards.dataset)
    original_df = pd.read_csv(snakemake.input[0])
    brenda_dict = defaultdict(lambda: defaultdict(list))

    if ec_nums:
        for ec_num in ec_nums:
            ec_dict = parse_proteins_for_ec(ec_num)
            original_df = pd.read_csv(snakemake.input[0])
            print(f"\nEC number is {ec_num}")
            print(f"\nBRENDA columns we are adding are {', '.join(brenda_cols)}")

            for prot_id, protein in sorted(ec_dict.items()):
                if protein.uniprot:
                    for refno, ref_dict in protein.references.items():
                        brenda_dict[protein.uniprot][f"BRENDA_REFERENCES_{refno}"].append(ref_dict["info"])
                        if "pubmed" in ref_dict.keys():
                            brenda_dict[protein.uniprot][f"BRENDA_REFERENCES_{refno}_PUBMED"].append(
                                (ref_dict["pubmed"]))

                    for bc in brenda_cols:
                        attrib_count = 0
                        attribs = getattr(protein, bc)
                        if attribs:
                            for attrib in attribs:
                                attrib_count += 1
                                add_val(brenda_dict, protein, attrib, attrib_count, bc)

        for entry_id, bd in brenda_dict.items():
            brenda_df = add_col_from_brenda_dict(original_df, entry_id, bd.keys(), bd)

        print(f"Writing out the BRENDA annotations to {snakemake.output[0]}")
        brenda_df.to_csv(snakemake.output[0], index=False)

    else:
        print("This dataset doesn't have an EC number associated with it")
        original_df.to_csv(snakemake.output[0], index=False)


if __name__ == "__main__":
    main()
