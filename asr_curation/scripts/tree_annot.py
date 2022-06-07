import argparse
import tree_code as tc
import sys
import random
import os
os.environ['QT_QPA_PLATFORM']='offscreen'


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", help="Path to tree", required=True)
    parser.add_argument("-a", "--align", help="Path to alignment", required=True)
    parser.add_argument("-c", "--csv", help="Path to csv", required=True)
    parser.add_argument(
        "-col", "--col", help="Column to annotate based on", required=True
    )
    parser.add_argument(
        "-s",
        "--seed",
        help="Set the random seed to an integer to make colour selection "
        "reproducible",
    )

    parser.add_argument("-o", "--outpath", help="Outpath", default="tree_annot.png")
    parser.add_argument(
        "-m",
        "--match_from",
        help="Match from a specific column in the annotations file",
        default="Name",
    )
    parser.add_argument(
        "-r",
        "--rotate",
        help="Internal nodes to rotate for visualisation",
        nargs="+",
        default=None,
    )

    parser.add_argument(
        "-w",
        "--whitespace_split",
        help="Split the column based on whitespace",
        action="store_true",
    )

    return parser.parse_args(args)


if __name__ == "__main__":

    print("\nRunning tree_annot")
    print("python version is ")
    print(sys.version)
    import pandas as pd

    # Parse the arguments
    parser = parse_args(sys.argv[1:])

    df = pd.read_csv(parser.csv)

    annot_dict = dict(zip(df[parser.match_from], df[parser.col]))

    for key, val in annot_dict.items():
        if pd.isnull(val):
            annot_dict[key] = None
        elif parser.whitespace_split:
            annot_dict[key] = annot_dict[key].split(" ")[0]

    # else:
    # annot_dict = df.to_dict(orient='index')

    print(annot_dict)

    col = parser.col.strip() if parser.col else df.columns[0]

    random_seed = parser.seed if parser.seed else random.randint(0, 999)

    print(
        "\nTo reproduce this colour scheme set the random seed as " + str(random_seed)
    )

    color_dict = tc.get_color_dict(annot_dict, random_seed)

    print(color_dict)

    color_dict[None] = "white"

    # Load tree
    tree = tc.load_tree(parser.tree, parser.align)

    print("parser rotate is ")
    print(parser.rotate)

    # Manually set the Clade (Mischko) colours to match Mischko et al. 2018
    if col == "Clade (Mischko)":

        color_dict[1] = "#f88485"
        color_dict[2] = "#96b9da"
        color_dict[3] = "#c4e0a4"
        color_dict[4] = "#ffdf80"
        color_dict[5] = "#d17de8"
        color_dict[6] = "#FFA533"
        color_dict[7] = "#905E22"

    # color_dict = {
    # 'Acidobacteria': '#e5daa3',
    # 'Cyanobacteria/Melainabacteria group': '#907bdd',
    # 'PVC group': '#79f6f0',
    # 'Nitrospirae': '#8eb979',
    # 'Proteobacteria': '#fd88f4',
    # 'Firmicutes': '#c87a86',
    # 'Chloroflexi': '#eef8fa',
    # 'FCB group': '#b6b1c8',
    # 'Actinobacteria': '#8bf1a7',
    # 'Type1': 'dodgerblue',
    # 'Type2b': 'gold',
    # 'Type2a': 'green',
    # 'Type3': 'purple'
    # }

    tree, ts = tc.get_example_tree(
        tree=tree,
        color_dict=color_dict,
        annot_dict=annot_dict,
        col=col,
        rotate=parser.rotate,
    )

    # Write to out path

    print("\nTree image has been written to " + parser.outpath)
    tree.render(parser.outpath, tree_style=ts, dpi=20)
