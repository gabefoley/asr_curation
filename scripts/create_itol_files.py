import pandas as pd
import distinctipy
import configs.itol_text as itol_text
import os

def get_color_dict_and_color_list(df, col):
    df = df.fillna("None")

    # Get unique values from the 'funfam_specific_summary' column
    unique_values = df[col].unique()

    print (unique_values)

    # Select an appropriate ColorBrewer color set based on the number of unique values
    color_set = [
        distinctipy.get_hex(x) for x in distinctipy.get_colors(len(unique_values))
    ]

    print (color_set)

    # Create a dictionary that maps unique keys to ColorBrewer colors
    color_dict = {
        value: color_set[i % len(color_set)] for i, value in enumerate(unique_values)
    }

    print (color_dict)

    colors = distinctipy.get_colors(len(unique_values))

    # Create a list that combines the 'info' value with the corresponding color
    color_list = [(df.at[i, "info"], df.at[i, col]) for i in df.index]

    return color_dict, color_list


def generate_itol_colorstrip(col, color_dict, color_list, output_filename):
    with open(output_filename, "w") as f:
        f.write(itol_text.dataset_colorstrip_text)
        for info, label in color_list:
            print (info)
            print (label)
            f.write(f"{info} {color_dict[label]} {label} \n")

    print(f"File '{output_filename}' has been created.")


def generate_itol_ranges(col, color_dict, color_list, output_filename):
    # START_NODE_ID,END_NODE_ID,FILL_COLOR,GRADIENT_COLOR,LINE_COLOR,LINE_STYLE,LINE_WIDTH,LABEL_TEXT,LABEL_COLOR,LABEL_SIZE_FACTOR,LABEL_STYLE
    # 9606,184922,#ffffff,#ff0000,#000000,dashed,2,Example range,#0000ff,1,italic

    with open(output_filename, "w") as f:
        f.write(itol_text.dataset_ranges_text)
        for info, label in color_list:
            f.write(
                f"{info},{info},{color_dict[label]},{color_dict[label]},{color_dict[label]},dashed,2,{label},black,italic\n"
            )

    print(f"File '{output_filename}' has been created.")


def generate_itol_shape(df):
    output_filename = "shapes.txt"

    selected_columns = df.columns[
        df.columns.get_loc("Activity_Zn_PnP") : df.columns.get_loc("Activity_Mn_4NPS")
        + 1
    ]
    selected_columns = selected_columns.insert(0, "UniProt")

    with open(output_filename, "w") as f:
        f.write(itol_text.dataset_shape_text)

        # Select columns from 'Activity_Zn_PnP' to 'Activity_Mn_4NPS' (inclusive)

        # Iterate through the rows of the DataFrame
        for index, row in df[selected_columns].iterrows():
            # Create a list of values where None is replaced with "None"
            row_values = ["-" if value is None else str(value) for value in row]

            if row_values.count("nan") < 5:
                # Join the row values into a single line with commas
                row_line = ", ".join(row_values)

                # Print the row line
                f.write(f"{row_line}\n")

    print(f"File '{output_filename}' has been created.")


def generate_dataset_style(col, color_dict, color_list, output_filename):
    output_filename = f"{col}_dataset_style.txt"
    # START_NODE_ID,END_NODE_ID,FILL_COLOR,GRADIENT_COLOR,LINE_COLOR,LINE_STYLE,LINE_WIDTH,LABEL_TEXT,LABEL_COLOR,LABEL_SIZE_FACTOR,LABEL_STYLE
    # 9606,184922,#ffffff,#ff0000,#000000,dashed,2,Example range,#0000ff,1,italic

    with open(output_filename, "w") as f:
        f.write(itol_text.dataset_style_text)
        for info, label in color_list:
            f.write(f"{info},branch,node,{color_dict[label]},1,normal\n")

    print(f"File '{output_filename}' has been created.")


if __name__ == "__main__":
    df = pd.read_csv(snakemake.input.csv)
    annotation_cols = snakemake.params.annotation_cols

    print ('Creating ITOL files')
    print (snakemake.output.tsv)

    outpath = os.path.dirname(snakemake.output.tsv)

    print (outpath)

    with open (snakemake.output.tsv, "w+") as output_text:
        for col in annotation_cols:

            # Skip the info column, which won't be informative
            if col != 'info':

                print (col)

                color_dict, color_list = get_color_dict_and_color_list(df, col)


                print ('color dict')
                print (color_dict)

                print ('color list')
                print (color_list)

                generate_itol_colorstrip(
                    col, color_dict, color_list, f"{outpath}/{col}itol_colorstrip.txt"
                )


                output_text.write(f"Colorstrip generated for {col} at {outpath}/{col}itol_colorstrip.txt ")

                generate_itol_ranges(
                    col, color_dict, color_list, f"{outpath}/{col}_itol_ranges.txt"
                )

                output_text.write(f"Color list generated for {col} at {outpath}/{col}_itol_ranges.txt ")


                generate_dataset_style(
                    col, color_dict, color_list, f"{outpath}/{col}itol_dataset_style.txt"
                )

                output_text.write(f"Dataset style generated for {col} at {outpath}/{col}itol_dataset_style.txt ")



