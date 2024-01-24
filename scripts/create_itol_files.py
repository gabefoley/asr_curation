import pandas as pd
import distinctipy
import configs.itol_text as itol_text
import os

def get_color_dict_and_info_list(df, col):
    df = df.fillna("None")

    # Get unique values from the 'funfam_specific_summary' column
    unique_values = df[col].unique()

    print (unique_values)

    # Select an appropriate ColorBrewer color set based on the number of unique values
    color_set = [
        distinctipy.get_hex(x) for x in distinctipy.get_colors(len(unique_values))
    ]

    # Create a dictionary that maps unique keys to ColorBrewer colors
    color_dict = {
        value: color_set[i % len(color_set)] for i, value in enumerate(unique_values)
    }

    print (color_dict)

    colors = distinctipy.get_colors(len(unique_values))

    # Create a list that combines the 'info' value with the corresponding color
    info_list = [(df.at[i, "info"], df.at[i, col]) for i in df.index]

    return color_dict, info_list


def get_single_info_dict(df, col):
    df = df.fillna("None")

    # Get unique values from the 'funfam_specific_summary' column
    unique_values = df[(df[col] != "None") & (df[col] != False)][col].unique()

    info_dict = {info: value for value in unique_values for info in df[df[col] == value]["info"].tolist()}

    return info_dict

def generate_itol_colorstrip(col, color_dict, info_list, output_filename):
    with open(output_filename, "w") as f:
        f.write(itol_text.dataset_colorstrip_text.replace("<custom_dataset_label>", col))
        for info, label in info_list:
            print (info)
            print (label)
            f.write(f"{info} {color_dict[label]} {label} \n")

    print(f"File '{output_filename}' has been created.")


def generate_itol_ranges(col, color_dict, info_list, output_filename):
    # START_NODE_ID,END_NODE_ID,FILL_COLOR,GRADIENT_COLOR,LINE_COLOR,LINE_STYLE,LINE_WIDTH,LABEL_TEXT,LABEL_COLOR,LABEL_SIZE_FACTOR,LABEL_STYLE
    # 9606,184922,#ffffff,#ff0000,#000000,dashed,2,Example range,#0000ff,1,italic

    with open(output_filename, "w") as f:
        f.write(itol_text.dataset_ranges_text.replace("<custom_dataset_label>", col))
        for info, label in info_list:
            f.write(
                f"{info},{info},{color_dict[label]},{color_dict[label]},{color_dict[label]},dashed,2,{label},black,italic\n"
            )

    print(f"File '{output_filename}' has been created.")

def get_lab_percentages(df):
    # Filter the DataFrame to include only rows with "Extant" in the "Enzyme Type" column
    extant_df = df[df['Enzyme Type'] == 'Extant']

    # Get the columns that start with "Activity"
    activity_columns = [col for col in extant_df.columns if col.startswith('Activity')]

    # Calculate the number of filled cells that start with "Activity" (excluding NaN values)
    filled_cells = extant_df[activity_columns].count().sum()

    # Calculate the total count of all cells (filled and unfilled) starting with "Activity"
    total_cells = 0

    for col in activity_columns:
        total_cells += extant_df[col].shape[0]

    percentage_filled = (filled_cells / total_cells) * 100 if total_cells > 0 else 0
    formatted_percentage = f"{percentage_filled:.2f}%"

    return filled_cells, total_cells, formatted_percentage



def generate_lab_assays(df, outpath):

    completed, total, percentage = get_lab_percentages(df)

    selected_columns = df.columns[
        df.columns.get_loc("Activity_Zn_PnP") : df.columns.get_loc("Activity_Mn_SLG")
        + 1
    ]
    selected_columns = selected_columns.insert(0, "UniProt")

    with open(outpath, "w+") as f:

        text_to_write = itol_text.lab_assays_text.replace("<completed>", str(completed))
        text_to_write = text_to_write.replace("<total>", str(total))
        f.write(text_to_write.replace("<percentage>", percentage))

        # Select columns from 'Activity_Zn_PnP' to 'Activity_Mn_4NPS' (inclusive)

        print (df[selected_columns])

        # Iterate through the rows of the DataFrame
        for index, row in df[selected_columns].iterrows():
            # Create a list of values where None is replaced with "None"
            row_values = ["-" if value is None else str(value).strip() for value in row]

            # if row_values.count("nan") < 5:
            #     # Join the row values into a single line with commas
            row_line = ", ".join(row_values)

                # Print the row line
            f.write(f"{row_line}\n")

    print(f"File '{outpath}' has been created.")


def generate_dataset_style(col, color_dict, info_list, output_filename):
    # output_filename = f"{col}_dataset_style.txt"
    # START_NODE_ID,END_NODE_ID,FILL_COLOR,GRADIENT_COLOR,LINE_COLOR,LINE_STYLE,LINE_WIDTH,LABEL_TEXT,LABEL_COLOR,LABEL_SIZE_FACTOR,LABEL_STYLE
    # 9606,184922,#ffffff,#ff0000,#000000,dashed,2,Example range,#0000ff,1,italic

    with open(output_filename, "w") as f:
        f.write(itol_text.dataset_style_text.replace("<custom_dataset_label>", col))
        for info, label in info_list:
            if label in color_dict:
                f.write(f"{info},branch,node,{color_dict[label]},1,normal\n")

    print(f"File '{output_filename}' has been created.")


def generate_shape_style(col, colour, info_dict, output_filename):

    with open(output_filename, "w") as f:

        text_to_write = itol_text.shape_text.replace("<custom_dataset_label>", col)
        text_to_write = text_to_write.replace("<col>", col)
        f.write(text_to_write.replace("<custom_color>", colour))


        for info, val in info_dict.items():
            f.write(f"{info},1\n")

    print(f"File '{output_filename}' has been created.")


if __name__ == "__main__":
    df = pd.read_csv(snakemake.input.csv)
    annotation_cols = snakemake.params.annotation_cols
    single_colour_annotation_cols = snakemake.params.single_colour_annotation_cols

    print ('Creating ITOL files')
    print (snakemake.output.tsv)

    outpath = os.path.dirname(snakemake.output.tsv)

    print (outpath)

    with open (snakemake.output.tsv, "w+") as output_text:

        output_text.write(f"Lab Assays written to {outpath}/lab_assays.txt\n")

        for col in annotation_cols:

            # Skip the info column, which won't be informative
            if col != 'info' and col in df:

                print (col)

                color_dict, info_list = get_color_dict_and_info_list(df, col)
                print ('color dict')
                print (color_dict)

                print ('color list')
                print (info_list)

                generate_itol_colorstrip(
                    col, color_dict, info_list, f"{outpath}/{col}itol_colorstrip.txt"
                )


                output_text.write(f"Colorstrip generated for {col} at {outpath}/{col}itol_colorstrip.txt ")

                generate_itol_ranges(
                    col, color_dict, info_list, f"{outpath}/{col}_itol_ranges.txt"
                )

                output_text.write(f"Color list generated for {col} at {outpath}/{col}_itol_ranges.txt\n")


        # single_colours = distinctipy.get_hex(distinctipy.get_colors(len(single_colour_annotation_cols), pastel_factor=0.7))
        single_colours = [
            distinctipy.get_hex(x) for x in distinctipy.get_colors(len(single_colour_annotation_cols), pastel_factor=0.7)
        ]
        print ('rasperyy')

        print (single_colours)
        idx = 0
        for single_col in single_colour_annotation_cols:
            print ('here is single col')
            print (single_col)

            if single_col in df.columns:
                info_dict = get_single_info_dict(df, single_col)


                print (single_colours[idx])

                generate_shape_style(
                    single_col, single_colours[idx], info_dict, f"{outpath}/{single_col}_itol_dataset_style.txt"
                )

                idx +=1

                output_text.write(f"Dataset style generated for {single_col} at {outpath}/{single_col}_itol_dataset_style.txt\n")





