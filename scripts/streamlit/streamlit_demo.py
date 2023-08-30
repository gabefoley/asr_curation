import streamlit as st
import pandas as pd
import streamlit_pandas as sp
import os

@st.cache_data
def load_data():
    df = pd.read_csv(file)
    return df
# file = "./workflows/example_workflow/datasets/als_example_ec_2_2_1_6/csv/custom/als_example_ec_2_2_1_6_generic_annotated.csv"
file = "./scripts/streamlit/streamlit.csv"

df = load_data()
create_data = {"info": "text",
                "gene_names": "multiselect",
                "xref_interpro": "multiselect"}

# all_widgets = sp.create_widgets(df, create_data, ignore_columns=["PassengerId"])
all_widgets = sp.create_widgets(df, create_data)

res = sp.filter_df(df, all_widgets)
st.title("Streamlit AutoPandas")
st.header("Original DataFrame")
st.write(df)

st.header("Result DataFrame")
st.write(res)

