import torch
from transformers import T5EncoderModel, T5Tokenizer
import re
import gc
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset
import faiss
from tm_vec.embed_structure_model import trans_basic_block, trans_basic_block_Config
from tm_vec.tm_vec_utils import featurize_prottrans, embed_tm_vec, encode
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import distinctipy as dp
from tqdm import tqdm
import os

def calculate_embeddings(sequences, model_name, model_checkpoint, config):

    tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False )
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
    gc.collect()

    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    model = model.eval()

    #TM-Vec model paths
    tm_vec_model_cpnt = model_checkpoint
    tm_vec_model_config = config


    #Load the TM-Vec model
    tm_vec_model_config = trans_basic_block_Config.from_json(tm_vec_model_config)
    model_deep = trans_basic_block.load_from_checkpoint(tm_vec_model_cpnt, config=tm_vec_model_config)
    model_deep = model_deep.to(device)
    model_deep = model_deep.eval()

    encoded_sequences = encode(sequences, model_deep, model, tokenizer, device)
    encoded_sequences = torch.tensor(encoded_sequences)

    mean_embedding = encoded_sequences.mean(dim=1)
    cls_embedding = encoded_sequences[:, 0]

    
    embeddings = {
        f"{model_name} Embed Encoded": encoded_sequences.tolist(),
        f"{model_name} Embed Mean": mean_embedding.tolist(),
        f"{model_name} Embed CLS": cls_embedding.tolist(),
    }

    return embeddings

def process_and_store_embeddings(df, model_name, model_checkpoint, model_config, embedding_df_path="embeddings.pkl"):
    if os.path.exists(embedding_df_path):
        embedding_df = pd.read_pickle(embedding_df_path)
    else:
        embedding_df = pd.DataFrame(columns=["info", "sequence", "model_name"])

    new_sequences = []

    for idx, row in df.iterrows():
        info = row['info']
        sequence = row["sequence"]

        existing_row = embedding_df[
            (embedding_df["sequence"] == sequence)
            & (embedding_df["model_name"] == model_name)
        ]

        if not existing_row.empty:
            embed_data = existing_row.iloc[0].to_dict()
            for key, value in embed_data.items():
                if key not in df.columns:
                    df[key] = [None] * len(df)
                df.at[idx, key] = value
        else:
            new_sequences.append((sequence, info)) 

    if new_sequences:
        embeddings = calculate_embeddings(new_sequences, model_name, model_checkpoint, model_config)
        
        new_embedding_data = []

        for idx, sequence_info in enumerate(new_sequences):
            sequence, info = sequence_info
            row_data = {"info": info, "sequence": sequence, "model_name": model_name}
            for embedding_type, embedding_vals in embeddings.items():
                embedding_val = embedding_vals[idx] if idx < len(embedding_vals) else None
                row_data[embedding_type] = embedding_val
            new_embedding_data.append(row_data)

        new_embedding_df = pd.DataFrame(new_embedding_data)

        # Merge the original DataFrame with the new embedding DataFrame based on 'sequence' column
        df = df.merge(new_embedding_df, on='sequence', how='left')

        # Concatenate the original embedding DataFrame with the new embedding DataFrame
        embedding_df = pd.concat([embedding_df, new_embedding_df], ignore_index=True)

        # Saving the updated embedding dataframe
        embedding_df.to_pickle(embedding_df_path)

    return df





