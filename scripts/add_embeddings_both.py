import pandas as pd
import os
from transformers import BertModel, AutoTokenizer, T5EncoderModel, T5Tokenizer
import torch
import gc
from tm_vec.embed_structure_model import trans_basic_block, trans_basic_block_Config
from tm_vec.tm_vec_utils import featurize_prottrans, embed_tm_vec, encode
from sklearn.manifold import TSNE

def calculate_t5_embeddings(sequences, model_deep, model, tokenizer, device):
    encoded_sequences = encode(sequences, model_deep, model, tokenizer, device)
    # Process and return T5 embeddings as needed
    return encoded_sequences

def process_and_store_embeddings(df, model_name, embedding_df_path, calculate_embeddings_func):
    model = None
    tokenizer = None

    if "DistilProtBert" in model_name:
        model = BertModel.from_pretrained(model_name, output_hidden_states=True)
        tokenizer = AutoTokenizer.from_pretrained(model_name)
    elif "prot_t5_xl_uniref50" in model_name:
        model = T5EncoderModel.from_pretrained(model_name)
        tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)
    
    if model is None or tokenizer is None:
        print("Unsupported model name.")
        return df

    if os.path.exists(embedding_df_path):
        embedding_df = pd.read_pickle(embedding_df_path)
    else:
        embedding_df = pd.DataFrame(columns=["sequence", "model_name"])

    for idx, row in df.iterrows():
        print(idx)
        sequence = row["sequence"]

        existing_row = embedding_df[
            (embedding_df["sequence"] == sequence)
            & (embedding_df["model_name"] == model_name)
        ]

        if not existing_row.empty:
            print("not empty")
            embed_data = existing_row.iloc[0].to_dict()
            for key, value in embed_data.items():
                if key not in df.columns:
                    df[key] = [None] * len(df)
                df.at[idx, key] = value
        else:
            embeddings = calculate_embeddings_func(sequence, model, tokenizer)
            
            new_row = {"sequence": sequence, "model_name": model_name, **embeddings}
            embedding_df = pd.concat([embedding_df, pd.DataFrame([new_row])], ignore_index=True)

    embedding_df.to_pickle(embedding_df_path)
    return df

# Example usage for BERT:
bert_model_name = "yarongef/DistilProtBert"
bert_embedding_df_path = "bert_embeddings.pkl"
process_and_store_embeddings(annot_df, bert_model_name, bert_embedding_df_path, calculate_bert_embeddings)

# Example usage for T5:
t5_model_name = "Rostlab/prot_t5_xl_uniref50"
t5_embedding_df_path = "t5_embeddings.pkl"
process_and_store_embeddings(annot_df, t5_model_name, t5_embedding_df_path, calculate_t5_embeddings)
