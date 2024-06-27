import pandas as pd
import os
from transformers import AutoModel, AutoTokenizer
import torch


def calculate_embeddings(sequence, model, tokenizer, model_type):
    """Calculate various embeddings for a given sequence."""
    inputs = tokenizer(" ".join(sequence), return_tensors='pt', padding=True, truncation=True)
    with torch.no_grad():
        if model_type == 'protbert':
            outputs = model(**inputs)
        elif model_type == 't5':
            outputs = model(**inputs.input_ids)
        else:
            raise ValueError(f"Unsupported model type: {model_type}")

    embeddings = outputs.last_hidden_state

    # Mean pooling
    mean_embedding = embeddings.mean(dim=1).squeeze().numpy()

    # CLS token pooling (if applicable, using first token)
    cls_embedding = embeddings[:, 0].squeeze().numpy()

    # Max pooling
    max_embedding = embeddings.max(dim=1).values.squeeze().numpy()

    # Weighted pooling
    weights = torch.linspace(0.1, 1.0, embeddings.size(1), device=embeddings.device)
    weights = weights.unsqueeze(0).unsqueeze(-1)  # Add extra dimensions for broadcasting
    weighted_embedding = (embeddings * weights).mean(dim=1).squeeze().numpy()

    return {
        f'{model_type}_mean_embedding': mean_embedding,
        f'{model_type}_cls_embedding': cls_embedding,
        f'{model_type}_max_embedding': max_embedding,
        f'{model_type}_weighted_embedding': weighted_embedding
    }


def process_and_store_embeddings(df, model_name, embedding_df_path, model_type):
    """Process and store multiple types of embeddings for sequences in the DataFrame."""
    model = AutoModel.from_pretrained(model_name, output_hidden_states=True)
    tokenizer = AutoTokenizer.from_pretrained(model_name)

    # Load existing embeddings if they exist
    if os.path.exists(embedding_df_path):
        embedding_df = pd.read_pickle(embedding_df_path)
    else:
        embedding_df = pd.DataFrame(columns=["info", "sequence", "model_name"])

    for idx, row in df.iterrows():
        info = row['info']
        sequence = row["sequence"]

        existing_row = embedding_df[
            (embedding_df["sequence"] == sequence)
            & (embedding_df["model_name"] == model_name)
        ]

        if not existing_row.empty:
            print(f"Embeddings for sequence {sequence} already exist.")
            continue  # Skip if embeddings for this sequence already exist

        try:
            embeddings = calculate_embeddings(sequence, model, tokenizer, model_type)
            new_row = {"info": info, "sequence": sequence, "model_name": model_name, **embeddings}
            embedding_df = pd.concat([embedding_df, pd.DataFrame([new_row])], ignore_index=True)

        except Exception as e:
            print(f"Failed to process sequence {sequence} with error: {e}")

    # Save embedding_df with full embeddings
    embedding_df.to_pickle(embedding_df_path)

    return embedding_df

