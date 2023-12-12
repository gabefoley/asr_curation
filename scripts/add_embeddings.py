import pandas as pd
import os
from transformers import BertModel, AutoTokenizer
import pandas as pd
import torch
import os


def calculate_embeddings(sequence, model, model_name, tokenizer):
    inputs = tokenizer(
        " ".join(sequence),
        return_tensors="pt",
        padding=True,
        truncation=True,
        max_length=512,
    )
    outputs = model(**inputs)

    last_hidden_state = outputs[0]
    embed_mean = last_hidden_state.mean(dim=1).detach().numpy().flatten().tolist()
    embed_cls = last_hidden_state[:, 0, :].detach().numpy().flatten().tolist()
    embed_max_pool = (
        last_hidden_state.max(dim=1).values.detach().numpy().flatten().tolist()
    )
    hidden_states = outputs[2]
    last_two_layers = torch.cat((hidden_states[-1], hidden_states[-2]), dim=-1)
    embed_last_two_concat = (
        last_two_layers.mean(dim=1).detach().numpy().flatten().tolist()
    )

    embeddings = {
        f"{model_name} Embed Mean": embed_mean,
        f"{model_name} Embed CLS": embed_cls,
        f"{model_name} Embed Max Pool": embed_max_pool,
        f"{model_name} Embed Last Two Concat": embed_last_two_concat,
    }

    # print (sequence)
    # print (embed_mean[0:5])
    # print (embed_cls[0:5])
    return embeddings


def process_and_store_embeddings(df, model_name, embedding_df_path="embeddings.pkl"):
    model = BertModel.from_pretrained(model_name, output_hidden_states=True)
    tokenizer = AutoTokenizer.from_pretrained(model_name)

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

        # print (existing_row)

        if not existing_row.empty:
            print("not empty")

            embed_data = existing_row.iloc[0].to_dict()
            for key, value in embed_data.items():
                # print (key)
                if key not in df.columns:
                    df[key] = [None] * len(df)
                df.at[idx, key] = value
        else:
            embeddings = calculate_embeddings(sequence, model, model_name, tokenizer)
            new_row = {"sequence": sequence, "model_name": model_name, **embeddings}
            embedding_df = pd.concat([embedding_df, pd.DataFrame([new_row])], ignore_index=True)

            # for key, value in embeddings.items():
            #     if key not in df.columns:
            #         df[key] = [None] * len(df)
            #     df.at[idx, key] = value

    # Saving the updated embedding dataframe
    embedding_df.to_pickle(embedding_df_path)

    return df
