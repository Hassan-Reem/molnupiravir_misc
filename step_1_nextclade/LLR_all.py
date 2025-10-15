import pandas as pd
import numpy as np
from collections import Counter
from scipy.stats import multinomial
import subprocess
import tqdm
import csv


probs_df=pd.read_csv("/Users/reem/Downloads/estimated_mutation_distribution.tsv", delimiter="\t")
pM=probs_df["Molnupiravir"].to_numpy(dtype=float)
pN=probs_df["Normal"].to_numpy(dtype=float)
mut_types=probs_df["MutationType"].str.replace("→",">").tolist()


def get_mut_type(mut_string):
    return mut_string[0] + '>' + mut_string[-1]


def get_likelihood_ratio(counts,pM,pN):
    counts=np.array(counts,dtype=float)
    llM= float(multinomial.logpmf(counts, n=np.sum(counts), p=pM))
    llN = float(multinomial.logpmf(counts, n=np.sum(counts), p=pN))
    llr=llM-llN
    return llr

def process_chunk(chunk):
    chunk["subs"]=chunk["privateNucMutations.unlabeledSubstitutions"].apply(
    lambda x: ','.join([get_mut_type(item) for item in x.split(',')]) 
    if pd.notna(x) and x != '' 
    else '')
    chunk["Counts"] =chunk["subs"].apply(lambda x: Counter(x.split(",") if x else {}))
    llr_list = []
    for counts_dict in chunk["Counts"]:
        counts = [counts_dict.get(mt, 0) for mt in mut_types]
        llr = get_likelihood_ratio(counts,pM,pN)
        llr_list.append(llr)
    chunk["LLR"] = llr_list
    return chunk


input_tsv = "/Users/reem/Mov/nextclade_results/final_results.tsv"
output_tsv="/Users/reem/Mov/nextclade_results/ll_final_results.tsv"

first_chunk=True
chunk_size = 10000 

for chunk in pd.read_csv(input_tsv, delimiter="\t",chunksize=chunk_size):
    processed = process_chunk(chunk)
    processed.to_csv(output_tsv,sep="\t",index=False,mode="w" if first_chunk else "a",
                     header=first_chunk)
    first_chunk=False

print("All done")























