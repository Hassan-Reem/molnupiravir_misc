"""This script calculates LLRs per substitution context (C>T, T>C, A>G) for each sequence based on their mutational spectra.
   It uses a threshold LLR of 6 to separate sequences into two groups (likely Molnupiravir-treated vs. likely non-treated), 
   computes mean proportions for each group, and then calculates LLRs for each sequence using these means.
   (Similar to All_llrs.ipynb but as a standalone script)."""


import pandas as pd
import os
import tqdm
from collections import Counter
import numpy as np
from scipy.stats import multinomial
import numpy as np
import ast


# Load the data iteratively to avoid crashing the kernel
iterator = pd.read_csv("/Users/reem/Mov/nextclade_results/LLR_final_results.tsv",sep="\t", usecols=['seqName','LLR','privateNucMutations.unlabeledSubstitutions'], dtype={"seqName":str, "LLR": np.float64, 'privateNucMutations.unlabeledSubstitutions':str } , chunksize=1000)
df = pd.concat([chunk for chunk in tqdm.tqdm(iterator, desc='Loading data')])


#Get Counts per substitution context
def count_GtoA(spectrum):
    counts = Counter()
    muts = spectrum.split(",")
    for mut in muts:
        if mut[2:5] == 'G>A':
            counts[mut]+=1
    return counts

def count_CtoT(spectrum):
    counts = Counter()
    muts = spectrum.split(",")
    for mut in muts:
        if mut[2:5] == 'C>T':
            counts[mut]+=1
    return counts

def count_AtoG(spectrum):
    counts = Counter()
    muts = spectrum.split(",")
    for mut in muts:
        if mut[2:5] == 'A>G':
            counts[mut]+=1
    return counts

df["G>A_counts"] = df["spectrum"].apply(count_GtoA)
df["C>T_counts"] = df["spectrum"].apply(count_CtoT)
df["A>G_counts"] = df["spectrum"].apply(count_AtoG)

# Get proportions per substitution context
def get_proportion(df):
    dict = {}
    total = sum(df.values())
    for key, value in df.items():
        dict[key] = value/total
       
    return dict

df["G>A_proportions"] = df.apply(lambda row: get_proportion(row["G>A_counts"]), axis=1)

df["C>T_proportions"] = df.apply(lambda row: get_proportion(row["C>T_counts"]), axis=1)

df["A>G_proportions"] = df.apply(lambda row: get_proportion(row["A>G_counts"]), axis=1)

#Load Molnupiravir_df (LLR>6)
df_high_llrs = df[df["LLR"]>6]
df_high_llrs.tail()
#Load Normal_df (LLR<6)
df_low_llrs = df[df["LLR"]>6]
df_low_llrs.tail()

#Create a pivot table for each mutational context per substitution:
proportions_list_high = pd.json_normalize(df_high_llrs['C>T_proportions'])
proportions_list_high.index = df_high_llrs.index 
proportions_list_high.head()
#Concatenate seqName and LLR from old df to pivot table
df_Mov = pd.concat([df_high_llrs[["seqName"]],proportions_list_high,df_high_llrs[["LLR"]]], axis=1)
df_Mov.tail()

#Calculate mean for every mutational context:
# calculate mean values for each G to A context and append the mean row to df_Mov.
context_cols = [col for col in df_Mov.columns if col not in ["seqName", "LLR"]]
mean_row = df_Mov[context_cols].fillna(0).mean()
mean_row['seqName'] = 'Mean'
#print(mean_row)
probs_df = pd.melt(df_Mov, id_vars=context_cols,
                       var_name = 'Molnupiravir',
                       value_name = mean_row)




df_Mov = pd.concat([df_Mov, pd.DataFrame([mean_row])], ignore_index=True)
print(df_Mov.tail())

df_Mov.to_csv("Mov_means.tsv",sep = "\t")

# Repeat the same process for LLRs < 6 (Likely Non-Mov seqs)
proportions_list_low = pd.json_normalize(df_low_llrs['C>T_proportions'])
proportions_list_low.index = df_low_llrs.index 
proportions_list_low.head()

df_Normal = pd.concat([df_low_llrs[["seqName"]],proportions_list_low,df_low_llrs[["LLR"]]], axis=1)
#df_Normal.head()
context_cols = [col for col in df_Normal.columns if col not in ["seqName", "LLR"]]
df_Normal[context_cols].dtypes
mean_row = df_Normal[context_cols].fillna(0).mean()
mean_row['seqName'] = 'Mean'
probs_df['Normal'] = df_Normal.melt(mean_row)

df_Normal = pd.concat([df_Normal, pd.DataFrame([mean_row])], ignore_index=True)
print(df_Normal.head())


df_Normal.to_csv("Normal_means.tsv",sep = "\t")

#Calculate LLRs using means:

counts=df[context_cols].values
pM=probs_df["Molnupiravir"].values.tolist()
pN=probs_df["Normal"].values.tolist()


llM= float(multinomial.logpmf(counts, n=np.sum(counts), p=pM))
llN = float(multinomial.logpmf(counts, n=np.sum(counts), p=pN))
llr=llM-llN
print(llr)







