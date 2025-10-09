import pandas as pd
import os
file = "/Users/reem/Mov/nextclade_results/final_results.tsv"
#df =pd.read_csv(file, delimiter="\t", low_memory=False)
#print(df["seqName"])

def GtoA_percent(muts):
    mut_list=muts.split(",")
    total=len(mut_list)
    total_muts=0
    G_to_A=0
    for m in mut_list:
        total_muts+=1
        if m.startswith("G") and m.endswith("A"):
            G_to_A+=1
    if total > 0:
        return((G_to_A/total_muts) * 100)

chunk_size=1000
chunk_reader=pd.read_csv(file, delimiter="\t", chunksize=chunk_size, low_memory=False)

for chunk in chunk_reader:
    df=chunk
    chunk= pd.DataFrame(chunk)
    chunk["G_to_A_%"] = chunk["privateNucMutations.unlabeledSubstitutions"].fillna("").apply(GtoA_percent)
    #print(chunk[["seqName", "privateNucMutations.unlabeledSubstitutions", "G_to_A_%"]])
    #break
os.makedirs('~/Documents', exist_ok=True)
chunk.to_csv('~/Documents/final.csv')