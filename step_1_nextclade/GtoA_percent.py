import pandas as pd
import os
file = "/Users/reem/Mov/final_llrs_with_sum.tsv"
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

processed_chunks = []
for chunk in chunk_reader:
    chunk= pd.DataFrame(chunk)
    chunk["G_to_A_%"] = chunk["privateNucMutations.unlabeledSubstitutions"].fillna("").apply(GtoA_percent)
    processed_chunks.append(chunk)

final_dataframe = pd.concat(processed_chunks, ignore_index=True)
final_dataframe.to_csv('/Users/reem/Documents/gtoa.csv')