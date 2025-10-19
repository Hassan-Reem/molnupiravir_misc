import pandas as pd
import os
from Bio import Entrez, SeqIO
from collections import Counter
import numpy as np
from scipy.stats import multinomial
import numpy as np

file="/Users/reem/Mov/nextclade_results/first_100.tsv"
df=pd.read_csv(file,delimiter="\t")
countries_list=df["seqName"].values.tolist()


def calculate_country_counts(group):
    """Calculate mutation counts for a group (country)"""
    total_counts = Counter()
    for counts_dict in group["Counts"]:
        if isinstance(counts_dict, dict):
            for mut, count in counts_dict.items():
                total_counts[mut] += count
    return dict(total_counts)

# Get mutation types
probs_df = pd.read_csv("/Users/reem/Downloads/estimated_mutation_distribution.tsv", delimiter="\t")
mutation_types = probs_df["MutationType"].apply(lambda x: x.replace('→', '>')).tolist()

# Group by country and calculate
country_results = []

for country, group in df.groupby("country"):
    # Get counts for this country
    country_counts = calculate_country_counts(group)
    
    # Convert to ordered list
    counts_list = [country_counts.get(mut, 0) for mut in mutation_types]
    
    # Calculate LLR
    def get_likelihood_ratio(counts,pM,pN):
        counts = np.array(counts)
        llM= float(multinomial.logpmf(counts, n=np.sum(counts), p=pM))
        llN = float(multinomial.logpmf(counts, n=np.sum(counts), p=pN))
        llr=llM-llN
        return llr
    pM=probs_df["Molnupiravir"].values.tolist()
    pN=probs_df["Normal"].values.tolist()

    llr = get_likelihood_ratio(counts_list, pM, pN)
    
    country_results.append({
        'country': country,
        'llr': llr,
        'total_mutations': sum(counts_list),
        'num_sequences': len(group),
        'mean_ga_percentage': group['GA_percentage'].mean() if 'GA_percentage' in group.columns else np.nan
    })

# Create results dataframe
results_df = pd.DataFrame(country_results)
results_df = results_df.sort_values('llr', ascending=False)

print(results_df.head(20))
        
