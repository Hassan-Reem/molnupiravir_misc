from Bio import SeqIO
import gzip
import subprocess
import tqdm
import pandas as pd
import csv

"""
Extracts FASTA sequences not present in the UShER phylogenetic tree based on
a list of unmatched sequence identifiers (unmatched.txt), and outputs them in batches for
subsequent addition using UShER.
"""


def read_fasta_batches(fasta_file, prefixes, batch_size=1000):
    """Read sequences from a compressed FASTA file in batches, matching those in the unmatched set."""
    batch = []
    cmd = 'xz -d -T0 -c "/Users/reem/Downloads/sequences_fasta_2025_09_21.tar.xz" | tar -xOf - sequences.fasta'

    stream =  subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True 
    ) 
    for record in SeqIO.parse(stream.stdout, 'fasta'):
            seq_id = record.name.split('|')[0]
            if seq_id.startswith("hCoV-19/"):
                seq_id =  seq_id[len("hCoV-19/"):]
                if seq_id in prefixes:
                    batch.append(record)
                    #print(f"{seq_id}")
                if len(batch) >= batch_size:
                    yield batch
                    batch = []
            
    if batch:
        yield batch


fasta_file = "/Users/reem/Downloads/sequences_fasta_2025_09_21.tar.xz"
output_file = "/Users/reem/missing_seqs618.fasta"
batch_size = 100


prefixes = set()
with open("/Users/reem/mov_not_in_usher.txt", 'r') as f:
    for line in f:
        name = line.strip()
        if name.startswith("hCoV-19/"):
            name = name[len("hCoV-19/"):]
        prefixes.add(name.split('|')[0])

    
first = True
iterator = read_fasta_batches(fasta_file, prefixes, batch_size=batch_size)
for batch in tqdm.tqdm(iterator):
    mode = 'w' if first else 'a'
    with open(output_file, mode) as f:
        SeqIO.write(batch, f, 'fasta')
    first = False