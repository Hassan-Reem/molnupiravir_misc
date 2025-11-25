from Bio import SeqIO
import gzip
import subprocess
import tqdm
import pandas as pd
import csv


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
            seq_id = record.name.split('|', 1)[0]
            if any(prefix in seq_id for prefix in prefixes):
                batch.append(record)
                if len(batch) >= batch_size:
                    yield batch
                    batch = []
    if batch:
        yield batch


fasta_file = "/Users/reem/Downloads/sequences_fasta_2025_09_21.tar.xz"
output_file = "/Users/reem/unmatched_sequences.fasta"
batch_size = 100



with open("/Users/reem/unmatched_names.txt", 'r') as f:
    prefixes = set(line.strip() for line in f)

iterator = read_fasta_batches(fasta_file, prefixes, batch_size=batch_size)
for batch in tqdm.tqdm(iterator):
    with open(output_file, 'a') as f:
        SeqIO.write(batch, f, 'fasta')




