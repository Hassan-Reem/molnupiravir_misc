from Bio import SeqIO
import gzip
import subprocess
import tqdm
import pandas as pd

def read_fasta_batches(fasta_file, exclude_if_name_is_in, batch_size=1000):
    batch = []
    with open(fasta_file, 'rt') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.name in exclude_if_name_is_in:
                continue
            batch.append(record)
            if len(batch) >= batch_size:
                yield batch
                batch = []
    if batch:
        yield batch

def write_batch_to_fasta(batch, output_file):
    with open(output_file, 'w') as f:
        SeqIO.write(batch, f, 'fasta')


def append_results_from_a_to_b(a,b, starting_from_scratch):
    with open(a, 'r') as a_file:
        with open(b, 'a') as b_file:
                lines = a_file.readlines()
                if not starting_from_scratch:
                    lines = lines[1:]
                b_file.write("".join(lines))

target_filename = "final_results.tsv"
starting_from_scratch = False
seqNames = []
try:
    results_so_far = pd.read_csv(target_filename,delimiter="\t")
    seqNames = results_so_far['seqName'].tolist() 
except FileNotFoundError:
    starting_from_scratch = True

seqNames = set(seqNames)


filename = "/Users/Reem/Downloads/sequences.fasta"
batch_size = 10e3
very_total = 17e6
remaining = very_total - len(seqNames)
iterator = read_fasta_batches(filename, seqNames, batch_size= batch_size)
for batch in tqdm.tqdm(iterator, total=remaining / batch_size):
    write_batch_to_fasta(batch, "batch.fa")
    subprocess.check_call("nextclade run batch.fa --output-tsv temp.tsv --dataset-name sars-cov-2", shell=True)
    append_results_from_a_to_b("temp.tsv", target_filename, starting_from_scratch)
    starting_from_scratch = False