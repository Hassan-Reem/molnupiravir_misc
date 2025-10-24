from Bio import SeqIO
import gzip
import subprocess
import tqdm
import pandas as pd
import csv


def read_fasta_batches(fasta_file, exclude_if_name_is_in, batch_size=1000):
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

target_filename = "/Users/reem/Mov/final_results.tsv"
starting_from_scratch = False

seqNamesAlreadyProcessed=set()

try:
    with open(target_filename, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seqNamesAlreadyProcessed.add(row['seqName'])
except FileNotFoundError:
    starting_from_scratch = True


filename = "/Users/Reem/Downloads/sequences.fasta"
batch_size = 1000
very_total = 17e6
remaining = very_total - len(seqNamesAlreadyProcessed)

iterator = read_fasta_batches(filename, seqNamesAlreadyProcessed, batch_size= batch_size)
for batch in tqdm.tqdm(iterator, total= remaining / batch_size):
    write_batch_to_fasta(batch, "batch.fa")
    subprocess.check_call("nextclade run batch.fa --output-tsv temp.tsv --dataset-name sars-cov-2", shell=True)
    append_results_from_a_to_b("temp.tsv", target_filename, starting_from_scratch)
    for record in batch:
        seqNamesAlreadyProcessed.add(record.name)
    starting_from_scratch = False

    