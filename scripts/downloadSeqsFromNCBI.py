import argparse
from Bio import Entrez, SeqIO

def read_id_list_from_file(file_path):
    with open(file_path, "r") as f:
        id_list = [line.strip() for line in f]
    return id_list

def download_sequences(email, api_key, id_list):
    Entrez.email = email
    Entrez.api_key = api_key

    for id in id_list:
        # Fetch the sequence
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()

        # Write the sequence to a .fasta file
        SeqIO.write(record, f"{id}.fasta", "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download sequences from NCBI')
    parser.add_argument('--email', required=True, help='Email address')
    parser.add_argument('--api_key', required=True, help='API key for NCBI')
    parser.add_argument('--input_file', required=True, help='Input file containing sequence IDs')

    args = parser.parse_args()

    id_list = read_id_list_from_file(args.input_file)
    download_sequences(args.email, args.api_key, id_list)
