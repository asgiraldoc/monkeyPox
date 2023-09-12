import argparse
import json
from Bio import Entrez

def read_id_list_from_file(file_path):
    with open(file_path, "r") as f:
        id_list = [line.strip() for line in f]
    return id_list

def download_metadata(email, api_key, id_list):
    Entrez.email = email
    Entrez.api_key = api_key

    for id in id_list:
        # Fetch the metadata
        handle = Entrez.esummary(db="nucleotide", id=id)
        record = Entrez.read(handle)
        handle.close()

        # Write the metadata to a .json file
        with open(f"{id}_metadata.json", "w") as f:
            json.dump(record[0], f, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download metadata from NCBI')
    parser.add_argument('--email', required=True, help='Email address')
    parser.add_argument('--api_key', required=True, help='API key for NCBI')
    parser.add_argument('--input_file', required=True, help='Input file containing sequence IDs')

    args = parser.parse_args()

    id_list = read_id_list_from_file(args.input_file)
    download_metadata(args.email, args.api_key, id_list)
