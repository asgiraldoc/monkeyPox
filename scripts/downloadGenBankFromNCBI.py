import argparse
from Bio import Entrez, SeqIO


def read_id_list_from_file(file_path):
    with open(file_path, "r") as f:
        id_list = [line.strip() for line in f]
    return id_list


def download_genbank_files(email, api_key, id_list):
    Entrez.email = email
    Entrez.api_key = api_key

    for seq_id in id_list:
        try:
            # Fetch the GenBank file
            handle = Entrez.efetch(
                db="nucleotide", id=seq_id, rettype="gb", retmode="text"
            )
            record = SeqIO.read(handle, "genbank")
            handle.close()

            # Write the GenBank file
            with open(f"{seq_id}.gb", "w") as f:
                SeqIO.write(record, f, "genbank")
        except Exception as e:
            print(f"An error occurred with ID {seq_id}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download GenBank files from NCBI")
    parser.add_argument("--email", required=True, help="Email address")
    parser.add_argument("--api_key", required=True, help="API key for NCBI")
    parser.add_argument(
        "--input_file", required=True, help="Input file containing sequence IDs"
    )

    args = parser.parse_args()

    id_list = read_id_list_from_file(args.input_file)
    download_genbank_files(args.email, args.api_key, id_list)
