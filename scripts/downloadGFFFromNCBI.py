import argparse
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm


def read_id_list_from_file(file_path):
    with open(file_path, "r") as f:
        id_list = [line.strip() for line in f]
    return id_list


def download_gff_file(seq_id, email, api_key, pbar):
    try:
        Entrez.email = email
        Entrez.api_key = api_key
        handle = Entrez.efetch(
            db="nucleotide", id=seq_id, rettype="gff3", retmode="text"
        )
        gff_data = handle.read()
        handle.close()

        with open(f"{seq_id}.gff", "w") as f:
            f.write(gff_data)

        pbar.update(1)
    except Exception as e:
        print(f"An error occurred with ID {seq_id}: {e}")


def download_gff_files(email, api_key, id_list, num_threads=10):
    with ThreadPoolExecutor(max_workers=num_threads) as executor, tqdm(
        total=len(id_list)
    ) as pbar:
        futures = {
            executor.submit(download_gff_file, seq_id, email, api_key, pbar): seq_id
            for seq_id in id_list
        }
        for future in futures:
            future.result()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download GFF files from NCBI")
    parser.add_argument("--email", required=True, help="Email address")
    parser.add_argument("--api_key", required=True, help="API key for NCBI")
    parser.add_argument(
        "--input_file", required=True, help="Input file containing sequence IDs"
    )
    parser.add_argument(
        "--num_threads", type=int, default=10, help="Number of threads to use"
    )

    args = parser.parse_args()

    id_list = read_id_list_from_file(args.input_file)
    download_gff_files(args.email, args.api_key, id_list, args.num_threads)
