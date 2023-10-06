import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def extract_sequences(msa_file, id_file, output_file):
    # Read the MSA file
    alignment = AlignIO.read(msa_file, "fasta")

    # Read the IDs of the sequences to be extracted
    with open(id_file, "r") as f:
        ids = set([line.strip() for line in f])

    # Filter the MSA to only contain sequences with IDs in the list
    new_alignment = MultipleSeqAlignment(
        [record for record in alignment if record.id in ids]
    )

    # Write the new MSA to the output file
    AlignIO.write(new_alignment, output_file, "fasta")
    print(f"New MSA file created with {len(new_alignment)} sequences: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract sequences from an MSA file based on a list of IDs and create a new MSA file."
    )

    parser.add_argument(
        "-m",
        "--msa_file",
        type=str,
        required=True,
        help="Input MSA file in FASTA format.",
    )
    parser.add_argument(
        "-i",
        "--id_file",
        type=str,
        required=True,
        help="Text file with the IDs of the sequences to extract.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="Output MSA file in FASTA format.",
    )

    args = parser.parse_args()

    extract_sequences(args.msa_file, args.id_file, args.output_file)


if __name__ == "__main__":
    main()
