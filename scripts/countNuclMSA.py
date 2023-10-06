import argparse
from Bio import AlignIO
from Bio.Data import IUPACData
import csv


def count_nucleotides_and_gaps(msa_file, output_file):
    alignment = AlignIO.read(msa_file, "fasta")

    # Include all nucleotides according to the IUPAC definition and the gap "-"
    nucleotide_codes = set(IUPACData.ambiguous_dna_values.keys())
    nucleotide_codes.add("-")

    with open(output_file, "wt") as out:
        writer = csv.writer(out, delimiter="\t")

        # Write header
        header = ["sequence_id"] + list(nucleotide_codes)
        writer.writerow(header)

        for record in alignment:
            counts = {nucleotide: 0 for nucleotide in nucleotide_codes}

            # Count occurrences of nucleotides and gaps
            for nucleotide in record.seq:
                if nucleotide in counts:
                    counts[nucleotide] += 1

            # Write results
            row = [record.id] + [counts[nucleotide] for nucleotide in nucleotide_codes]
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Count nucleotides and gaps from an MSA file."
    )
    parser.add_argument(
        "--msa_file", required=True, help="Path to the MSA file in FASTA format."
    )
    parser.add_argument(
        "--output_file", required=True, help="Path to the output TSV file."
    )

    args = parser.parse_args()

    count_nucleotides_and_gaps(args.msa_file, args.output_file)


if __name__ == "__main__":
    main()
