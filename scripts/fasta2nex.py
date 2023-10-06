from Bio import SeqIO
from collections import defaultdict
import re
import argparse


def parse_fasta(fasta_file):
    """Parse FASTA file."""
    return SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))


def parse_gff(gff_file):
    """Parse GFF file and return gene annotations."""
    regions = []
    with open(gff_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                if parts[2] == "gene":
                    start, end, attributes = int(parts[3]), int(parts[4]), parts[8]
                    match = re.search(r'gene_id "([^"]+)"', attributes)
                    gene_id = match.group(1) if match else "unknown_id"
                    regions.append((start, end, gene_id))
    return regions


def determine_intergenic_regions(regions):
    """Determine intergenic regions based on gene annotations."""
    previous_end = 0
    intergenic_count = 1
    final_regions = []
    for start, end, gene_id in regions:
        if start - previous_end > 1:
            final_regions.append(
                (previous_end + 1, start - 1, f"intergenic_region_{intergenic_count}")
            )
            intergenic_count += 1
        final_regions.append((start, end, gene_id))
        previous_end = end
    return final_regions


def write_nexus_block(output_file, name, sequences, start, end):
    """Write blocks to the NEXUS file."""
    output_file.write("#NEXUS\n\n")
    output_file.write(f"[Block: {name}]\n")
    output_file.write(f"begin taxa;\n")
    output_file.write(f"  dimensions ntax={len(sequences)};\n")
    output_file.write(f"  taxlabels {' '.join(sequences.keys())};\n")
    output_file.write(f"end;\n\n")
    output_file.write(f"begin characters;\n")
    output_file.write(f"  dimensions nchar={end - start + 1};\n")
    output_file.write(f"  format datatype=dna missing=? gap=- matchchar=.;\n")
    output_file.write(f"  matrix\n")
    for seq_name, seq in sequences.items():
        output_file.write(f"{seq_name} {str(seq.seq)[start:end+1]}\n")
    output_file.write(f";\n")
    output_file.write(f"end;\n\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a NEXUS file based on gene annotations."
    )
    parser.add_argument(
        "--input_fasta", required=True, help="Path to the input FASTA file."
    )
    parser.add_argument(
        "--input_gff", required=True, help="Path to the input GFF file."
    )
    parser.add_argument(
        "--output_nexus", required=True, help="Path to the output NEXUS file."
    )
    args = parser.parse_args()

    alignment = parse_fasta(args.input_fasta)
    regions = parse_gff(args.input_gff)
    regions.sort(key=lambda x: x[0])
    final_regions = determine_intergenic_regions(regions)

    with open(args.output_nexus, "w") as nexus_file:
        for start, end, name in final_regions:
            write_nexus_block(nexus_file, name, alignment, start - 1, end - 1)


if __name__ == "__main__":
    main()
