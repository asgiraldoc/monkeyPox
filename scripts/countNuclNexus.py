import re
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(
    description="Count IUPAC nucleotides and gaps in a NEXUS file."
)
parser.add_argument(
    "--input", required=True, help="Path to the NEXUS file.", default="alignment.nexus"
)
parser.add_argument(
    "--output",
    required=True,
    help="Path to the output TSV file.",
    default="iupac_counts.tsv",
)
args = parser.parse_args()


# Function to count IUPAC nucleotides and gaps
def count_iupac_and_gaps(sequence):
    iupac_chars = set("YDXRWHNBKMSV")
    return sum(1 for char in sequence if char in iupac_chars or char == "-")


# List to store extracted data
data = {}

# Read NEXUS file and extract relevant data
with open(args.input, "r") as nexus_file:
    content = nexus_file.read()
    blocks = re.findall(r"\[Block: (.*?)\]\n.*?matrix\n(.*?);\n", content, re.DOTALL)

    # Filter blocks that do not correspond to intergenic regions
    gene_blocks = [
        (gene_name, sequences)
        for gene_name, sequences in blocks
        if not gene_name.startswith("intergenic_region_")
    ]

    for gene_name, sequences in gene_blocks:
        sequences = sequences.strip().split("\n")
        for seq in sequences:
            parts = seq.split()
            seq_name = parts[0]
            seq_data = "".join(parts[1:])
            count = count_iupac_and_gaps(seq_data)

            if seq_name not in data:
                data[seq_name] = {}
            data[seq_name][gene_name] = count

# Write results to TSV file
with open(args.output, "w") as tsv_file:
    # Header
    header = ["Sequence_ID"] + [gene for gene, _ in gene_blocks]
    tsv_file.write("\t".join(header) + "\n")

    # Length of each gene
    lengths = ["Length"] + [
        str(len(gene_blocks[i][1].split()[1])) for i in range(len(gene_blocks))
    ]
    tsv_file.write("\t".join(lengths) + "\n")

    # Data for each sample
    for seq_name, gene_data in data.items():
        line_data = [seq_name] + [
            str(gene_data.get(gene, 0)) for gene, _ in gene_blocks
        ]
        tsv_file.write("\t".join(line_data) + "\n")
