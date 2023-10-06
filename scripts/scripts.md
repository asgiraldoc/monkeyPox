The `downloadSeqsFromNCBI.py` script downloads nucleotide sequences from NCBI's database using specified sequence IDs, and saves each sequence as a separate FASTA file.

How to run it: 


```
python downloadSeqsFromNCBI.py --email "tul54064@temple.edu" \
--api_key "44816979e53e34d64b1952d71518db87ab08" \
--input_file "../data/data_ncbi/ids.txt"
```

The `downloadGenBankFromNCBI.py` script downloads the metadata from NCBI's database using specified sequence IDs, and saves each record as a separate GenBank file.

How to run it: 


```
python downloadGenBankFromNCBI.py --email "tul54064@temple.edu" \
--api_key "44816979e53e34d64b1952d71518db87ab08" \
--input_file "../data/data_ncbi/ids.txt"
```

The `downloadGFFFromNCBI.py` script downloads the annotation genome information from NCBI's database using specified sequence IDs, and saves each record as a separate GFF file.

How to run it: 


```
python downloadGFFFromNCBI.py --email "tul54064@temple.edu" \
--api_key "44816979e53e34d64b1952d71518db87ab08" \
--input_file "../data/data_ncbi/ids.txt"
```


The `fasta2nex.py` script converts MSA fasta file to an annotated NEXUS file using the refSeq gff file.

How to run it: 
```
python fasta2nex.py --input_fasta ../data/data_ncbi/MSA.fasta \
--input_gff ../data/refSeq/refSeq.gff --output_nexus output_file.nex\
```

The `countNuclMSA.py` script counts the number of nucleotides (including IUPAC nucleotides and gaps) in each position of the MSA fasta file for each sequence and generate a TSV file with the result.

How to run it: 
```
python countNuclMSA.py --msa_file ../data/data_ncbi/MSA.fasta --output_file output_file.tsv
```

The `countNuclNexus.py` script counts the number of nucleotides (including IUPAC nucleotides and gaps) in each position of the nexus file for each gene block and generate a TSV file with the result.

How to run it: 
```
python countNuclNexus.py --input ../data/data_ncbi/MSA.nex --output output_file.tsv
```