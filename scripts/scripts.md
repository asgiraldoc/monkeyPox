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