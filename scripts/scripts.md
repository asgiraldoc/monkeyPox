The `downloadSeqsFromNCBI.py` script downloads nucleotide sequences from NCBI's database using specified sequence IDs, and saves each sequence as a separate FASTA file.

How to run it: 


```
python downloadSeqsFromNCBI.py --email "tul54064@temple.edu" \
--api_key "44816979e53e34d64b1952d71518db87ab08" \
--input_file "../data/data_ncbi/ids.txt"
```

The `downloadMetadataFromNCBI.py` script downloads the metadata from NCBI's database using specified sequence IDs, and saves each record as a separate JSON file.

How to run it: 


```
python downloadMetadataFromNCBI.py --email "tul54064@temple.edu" \
--api_key "44816979e53e34d64b1952d71518db87ab08" \
--input_file "../data/data_ncbi/ids.txt"
```