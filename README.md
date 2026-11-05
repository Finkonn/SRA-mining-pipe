# SRA-mining-pipe

Snakemake pipeline for viral sequence discovery from SRA metagenomes. Performs read classification, host/contaminant removal, assembly, HMM-based profile search, and DIAMOND annotation.

## Requirements

- Snakemake
- Conda
- Kraken2 database
- DIAMOND database

## Usage

```bash
python launch_pipe.py -p {paired|single} -s <suffix> -i <input_dir> -o <output_dir> \
    -f <hmm_dir> -n <hmm_info.csv> -k <kraken_db> -d <diamond_db> \
    -m <metadata.tsv> -t <threads>
```