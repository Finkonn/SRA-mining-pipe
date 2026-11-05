import pandas as pd
import numpy as np
import argparse
import os
import sys
from Bio import SeqIO
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def parse_arguments():
    p = argparse.ArgumentParser(description="Build HMM report from analyze_seq_hmm.py output.")
    p.add_argument('-i',  '--input',          required=True)
    p.add_argument('-i2', '--input2')
    p.add_argument('-c',  '--contigs_file',   required=True)
    p.add_argument('-o',  '--output',         required=True)
    p.add_argument('-m',  '--hmm_info',       required=True)
    p.add_argument('-t',  '--hmm_out')
    p.add_argument('-n',  '--output_contigs', required=True)
    p.add_argument('-r',  '--score_ratio',    type=float, default=0.5)
    return p.parse_args()


def check_file_exists(path, msg):
    if not os.path.isfile(path):
        logging.error(msg)
        sys.exit(1)


def load_csv(path):
    """Read CSV without consuming the first column as the index."""
    return pd.read_csv(path)


def _strip_bytes_repr(s: pd.Series) -> pd.Series:
    """Remove b'...' wrapper if present (defense-in-depth; scanner should already decode)."""
    return s.astype(str).str.replace(r"^b'|'$", "", regex=True)


def parse_name_column(names: pd.Series) -> pd.DataFrame:
    """Parse 'k141_207|-|frame0|357-582' into columns."""
    cleaned = _strip_bytes_repr(names)
    parts = cleaned.str.split('|', expand=True)
    if parts.shape[1] != 4:
        raise ValueError(f"Expected 4 fields in Name, got {parts.shape[1]}. "
                         f"First value: {cleaned.iloc[0]!r}")
    coords = parts[3].str.split('-', expand=True)
    return pd.DataFrame({
        'Name':      parts[0],
        'Strand':    parts[1],
        'Frame':     parts[2].str.replace('frame', '', regex=False).astype('int32'),
        'Orf_start': coords[0].astype('int32'),
        'Orf_end':   coords[1].astype('int32'),
    })


def map_coordinates(df: pd.DataFrame) -> pd.DataFrame:
    """Map protein-local (From, To) back to contig nucleotide coords."""
    fwd = df['Strand'] == '+'
    aa_from   = df['From'].astype('int32')
    aa_to     = df['To'].astype('int32')
    orf_start = df['Orf_start'].astype('int32')
    orf_end   = df['Orf_end'].astype('int32')
    df['From_nt'] = np.where(fwd, orf_start + aa_from * 3, orf_end - aa_to   * 3)
    df['To_nt']   = np.where(fwd, orf_start + aa_to   * 3, orf_end - aa_from * 3)
    return df.drop(columns=['From', 'To', 'Orf_start', 'Orf_end'])


def create_report(hits: pd.DataFrame, hmm_info: pd.DataFrame) -> pd.DataFrame:
    hits = hits.copy()
    hits['Query'] = _strip_bytes_repr(hits['Query'])

    report = hits.merge(hmm_info, left_on='Query', right_on='Model_ID', how='left')
    if 'Model_ID' in report.columns:
        report = report.drop(columns=['Model_ID'])

    if 'Threshold' not in report.columns:
        raise ValueError("HMM info CSV must contain a 'Threshold' column")
    if report['Threshold'].isna().all():
        raise ValueError("No hits matched any Model_ID — check Query vs Model_ID formatting")

    report['Score_ratio'] = report['Score'] / report['Threshold']

    parsed = parse_name_column(report['Name'])
    report = report.drop(columns=['Name']).join(parsed)
    report = map_coordinates(report)
    return report.sort_values('Score_ratio', ascending=False)


def main():
    args = parse_arguments()
    check_file_exists(args.input,        "Input hits CSV does not exist")
    check_file_exists(args.contigs_file, "Contigs FASTA does not exist")
    check_file_exists(args.hmm_info,     "HMM info CSV does not exist")

    hits     = load_csv(args.input)
    hmm_info = load_csv(args.hmm_info)

    if hits.empty:
        logging.warning("The hits file is empty (no HMM matches found). Creating empty outputs.")
        pd.DataFrame(columns=['Query','Name','Score','Score_ratio']).to_csv(args.output, index=False)
        if args.hmm_out:
            open(args.hmm_out, 'w').close()
        open(args.output_contigs, 'w').close()
        return
    report   = create_report(hits, hmm_info)

    if args.input2:
        check_file_exists(args.input2, "Second input hits CSV does not exist")
        hits2  = load_csv(args.input2)
        report = pd.concat([report, create_report(hits2, hmm_info)], ignore_index=True)

    report = report.query('Score_ratio > @args.score_ratio')

    if args.hmm_out:
        report[['Query']].drop_duplicates().to_csv(args.hmm_out, header=False, index=False)

    report.reset_index(drop=True).to_csv(args.output, index=False)

    wanted = set(report['Name'].unique())
    kept = []
    with open(args.contigs_file) as handle:
        for rec in SeqIO.parse(handle, 'fasta'):
            if rec.id in wanted:
                kept.append(rec)
    SeqIO.write(kept, args.output_contigs, 'fasta')
    logging.info(f"Wrote {len(kept)} contigs with hits; {len(report)} hit rows in report.")


if __name__ == '__main__':
    main()