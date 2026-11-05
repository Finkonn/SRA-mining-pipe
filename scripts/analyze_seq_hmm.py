import os, glob, argparse
from typing import List, Optional, Set
import numpy as np
import pandas as pd
import pyhmmer
from pyhmmer.easel import SequenceFile, TextSequence
from Bio.Seq import Seq

DEFAULT_THREADS = 1
DEFAULT_BATCH_SIZE = 10000
DEFAULT_MIN_AA_LEN = 30


def parse_arguments():
    p = argparse.ArgumentParser()
    p.add_argument('-f', '--file', required=True)
    p.add_argument('-o', '--output', required=True)
    p.add_argument('-m', '--hmm', required=True)
    p.add_argument('-t', '--threads', type=int, default=DEFAULT_THREADS)
    p.add_argument('-b', '--batch', type=int, default=DEFAULT_BATCH_SIZE)
    p.add_argument('-r', '--min_aa_len', type=int, default=DEFAULT_MIN_AA_LEN,
                   help='Minimum ORF length in amino acids (filters short frames)')
    p.add_argument('-l', '--positive_list', type=str)
    return p.parse_args()


def _to_str(x):
    """Decode pyhmmer bytes to str so pandas writes clean CSV (no b'...')."""
    return x.decode() if isinstance(x, (bytes, bytearray)) else x


def load_hmms(path, positive_list=None):
    hmms = []
    files = [path] if os.path.isfile(path) else glob.glob(os.path.join(path, '*'))
    for f in files:
        try:
            with pyhmmer.plan7.HMMFile(f) as hf:
                for hmm in hf:
                    hmms.append(hmm)
        except Exception as e:
            print(f"skip {f}: {e}")
    if positive_list:
        hmms = [h for h in hmms if _to_str(h.name) in positive_list]
    if not hmms:
        raise ValueError("No HMMs loaded")
    return hmms


def six_frame_orfs(seq: str, min_aa_len: int):
    """Yield (protein, start_nt, end_nt, strand, frame) for all 6 frames.
    start/end are 0-based half-open coords on the original contig."""
    L = len(seq)
    rc = str(Seq(seq).reverse_complement())
    for strand, s in (('+', seq), ('-', rc)):
        for frame in range(3):
            prot = str(Seq(s[frame:]).translate(to_stop=True))
            if len(prot) < min_aa_len:
                continue
            if strand == '+':
                start_nt = frame
                end_nt   = frame + len(prot) * 3
            else:
                end_nt   = L - frame
                start_nt = L - (frame + len(prot) * 3)
            yield prot, start_nt, end_nt, strand, frame


def analyze_file(seq_file, hmms, batch_size, threads, min_aa_len):
    alphabet = hmms[0].alphabet
    rows = []

    def flush(seqs):
        if not seqs:
            return
        for hit in pyhmmer.hmmer.hmmsearch(hmms, seqs, cpus=threads):
            q = _to_str(hit.query_name)
            for domain in hit:
                d_name = _to_str(domain.name)
                for d in domain.domains:
                    rows.append((q, d_name, d.score,
                                 d.alignment.target_from, d.alignment.target_to))

    seqs = []
    with SequenceFile(seq_file, digital=True, alphabet=alphabet) as sf:
        for dna in sf:
            name = dna.name.decode()
            seq = dna.textize().sequence
            for prot, start, end, strand, frame in six_frame_orfs(seq, min_aa_len):
                meta = f"{name}|{strand}|frame{frame}|{start}-{end}"
                seqs.append(TextSequence(sequence=prot, name=meta.encode()).digitize(alphabet))
                if len(seqs) >= batch_size:
                    flush(seqs); seqs = []
    flush(seqs)

    if not rows:
        return pd.DataFrame(columns=['Query', 'Name', 'Score', 'From', 'To'])
    df = pd.DataFrame(rows, columns=['Query', 'Name', 'Score', 'From', 'To'])
    df['Score'] = df['Score'].astype('float32')
    df['From']  = df['From'].astype('int32')
    df['To']    = df['To'].astype('int32')
    return df


def main():
    args = parse_arguments()
    positive = set(pd.read_csv(args.positive_list, header=None)[0]) if args.positive_list else None
    hmms = load_hmms(args.hmm, positive)

    files = glob.glob(args.file) if os.path.isdir(args.file) or '*' in args.file else [args.file]
    for f in files:
        name = os.path.basename(f).split('.')[0]
        print(f"Analyzing {name}")
        df = analyze_file(f, hmms, args.batch, args.threads, args.min_aa_len)
        out = args.output if args.output.endswith(('.csv', '.tsv')) else os.path.join(args.output, f"{name}.csv")
        os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
        df.to_csv(out, index=False)


if __name__ == '__main__':
    main()