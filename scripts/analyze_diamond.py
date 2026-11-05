#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd

DIAMOND_COLS = [
    "qseqid", "sseqid", "pident", "align_len", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qcovhsp", "scovhsp", "stitle",
]

OUT_COLS = [
    "qseqid", "sseqid", "pident", "align_len", "evalue", "bitscore",
    "qlen", "stitle", "top3_blast", "top5_profiles",
]

ANNOT_COLS = [
    "Sample", "Country", "Year", "Host", "Contig_length_nt", "Profile",
    "Percent_identity", "Alignment_length_aa", "Mismatches_aa",
    "Contig_start", "Contig_end",
    "Reference_start", "Reference_end", "Reference_length_aa",
    "Contig_coverage_pct", "Reference_coverage_pct",
    "E_value", "Bit_score",
    "Blast_hit", "Info", "Top5_profiles",
]


def parse_args():
    p = argparse.ArgumentParser(description="Post-process DIAMOND results per sample.")
    p.add_argument("-r",  "--raw",        required=True, help="Raw DIAMOND TSV")
    p.add_argument("-c",  "--report",     required=True, help="HMM report CSV for this sample")
    p.add_argument("-m",  "--metadata",   required=True, help="Metadata table (quoted TSV)")
    p.add_argument("-s",  "--sample",     required=True, help="Sample name (for qseqid prefix)")
    p.add_argument("-ob", "--out-best",   required=True, help="Output best_blast_hits TSV")
    p.add_argument("-oa", "--out-annot",  required=True, help="Output annotated_blast TSV")
    return p.parse_args()


def _clean(s):
    """Sanitize a field for embedding in a FASTA-style qseqid."""
    return (str(s)
            .replace(";", "_").replace(",", "_").replace("|", "_")
            .replace("/", "_").replace(" ", "_").strip())


def load_profiles_and_best(path):
    """From the HMM report CSV return (profiles, best) keyed by contig name.
    profiles[contig] -> 'top5 profiles' string
    best[contig]     -> {taxon, positive_terms, score_ratio} for the top-1 hit."""
    df = pd.read_csv(path)
    df = df.loc[:, ~df.columns.str.startswith("Unnamed:")]
    if df.empty:
        return {}, {}

    df = df.sort_values(["Name", "Score_ratio"], ascending=[True, False])
    profiles, best = {}, {}
    for contig, group in df.groupby("Name", sort=False):
        top = []
        for _, row in group.head(5).iterrows():
            taxon = str(row.get("Taxon", "Unknown"))
            term  = str(row.get("Positive terms", "") or "")
            score = float(row.get("Score_ratio", 0.0))
            if term and term.lower() != "nan":
                top.append(f"{taxon} *{term}* ({score:.2f})")
            else:
                top.append(f"{taxon} ({score:.2f})")
        profiles[contig] = "; ".join(top)

        r0 = group.iloc[0]
        best[contig] = {
            "taxon":          str(r0.get("Taxon", "Unknown")),
            "positive_terms": str(r0.get("Positive terms", "") or ""),
            "score_ratio":    float(r0.get("Score_ratio", 0.0)),
        }
    return profiles, best


def format_top3(rows):
    out = []
    for _, r in rows.head(3).iterrows():
        stitle = str(r.get("stitle", "") or "").strip()
        sid    = str(r.get("sseqid", "") or "").strip()
        bs     = float(r.get("bitscore", 0.0))
        out.append(f"{stitle or sid} ({bs:.1f})")
    return "; ".join(out)


def load_metadata(path):
    """Read a TSV whose header and data rows are wrapped in double quotes.
    Expected columns: Run, Country, Host, Year, BioProject.
    Returns {Run: {Country, Host, Year, BioProject}}."""
    rows = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith('"') and line.endswith('"'):
                line = line[1:-1]
            rows.append(line.split("\t"))

    if not rows:
        return {}
    header = [h.strip() for h in rows[0]]
    out = {}
    for r in rows[1:]:
        if len(r) < len(header):
            continue
        rec = dict(zip(header, [x.strip() for x in r]))
        run = rec.get("Run", "").strip()
        if not run:
            continue
        out[run] = {
            "Country":    rec.get("Country", ""),
            "Host":       rec.get("Host", ""),
            "Year":       rec.get("Year", ""),
            "BioProject": rec.get("BioProject", ""),
        }
    return out


def parse_qseqid(q):
    """Parse 'sample|contig/taxon/len/positive_terms'."""
    if not q:
        return {"sample": "", "contig": "", "taxon": "",
                "contig_len": "", "positive_terms": ""}
    if "|" in q:
        sample, rest = q.split("|", 1)
    else:
        sample, rest = "", q
    fields = rest.split("/")
    return {
        "sample":         sample,
        "contig":         fields[0] if len(fields) > 0 else "",
        "taxon":          fields[1] if len(fields) > 1 else "",
        "contig_len":     fields[2] if len(fields) > 2 else "",
        "positive_terms": "/".join(fields[3:]) if len(fields) > 3 else "",
    }


def build_qseqid(sample, contig, best_info, contig_len):
    """SRR28250302|k141_14479/Hepatovirus/9793/polyprotein"""
    taxon = _clean(best_info.get("taxon", "Unknown"))
    term  = _clean(best_info.get("positive_terms", "") or "")
    parts = [sample, contig, taxon, str(contig_len)]
    if term and term.lower() != "nan":
        parts.append(term)
    return "|".join([parts[0], "/".join(parts[1:])])


def main():
    args = parse_args()
    for p in (args.raw, args.report, args.metadata):
        if not os.path.isfile(p):
            print(f"ERROR: missing input: {p}", file=sys.stderr)
            sys.exit(1)

    profiles, best_info = load_profiles_and_best(args.report)
    meta_dict = load_metadata(args.metadata)

    if os.path.getsize(args.raw) == 0:
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_best, sep="\t", index=False)
        pd.DataFrame(columns=ANNOT_COLS).to_csv(args.out_annot, sep="\t", index=False)
        print(f"Empty DIAMOND output for {args.sample}; wrote empty tables", file=sys.stderr)
        return

    df = pd.read_csv(args.raw, sep="\t", header=None, names=DIAMOND_COLS)
    df_sorted = df.sort_values("bitscore", ascending=False)

    top3 = df_sorted.groupby("qseqid", sort=False).head(3)
    top3_map = (
        top3.groupby("qseqid", sort=False)
            .apply(format_top3, include_groups=False)
            .to_dict()
    )

    best = df_sorted.drop_duplicates("qseqid").copy()
    best["top3_blast"]    = best["qseqid"].map(top3_map).fillna("")
    best["top5_profiles"] = best["qseqid"].map(profiles).fillna("")

    qlen_by_qseqid = dict(zip(best["qseqid"], best["qlen"].astype(int)))

    def make_qseqid(q):
        bi = best_info.get(q, {})
        return build_qseqid(args.sample, q, bi, qlen_by_qseqid[q])

    best["qseqid"] = best["qseqid"].map(make_qseqid)
    best[OUT_COLS].to_csv(args.out_best, sep="\t", index=False)
    print(f"Wrote {len(best)} best hits to {args.out_best}")

    meta_info = meta_dict.get(args.sample, {})
    rows = []
    for _, r in best.iterrows():
        parsed = parse_qseqid(r["qseqid"])
        contig     = parsed["contig"]
        taxon      = parsed["taxon"]
        term       = parsed["positive_terms"]
        contig_len = parsed["contig_len"]

        profile = f"{taxon}|{term}" if term else taxon

        blast_hit = str(r["stitle"]).replace(";", " ") if pd.notna(r["stitle"]) else ""

        info_parts = [f"Contig_ID: {contig}"]
        if meta_info.get("BioProject"):
            info_parts.append(f"BioProject: {meta_info['BioProject']}")
        info = " | ".join(info_parts)

        rows.append({
            "Sample":                args.sample,
            "Country":               meta_info.get("Country", ""),
            "Year":                  meta_info.get("Year", ""),
            "Host":                  meta_info.get("Host", ""),
            "Contig_length_nt":      contig_len,
            "Profile":               profile,
            "Percent_identity":      r["pident"],
            "Alignment_length_aa":   r["align_len"],
            "Mismatches_aa":         r["mismatch"],
            "Contig_start":          r["qstart"],
            "Contig_end":            r["qend"],
            "Reference_start":       r["sstart"],
            "Reference_end":         r["send"],
            "Reference_length_aa":   r["slen"],
            "Contig_coverage_pct":   r["qcovhsp"],
            "Reference_coverage_pct": r["scovhsp"],
            "E_value":               r["evalue"],
            "Bit_score":             r["bitscore"],
            "Blast_hit":             blast_hit,
            "Info":                  info,
            "Top5_profiles":         r["top5_profiles"],
        })

    annot = pd.DataFrame(rows, columns=ANNOT_COLS)
    if not annot.empty:
        annot["E_value"] = annot["E_value"].apply(
            lambda x: f"{x:.2e}" if pd.notna(x) else "N/A"
        )
    annot.to_csv(args.out_annot, sep="\t", index=False)
    print(f"Wrote {len(annot)} annotated rows to {args.out_annot}")


if __name__ == "__main__":
    main()