import os

# Config
WORKDIR = config["workdir"]
OUTPUT_DIR = config["output_directory"]
INPUT_DIR = config["input_directory"]  
HMM_DIR = config["hmm_folder"]
HMM_INFO = config["hmm_info"]
KRAKEN_DB = config["kraken_database"] # Should point to the directory with hash.k2d file etc
DIAMOND_DB = config['diamond_database'] # Should point to .dmnd file
METADATA = config["metadata"]

MODE = config.get("mode", "paired").lower()
if MODE not in ("paired", "single"):
    raise ValueError(f"config['mode'] must be 'paired' or 'single', got {MODE!r}")

READ_SUFFIXES = [config["suffix"].replace('#', str(i)) for i in (1, 2)]

# INPUT_DIR/ERR6090673_1.fastq.gz
SAMPLES = glob_wildcards(os.path.join(INPUT_DIR, "{sample}" + READ_SUFFIXES[0])).sample
SAMPLES = [s for s in SAMPLES if not s.startswith(('report', 'general', 'summary'))]


rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "diamond/{sample}_annotated_blast.tsv"), sample=SAMPLES),


rule kraken2_classification:
    input:
        r1 = os.path.join(INPUT_DIR, "{sample}" + READ_SUFFIXES[0]),
        r2 = os.path.join(INPUT_DIR, "{sample}" + READ_SUFFIXES[1]),
    output:
        report = os.path.join(OUTPUT_DIR, "kraken_results/{sample}.report"),
        out    = os.path.join(OUTPUT_DIR, "kraken_results/{sample}.out"),
        uc_r1  = os.path.join(OUTPUT_DIR, "not_classified_fastq/{sample}_1.fastq.gz"),
        uc_r2  = os.path.join(OUTPUT_DIR, "not_classified_fastq/{sample}_2.fastq.gz"),
    params:
        db = KRAKEN_DB,
        uc_tmp = os.path.join(OUTPUT_DIR, "not_classified_fastq", "{sample}#.fastq.gz.tmp"),
    threads: 8
    priority: 1
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/kraken_results {OUTPUT_DIR}/not_classified_fastq
        kraken2 --threads {threads} --memory-mapping --confidence 0.7 --db {params.db} \
            {input.r1} {input.r2} --use-names --report {output.report} \
            --output {output.out} \
            --unclassified-out {params.uc_tmp} \
            --paired
        gzip -c {output.uc_r1}.tmp > {output.uc_r1}
        gzip -c {output.uc_r2}.tmp > {output.uc_r2}
        rm {output.uc_r1}.tmp {output.uc_r2}.tmp
        """


rule deduplicate_reads:
    input:
        r1 = os.path.join(OUTPUT_DIR, "not_classified_fastq/{sample}_1.fastq.gz"),
        r2 = os.path.join(OUTPUT_DIR, "not_classified_fastq/{sample}_2.fastq.gz"),
    output:
        r1 = os.path.join(OUTPUT_DIR, "deduplicated_fastq/{sample}_1.fastq.gz"),
        r2 = os.path.join(OUTPUT_DIR, "deduplicated_fastq/{sample}_2.fastq.gz"),
    threads: 8
    priority: 2
    conda:
        "envs/fastp.yaml"
    shell:
        """
        fastp -w {threads} -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --dedup \
            --detect_adapter_for_pe \
            -l 50 \
            -h {OUTPUT_DIR}/fastp_reports/{wildcards.sample}.html \
            -j {OUTPUT_DIR}/fastp_reports/{wildcards.sample}.json
        """


rule megahit_assembly:
    input:
        r1 = rules.deduplicate_reads.output.r1,
        r2 = rules.deduplicate_reads.output.r2,
    output:
        contigs = os.path.join(OUTPUT_DIR, "megahit/{sample}/final.contigs.fa"),
    params:
        memory = 0.8
    threads: 16
    priority: 3
    conda:
        "envs/megahit.yaml"
    shell:
        """
        rm -rf {OUTPUT_DIR}/megahit/{wildcards.sample}
        megahit -1 {input.r1} -2 {input.r2} \
                -o {OUTPUT_DIR}/megahit/{wildcards.sample} \
                -t {threads} --memory {params.memory}
        """


rule hmm_scan:
    input:
        contigs = rules.megahit_assembly.output.contigs,
    output:
        hits = os.path.join(OUTPUT_DIR, "hmm_results/{sample}.csv"),
    params:
        hmm_dir = HMM_DIR,
        min_aa_len = 50,
    threads: 16
    priority: 4
    conda:
        "envs/hmmer.yaml"
    shell:
        """
        python scripts/analyze_seq_hmm.py \
            -f {input.contigs} \
            -o {output.hits} \
            -m "{params.hmm_dir}" \
            -t {threads} \
            --batch 10000 \
            --min_aa_len {params.min_aa_len}
        """


rule generate_hmm_report:
    input:
        hits = rules.hmm_scan.output.hits,
        contigs = rules.megahit_assembly.output.contigs,
    output:
        report = os.path.join(OUTPUT_DIR, "hmm_reports/{sample}.csv"),
        matched_contigs = os.path.join(OUTPUT_DIR, "hmm_reports/{sample}_matched_contigs.fasta"),
        hmm_list = os.path.join(OUTPUT_DIR, "hmm_reports/{sample}_hmm_list.csv"),
    params:
        hmm_info = HMM_INFO,
        score_ratio = 0.5,
    priority: 5
    conda:
        "envs/hmmer.yaml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/hmm_reports
        python scripts/generate_hmm_report.py \
            -i "{input.hits}" \
            -c "{input.contigs}" \
            -o "{output.report}" \
            -m "{params.hmm_info}" \
            -t "{output.hmm_list}" \
            -n "{output.matched_contigs}" \
            -r {params.score_ratio}
        """


rule diamond_blast:
    input:
        contigs = rules.generate_hmm_report.output.matched_contigs,
    output:
        raw = os.path.join(OUTPUT_DIR, "diamond/{sample}_raw_diamond_results.tsv"),
    params:
        db = DIAMOND_DB,
        fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend "
               "sstart send evalue bitscore qlen slen qcovhsp scovhsp stitle",
    threads: 16
    priority: 6
    conda:
        "envs/diamond.yaml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/diamond
        diamond blastx \
            --ultra-sensitive \
            -d "{params.db}" \
            -q "{input.contigs}" \
            -o "{output.raw}" \
            -f {params.fmt} \
            -k 5 \
            -p {threads}
        """


rule analyze_diamond:
    input:
        raw      = rules.diamond_blast.output.raw,
        report   = rules.generate_hmm_report.output.report,
        metadata = METADATA,
    output:
        best  = os.path.join(OUTPUT_DIR, "diamond/{sample}_best_blast_hits.tsv"),
        annot = os.path.join(OUTPUT_DIR, "diamond/{sample}_annotated_blast.tsv"),
    priority: 5
    conda:
        "envs/hmmer.yaml"
    shell:
        """
        python scripts/analyze_diamond.py \
            -r "{input.raw}" \
            -c "{input.report}" \
            -m "{input.metadata}" \
            -s "{wildcards.sample}" \
            -ob "{output.best}" \
            -oa "{output.annot}"
        """