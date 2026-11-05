import subprocess
import argparse
import os
import sys

description = """pipe wrapper"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-p', '--paired', type=str, required=True, choices=['paired', 'single'],
                    help="Mode: 'paired' or 'single'")
parser.add_argument('-s', '--suffix', type=str, default='_#.fastq.gz')
parser.add_argument('-i', '--input_directory', type=str, default="test_input/")
parser.add_argument('-o', '--output_directory', type=str, default="test_output/")
parser.add_argument('-f', '--hmm_folder', type=str, default='test_profiles/')
parser.add_argument('-n', '--hmm_info', type=str, default='hmm_info.csv')
parser.add_argument('-k', '--kraken_database', type=str, required=True)
parser.add_argument('-d', '--diamond_database', type=str, required=True)
parser.add_argument('-m', '--metadata', type=str, required=True)
parser.add_argument('-t', '--threads', type=int, default=32)

args = parser.parse_args()
conf = vars(args)

paired = conf.pop('paired')
threads = conf.pop('threads')
current_dir = os.getcwd()

if paired == "paired":
    snakefile = "scripts/pipe_paired.py"
    if '#' not in conf['suffix']:
        print(f"ERROR: paired mode requires '#' in --suffix, got {conf['suffix']!r}",
              file=sys.stderr)
        sys.exit(1)
elif paired == "single":
    snakefile = "scripts/pipe_single.py"
    conf['suffix'] = conf['suffix'].replace('_#', '').replace('#', '')
    if '#' in conf['suffix']:
        print(f"ERROR: could not strip '#' from suffix, got {conf['suffix']!r}",
              file=sys.stderr)
        sys.exit(1)
else:
    print(f"ERROR: unknown mode: {paired}", file=sys.stderr)
    sys.exit(1)

conf['workdir'] = current_dir

config_args = [f'{k}={v}' for k, v in conf.items()]

cmd = [
    'snakemake',
    '-s', os.path.abspath(snakefile),
    '--use-conda', '-k', '-p', '--rerun-incomplete',
    '-j', str(threads),
    '--config', *config_args,
]

print('Running:', ' '.join(cmd), file=sys.stderr)
subprocess.run(cmd, check=True, text=True)