#!/usr/bin/env python3
"""

This script generates fastq.gz mixtures by combining single-source fastq.gz profiles based on user-defined mixture parameters.

Usage:
  MixtureSimulator1.0.py -i <input_dir> -o <output_dir> -b <base_line> -r <repetitions> -n <number_ind> -ra <ratios> [--paired]
  MixtureSimulator1.0.py --input_dir <input_dir> --output_dir <output_dir> --base_line <base_line> --repetitions <repetitions> --number_ind <number_ind> --ratios <ratios> [--paired]
  MixtureSimulator1.0.py (-h | --help)
  MixtureSimulator1.0.py --version

Options:
  -i --input_dir <input_dir>     Path to the folder containing single-source fastq/fastq.gz files.
  -o --output_dir <output_dir>   Path to the output folder for writing mixtures and logs.
  -b --base_line <base_line> Number of reads in simulated mixtures.
  -re --repetitions <repetitions>  
                             Number of repetitions per mixture.
  -n --number_ind <number_ind>  
                             Number of individuals per mixture (e.g., 2, 3, 4).
  -ra --ratios <ratios>      Comma-separated list of mixture ratios (e.g., "1_1,2_1,3_2_1").
  --paired                   Flag activating script for paired-end sequencing data.
  -h --help                  Show this help message and exit.
  --version                  Show script version.

Example Usage:
  MixtureSimulator1.0.py -i ./data -o ./mixtures  -b 220000 -re 3 -n 3 -ra 1_1_1,2_1_1 --paired
  MixtureSimulator1.0.py --input_dir ./samples --output_dir ./results --base_line 100000 --repetitions 5 --number_ind 2 --ratios 1_1,3_1

Notes:
  - The script will automatically verify that the number of available fastq/fastq.gz files is sufficient for the requested mixtures.
  - Each sample in a mixture is **shuffled independently** before subsampling.
  - If a file has **fewer reads than required**, the script will warn the user.
  - The script will generate a log file containing the information about the dependencies used, timestamps, warnings, and errors if any.

"""

import argparse
import platform
import Bio
import sys
import subprocess
import os
import gzip
import shutil
import random
from datetime import datetime
from pathlib import Path
from itertools import combinations
from Bio import SeqIO


log_lines = []
log_file_path = None

def init_logger(output_dir, args):
    global log_file_path
    timestamp = datetime.now().strftime('%Y%b%d_%H%M%S')
    log_file_path = output_dir / f"log_{timestamp}.log"
    
    # Header
    log_event(f"Analysis date: {timestamp}", timestamp=True)

    # Parameters
    log_event("\nCommandline parameters:", timestamp=False)
    log_event(f"Path to folder: {Path(args.input_dir).resolve()}", timestamp=False)
    log_event(f"Baseline read count: {args.base_line}", timestamp=False)
    log_event(f"Repetitions: {args.repetitions}", timestamp=False)
    log_event(f"Number of individuals per mixture: {args.number_ind}", timestamp=False)
    log_event(f"Mixture ratios: {args.ratios}", timestamp=False)
    log_event(f"Paired-end mode: {'Enabled' if args.paired else 'Disabled'}", timestamp=False)

    # Dependencies
    log_event("\nDependency versions:", timestamp=False)
    log_event(f"Python version: {platform.python_version()}", timestamp=False)
    log_event(f"Bio version: {Bio.__version__}\n", timestamp=False)



def log_event(message, timestamp=True):
    """Append a log entry with optional timestamp."""
    entry = f"> timestamp {datetime.now().strftime('%H:%M:%S')}\n{message}" if timestamp else message
    log_lines.append(entry)
    if log_file_path:
        with open(log_file_path, 'a') as f:
            f.write(entry + '\n')


def convert_fastq_to_gz(fastq_file):
    """Convert .fastq to .fastq.gz"""
    gz_file = fastq_file.with_suffix('.fastq.gz')
    with open(fastq_file, 'rb') as f_in, gzip.open(gz_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(fastq_file)  # Remove the original .fastq file
    return gz_file


def parse_ratios(ratio_str):
    """Convert ratio string (e.g., '1_1,2_1' -> [[1,1], [2,1]])"""
    return [list(map(int, r.split("_"))) for r in ratio_str.split(",")]

def compute_mixture_counts(ratios, base_line):
    """Compute read counts per individual based on the provided ratios."""
    mixture_counts = []
    for ratio in ratios:
        total_parts = sum(ratio)  # Sum of all ratio parts
        counts = [int((part / total_parts) * base_line) for part in ratio]  # Scale to base_line
        mixture_counts.append(counts)
    return mixture_counts

def index_fastq_reads(file_path):
    """Read the fastq files and save index and number of reads"""
    with gzip.open(file_path, "rt") as f:
        reads = list(SeqIO.parse(f, "fastq"))
        return reads, len(reads)

def sample_reads(reads, sample_size):
    """Randomly shuffles reads and select based on sample_size"""
    random.shuffle(reads)
    return random.sample(reads, sample_size)

def main():
    parser = argparse.ArgumentParser(description="Generate FASTQ mixtures.")
    parser.add_argument("-i", "--input_dir", type=str, required=True, help="Path to folder with input fastq/fastq.gz files")
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="Path to folder where mixtures and log will be saved")
    parser.add_argument("-b", "--base_line", type=int, required=True, help="Baseline read count for normalization")
    parser.add_argument("-re", "--repetitions", type=int, required=True, help="Number of repetitions")
    parser.add_argument("-n", "--number_ind", type=int, required=True, help="Number of individuals per mixture")
    parser.add_argument("-ra", "--ratios", type=str, required=True, help="Comma-separated mixture ratios")
    parser.add_argument("--paired", action="store_true", help="Enable paired-end mode")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    if not input_dir.exists():
        print(f"Error: The specified path '{args.input_dir}' does not exist.")
        return
    elif not input_dir.is_dir():
        print(f"Error: The specified path '{args.input_dir}' is not a directory.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Convert .fastq files to .fastq.gz files if necessary
    fastq_files = list(input_dir.glob("*.fastq"))
    if fastq_files:
        print("Converting fastq files to fastq.gz...")
        for fastq_file in fastq_files:
            convert_fastq_to_gz(fastq_file)
        print("Converted .fastq files to .fastq.gz.")

    # Find all .fastq.gz files in the folder
    fastq_files = list(input_dir.glob("*.fastq.gz"))
    if not fastq_files:
        print("Error: No .fastq.gz files found.")
        return

    # Write in the log_file
    init_logger(output_dir, args)
    log_event("fastq.gz files found: " + str([f.name for f in fastq_files]))
    
    # Section for R2 files
    if args.paired:
        r1_files = sorted([f for f in fastq_files if "_R1" in f.name])
        r2_files = sorted([f for f in fastq_files if "_R2" in f.name])
        r1_to_r2 = {}

        # Check if any R1 or R2 files were found
        if not r1_files or not r2_files:
            print(
                "Error: There is no _R1 and/or _R2 file found. Possible reasons: wrong file naming or data is not paired-end.")
            return

        # Checking all R1 files have an R2
        for r1 in r1_files:
            base = r1.name.replace("_R1", "")
            r2_match = next((r2 for r2 in r2_files if r2.name.replace("_R2", "") == base), None)
            if not r2_match:
                print(f"Error: Could not find R2 pair for {r1.name}")
                return
            r1_to_r2[r1] = r2_match
        fastq_files = list(r1_to_r2.keys())  # Only use R1s for combinations

    if len(fastq_files) < args.number_ind:
        print(f"Error: You requested {args.number_ind} individuals per mixture, but only {len(fastq_files)} fastq.gz files were found.")
        return

    # Check if the number of ratios matches the expected number of individuals
    ratios = parse_ratios(args.ratios)
    for ratio in ratios:
        if len(ratio) != args.number_ind:
            print(f"Error: Ratio {ratio} has {len(ratio)} parts, expected {args.number_ind}.")
            return

    print("fastq.gz files found:", [f.name for f in fastq_files])

    # Indexing and counting R1 and R2 files
    print("Indexing FASTQ files...")
    read_counts = {}
    reads_index = {}

    # Check for empty files and match R1&R2
    for f in fastq_files:
        try:
            reads, count = index_fastq_reads(f)
            if count == 0:
                print(f"Warning: The file '{f.name}' is empty and will be skipped.")
                continue
            read_counts[f] = count
            reads_index[f] = reads
        except Exception as e:
            print(f"Error reading file {f.name}: {e}")
            continue

    if args.paired:
        for r1, r2 in r1_to_r2.items():
            try:
                reads, count = index_fastq_reads(r2)
                if count == 0:
                    print(f"Warning: The file '{r2.name}' is empty and will be skipped.")
                    continue
                read_counts[r2] = count
                reads_index[r2] = reads
                if read_counts[r1] != read_counts[r2]:
                    print(f"Warning: total read count mismatch for {r1.name} and {r2.name}")
            except Exception as e:
                print(f"Error reading paired file {r2.name}: {e}")
                continue


    # Calculating mixture proportion and read counts of the single source samples
    base_line_mixtures_list = compute_mixture_counts(ratios, args.base_line)
    unique_combi = list(combinations(fastq_files, args.number_ind))

    # Main loop
    for pair in unique_combi:
        print(f"Processing pair: {[p.name for p in pair]}")

        for idx, mix in enumerate(base_line_mixtures_list):
            mix_label = args.ratios.split(",")[idx]
            for rep in range(args.repetitions):

                # Output file names
                base_names = "_".join([p.stem for p in pair])
                if args.paired:
                    final_output_file_r1 = output_dir / f"{base_names}_{mix_label}_R{rep}_R1.fastq.gz"
                    final_output_file_r2 = output_dir / f"{base_names}_{mix_label}_R{rep}_R2.fastq.gz"
                else:
                    final_output_file_r1 = output_dir / f"{base_names}_{mix_label}_R{rep}.fastq.gz"

                merged_r1 = []
                merged_r2 = []

                log_event(f"Processing pair: {[p.name for p in pair]} with ratio {mix_label} (rep {rep})")

                for i in range(args.number_ind):
                    # Read R1, shuffle, sample
                    required_reads = mix[i]
                    available_reads = read_counts[pair[i]]
                    if available_reads < required_reads:
                        print(f"Error: Not enough reads in '{pair[i].name}': required {required_reads}, found {available_reads}")
                        sys.exit(1)

                    sampled_r1 = sample_reads(reads_index[pair[i]], required_reads)
                    merged_r1.extend(sampled_r1)


                    if args.paired:
                        # Extract IDs safely from R1 reads
                        ids_to_keep = {r.id.split()[0] for r in sampled_r1}
                        
                        r2_file = r1_to_r2.get(pair[i])

                        # Match R2 reads based on the sampled R1 IDs
                        matched_r2 = [r for r in reads_index[r2_file] if r.id.split()[0] in ids_to_keep]
                        if len(matched_r2) != len(sampled_r1):
                            print(f"Warning: Paired read count mismatch for {pair[i]}: R1={len(sampled_r1)}, R2={len(matched_r2)}")
                        merged_r2.extend(matched_r2)

                with gzip.open(final_output_file_r1, "wt") as fout:
                    SeqIO.write(merged_r1, fout, "fastq")

                if args.paired:
                    with gzip.open(final_output_file_r2, "wt") as fout2:
                        SeqIO.write(merged_r2, fout2, "fastq")
    
    log_event("> END " + datetime.now().strftime('%H:%M:%S'))
    print(f"Log written to: {log_file_path}")
    print("DONE")    

if __name__ == "__main__":
    main()
