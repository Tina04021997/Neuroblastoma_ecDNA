#!/usr/bin/env python3
"""Reformat the PASS-filtered DELLY calls (02_filter_delly_vcf.sh, <sample>_sv.txt) into the 6-column SV table
ShatterSeek expects (chr, start, chr2, end, svtype, ct), written as sv_<sample>.txt.

Usage:
  python 03_format_sv_for_shatterseek.py --in <dir with <sample>_sv.txt> --out <dir for sv_<sample>.txt>
Exits non-zero if no input file is found or a record lacks END/SVTYPE/CT.
"""
import argparse
import glob
import os
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--in",  dest="in_dir",  required=True, help="folder with <sample>_sv.txt from 02_filter_delly_vcf.sh")
ap.add_argument("--out", dest="out_dir", required=True, help="output folder for sv_<sample>.txt (05_run_shatterseek.R input)")
a = ap.parse_args()

input_files = sorted(glob.glob(os.path.join(a.in_dir, "*_sv.txt")))
if not input_files:
    sys.exit(f"ERROR: no *_sv.txt files in {a.in_dir}")
os.makedirs(a.out_dir, exist_ok=True)
output_directory = a.out_dir

# Process each input file
for input_file_path in input_files:
    # Extract the input file name without the "_sv" suffix
    input_file_name = os.path.basename(input_file_path).replace("_sv", "")
    # Generate the output file path based on the modified input file name
    output_file_path = os.path.join(output_directory, f"sv_{input_file_name}")

    print(f"Processing {input_file_name}")

    # Open the input file for reading and the output file for writing
    with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
        # Write the header to the output file
        output_file.write("chr\tstart\tchr2\tend\tsvtype\tct\n")
        # Process each line in the input file
        for input_line in input_file:
            # Split the input line by tabs
            if not input_line.strip():
                continue
            fields = input_line.strip().split("\t")
            if len(fields) < 8:
                sys.exit(f"ERROR: {input_file_path}: fewer than 8 tab-separated columns (not a VCF record): {input_line.strip()[:80]}")
            chr = fields[0].replace("chr", "")
            start = fields[1]
            info_field = fields[7].split(";")
            try:
                # DELLY 1.x writes CHR2 only for inter-chromosomal records (BND/TRA); intra-chromosomal DEL/DUP/INV
                # carry END on the same chromosome, so chr2 defaults to chr (as in the executed 2026 run).
                chr2_list = [item.split("=")[1] for item in info_field if item.startswith("CHR2=")]
                chr2 = chr2_list[0].replace("chr", "") if chr2_list else chr
                end = [item.split("=")[1] for item in info_field if item.startswith("END=")][0]
                svtype = [item.split("=")[1] for item in info_field if item.startswith("SVTYPE=")][0]
                ct = [item.split("=")[1] for item in info_field if item.startswith("CT=")][0]
            except IndexError:
                sys.exit(f"ERROR: {input_file_path}: record without END/SVTYPE/CT: {fields[0]}:{fields[1]}")

            # Format the extracted information
            output_line = f"{chr}\t{start}\t{chr2}\t{end}\t{svtype}\t{ct}"

            # Write the formatted line to the output file
            output_file.write(output_line + '\n')

    print(f"Processing complete. Results saved to {output_file_path}")
