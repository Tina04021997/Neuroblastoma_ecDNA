import glob
import os

# ── Paths: edit before running ──────────────────────────────────────────────
# Input: <sample>_sv.txt files from 02_filter_delly_vcf.sh
input_files_pattern = '/path/to/shatterseek/sv_files/*.txt'

# Specify the output directory
output_directory = '/path/to/shatterseek/sv_files/SVs'

# Get a list of input files that match the pattern
input_files = glob.glob(input_files_pattern)

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
            fields = input_line.strip().split("\t")

            chr = fields[0].replace("chr", "")
            start = fields[1]
            info_field = fields[7].split(";")
            chr2 = [item.split("=")[1] for item in info_field if item.startswith("CHR2=")][0].replace("chr", "")
            end = [item.split("=")[1] for item in info_field if item.startswith("END=")][0]
            svtype = [item.split("=")[1] for item in info_field if item.startswith("SVTYPE=")][0]
            ct = [item.split("=")[1] for item in info_field if item.startswith("CT=")][0]

            # Format the extracted information
            output_line = f"{chr}\t{start}\t{chr2}\t{end}\t{svtype}\t{ct}"

            # Write the formatted line to the output file
            output_file.write(output_line + '\n')

    print(f"Processing complete. Results saved to {output_file_path}")
