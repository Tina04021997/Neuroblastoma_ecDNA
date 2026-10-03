import os

# ── Paths: edit before running ──────────────────────────────────────────────
# Input: <sample>.cs.rmdup_CNV_CALLS.bed files from AmpliconSuite (CNVkit output)
input_directory = '/path/to/shatterseek/cn_files'
output_directory = '/path/to/shatterseek/cn_files/CNs'

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

# Loop through all files in the input directory
for filename in os.listdir(input_directory):
    if filename.endswith('.bed'):
        input_path = os.path.join(input_directory, filename)
        
        # Modify output file name
        output_filename = filename.replace('.cs.rmdup_CNV_CALLS.bed', '_cn.txt')
        output_path = os.path.join(output_directory, output_filename)
        print(f"Processing {output_filename}")

        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            # Write the header to the output file
            outfile.write("chromosome\tstart\tend\tCN\n")
            for line in infile:
                # Split the line into columns
                columns = line.strip().split('\t')

                # Extract chromosome number
                chromosome = columns[0].replace('chr', '')

                # Check if the chromosome is 1-22 or X
                if chromosome.isdigit() and 1 <= int(chromosome) <= 22 or chromosome == 'X':
                    # Remove 'chr' and column 4, round the value in column 5
                    rounded_value = round(float(columns[4]))
                    new_line = '\t'.join([chromosome] + columns[1:3] + [str(rounded_value)]) + '\n'
                    outfile.write(new_line)

print("Processing complete. Check the output directory for the results.")
