import argparse
#Example use is 
# python3 parseFastq.py --fastq /home/rbif/week6/hawkins_pooled_sequences.fastq


################################################
# Pre-made script
class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
 
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
         
    def __iter__(self):
        return self
     
    def __next__(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
         
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
         
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)
##########################################################################
############################## STEP 1 ##############################

import argparse
import os

# Trim reads based on quality scores
def trim_reads(reads, quality_scores):
    for i in range(len(quality_scores)):
        qs = quality_scores[i]  
        end_index = len(qs) - 1 
        # Start iterating from the end until it runs into D or F
        while end_index >= 0 and (qs[end_index] == 'D' or qs[end_index] == 'F'):
            end_index -= 1
        # Check if 'end_index' moved from its initial position
        if end_index < len(qs) - 1:
            reads[i] = reads[i][:end_index+1]
            # Trim the quality scores
            quality_scores[i] = quality_scores[i][:end_index+1]
    return reads, quality_scores

# Check if script is being executed alone or as a module in another script
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastq", required=True, help="Place fastq inside here")
    args = parser.parse_args()
    
    # Read the barcode information from the text file
    barcodes_file = 'harrington_clinical_data.txt'  
    barcodes_dict = {}
    with open(barcodes_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            name = cols[0]
            barcode = cols[2]
            barcodes_dict[barcode] = name

    # Process the fastq file
    fastq_file = args.fastq
    output_dir = 'fastqs'  # Output directory

    # Create the output directory if it doesn't exist (prevents duplicates while testing)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    with open(fastq_file, 'r') as f:
        for fastq_obj in ParseFastQ(fastq_file):
            barcode_len = 5
            sequence = fastq_obj[1][barcode_len:]
            quality_score = fastq_obj[3][barcode_len:]
            barcode = fastq_obj[1][:barcode_len]  
            if barcode in barcodes_dict:
                name = barcodes_dict[barcode]
                output_file = os.path.join(output_dir, '{}_trimmed.fastq'.format(name))
                sequence, quality_score = trim_reads([sequence], [quality_score])
                with open(output_file, 'a') as out_f:
                    out_f.write(f"{fastq_obj[0]}\n{sequence[0]}\n{fastq_obj[2]}\n{quality_score[0]}\n")


############################## STEP 2 ##############################

import os

# create output directory for the sam files
output_dir = "sam_files"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# Create output directory for BWA index files
index_dir = "bwa_index"
if not os.path.exists(index_dir):
    os.mkdir(index_dir)

# Index the reference genome and generate index files within the index directory
reference_genome = "dgorgon_reference.fa"
os.system(f"bwa index -p {index_dir}/dgorgon_reference {reference_genome}")

# loop through all .fastq files in the fastqs directory
for file in os.listdir("fastqs"):
    if file.endswith("_trimmed.fastq"):
        # run bwa mem command
        name = file.split("_")[0]
        input_file = os.path.join("fastqs", file)
        output_file = os.path.join(output_dir, f"{name}.sam")
        os.system(f"bwa mem {index_dir}/dgorgon_reference {input_file} > {output_file}")


############################## STEP 3 ##############################

import os

# Create output directory for sorted BAM files if it doesn't exist
output_dir = "sorted_bam"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# Loop through all .sam files in the directory
for file in os.listdir("sam_files"):
    if file.endswith(".sam"):
        # Construct paths
        name = file[:-4]
        input_file = os.path.join("sam_files", file)
        bam_file = os.path.join(output_dir, f"{name}.bam")  # Store BAM files in the sorted_bam directory

        # Convert SAM to BAM
        os.system(f"samtools view -bS {input_file} > {bam_file}")

        # Sort BAM file
        sorted_file = os.path.join(output_dir, f"{name}.sorted.bam")
        os.system(f"samtools sort -m 100M -o {sorted_file} {bam_file}")

        # Index sorted BAM file
        os.system(f"samtools index {sorted_file}")

        # Remove SAM and unsorted BAM files
        os.remove(input_file)
        os.remove(bam_file)

# Remove all remaining BAM files in sam_files directory
for file in os.listdir("sam_files"):
    if file.endswith(".bam"):
        os.remove(os.path.join("sam_files", file))

# Remove the empty "sam_files" directory
os.rmdir("sam_files")


############################## STEP 4 ##############################

import pysam
import os

def load_reference_sequence(reference_fa_path):
    # Load the reference sequence from a FASTA file
    with open(reference_fa_path, "r") as reference_file:
        reference_sequence = ""
        for line in reference_file:
            if not line.startswith(">"):  # Skip header lines
                reference_sequence += line.strip()
    return reference_sequence

def pileup():
    # Define the threshold for SNP detection
    threshold = 0.0

    # Load the reference sequence
    reference_fa_path = "dgorgon_reference.fa"  # Update with the path to your reference sequence
    reference_sequence = load_reference_sequence(reference_fa_path)

    # Create a file to write the output to
    output_file = open("mutations-report.txt", "w")

    # Loop through all files in the 'sorted_bam' directory
    for filename in os.listdir('sorted_bam'):
        if filename.endswith('.sorted.bam'):
            filepath = os.path.join('sorted_bam', filename)
            output_file.write(f"Processing file: {filepath}\n")
            print(f"Processing file: {filepath}")
            # Open the BAM file
            samfile = pysam.AlignmentFile(filepath, "rb")

            for pileupcolumn in samfile.pileup():
                output_file.write(f"coverage at base {pileupcolumn.pos} = {pileupcolumn.n}\n")

                # Reset the dictionary and total reads for each pileup column
                ntdict = {}
                total_reads = 0

                # Compare bases to the reference sequence
                reference_base = reference_sequence[pileupcolumn.pos].upper()  # Convert to uppercase for consistency
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        base = pileupread.alignment.query_sequence[pileupread.query_position].upper()  # Convert to uppercase

                        # Populate the ntdict with the counts of each base
                        if base in ntdict:
                            ntdict[base] += 1
                        else:
                            ntdict[base] = 1

                        total_reads += 1

                # Determine if it's a SNP based on the reference sequence
                for base, count in ntdict.items():
                    if base != "N" and base != reference_base:
                        frequency = count / total_reads
                        if frequency > threshold:
                            output_file.write(f"SNP at position {pileupcolumn.pos + 1}: {base} (Frequency: {frequency:.2f})\n")

            samfile.close()

    # Close the output file
    output_file.close()

if __name__ == "__main__":
    pileup()


############################## STEP 5 ##############################

import re
from collections import defaultdict

# read clinical data file
clinical_data = {}
with open('harrington_clinical_data.txt') as f:
    for line in f:
        name, color, barcode = line.strip().split('\t')
        clinical_data[name] = (color, barcode)

# read fastq files and count reads
read_counts = defaultdict(int)
import glob
for filename in glob.glob('fastqs/*_trimmed.fastq'):
    name = filename.split('/')[1].split('_')[0]
    with open(filename) as f:
        for line in f:
            if line.startswith('@'):
                read_counts[name] += 1

# process mutations report file and write output to report.txt
with open('mutations-report.txt') as f:
    with open('report.txt', 'w') as out:
        current_name = ""
        sample_info = ""
        for line in f:
            if line.startswith('Processing file'):
                name = line.strip().split('/')[1].split('.')[0]
                color, barcode = clinical_data.get(name, ("Unknown", "Unknown"))
                read_count = read_counts.get(name, 0)
                sample_info = f'Sample {name} had a {color} mold, {read_count} reads'
                current_name = name
            elif line.startswith('SNP at position'):
                match = re.search(r'SNP at position (\d+): (\w) \(Frequency: (\d+\.\d+)', line)
                if match:
                    position = int(match.group(1))
                    base = match.group(2)
                    frequency = float(match.group(3))
                    mutation = f'{base} at position {position}'
                    sample_info += f', and had {frequency:.2%} of the reads at {mutation}.'
                    out.write(sample_info + '\n')

print("All done! Please review the following contents: directory 'fastqs', directory 'sorted_bam', and file 'report.txt'.")
