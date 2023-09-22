import argparse
import requests
from Bio.Seq import Seq

# Define the MyGene API endpoint
mygene_url = "https://mygene.info/v3/query"

# Define the command-line argument parser
parser = argparse.ArgumentParser(description='Retrieve the DNA sequence for a gene')
parser.add_argument('gene_symbol', type=str, help='the symbol of the gene to retrieve the sequence for')

# Parse the command-line arguments
args = parser.parse_args()

# Set the gene symbol and parameters for the MyGene API
gene_symbol = args.gene_symbol
params = {"q": f"symbol:{gene_symbol}", "fields": "entrezgene,ensembl.gene"}

# Send a GET request to the MyGene API
response = requests.get(mygene_url, params=params)

# Check if the response was successful
if response.status_code == 200:
    # Extract Ensembl ID from the response
    hits = response.json()["hits"]
    if len(hits) > 0:
        ensembl_id = hits[0]["ensembl"]["gene"]
        print(f"Ensembl ID for {gene_symbol} is {ensembl_id}")
    else:
        print(f"No hits found for the specified gene symbol {gene_symbol}.")
else:
    print("Error: MyGene API request failed.")

# Define the Ensembl gene ID
ID = ensembl_id
str(ensembl_id)

server = "https://rest.ensembl.org"
ext = "/sequence/id/" + ensembl_id + "?content-type=text/plain"
r = requests.get(server + ext, headers={"Content-Type": "application/json"})

# Remove the extra characters from the sequence
dna_sequence = ''.join(filter(lambda char: char in ['T', 'G', 'A', 'C'], r.text))

# Find longest ORF
def find_longest_orf_and_translate(sequence):
    orfs = []
    for frame in [0, 1, 2]:
        translate_seq = sequence[frame:].translate(table=1, to_stop=True)
        orfs.append(translate_seq)

    longest_orf = max(orfs, key=lambda x: len(x))
    return longest_orf

# Use Bio.Seq to convert the DNA sequence to an amino acid sequence
aa_sequence = find_longest_orf_and_translate(Seq(dna_sequence))

# Write the sequences to a FASTA file
filename = f"{gene_symbol}_sequence.fasta"
with open(filename, 'w') as file:
    file.write(f'>{gene_symbol} nucleotide sequence\n')
    file.write(dna_sequence + '\n')
    file.write(f'>{gene_symbol} amino acid sequence\n')
    file.write(str(aa_sequence))

print(f"FASTA file can be found in {filename}")

# Define the search parameters
gene_id = ensembl_id
homology_type = 'orthologues'
content_type = 'application/json'

# Create the URL for the RESTful API call
ensembl_url = f'https://rest.ensembl.org/homology/id/{gene_id}?' \
      f'type={homology_type}&content-type={content_type}'

# Make the API call and get the response
try:
    response = requests.get(ensembl_url)
    response.raise_for_status()
except requests.exceptions.RequestException as e:
    print(f"Failed to retrieve homologs. {e}")
    exit()

# Parse the response
try:
    data = response.json()
    homologs = data['data'][0]['homologies']
except (KeyError, IndexError):
    print(f"No homologs found for gene {gene_id}")
    exit()

# Filter out duplicates by species
unique_species = set()
filtered_homologs = []
for homology in homologs:
    target_species = homology['target']['species']
    if target_species not in unique_species:
        filtered_homologs.append(homology)
        unique_species.add(target_species)

# Write the results to a file
with open(f"{gene_symbol}_homologs_list.txt", 'w') as f:
    f.write(f"Homologs for gene {gene_id}:\n")
    for homology in filtered_homologs:
        target_gene = homology['target'].get('gene', {}).get('display_label')
        target_species = homology['target']['species'].split()[0]
        f.write(f"\t {target_species}\n")
print(f"Results written to {gene_symbol}_homologs_list.txt")
