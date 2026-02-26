from Bio import Entrez, SeqIO
import requests
import os
from io import StringIO 
from Bio.SeqRecord import SeqRecord


Entrez.email = "steph@gmail.com"

def fetch_from_ncbi(accession_list, database="protein", output_format="fasta"):
    sequences = {}

    for acc in accession_list:
        try:
            handle = Entrez.efetch(db=database, id=acc, rettype=output_format, retmode="text")
            records = list(SeqIO.parse(handle, output_format))  # Convert iterator to list
            sequences[acc] = records if len(records) == 1 else records[0]  # Store single record or list
            handle.close()
        except Exception as e:
            print(f"Encountered an error with {acc}: {str(e)}") 
    
    return sequences 

def fetch_from_uniprot(uni_ids):
    sequences = {}
    for up_id in uni_ids:
        try:
            # Updated to use modern UniProt REST API
            url = f"https://www.uniprot.org/uniprot/{up_id}.fasta"
            response = requests.get(url)
            if response.status_code == 200:
                handle = StringIO(response.text)
                record = SeqIO.read(handle, "fasta")
                sequences[up_id] = record 
            else:
                print(f"Failed to fetch {up_id}: HTTP {response.status_code}")
        except Exception as e:
            print(f"There was an error with {up_id}: {str(e)}")
    
    return sequences  

def download_sample_dataset():

    human_alpha = fetch_from_ncbi(['NP_000558'])
    human_beta = fetch_from_ncbi(['NP_000509'])

    other_species = fetch_from_ncbi([
        'NP_032243',
        'NP_032244',
        'NP_001107728',
        'XP_015733388',
        'NP_001094389'
    ])

    # print(human_alpha)
    # print("Next")
    # print(human_beta)
    # print("Next")
    # print(other_species)
    all_sequences = {**human_alpha, **human_beta, **other_species}
    # print("Next")
    # print(all_sequences)

    output_file = "haemoglobin_sequences.fasta"

    with open(output_file, "w") as f:
        for acc, record in all_sequences.items():
            print(record)
            SeqIO.write(record, f, "fasta")

    print(f"\n✓ Saved {len(all_sequences)} sequences to {output_file}")
    
    return all_sequences


if __name__ == "__main__":
    sequences = download_sample_dataset()

