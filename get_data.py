import requests
import json

def fetch_premrna_data(gene_ids):
    base_url = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json"}
    premrna_data = []

    for gene_id in gene_ids:
        # Obter a sequência de DNA do gene
        endpoint = f"/sequence/id/{gene_id}?type=genomic"
        response = requests.get(base_url + endpoint, headers=headers)
        if response.status_code == 200:
            data = response.json()
            dna_sequence = data["seq"]
            rna_sequence = transcribe_dna_to_rna(dna_sequence)
            
            # Obter informações de transcritos, éxons e íntrons
            transcript_info = fetch_transcript_info(gene_id)
            
            premrna_data.append({
                "gene_id": gene_id,
                "dna_sequence": dna_sequence,
                "rna_sequence": rna_sequence,
                "transcript_info": transcript_info
            })
        else:
            print(f"Failed to fetch data for gene ID: {gene_id}")

    return premrna_data

def transcribe_dna_to_rna(dna_sequence):
    return dna_sequence.replace('T', 'U')

def fetch_transcript_info(gene_id):
    base_url = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json"}
    endpoint = f"/overlap/id/{gene_id}?feature=transcript"
    response = requests.get(base_url + endpoint, headers=headers)
    
    if response.status_code == 200:
        transcripts = response.json()
        transcript_data = []
        
        for transcript in transcripts:
            transcript_id = transcript["id"]
            exons_introns_info = fetch_exons_introns_info(transcript_id)
            
            transcript_data.append({
                "transcript_id": transcript_id,
                "exons_count": exons_introns_info["exons_count"],
                "introns_count": exons_introns_info["introns_count"],
                "exons": exons_introns_info["exons"],
                "introns": exons_introns_info["introns"]
            })
        
        return transcript_data
    else:
        print(f"Failed to fetch transcripts for gene ID: {gene_id}")
        return []

def fetch_exons_introns_info(transcript_id):
    base_url = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json"}
    endpoint = f"/lookup/id/{transcript_id}?expand=1"
    response = requests.get(base_url + endpoint, headers=headers)
    
    if response.status_code == 200:
        data = response.json()
        exons = data.get("Exon", [])
        introns = calculate_introns(exons)
        return {
            "exons_count": len(exons),
            "introns_count": len(introns),
            "exons": exons,
            "introns": introns
        }
    else:
        print(f"Failed to fetch exons and introns info for transcript ID: {transcript_id}")
        return None

def calculate_introns(exons):
    introns = []
    for i in range(len(exons) - 1):
        intron = {
            "start": exons[i]["end"] + 1,
            "end": exons[i+1]["start"] - 1
        }
        introns.append(intron)
    return introns

def save_to_json(data, filename):
    with open(filename, 'w') as json_file:
        json.dump(data, json_file, indent=4)

gene_ids = [
    "ENSG00000139618", "ENSG00000157764", "ENSG00000141510",
    "ENSG00000121410", "ENSG00000148584", "ENSG00000198727",
    "ENSG00000284733", "ENSG00000284734", "ENSG00000284735",
    "ENSG00000142611", "ENSG00000157911", "ENSG00000260972", 
    "ENSG00000284616", "ENSG00000142655", "ENSG00000226374", 
    "ENSG00000232596", "ENSG00000228037", "ENSG00000231510",
    "ENSG00000235054", "ENSG00000149527", "ENSG00000157916",
    "ENSG00000287586", "ENSG00000159423", "ENSG00000227169",
    "ENSG00000283356", "ENSG00000131697", "ENSG00000162444",
    "ENSG00000184677", "ENSG00000157881", "ENSG00000180758",
    "ENSG00000048707", "ENSG00000197921", "ENSG00000055070",
    "ENSG00000204138", "ENSG00000280222", "ENSG00000284745",
    "ENSG00000285833", "ENSG00000284692", "ENSG00000288398",
    "ENSG00000225196"
]

premrna_data = fetch_premrna_data(gene_ids)

save_to_json(premrna_data, "premrna_data.json")

print("Dados de pré-mRNA salvos em premrna_data.json")
