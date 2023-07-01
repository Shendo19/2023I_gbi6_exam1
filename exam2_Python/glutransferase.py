import requests
from Bio import Entrez, SeqIO
from Bio.SeqUtils import ProtParam
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter

def get_source_info(accession):
    Entrez.email = "shendo.carbonelly@est.ikiam.edu.ec"
    handle = Entrez.efetch(db="nucleotide", id=accession, retmode="xml")
    record = Entrez.read(handle)
    handle.close()

    source_info = record[0]["GBSeq_organism"]
    return source_info

def Source(accessions):
    source_info_list = []
    for accession in accessions:
        source_info = get_source_info(accession)
        source_info_list.append(source_info)

    species_count = Counter(source_info_list)
    species_df = pd.DataFrame({"Species": species_count.keys(), "Frequency": species_count.values()})
    species_df.to_csv("results/source.csv", index=False)

def get_sequence(accession):
    Entrez.email = "shendo.carbonelly@est.ikiam.edu.ec"  
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record.seq

def translate_sequence(dna_sequence):
    coding_dna = dna_sequence.translate(to_stop=True)
    return coding_dna

def get_peptides(dna_sequence):
    coding_dna = translate_sequence(dna_sequence)
    peptides = []
    for i in range(len(coding_dna)):
        peptide = coding_dna[i:]
        if peptide.startswith("M"):
            peptides.append(peptide)
    return peptides

def calculate_molecular_weight(peptide):
    return ProtParam.ProteinAnalysis(str(peptide)).molecular_weight()

def calculate_instability_index(peptide):
    return ProtParam.ProteinAnalysis(str(peptide)).instability_index()

def generate_joinplot(peptides):
    data = pd.DataFrame({
        "Molecular Weight": [calculate_molecular_weight(peptide) for peptide in peptides],
        "Instability Index": [calculate_instability_index(peptide) for peptide in peptides]
    })

    return data

def Sequences(accessions):
    all_data = pd.DataFrame()

    for accession in accessions:
        dna_sequence = get_sequence(accession)
        peptides = get_peptides(dna_sequence)
        data = generate_joinplot(peptides)
        all_data = pd.concat([all_data, data])

    sns.set(style="ticks")
    g = sns.jointplot(data=all_data, x="Molecular Weight", y="Instability Index", kind="scatter", color="red", height=6)
    g.ax_joint.set_xlabel("Molecular Weight")
    g.ax_joint.set_ylabel("Instability Index")
    g.fig.suptitle("Peptide Molecular Weight vs Instability Index")
    plt.savefig("results/glupeptides.png")
    plt.close()

    all_data.to_csv("results/glupeptides.csv", index=False)