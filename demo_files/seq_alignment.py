# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 09:20:24 2026

@author: Marina Moro López
"""
    
from tkinter.filedialog import askopenfile
import pandas as pd

def read_sequence_from_file(prompt):
    print(prompt)
    file = askopenfile(mode='r')
    sequence = file.readlines()[1:]
    sequence = ''.join(sequence).replace('\n', '')
    return sequence

def align_sequences(gene_seq, patient_seq):
    alignment = [gene_seq[position] == patient_seq[position] for position in range(len(gene_seq))]
    return alignment

def get_mutation_positions(alignment):
    mutation_positions = [position for position, match in enumerate(alignment) if not match]
    return mutation_positions

def write_vcf(mutation_positions, gene_seq, patient_seq, filename="variants.vcf"):
    snps_detected = []
    vcf = open(filename, "w")
    vcf.write("##fileformat=VCFv4.2\n")
    vcf.write("##source=T3chFestsequencer\n")
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for pos in mutation_positions:
        ref = gene_seq[pos]
        alt = patient_seq[pos]
        if len(ref) == 1 and len(alt) == 1:
            snp_entry = {
                "chrom": "1",
                "pos": pos + 1,
                "ref": ref,
                "alt": alt,
            }
            snps_detected.append(snp_entry)
            vcf.write(f"1\t{pos+1}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")
    vcf.close()
    return snps_detected

def annotate_snps(snps_detected):
    annotated_results = []
    snp_database = {
        ("1", 2380): {"gene": "TPMT", "snpref": "rs1142345", "ref": "G", "alt": "A", "effect": "Alta toxicidad en tiopurinas"},
        ("1", 75203): {"gene": "TCF7L2", "snpref": "rs7903146", "ref": "C", "alt": "T", "effect": "Aumenta riesgo de diabetes tipo 2"},
        ("1", 245761): {"gene": "HBB", "snpref": "rs334", "ref": "A", "alt": "T", "effect": "Causa anemia falciforme"}
        }
    for snp in snps_detected:
        key = (snp["chrom"], snp["pos"])
        if key in snp_database:
            snp_info = snp_database[key]
            if snp["ref"] == snp_info["ref"] and snp["alt"] == snp_info["alt"]:
                gene = snp_info["gene"]
                effect = snp_info["effect"]
                snpref = snp_info["snpref"]
            else:
                gene = "Desconocido"
                effect = "Sin anotación"
                snpref = "-"
        else:
            gene = "Desconocido"
            effect = "Sin anotación"
            snpref = "-"
        annotated_results.append({
            "Cromosoma": snp["chrom"],
            "Posición": snp["pos"],
            "Base ref": snp["ref"],
            "Base alt": snp["alt"],
            "Gen": gene,
            "SNP": snpref,
            "Efecto": effect
        })
    return pd.DataFrame(annotated_results)


def main():

    gene_seq = read_sequence_from_file('Please select the file with the gene of reference')
    patient_seq = read_sequence_from_file('Please select the file with the patient sequence')
    alignment = align_sequences(gene_seq, patient_seq)
    mutation_positions = get_mutation_positions(alignment)
    snps_detected = write_vcf(mutation_positions, gene_seq, patient_seq)
    annotated_snps = annotate_snps(snps_detected)
    print("\nAnnotated SNPs:")
    print(annotated_snps.to_string(index=False))

if __name__ == "__main__":
    main()
