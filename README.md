# SNP snalysis for personalized medicine with Python
Python application developed for the T3chFest 2026 talk "Personalized medicine: how 0.1% of your DNA can change everything".
This project demonstrates how Python can be used to identify single nucleotide polymorphisms (SNPs), evaluate their potential biological impact and generate Variant Call Format (VCF) files.

## Features
- Load reference and patient DNA sequences from FASTA files
- Detect single nucleotide polymorphisms (SNPs)
- Annotate detected variants and their potential effects
- Generate VCF files for downstream genomic analyses
- Produce a simple mutation report

## Installation
You simply have to download the .py and fasta files (or get your own fasta files from a database as [NCBI](https://www.ncbi.nlm.nih.gov/)).

## Usage
Run the script and select:
1. The reference DNA sequence.
2. The patient's DNA sequence.
The program compares both sequences and generates:
- A list of detected SNPs
- A summary of their predicted effects
- A VCF file containing all identified variants

## Example
Input
```
Reference sequence: reference_sequence.fasta
Patient sequence: patient_sequence.fasta
```
Output
```
Detected SNPs: 12
Position: 324
Reference: C
Patient: T
Effect: Increased risk of diabetes
Position: 728
Reference: A
Patient: G
Effect: Unknown
Output file:
variants.vcf
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)
