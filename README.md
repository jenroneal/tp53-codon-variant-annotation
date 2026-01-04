# TP53 Codon-Level Functional Annotation of Somatic Variants
This project annotates the effects of somatic mutations in the TP53 gene using codon-level translation and sequence comparison. Each variant is functionally classified as synonymous, missense, or nonsense. The output also includes a longest common subsequence (LCS) comparison between the wild-type and mutated protein sequences to estimate functional disruption.

## Description
- The script downloads the standard TP53 coding sequence from Ensembl using its transcript ID.
- It parses a list of SNVs from a tab-separated input file (variants.tsv).
- Each variant is checked for a reference match and mapped to its corresponding codon.
- A longest common subsequence (LCS) is calculated between the wild-type and mutated proteins.
- The final output summarizes all annotations in a .tsv file.

## File Input
- variants.tsv - A tab-separated file listing the position, reference base, and alternate base for each mutation to be analyzed.
- tp53_cds.fasta — The TP53 coding sequence used to identify affected codons. This file is automatically downloaded from Ensembl using the standard transcript (ENST00000269305), which is commonly used in bioinformatics analyses and literature.

## File Output
- tp53_cds.fasta – the downloaded TP53 coding sequence
- variant_annotation_results.tsv – tab-separated summary of mutation effects
- LCS summary – printed to terminal showing percent similarity between WT and mutated protein


## How to Run
1. Ensure Python 3 is installed.
2. Install dependencies:
   pip install requests
3. Prepare a variants file (see variants_example.tsv).
4. Run:
   python tp53_variant_annotation.py
