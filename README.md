# TP53 Codon-Level Functional Annotation of Somatic Variants

This project performs codon-level functional annotation of somatic variants in the TP53 gene using sequence translation and comparative protein analysis.

Single-nucleotide variants (SNVs) are classified as synonymous, missense, or nonsense mutations based on codon-level changes in the translated protein sequence. The workflow also performs longest common subsequence (LCS) analysis between wild-type and mutated protein sequences to estimate sequence-level functional disruption.

## Workflow Overview

The workflow performs the following steps:

1. Download the canonical TP53 coding sequence from Ensembl
2. Parse somatic variants from a tab-separated input file
3. Validate reference nucleotide matches
4. Map variants to affected codons
5. Translate mutated coding sequences into protein sequences
6. Classify variants as synonymous, missense, or nonsense mutations
7. Compute longest common subsequence (LCS) similarity between wild-type and mutated proteins
8. Export annotated variant results to a tab-separated output file

## Input Files

### `variants.tsv`

Tab-separated file containing:

- Variant position
- Reference nucleotide
- Alternate nucleotide

### `tp53_cds.fasta`

Canonical TP53 coding sequence downloaded automatically from Ensembl using transcript:

`ENST00000269305`

## Output Files

### `variant_annotation_results.tsv`

Tab-separated summary of:

- Codon changes
- Amino acid changes
- Functional mutation classification
- Protein similarity metrics

### LCS Summary

Terminal output displaying longest common subsequence similarity between wild-type and mutated protein sequences.

## Tools and Technologies

- Python
- Ensembl REST API
- Sequence translation logic
- Longest Common Subsequence (LCS) analysis

## Biological Relevance

TP53 is one of the most commonly mutated tumor suppressor genes in human cancer. This workflow demonstrates foundational concepts in computational cancer genomics, including codon-level variant interpretation, protein consequence analysis, and sequence-based functional annotation.

## Usage

Ensure Python 3 is installed.

Install required dependency:

```bash
pip install requests
