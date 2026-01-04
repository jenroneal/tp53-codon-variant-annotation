import requests

def fetch_tp53cds():
    url = "https://rest.ensembl.org/sequence/id/ENST00000269305?type=cds"
    headers = { "Content-Type" : "text/x-fasta" }
    response = requests.get(url, headers=headers)
    if not response.ok:
        raise Exception("Failed to obtain data from Ensembl")
    with open("tp53_cds.fasta", "w") as f:
        f.write(response.text)

print("TP53 CDS saved to the tp53_cds.fasta")
fetch_tp53cds()

# Load FASTA file
fasta_path = "tp53_cds.fasta"
sequence = ''
with open(fasta_path, 'r') as f:
    for line in f:
        if not line.startswith('>'):
            sequence += line.strip().upper()
print(" Loaded TP53 CDS:")
print(sequence[:60] + "...")
print(f"Total length: {len(sequence)} bp")

def load_variants(file_path):
    variants = []
    with open(file_path, 'r') as f:
        next(f)
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) != 3:
                raise ValueError(f"Malformed line in TSV: {repr(line)}")
            pos = int(parts[0]) - 1
            ref = parts[1].upper()
            alt = parts[2].upper()
            variants.append((pos, ref, alt))
    return variants

codon_table = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def translate_codon(codon):
    return codon_table.get(codon.upper(), '?')

def translate_full_cds(dna_seq):
    protein_seq = ''
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        protein_seq += codon_table.get(codon, '?')
    return protein_seq

def longest_subsequence(seq1, seq2):
    m, n = len(seq1), len(seq2)
    dp = [["" for _ in range(n + 1)] for _ in range(m + 1)]
    for i in range(m):
        for j in range(n):
            if seq1[i] == seq2[j]:
                dp[i + 1][j + 1] = dp[i][j] + seq1[i]
            else:
                dp[i + 1][j + 1] = max(dp[i][j + 1], dp[i + 1][j], key=len)
    return dp[m][n]

def annotate_variant(sequence, position, ref, alt):
    if sequence[position] != ref:
        return "Reference mismatch", None
    codon_start = position - (position % 3)
    wt_codon = sequence[codon_start:codon_start+3]
    codon_list = list(wt_codon)
    codon_list[position % 3] = alt
    mut_codon = ''.join(codon_list)
    wt_aa = translate_codon(wt_codon)
    mut_aa = translate_codon(mut_codon)
    if wt_aa == mut_aa:
        effect = "Synonymous"
    elif mut_aa == '*':
        effect = "Nonsense"
    else:
        effect = "Missense"
    return effect, (wt_codon, mut_codon, wt_aa, mut_aa)

with open("variant_annotation_results.tsv", "w") as out:
    print("Variant annotation complete.")
    out.write("Position\tRef\tAlt\tWT_Codon\tMut_Codon\tWT_AA\tMut_AA\tEffect\n")
    for pos, ref, alt in load_variants("variants.tsv"):
        effect, details = annotate_variant(sequence, pos, ref, alt)
        if details:
            wt_codon, mut_codon, wt_aa, mut_aa = details
            out.write(f"{pos+1}\t{ref}\t{alt}\t{wt_codon}\t{mut_codon}\t{wt_aa}\t{mut_aa}\t{effect}\n")
        else:
            out.write(f"{pos+1}\t{ref}\t{alt}\tN/A\tN/A\tN/A\tN/A\t{effect}\n")


mutated_sequence = list(sequence)
for pos, ref, alt in load_variants("variants.tsv"):
    if sequence[pos] == ref:
        mutated_sequence[pos] = alt
mutated_sequence = ''.join(mutated_sequence)

wt_protein = translate_full_cds(sequence)
mut_protein = translate_full_cds(mutated_sequence)
lcs = longest_subsequence(wt_protein, mut_protein)
lcs_len = len(lcs)
protein_len = len(wt_protein)
lcs_percent = round((lcs_len / protein_len) * 100, 2)

print("\nLongest Common Subsequence (LCS):")
print(f"LCS ({lcs_len} aa): {lcs}")
print(f"Protein Length: {protein_len}")
print(f"LCS % Similarity: {lcs_percent}%")
