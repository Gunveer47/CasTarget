from Bio.Seq import Seq
import csv
def find_crispr_guides(dna_sequence):
    """
    Scans a DNA sequence for SpCas9 target sites (N20 + NGG)
    on both forward and reverse strands.
    """
    seq = Seq(dna_sequence)
    results = []
    
    strands = {"Forward": seq, "Reverse": seq.reverse_complement()}
    
    for strand_name, s in strands.items():
        for i in range(len(s) - 3):
            pam = s[i+20:i+23]
            if pam.endswith("GG") and len(pam) == 3:
                spacer = s[i:i+20]
                
                gc_content = (spacer.count("G") + spacer.count("C")) / len(spacer) * 100
                
                if 40 <= gc_content <= 60:
                    results.append({
                        "Strand": strand_name,
                        "Position": i,
                        "Spacer": str(spacer),
                        "PAM": str(pam),
                        "GC_Content": round(gc_content, 2)
                    })
                    
    return results

#test_gene = "ATGGCCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCCGGATAGCTAGCTAGCTAG"
#guides = find_crispr_guides(test_gene)
print("--- Welcome to CasTarget: CRISPR Guide RNA Finder ---")
user_sequence = input("Enter or paste your DNA sequence: ").strip().upper()

guides = find_crispr_guides(user_sequence)
for g in guides:
    print(f"Found on {g['Strand']} strand at pos {g['Position']}: Spacer = {g['Spacer']}, PAM = {g['PAM']}, GC% = {g['GC_Content']}")
csv_filename = "crispr_results.csv"
keys = ["Strand", "Position", "Spacer", "PAM", "GC_Content"]

with open(csv_filename, mode="w", newline="") as output_file:
    dict_writer = csv.DictWriter(output_file, fieldnames=keys)
    dict_writer.writeheader()
    dict_writer.writerows(guides)

print(f"\nSuccessfully saved {len(guides)} guide targets to {csv_filename}!")