import streamlit as st
from Bio.Seq import Seq
import pandas as pd

st.set_page_config(page_title="CasTarget: CRISPR Guide Finder", page_icon="🧬", layout="centered")

st.title("🧬 CasTarget: CRISPR Guide RNA Finder")
st.write("Scan DNA sequences for SpCas9 targets, validate NGG PAM sites, and filter optimal spacers.")

dna_input = st.text_area(
    "Paste your target DNA sequence here:", 
    placeholder="e.g., ATGGCCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCCGGATAGCTAGCTAGCTAG"
)

def find_crispr_guides(dna_sequence):
    seq = Seq(dna_sequence.strip().upper())
    strands = {"Forward": seq, "Reverse": seq.reverse_complement()}
    guides = []
    
    for strand_name, s in strands.items():
        for i in range(len(s) - 23):
            spacer = s[i:i+20]
            pam = s[i+20:i+23]
            
            if len(pam) == 3 and pam.endswith("GG"):
                gc_content = (spacer.count("G") + spacer.count("C")) / len(spacer) * 100
                if 40 <= gc_content <= 60:
                    guides.append({
                        "Strand": strand_name,
                        "Position": i + 1,
                        "Spacer": str(spacer),
                        "PAM": str(pam),
                        "GC_Content": round(gc_content, 1)
                    })
    return guides

if st.button("Run CRISPR Analysis", key="run_analysis_button"):
    if not dna_input.strip():
        st.warning("Please enter a DNA sequence first.")
    elif len(dna_input.strip()) < 23:
        st.error("Sequence is too short! Please enter at least 23 bases.")
    else:
        guides = find_crispr_guides(dna_input)
        if len(guides) > 0:
            st.success(f"Successfully found {len(guides)} optimal guide targets!")
            df = pd.DataFrame(guides)
            st.dataframe(df, use_container_width=True)
            
            csv_data = df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 Download Results as CSV",
                data=csv_data,
                file_name="cas_target_results.csv",
                mime="text/csv",
                key="download_csv_button"
            )
        else:
            st.warning("No valid guides found matching the NGG and GC content criteria.")
