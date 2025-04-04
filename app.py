import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
from collections import defaultdict
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import biotite.structure as bs

# Set page config
st.set_page_config(layout='wide')

# Session state initialization
if 'pdb_string' not in st.session_state:
    st.session_state.pdb_string = None
if 'b_value' not in st.session_state:
    st.session_state.b_value = None
if 'protein_name' not in st.session_state:
    st.session_state.protein_name = "Predicted Protein"

# Shared function to predict structure
def update(sequence):
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', 
                           headers=headers, 
                           data=sequence)
    st.session_state.pdb_string = response.content.decode('utf-8')
    with open('predicted.pdb', 'w') as f:
        f.write(st.session_state.pdb_string)
    struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
    st.session_state.b_value = round(struct.b_factor.mean(), 4)

# Predictor tab
def emsfold_app():
    st.sidebar.title('ProtoAnalyzer')
    st.sidebar.write("Predict protein structures from sequence")

    DEFAULT_SEQ = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
    txt = st.sidebar.text_area('Input sequence', DEFAULT_SEQ, height=275)
    
    custom_name = st.sidebar.text_input("Protein Name", st.session_state.protein_name)
    st.session_state.protein_name = custom_name

    if st.sidebar.button('⏳ Predict Structure'):
        with st.spinner('Predicting structure...'):
            update(txt)
            st.success("Prediction complete! Switch to the Analyzer tab to view details.")

    st.sidebar.title('Display Options')
    background_color = st.sidebar.color_picker("Background", "#000000")
    show_labels = st.sidebar.checkbox("Show Residue Labels", False)

    if st.session_state.pdb_string:
        col1, col2 = st.columns([2, 1])
        with col1:
            st.subheader(f'🧬 {st.session_state.protein_name} Structure')
            st.caption(f"Confidence score: {st.session_state.b_value}")
            
            pdbview = py3Dmol.view()
            pdbview.addModel(st.session_state.pdb_string, 'pdb')
            pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
            pdbview.setBackgroundColor(background_color.lower())
            if show_labels:
                pdbview.addResLabels()
            pdbview.zoomTo()
            pdbview.spin(True)
            showmol(pdbview, height=500, width=800)

            st.download_button("📥 Download PDB",
                data=st.session_state.pdb_string,
                file_name=f'{st.session_state.protein_name.replace(" ", "_")}.pdb',
                mime='text/plain')
        
        with col2:
            st.subheader('📊 Confidence Scores')
            st.markdown("""
            | Color | plDDT Score | Confidence Level |
            |-------|-------------|------------------|
            | 🔵    | 90-100      | Very High        |
            | 🟢    | 70-90       | High             |
            | 🟡    | 50-70       | Medium           |
            | 🔴    | <50         | Low              |
            """)

            st.subheader('🧪 Protein Properties')
            protein_seq = ProteinAnalysis(txt)
            hydrophobic = sum(txt.count(res) for res in 'AILMFWYV')
            data = {
                "Property": ["Length", "MW (Da)", "Hydrophobicity", "Net Charge", "Avg Confidence"],
                "Value": [
                    len(txt),
                    f"{protein_seq.molecular_weight()/1000:.1f} kDa",
                    f"{hydrophobic/len(txt)*100:.1f}%",
                    sum(txt.count(res) for res in 'KRH') - sum(txt.count(res) for res in 'DE'),
                    st.session_state.b_value
                ]
            }
            st.table(data)
    else:
        st.info("💡 Enter a protein sequence and click 'Predict Structure'")

# Analyzer tab
def ranaatom_app():
    st.title("🔍 PDB Analysis Toolkit")

    if not st.session_state.pdb_string:
        st.warning("No structure available. Please predict a structure first in the Predictor tab.")
        return

    st.write(f"Analyzing: {st.session_state.protein_name}")

    st.sidebar.title('Analysis Settings')
    style = st.sidebar.selectbox("Style", ["cartoon", "sphere", "stick", "surface"], index=0)
    color_scheme = st.sidebar.selectbox("Color Scheme", ["spectrum", "chain", "residue"], index=0)
    show_labels = st.sidebar.checkbox("Show Atom Labels", False)

    col1, col2 = st.columns(2)
    with col1:
        st.subheader(f"🔬 {st.session_state.protein_name} Structure")
        atom_lines = [line for line in st.session_state.pdb_string.split('\n') if line.startswith("ATOM")]
        num_atoms = len(atom_lines)
        st.caption(f"{num_atoms:,} atoms | {len(set(line[21] for line in atom_lines))} chains")

        view = py3Dmol.view(width=600, height=400)
        view.addModel(st.session_state.pdb_string, "pdb")

        if style == "cartoon":
            view.setStyle({'cartoon': {'color': color_scheme}})
        elif style == "sphere":
            view.setStyle({'sphere': {'colorscheme': color_scheme}})
        elif style == "stick":
            view.setStyle({'stick': {'colorscheme': color_scheme}})
        elif style == "surface":
            view.addSurface(py3Dmol.VDW, {'opacity': 0.7, 'color': 'white'})

        if show_labels:
            view.addResLabels()

        view.zoomTo()
        showmol(view, height=400)

        st.download_button("⬇️ Download PDB",
            data=st.session_state.pdb_string,
            file_name=f'{st.session_state.protein_name.replace(" ", "_")}.pdb',
            mime='text/plain')

    with col2:
        st.subheader("Residue Distribution")
        res_counts = defaultdict(int)
        for line in atom_lines:
            res_name = line[17:20].strip()
            res_counts[res_name] += 1
        res_df = pd.DataFrame.from_dict(res_counts, orient='index', columns=['Count'])
        st.bar_chart(res_df)

        # Heatmap of residue-residue distances (C-alpha atoms)
        st.subheader("Residue-Residue Distance Heatmap")
        struct = bsio.load_structure("predicted.pdb")
        ca_atoms = struct[struct.atom_name == "CA"]
        coords = ca_atoms.coord
        dist_matrix = np.linalg.norm(coords[:, np.newaxis, :] - coords[np.newaxis, :, :], axis=-1)

        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(dist_matrix, cmap='coolwarm', ax=ax)
        ax.set_title("Distance Map (Cα-Cα)")
        ax.set_xlabel("Residue Index")
        ax.set_ylabel("Residue Index")
        st.pyplot(fig)

        # Additional analysis: Residue Hydrophobicity Profile
        st.subheader("Hydrophobicity Profile")
        sequence = ''.join([line[13] for line in atom_lines if line[12:16].strip() == "CA"])
        hydropathy_index = {
            'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
            'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
            'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
            'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
        }
        hydrophobicity = [hydropathy_index.get(aa, 0.0) for aa in sequence]
        fig2, ax2 = plt.subplots(figsize=(8, 3))
        ax2.plot(hydrophobicity, label="Hydropathy")
        ax2.set_xlabel("Residue Index")
        ax2.set_ylabel("Hydropathy")
        ax2.set_title("Hydrophobicity Profile")
        ax2.legend()
        st.pyplot(fig2)

# Main layout
st.markdown(''' 
    <div style="text-align: center;">
        <h1 style="font-family: 'Dustosmo Roman', 'Times New Roman', serif; color:lightyellow;">
            ProtoAnalyzer
        </h1>
        <p>
            <em>A Streamlit <strong>Component</strong> for creating and analyzing protein structures.</em>
        </p>
    </div>
    ''', unsafe_allow_html=True) 

tab1, tab2 = st.tabs(["🧬 Predictor", "🔍 Analyzer"])
with tab1:
    emsfold_app()
with tab2:
    ranaatom_app()
