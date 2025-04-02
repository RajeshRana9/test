import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
from collections import defaultdict
import re

# Set page config
st.set_page_config(layout='wide')

# Initialize session state
if 'pdb_string' not in st.session_state:
    st.session_state.pdb_string = None
if 'b_value' not in st.session_state:
    st.session_state.b_value = None
if 'protein_name' not in st.session_state:
    st.session_state.protein_name = "Predicted Protein"

# ========================
# SHARED FUNCTIONS
# ========================
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

# ========================
# EMSFold APP (Prediction)
# ========================
def emsfold_app():
    st.sidebar.title('ProtoAnalyzer')
    st.sidebar.write("Predict protein structures from sequence")
    
    # Sequence input
    DEFAULT_SEQ = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
    txt = st.sidebar.text_area('Input sequence', DEFAULT_SEQ, height=275)
    
    # Protein name input
    custom_name = st.sidebar.text_input("Protein Name", st.session_state.protein_name)
    st.session_state.protein_name = custom_name

    if st.sidebar.button('⏳ Predict Structure'):
        with st.spinner('Predicting structure...'):
            update(txt)
            st.success("Prediction complete! Switch to the Analyzer tab to view details.")
            st.session_state.pdb_string = st.session_state.pdb_string  # Ensure persistence
    
    # Visualization settings
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
            
            st.download_button(
                label="📥 Download PDB",
                data=st.session_state.pdb_string,
                file_name=f'{st.session_state.protein_name.replace(" ", "_")}.pdb',
                mime='text/plain'
            )
        
        with col2:
            st.subheader('📊 Confidence Scores')
            color_table = """
            | Color | plDDT Score | Confidence Level |
            |-------|------------|------------------|
            | 🔵  | 90-100 | Very High |
            | 🟢  | 70-90 | High |
            | 🟡  | 50-70 | Medium |
            | 🔴  | <50 | Low |
            """
            st.markdown(color_table)
    else:
        st.info("💡 Enter a protein sequence and click 'Predict Structure'")

# ========================
# ranaatom APP (Analysis)
# ========================
def ranaatom_app():
    st.title("🔍 PDB Analysis Toolkit")
    
    if not st.session_state.pdb_string:
        st.warning("No structure available. Please predict a structure first in the Predictor tab.")
        return
    
    st.write(f"Analyzing: {st.session_state.protein_name}")
    
    # Protein Properties
    st.subheader('🧪 Protein Properties')
    txt = st.session_state.get("input_sequence", "")
    if txt:
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

# ========================
# MAIN APP
# ========================
st.markdown(''' 
    <div style="text-align: center;">
        <h1 style="font-family: 'Dustosmo Roman', 'Times New Roman', serif; color:lightyellow;">
            ProtoAnalyzer
        </h1>
        <p>
            <em>A Streamlit <strong>Component</strong> for creating Speck molecular structures within Streamlit Web app.</em>
        </p>
    </div>
    ''', unsafe_allow_html=True) 
tab1, tab2 = st.tabs(["🧬 Predictor", "🔍 Analyzer"])
with tab1:
    emsfold_app()
with tab2:
    ranaatom_app()
