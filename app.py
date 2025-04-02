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
    txt = st.session_state.get("input_sequence", "")
    if txt:
    hydrophobic = sum(txt.count(res) for res in 'AILMFWYV')
    protein_seq = ProteinAnalysis(txt)  # Ensure txt is not empty before analysis
    
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
else:
    data = {
        "Property": ["Length", "MW (Da)", "Hydrophobicity", "Net Charge", "Avg Confidence"],
        "Value": ["N/A", "N/A", "N/A", "N/A", st.session_state.b_value]
    }

st.table(data)


    # Advanced visualization controls
    st.sidebar.title('Analysis Settings')
    style = st.sidebar.selectbox("Style", 
                               ["cartoon", "sphere", "stick", "surface"],
                               index=0)
    color_scheme = st.sidebar.selectbox("Color Scheme",
                                      ["spectrum", "chain", "residue"],
                                      index=0)
    show_labels = st.sidebar.checkbox("Show Atom Labels", False)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader(f"🔬 {st.session_state.protein_name} Structure")
        
        # Display basic info
        atom_lines = [line for line in st.session_state.pdb_string.split('\n') if line.startswith("ATOM")]
        num_atoms = len(atom_lines)
        st.caption(f"{num_atoms:,} atoms | {len(set(line[21] for line in atom_lines))} chains")
        
        # 3D Visualization
        view = py3Dmol.view(width=600, height=400)
        view.addModel(st.session_state.pdb_string, "pdb")
        
        if style == "cartoon":
            view.setStyle({'cartoon': {'color': color_scheme}})
        elif style == "sphere":
            view.setStyle({'sphere': {'colorscheme': color_scheme}})
        elif style == "stick":
            view.setStyle({'stick': {'colorscheme': color_scheme}})
        elif style == "surface":
            view.addSurface(py3Dmol.VDW, {'opacity':0.7, 'color':'white'})
        
        if show_labels:
            view.addResLabels()
            
        view.zoomTo()
        showmol(view, height=400)
        
        st.download_button(
            label="⬇️ Download PDB",
            data=st.session_state.pdb_string,
            file_name=f'{st.session_state.protein_name.replace(" ", "_")}.pdb',
            mime='text/plain'
        )
    
    with col2:
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
        st.subheader("Structure Analysis")
        
        # Secondary structure analysis
        helix_count = sum(1 for line in st.session_state.pdb_string.split('\n') if line.startswith("HELIX"))
        sheet_count = sum(1 for line in st.session_state.pdb_string.split('\n') if line.startswith("SHEET"))
        coil_count = sum(1 for line in st.session_state.pdb_string.split('\n') if line.startswith("coil"))
        ss_data = pd.DataFrame({
            'Type': ['Helix', 'Sheet', 'Coil'],
            'Count': [helix_count, sheet_count, num_atoms - (helix_count + sheet_count)]
        })
        
        st.subheader("Secondary Structure")
        st.bar_chart(ss_data.set_index('Type'))
        
        # Residue type distribution
        res_counts = defaultdict(int)
        for line in atom_lines:
            res_name = line[17:20].strip()
            res_counts[res_name] += 1
        
        res_df = pd.DataFrame.from_dict(res_counts, orient='index', columns=['Count'])
        st.subheader("Residue Distribution")
        st.bar_chart(res_df)

        
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
