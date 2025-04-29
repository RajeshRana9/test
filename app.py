import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
from collections import defaultdict

# Page setup
st.set_page_config(layout='wide')

# Session state initialization
if 'pdb_string' not in st.session_state:
    st.session_state.pdb_string = None
if 'b_value' not in st.session_state:
    st.session_state.b_value = None
if 'protein_name' not in st.session_state:
    st.session_state.protein_name = "Predicted Protein"
if 'sequence' not in st.session_state:
    st.session_state.sequence = ""

# Prediction function using ESMFold
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
# Predictor Tab
# ========================
def emsfold_app():
    st.sidebar.title('Input Options')

    DEFAULT_SEQ = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQASALSLSSSTSTWPEGLDATARAPPALVVTANIGQAGGSSSRQFRQRALGTSDSPVLFIHCPGAAGTAQGLEYRGRRVTTELVWEEVDSSPQPQGSESLPAQPPAQPAPQPEPQQAREPSPEVSCCGLWPRRPQRSQN"
    st.session_state.sequence = st.sidebar.text_area('Input sequence', DEFAULT_SEQ, height=275)

    custom_name = st.sidebar.text_input("Protein Name", st.session_state.protein_name)
    st.session_state.protein_name = custom_name

    if st.sidebar.button('⏳ Predict Structure'):
        with st.spinner('Predicting structure...'):
            update(st.session_state.sequence)
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
            |-------|-------------|------------------|
            | 🔵     | 90–100       | Very High        |
            | 🟢     | 70–90        | High             |
            | 🟡     | 50–70        | Medium           |
            | 🔴     | <50          | Low              |
            """
            st.markdown(color_table)

            st.subheader('🧪 Protein Properties')
            seq = st.session_state.sequence
            protein_seq = ProteinAnalysis(seq)
            hydrophobic = sum(seq.count(res) for res in 'AILMFWYV')
            data = {
                "Property": ["Length", "MW (Da)", "Hydrophobicity", "Net Charge", "Avg Confidence"],
                "Value": [
                    len(seq),
                    f"{protein_seq.molecular_weight()/1000:.1f} kDa",
                    f"{hydrophobic/len(seq)*100:.1f}%",
                    sum(seq.count(res) for res in 'KRH') - sum(seq.count(res) for res in 'DE'),
                    st.session_state.b_value
                ]
            }
            st.table(data)
    else:
        st.info("💡 Enter a protein sequence and click 'Predict Structure'")

# ========================
# Analyzer Tab
# ========================
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
        st.subheader("⚗️ Residue Property Distribution")
        seq = st.session_state.get("sequence", "").upper()

        categories = {
            "Hydrophobic": "AILMFWYV",
            "Polar": "STNQ",
            "Positively Charged": "KRH",
            "Negatively Charged": "DE",
            "Special Cases": "CGP"
        }

        prop_counts = {
            key: sum(seq.count(res) for res in residues)
            for key, residues in categories.items()
        }

        prop_df = pd.DataFrame.from_dict(prop_counts, orient='index', columns=['Count'])
        st.bar_chart(prop_df)

        st.subheader("Residue Type Distribution")
        res_counts = defaultdict(int)
        for line in atom_lines:
            res_name = line[17:20].strip()
            res_counts[res_name] += 1
        res_df = pd.DataFrame.from_dict(res_counts, orient='index', columns=['Count'])
        st.bar_chart(res_df)

# ========================
# Main App Layout
# ========================
st.markdown(''' 
    <div style="text-align: center;">
        <h1 style="font-family: 'Dustosmo Roman', 'Times New Roman', serif; color:purple;">
            PRo-EDICTOR
        </h1>
        <p>
            <em>Developed by </em> <strong>Rajesh Rana </strong> <em>  Department of BIOINFORMATICS, CPGS, OUAT</em>
        </p>
    </div>
    ''', unsafe_allow_html=True)

tab1, tab2 = st.tabs(["🧬 Predictor", "🔍 Analyzer"])
with tab1:
    emsfold_app()
with tab2:
    ranaatom_app()
