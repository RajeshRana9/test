import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Custom CSS to make sidebar slide from top to bottom
st.markdown("""
<style>
    [data-testid="stSidebar"] {
        top: 0;
        left: 0;
        width: 100%;
        height: auto;
        max-height: 0;
        overflow: hidden;
        transition: max-height 0.5s ease;
        position: relative;
        z-index: 999;
    }
    [data-testid="stSidebar"]:hover {
        max-height: 500px;
    }
    .sidebar-content {
        padding: 1rem;
        background-color: #f0f2f6;
        border-bottom: 1px solid #ddd;
    }
    .sidebar-title {
        text-align: center;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

st.set_page_config(layout='wide')

# Top sliding sidebar
with st.sidebar:
    st.markdown('<div class="sidebar-content">', unsafe_allow_html=True)
    st.markdown('<div class="sidebar-title">PRo-EDICTOR</div>', unsafe_allow_html=True)
    st.write("Built By Department Of BIOINFORMATICS")

    # Default protein sequence
    DEFAULT_SEQ = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
    txt = st.text_area('Input sequence', DEFAULT_SEQ, height=150)
    
    st.markdown('</div>', unsafe_allow_html=True)

# Rest of your existing code remains the same...
[keep all the existing code for prediction functions, session state, and main content display]
