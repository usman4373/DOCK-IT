import streamlit as st
import pandas as pd


def initialize_session_state():
    """Initialize all session state variables"""
    if 'stage' not in st.session_state:
        st.session_state.stage = 0 
    if 'workflow_type' not in st.session_state:
        st.session_state.workflow_type = "blind"  
    if 'site_specific_workflow' not in st.session_state:
        st.session_state.site_specific_workflow = "reference_ligand"
    if 'proteins_dir' not in st.session_state:
        st.session_state.proteins_dir = ""
    if 'ligands_input_type' not in st.session_state:
        st.session_state.ligands_input_type = "pdb directory"
    if 'ligands_folder' not in st.session_state:
        st.session_state.ligands_folder = ""
    if 'reference_ligands_dir' not in st.session_state:
        st.session_state.reference_ligands_dir = ""
    if 'smiles_file' not in st.session_state:
        st.session_state.smiles_file = None
    if 'smiles_column' not in st.session_state:
        st.session_state.smiles_column = ""
    if 'id_column' not in st.session_state:
        st.session_state.id_column = ""
    if 'ligands_output_dir' not in st.session_state:
        st.session_state.ligands_output_dir = ""
    if 'output_dir' not in st.session_state:
        st.session_state.output_dir = ""
    if 'scoring' not in st.session_state:
        st.session_state.scoring = "vina"
    if 'exhaustiveness' not in st.session_state:
        st.session_state.exhaustiveness = 64
    if 'cnn_scoring' not in st.session_state:
        st.session_state.cnn_scoring = "rescore"
    if 'num_modes' not in st.session_state:
        st.session_state.num_modes = 9
    if 'minimize' not in st.session_state:
        st.session_state.minimize = False
    if 'cpu' not in st.session_state:
        st.session_state.cpu = 6
    if 'gpu' not in st.session_state:
        st.session_state.gpu = False
    if 'device' not in st.session_state:
        st.session_state.device = 0
    if 'results_df' not in st.session_state:
        st.session_state.results_df = pd.DataFrame(columns=["protein", "ligand", "affinity", "pose"])
    if 'docking_completed' not in st.session_state:
        st.session_state.docking_completed = False
    if 'error_log' not in st.session_state:
        st.session_state.error_log = []
    if 'docking_log' not in st.session_state:
        st.session_state.docking_log = []
    if 'protein_ref_mapping' not in st.session_state:
        st.session_state.protein_ref_mapping = {}
    if 'extracted_ligands_dir' not in st.session_state:
        st.session_state.extracted_ligands_dir = ""
    if 'autobox_add' not in st.session_state:
        st.session_state.autobox_add = 4.0
    if 'residue_mapping' not in st.session_state:
        st.session_state.residue_mapping = {}
    if 'bbox_params' not in st.session_state:
        st.session_state.bbox_params = {}
    if 'padding' not in st.session_state:
        st.session_state.padding = 2.0
    if 'current_protein' not in st.session_state:
        st.session_state.current_protein = ""
    if 'padding_mode' not in st.session_state:
        st.session_state.padding_mode = "global"
    if 'global_padding' not in st.session_state:
        st.session_state.global_padding = 2.0
    if 'individual_padding' not in st.session_state:
        st.session_state.individual_padding = {}
    
    if 'results_df' not in st.session_state:
        st.session_state.results_df = pd.DataFrame(columns=[
            "protein", "ligand", "affinity", "pose", 
            "number of interactions", "types of interactions"
        ])
    
    # Constants
    st.session_state.EXCLUDED_IONS = {
        "HOH", "WAT", "Na+", "K+", "Cl-", "Ca2+", "Mg2+", "Zn2+", "Mn2+", 
        "Fe2+", "Fe3+", "Cu+", "Cu2+", "Co2+", "Ni2+", "Cd2+", "Hg2+", 
        "As3+", "Se4+", "SO4", "PO4", "NO3", "CO3"
    }
    st.session_state.standard_amino_acids = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }