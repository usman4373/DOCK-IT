import streamlit as st
import streamlit.components.v1 as components
from components.workflow_selection import render_workflow_selection, render_site_specific_workflow_selection
from components.blind_docking import render_blind_docking_inputs, render_blind_docking_parameters, render_blind_docking_execution
from components.site_specific_docking import render_site_specific_docking_inputs, render_site_specific_docking_parameters, render_site_specific_docking_execution
from components.residue_based_docking import render_residue_based_docking_inputs, render_residue_specification, render_bounding_box_adjustment, render_residue_based_docking_parameters, render_residue_based_docking_execution
from components.results_processing import render_results_processing
from utils.config import initialize_session_state


# Initialize session state
initialize_session_state()

# Set up page configuration
st.set_page_config(
    page_title="DOCK-IT",
    layout="wide"
)

# Create a two-column layout for the title
col1, col2 = st.columns([1, 8])

# Add image in the first column
with col1:
    # You can replace this with your actual image path
    st.image("icon/icon.png", width=150)

# Add app name in the second column
with col2:
    st.title("DOCK-IT")

# Add tagline beneath the title
st.markdown("GNINA-Powered Virtual Screening")

# Router based on stage
if st.session_state.stage == 0:
    render_workflow_selection()
elif st.session_state.stage == 1.1:
    render_site_specific_workflow_selection()
elif st.session_state.stage == 1 and st.session_state.workflow_type == "blind":
    render_blind_docking_inputs()
elif st.session_state.stage == 1.2 and st.session_state.workflow_type == "site_specific":
    render_site_specific_docking_inputs()
elif st.session_state.stage == 1.3 and st.session_state.workflow_type == "site_specific":
    render_residue_based_docking_inputs()
elif st.session_state.stage == 1.4 and st.session_state.workflow_type == "site_specific":
    render_residue_specification()
elif st.session_state.stage == 1.5 and st.session_state.workflow_type == "site_specific":
    render_bounding_box_adjustment()
elif st.session_state.stage == 2 and st.session_state.workflow_type == "blind":
    render_blind_docking_parameters()
elif st.session_state.stage == 2.1 and st.session_state.workflow_type == "site_specific":
    render_site_specific_docking_parameters()
elif st.session_state.stage == 2.2 and st.session_state.workflow_type == "site_specific":
    render_residue_based_docking_parameters()
elif st.session_state.stage == 3 and st.session_state.workflow_type == "blind":
    render_blind_docking_execution()
elif st.session_state.stage == 3.1 and st.session_state.workflow_type == "site_specific":
    render_site_specific_docking_execution()
elif st.session_state.stage == 3.2 and st.session_state.workflow_type == "site_specific":
    render_residue_based_docking_execution()
elif st.session_state.docking_completed:
    render_results_processing()

# Instructions sidebar
st.sidebar.header("Instructions")
st.sidebar.info("""
1. Select the docking workflow (Blind or Site-specific)
2. Provide input paths for proteins and ligands
3. Select docking parameters
4. Run the screening
5. View results

Note: This app requires GNINA to be installed and available in your PATH.
""")