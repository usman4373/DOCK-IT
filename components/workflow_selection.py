import streamlit as st


def render_workflow_selection():
    """Render the workflow selection interface"""
    st.header("Select Docking Workflow")
    
    st.session_state.workflow_type = st.radio(
        "Choose Docking Workflow:",
        ["Blind Docking (whole protein docking)", "Site-specific Docking"],
        index=0
    )
    
    if st.button("Continue"):
        if "Blind" in st.session_state.workflow_type:
            st.session_state.workflow_type = "blind"
            st.session_state.stage = 1
        else:
            st.session_state.workflow_type = "site_specific"
            st.session_state.stage = 1.1  # Special stage for site-specific workflow selection
        st.rerun()

def render_site_specific_workflow_selection():
    """Render the site-specific workflow selection interface"""
    st.header("Site-specific Docking Workflow")
    
    st.session_state.site_specific_workflow = st.radio(
        "Choose Site-specific Docking Method:",
        ["Dock using reference ligand", "Dock by providing specific protein residues"],
        index=0
    )
    
    if st.button("Continue"):
        if "reference" in st.session_state.site_specific_workflow:
            st.session_state.site_specific_workflow = "reference_ligand"
            st.session_state.stage = 1.2  # Input paths for reference ligand workflow
        else:
            st.session_state.site_specific_workflow = "residues"
            st.session_state.stage = 1.3  # Input paths for residue-based workflow
        st.rerun()