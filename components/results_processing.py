import streamlit as st
import pandas as pd
import os


def render_results_processing():
    """Render the results processing interface"""
    st.header("Docking Results")
    
    if not st.session_state.results_df.empty:
        # Show top 3 results for each protein
        st.subheader("Top 3 Docking Results per Protein")
        
        # Group by protein and get top 3 affinities for each
        top_results = pd.DataFrame()
        for protein in st.session_state.results_df['protein'].unique():
            protein_results = st.session_state.results_df[st.session_state.results_df['protein'] == protein]
            protein_top3 = protein_results.nsmallest(3, 'affinity')
            top_results = pd.concat([top_results, protein_top3])
        
        st.dataframe(top_results)
        
        # Download button for full results
        csv = st.session_state.results_df.to_csv(index=False)
        st.download_button(
            label="Download Full Results as CSV",
            data=csv,
            file_name="docking_results.csv",
            mime="text/csv"
        )
    else:
        st.warning("No successful docking results were obtained.")
    
    # Show error log if there are errors
    if st.session_state.error_log:
        st.subheader("Error Log")
        with st.expander("View Errors"):
            for error in st.session_state.error_log:
                st.error(error)
    
    # Show docking log
    st.subheader("Docking Log")
    with st.expander("View Detailed Log"):
        for log_entry in st.session_state.docking_log:
            st.text(log_entry)
    
    # Option to start a new screening
    if st.button("Start New Screening"):
        # Reset session state
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        st.session_state.stage = 0
        st.rerun()