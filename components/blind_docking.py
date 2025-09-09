import streamlit as st
import os
import pandas as pd
import glob
from utils.file_processing import detect_smiles_columns, detect_id_columns
from utils.smiles_processing import smiles_to_pdb
from utils.docking_functions import run_gnina_docking
from utils.file_processing import update_results_xlsx, process_docking_results


def render_blind_docking_inputs():
    """Render the blind docking input interface"""
    st.header("Step 1: Input Parameters - Blind Docking")
    
    # Proteins input
    st.session_state.proteins_dir = st.text_input(
        "Input Proteins Directory Path:",
        value=st.session_state.proteins_dir,
        help="Path to the directory containing protein PDB files"
    )
    
    # Ligands input type
    st.session_state.ligands_input_type = st.radio(
        "Ligands Input Type:",
        ["pdb directory", "SMILES CSV"],
        index=0 if st.session_state.ligands_input_type == "pdb directory" else 1
    )
    
    if st.session_state.ligands_input_type == "pdb directory":
        st.session_state.ligands_folder = st.text_input(
            "Ligands Directory Path:",
            value=st.session_state.ligands_folder,
            help="Path to the directory containing ligand PDB files"
        )
    else:
        uploaded_file = st.file_uploader("Upload SMILES CSV File:", type=["csv"])
        if uploaded_file is not None:
            # Read the CSV file
            df = pd.read_csv(uploaded_file)
            st.session_state.smiles_file = df
    
            # Display the dataframe
            st.write("Preview of uploaded file:")
            st.dataframe(df.head())
    
            # Ask user for SMILES column
            st.session_state.smiles_column = st.selectbox(
                "Select SMILES Column:",
                df.columns.tolist(),
                index=0
            )
    
            # Ask user for ID column
            st.session_state.id_column = st.selectbox(
                "Select Compound ID/Name Column:",
                df.columns.tolist(),
                index=0
            )
        
        st.session_state.ligands_output_dir = st.text_input(
            "Output Directory for Converted Ligands (PDB format):",
            value=st.session_state.ligands_output_dir,
            help="Path where converted ligand PDB files will be saved"
    )
    
    # Output directory
    st.session_state.output_dir = st.text_input(
        "Docking Output Directory Path:",
        value=st.session_state.output_dir,
        help="Path where docking results will be saved"
    )
    
    # Proceed button
    if st.button("Proceed to Parameters"):
        # Validate inputs
        errors = []
        if not os.path.isdir(st.session_state.proteins_dir):
            errors.append("Proteins directory does not exist.")
        
        if st.session_state.ligands_input_type == "pdb directory":
            if not os.path.isdir(st.session_state.ligands_folder):
                errors.append("Ligands directory does not exist.")
        else:
            if st.session_state.smiles_file is None:
                errors.append("Please upload a SMILES CSV file.")
            if not st.session_state.ligands_output_dir:
                errors.append("Please specify an output directory for converted ligands.")
            else:
                os.makedirs(st.session_state.ligands_output_dir, exist_ok=True)
        
        if not st.session_state.output_dir:
            errors.append("Please specify an output directory for docking results.")
        else:
            os.makedirs(st.session_state.output_dir, exist_ok=True)
        
        if errors:
            for error in errors:
                st.error(error)
        else:
            st.session_state.stage = 2
            st.rerun()

def render_blind_docking_parameters():
    """Render the blind docking parameters interface"""
    st.header("Step 2: Docking Parameters - Blind Docking")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.session_state.scoring = st.selectbox(
            "Scoring Function:",
            ["vina", "ad4_scoring", "dkoes_fast", "dkoes_scoring", "dkoes_scoring_old", "vinardo"],
            index=0
        )
        
        st.session_state.exhaustiveness = st.number_input(
            "Exhaustiveness:",
            min_value=4,
            max_value=200,
            value=64,
            step=4
        )
        
        st.session_state.cnn_scoring = st.selectbox(
            "CNN Scoring:",
            ["none", "rescore", "refinement", "metrorescore", "metrorefine", "all"],
            index=0
        )
        
        st.session_state.num_modes = st.number_input(
            "Number of Models:",
            min_value=3,
            max_value=20,
            value=9,
            step=1
        )
        
    with col2:
        st.session_state.cpu = st.number_input(
            "CPU Cores:",
            min_value=2,
            max_value=128,
            value=6,
            step=2
        )
        
        st.session_state.minimize = st.checkbox(
            "Minimize After Docking",
            value=False,
            help="When enabled, GNINA will generate only one minimized structure instead of multiple poses."
        )
        
        # Show note immediately if minimize is checked
        if st.session_state.minimize:
            st.info("Note: When minimization is enabled, GNINA will generate only one minimized structure instead of multiple poses.")
        
        st.session_state.gpu = st.checkbox(
            "Use GPU",
            value=False
        )
        
        if st.session_state.gpu:
            st.session_state.device = st.number_input(
                "GPU Device:",
                min_value=0,
                max_value=7,
                value=0
            )
            # Add note about GPU devices
            st.info("Note: GPU device numbers typically start from 0. If you have multiple GPUs, device 0 is usually the first GPU, device 1 is the second, etc.")
    
    col3, col4 = st.columns(2)
    with col3:
        if st.button("Back to Inputs"):
            st.session_state.stage = 1
            st.rerun()
    
    with col4:
        if st.button("Run Screening"):
            st.session_state.stage = 3
            st.rerun()

def render_blind_docking_execution():
    """Render the blind docking execution interface"""
    st.header("Step 3: Running Virtual Screening - Blind Docking")
    
    # Get protein files
    protein_files = glob.glob(os.path.join(st.session_state.proteins_dir, "*.pdb"))
    if not protein_files:
        st.error("No PDB files found in the proteins directory.")
        st.stop()
    
    # Get ligand files or convert SMILES
    if st.session_state.ligands_input_type == "pdb directory":
        ligand_files = glob.glob(os.path.join(st.session_state.ligands_folder, "*.pdb"))
        if not ligand_files:
            st.error("No PDB files found in the ligands directory.")
            st.stop()
        ligand_names = [os.path.basename(f).replace('.pdb', '') for f in ligand_files]
    else:
        # Convert SMILES to PDB
        st.info("Converting SMILES to PDB format...")
        df = st.session_state.smiles_file
        ligand_files = []
        ligand_names = []
        
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        for i, row in df.iterrows():
            smiles = row[st.session_state.smiles_column]
            if pd.isna(smiles):
                continue
                
            if st.session_state.id_column and st.session_state.id_column in row:
                ligand_name = str(row[st.session_state.id_column])
            else:
                ligand_name = f"compound_{i}"
                
            output_path = os.path.join(st.session_state.ligands_output_dir, f"{ligand_name}.pdb")
            
            status_text.text(f"Converting {ligand_name}...")
            if smiles_to_pdb(smiles, output_path, ligand_name):
                ligand_files.append(output_path)
                ligand_names.append(ligand_name)
            
            progress_bar.progress((i + 1) / len(df))
        
        status_text.text("Conversion completed!")
        progress_bar.empty()
        
        if not ligand_files:
            st.error("No ligands were successfully converted.")
            st.stop()
    
    # Create output directory structure
    for protein_file in protein_files:
        protein_name = os.path.basename(protein_file).replace('.pdb', '')
        protein_output_dir = os.path.join(st.session_state.output_dir, protein_name)
        os.makedirs(protein_output_dir, exist_ok=True)
    
    log_dir = os.path.join(st.session_state.output_dir, "log")
    os.makedirs(log_dir, exist_ok=True)
    
    # Initialize results dataframe
    st.session_state.results_df = pd.DataFrame(columns=["protein", "ligand", "affinity", "pose"])
    st.session_state.error_log = []
    st.session_state.docking_log = []
    
    # Run docking
    st.info("Starting docking process...")
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    total_operations = len(ligand_files) * len(protein_files)
    completed_operations = 0
    
    # Create a scrollable log box
    log_container = st.container()
    with log_container:
        st.markdown("**Docking Progress Log**")
        log_display = st.empty()
    
    for lig_idx, (ligand_file, ligand_name) in enumerate(zip(ligand_files, ligand_names)):
        for prot_idx, protein_file in enumerate(protein_files):
            protein_name = os.path.basename(protein_file).replace('.pdb', '')
            
            # Add to log
            log_entry = f"Docking {ligand_name} with {protein_name}"
            st.session_state.docking_log.append(log_entry)
            
            # Update log display (show last 20 entries)
            with log_display:
                st.text_area("", "\n".join(st.session_state.docking_log[-20:]), height=200)
            
            # Create output paths
            protein_output_dir = os.path.join(st.session_state.output_dir, protein_name)
            output_name = f"{protein_name}_{ligand_name}.pdbqt"
            output_path = os.path.join(protein_output_dir, output_name)
            log_path = os.path.join(log_dir, f"{protein_name}_{ligand_name}.log")
            
            # Run docking
            affinity, pose = run_gnina_docking(
                protein_file, 
                ligand_file, 
                output_path, 
                log_path,
                st.session_state.scoring,
                st.session_state.exhaustiveness,
                st.session_state.cnn_scoring,
                st.session_state.num_modes,
                st.session_state.minimize,
                st.session_state.cpu,
                st.session_state.gpu,
                st.session_state.device
            )
            
            # Update results if docking was successful
            if affinity is not None:
                update_results_xlsx(protein_name, ligand_name, affinity, pose, st.session_state.output_dir)
                log_entry = f"✓ Completed: {ligand_name} with {protein_name} ➡ Affinity: {affinity:.2f}, Pose: {pose}"
            else:
                log_entry = f"✗ Failed: {ligand_name} with {protein_name}"
            
            st.session_state.docking_log.append(log_entry)
            
            # Update log display
            with log_display:
                st.text_area("", "\n".join(st.session_state.docking_log[-20:]), height=200)
            
            # Update progress
            completed_operations += 1
            progress_bar.progress(completed_operations / total_operations)
    
    # Final status
    progress_bar.empty()
    status_text.text("Docking completed!")
    
    # Save error log
    if st.session_state.error_log:
        error_log_path = os.path.join(st.session_state.output_dir, "error_log.txt")
        with open(error_log_path, 'w') as f:
            for error in st.session_state.error_log:
                f.write(error + "\n")
        
        st.warning(f"Some errors occurred during docking. See {error_log_path} for details.")
    
    # Post-processing
    st.info("Post-processing docking results (splitting poses and creating complexes)...")
    if process_docking_results(st.session_state.output_dir, st.session_state.minimize):
        st.success("Post-processing completed successfully!")
    else:
        st.warning("Post-processing completed with some errors. Check error log for details.")
    
    # Show top 3 results for each protein
    st.subheader("Top 3 Docking Results per Protein")
    if not st.session_state.results_df.empty:
        # Group by protein and get top 3 affinities for each
        top_results = pd.DataFrame()
        for protein in st.session_state.results_df['protein'].unique():
            protein_results = st.session_state.results_df[st.session_state.results_df['protein'] == protein]
            protein_top3 = protein_results.nsmallest(3, 'affinity')
            top_results = pd.concat([top_results, protein_top3])
        
        st.dataframe(top_results)
    else:
        st.warning("No successful docking results were obtained.")
    
    st.session_state.docking_completed = True