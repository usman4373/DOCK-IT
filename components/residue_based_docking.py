import streamlit as st
import os
import pandas as pd
import glob
import numpy as np
from utils.pdb_utils import parse_residue_string, calculate_bounding_box, get_protein_structure
from utils.smiles_processing import smiles_to_pdb
from utils.docking_functions import run_gnina_custom_box
from utils.file_processing import update_results_xlsx, process_docking_results
from components.visualization import create_3dmol_view
import streamlit.components.v1 as components


def render_residue_based_docking_inputs():
    """Render the residue-based docking input interface"""
    st.header("Step 1: Input Parameters - Residue-based Docking")
    
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
    if st.button("Specify Residue Positions"):
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
            st.session_state.stage = 1.4
            st.rerun()

def render_residue_specification():
    """Render the residue specification interface"""
    st.header("Step 2: Specify Residue Positions")
    
    st.info("""
    **Instructions for specifying residue positions:**
    - For individual residues: Use the format `45+67+23`
    - For residue ranges: Use the format `79-100` or `40-45+23+32-39`
    - Example: `45+67+23` or `79-100` or `40-45+23+32-39`
    """)
    
    # Get protein files
    protein_files = glob.glob(os.path.join(st.session_state.proteins_dir, "*.pdb"))
    if not protein_files:
        st.error("No PDB files found in the proteins directory.")
        st.stop()
    
    # Initialize residue mapping if not exists
    if not st.session_state.residue_mapping:
        for protein_file in protein_files:
            st.session_state.residue_mapping[protein_file] = ""
    
    # Display input fields for each protein
    for protein_file in protein_files:
        protein_name = os.path.basename(protein_file).replace('.pdb', '')
        
        st.subheader(f"Protein: {protein_name}")
        st.session_state.residue_mapping[protein_file] = st.text_input(
            f"Residue positions for {protein_name}:",
            value=st.session_state.residue_mapping[protein_file],
            key=f"residues_{protein_file}"
        )
    
    if st.button("Adjust Bounding Box Size"):
        # Validate residue inputs
        valid = True
        for protein_file, res_str in st.session_state.residue_mapping.items():
            if not res_str.strip():
                st.error(f"Please specify residue positions for {os.path.basename(protein_file)}")
                valid = False
                break
            
            # Try to parse the residue string
            residues = parse_residue_string(res_str)
            if not residues:
                st.error(f"Invalid residue format for {os.path.basename(protein_file)}")
                valid = False
                break
        
        if valid:
            st.session_state.stage = 1.5
            st.rerun()

def render_bounding_box_adjustment():
    """Render the bounding box adjustment interface"""
    st.header("Step 3: Adjust Bounding Box Size")
    
    st.info("""
    **Note:** The bounding box defines the search space for docking. 
    You can adjust the padding around the selected residues to increase or decrease the box size.
    Proceed to parameters if you want to keep the default size.
    """)
    
    # Padding mode selection
    padding_options = ["Set padding globally for all proteins", "Set padding individually for each protein"]
    st.session_state.padding_mode = st.radio(
        "Padding Mode:",
        options=padding_options,
        index=0 if st.session_state.get('padding_mode', 'global') == "global" else 1
    )
    
    # Simplify padding mode to just 'global' or 'individual'
    is_global_mode = st.session_state.padding_mode == padding_options[0]
    st.session_state.padding_mode_simple = "global" if is_global_mode else "individual"
    
    # Get protein files
    protein_files = list(st.session_state.residue_mapping.keys())
    
    # Initialize individual padding if not exists
    for protein_file in protein_files:
        if protein_file not in st.session_state.individual_padding:
            st.session_state.individual_padding[protein_file] = 2.0
    
    # Apply global padding to all proteins if in global mode
    if is_global_mode:
        # If we just switched to global mode, set all to current global padding
        for protein_file in protein_files:
            st.session_state.individual_padding[protein_file] = st.session_state.global_padding
    
    # Create sidebar for protein selection
    st.sidebar.header("Select Protein")
    selected_protein = st.sidebar.selectbox(
        "Choose a protein to view:",
        [os.path.basename(f).replace('.pdb', '') for f in protein_files],
        index=0
    )
    
    # Find the selected protein file
    selected_file = None
    for protein_file in protein_files:
        if os.path.basename(protein_file).replace('.pdb', '') == selected_protein:
            selected_file = protein_file
            break
    
    if selected_file:
        st.session_state.current_protein = selected_file
        
        # Get the padding value for this protein
        padding_value = st.session_state.individual_padding[selected_file]
        
        # Parse residues
        residues = parse_residue_string(st.session_state.residue_mapping[selected_file])
        
        # Calculate bounding box
        center, size = calculate_bounding_box(selected_file, residues, padding_value)
        
        if center is not None and size is not None:
            # Store bounding box parameters
            st.session_state.bbox_params[selected_file] = {
                'center': center,
                'size': size,
                'residues': residues,
                'padding': padding_value
            }
            
            # Get protein structure
            pdb_data = get_protein_structure(selected_file)
            
            # Create 3D view
            html_content = create_3dmol_view(pdb_data, residues, center, size, padding_value)

            # Display the view and bounding box info side by side
            col1, col2 = st.columns([2, 1])
            
            with col1:
                st.subheader(f"3D View: {selected_protein}")
                components.html(html_content, height=500)
            
            with col2:
                st.subheader("Bounding Box Info")
                # Show bounding box dimensions
                st.write(f"**Bounding Box Dimensions:**")
                st.write(f"- X: {size[0]:.2f} Å")
                st.write(f"- Y: {size[1]:.2f} Å")
                st.write(f"- Z: {size[2]:.2f} Å")
                st.write(f"- Center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
                
                # Show padding information
                st.write(f"- Padding: {padding_value:.1f} Å")
                
                # Show the appropriate slider based on mode
                if is_global_mode:
                    new_padding = st.slider(
                        "Adjust Global Padding (Å):",
                        min_value=0.0,
                        max_value=50.0,
                        value=st.session_state.global_padding,
                        step=0.5,
                        key="global_padding_right_slider"
                    )
                    # Update the global padding value if changed
                    if new_padding != st.session_state.global_padding:
                        st.session_state.global_padding = new_padding
                        # Apply to all proteins
                        for protein_file in protein_files:
                            st.session_state.individual_padding[protein_file] = new_padding
                        st.rerun()
                else:
                    new_padding = st.slider(
                        f"Adjust Padding for {selected_protein} (Å):",
                        min_value=0.0,
                        max_value=50.0,
                        value=padding_value,
                        step=0.5,
                        key=f"individual_padding_right_{selected_file}"
                    )
                    # Update the individual padding value if changed
                    if new_padding != padding_value:
                        st.session_state.individual_padding[selected_file] = new_padding
                        st.rerun()
        else:
            st.error(f"Failed to calculate bounding box for {selected_protein}")
    
    # Move the button to the bottom
    if st.button("Proceed to Parameters"):
        # Calculate bounding box for any proteins that haven't been processed
        for protein_file in protein_files:
            if protein_file not in st.session_state.bbox_params:
                protein_name = os.path.basename(protein_file).replace('.pdb', '')
                
                # Get the padding value for this protein
                padding_value = st.session_state.individual_padding[protein_file]
                
                # Parse residues
                residues = parse_residue_string(st.session_state.residue_mapping[protein_file])

                # Calculate bounding box with current padding
                center, size = calculate_bounding_box(protein_file, residues, padding_value)

                if center is not None and size is not None:
                    # Store bounding box parameters
                    st.session_state.bbox_params[protein_file] = {
                        'center': center,
                        'size': size,
                        'residues': residues,
                        'padding': padding_value
                    }
                    st.success(f"Calculated bounding box for {protein_name} with padding {padding_value} Å")
                else:
                    st.error(f"Failed to calculate bounding box for {protein_name}")

        # Proceed to parameters
        st.session_state.stage = 2.2
        st.rerun()

def render_residue_based_docking_parameters():
    """Render the residue-based docking parameters interface"""
    st.header("Step 4: Docking Parameters - Residue-based Docking")
    
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
    if st.session_state.gpu:
        st.info("Note: GPU device numbers typically start from 0. If you have multiple GPUs, device 0 is usually the first GPU, device 1 is the second, etc.")
    
    col3, col4 = st.columns(2)
    with col3:
        if st.button("Back to Bounding Box Adjustment"):
            st.session_state.stage = 1.5
            st.rerun()
    
    with col4:
        if st.button("Run Screening"):
            st.session_state.stage = 3.2
            st.rerun()

def render_residue_based_docking_execution():
    """Render the residue-based docking execution interface"""
    st.header("Step 5: Running Virtual Screening - Residue-based Docking")
    
    # Get protein files
    protein_files = list(st.session_state.residue_mapping.keys())
    if not protein_files:
        st.error("No protein files found.")
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
            
            # Get bounding box parameters
            if protein_file not in st.session_state.bbox_params:
                st.session_state.error_log.append(f"No bounding box parameters found for {protein_name}")
                continue
                
            bbox_params = st.session_state.bbox_params[protein_file]
            center = bbox_params['center']
            size = bbox_params['size']
            
            # Add to log
            log_entry = f"Docking {ligand_name} with {protein_name} using custom box"
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
            affinity, pose = run_gnina_custom_box(
                protein_file, 
                ligand_file, 
                center, 
                size,
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