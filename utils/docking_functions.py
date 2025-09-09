import subprocess
import re
import time
import streamlit as st


def run_gnina_docking(protein_path, ligand_path, output_path, log_path, scoring, exhaustiveness, 
                     cnn_scoring, num_modes, minimize, cpu, gpu, device):
    """Run GNINA docking with the specified parameters"""
    try:
        # Build the command
        cmd = [
            "gnina",
            "-r", protein_path,
            "-l", ligand_path,
            "--autobox_ligand", protein_path,
            "-o", output_path,
            "--log", log_path,
            f"--scoring={scoring}",
            f"--exhaustiveness={exhaustiveness}",
            f"--cnn_scoring={cnn_scoring}",
            f"--num_modes={num_modes}",
            f"--cpu={cpu}"
        ]
        
        if minimize:
            cmd.append("--minimize")
            
        if gpu:
            cmd.append(f"--device={device}")
        else:
            cmd.append("--no_gpu")
            
        # Run the command
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            error_msg = f"GNINA error for {ligand_path} with {protein_path}: {result.stderr}"
            st.session_state.error_log.append(error_msg)
            return None, None
        
        # Wait a moment to ensure the log file is completely written
        time.sleep(0.5)
        
        # Parse the output to get affinity scores and pose number
        with open(log_path, 'r') as f:
            log_content = f.read()
        
        # Extract affinity scores and pose number
        affinity = None
        pose = 1  # Default to pose 1
        
        if minimize:
            # Format for minimization: "Affinity: -0.00377  0.15176 (kcal/mol)"
            affinity_pattern = r"Affinity:\s+([-\d.]+)"
            match = re.search(affinity_pattern, log_content)
            if match:
                affinity = float(match.group(1))
        else:
            # Try to find multiple poses first
            affinity_pattern = r"\s+(\d+)\s+([-\d.]+)\s+"
            matches = re.findall(affinity_pattern, log_content)
            
            if matches:
                # Format for regular docking: table with multiple modes
                # Convert to list of tuples (pose, affinity)
                pose_affinities = [(int(pose), float(affinity)) for pose, affinity in matches]
                # Find the best (lowest) affinity and its pose number
                best_pose, best_affinity = min(pose_affinities, key=lambda x: x[1])
                affinity = best_affinity
                pose = best_pose
            else:
                # Try to find single pose format if multiple poses not found
                affinity_pattern = r"Affinity:\s+([-\d.]+)"
                match = re.search(affinity_pattern, log_content)
                if match:
                    affinity = float(match.group(1))
        
        if affinity is not None:
            return affinity, pose
        else:
            error_msg = f"No affinity scores found in log for {ligand_path} with {protein_path}"
            st.session_state.error_log.append(error_msg)
            return None, None
            
    except Exception as e:
        error_msg = f"Exception during docking of {ligand_path} with {protein_path}: {str(e)}"
        st.session_state.error_log.append(error_msg)
        return None, None

def run_gnina_site_specific_docking(protein_path, ligand_path, ref_ligand_path, output_path, log_path, scoring, exhaustiveness, 
                                   cnn_scoring, num_modes, minimize, cpu, gpu, device, autobox_add):
    """Run GNINA docking with site-specific parameters"""
    try:
        # Build the command
        cmd = [
            "gnina",
            "-r", protein_path,
            "-l", ligand_path,
            "--autobox_ligand", ref_ligand_path,
            "--autobox_add", str(autobox_add),
            "-o", output_path,
            "--log", log_path,
            f"--scoring={scoring}",
            f"--exhaustiveness={exhaustiveness}",
            f"--cnn_scoring={cnn_scoring}",
            f"--num_modes={num_modes}",
            f"--cpu={cpu}"
        ]
        
        if minimize:
            cmd.append("--minimize")
            
        if gpu:
            cmd.append(f"--device={device}")
        else:
            cmd.append("--no_gpu")
            
        # Run the command
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            error_msg = f"GNINA error for {ligand_path} with {protein_path}: {result.stderr}"
            st.session_state.error_log.append(error_msg)
            return None, None
        
        # Wait a moment to ensure the log file is completely written
        time.sleep(0.5)
        
        # Parse the output to get affinity scores and pose number
        with open(log_path, 'r') as f:
            log_content = f.read()
        
        # Extract affinity scores and pose number
        affinity = None
        pose = 1  # Default to pose 1
        
        if minimize:
            # Format for minimization: "Affinity: -0.00377  0.15176 (kcal/mol)"
            affinity_pattern = r"Affinity:\s+([-\d.]+)"
            match = re.search(affinity_pattern, log_content)
            if match:
                affinity = float(match.group(1))
        else:
            # Try to find multiple poses first
            affinity_pattern = r"\s+(\d+)\s+([-\d.]+)\s+"
            matches = re.findall(affinity_pattern, log_content)
            
            if matches:
                # Format for regular docking: table with multiple modes
                # Convert to list of tuples (pose, affinity)
                pose_affinities = [(int(pose), float(affinity)) for pose, affinity in matches]
                # Find the best (lowest) affinity and its pose number
                best_pose, best_affinity = min(pose_affinities, key=lambda x: x[1])
                affinity = best_affinity
                pose = best_pose
            else:
                # Try to find single pose format if multiple poses not found
                affinity_pattern = r"Affinity:\s+([-\d.]+)"
                match = re.search(affinity_pattern, log_content)
                if match:
                    affinity = float(match.group(1))
        
        if affinity is not None:
            return affinity, pose
        else:
            error_msg = f"No affinity scores found in log for {ligand_path} with {protein_path}"
            st.session_state.error_log.append(error_msg)
            return None, None
            
    except Exception as e:
        error_msg = f"Exception during docking of {ligand_path} with {protein_path}: {str(e)}"
        st.session_state.error_log.append(error_msg)
        return None, None

def run_gnina_custom_box(protein_path, ligand_path, center, size, output_path, log_path, scoring, exhaustiveness, 
                        cnn_scoring, num_modes, minimize, cpu, gpu, device):
    """Run GNINA docking with custom box parameters"""
    try:
        # Build the command
        cmd = [
            "gnina",
            "-r", protein_path,
            "-l", ligand_path,
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "-o", output_path,
            "--log", log_path,
            f"--scoring={scoring}",
            f"--exhaustiveness={exhaustiveness}",
            f"--cnn_scoring={cnn_scoring}",
            f"--num_modes={num_modes}",
            f"--cpu={cpu}"
        ]
        
        if minimize:
            cmd.append("--minimize")
            
        if gpu:
            cmd.append(f"--device={device}")
        else:
            cmd.append("--no_gpu")
            
        # Run the command
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            error_msg = f"GNINA error for {ligand_path} with {protein_path}: {result.stderr}"
            st.session_state.error_log.append(error_msg)
            return None, None
        
        # Wait a moment to ensure the log file is completely written
        time.sleep(0.5)
        
        # Parse the output to get affinity scores and pose number
        with open(log_path, 'r') as f:
            log_content = f.read()
        
        # Extract affinity scores and pose number
        affinity = None
        pose = 1  # Default to pose 1
        
        if minimize:
            # Format for minimization: "Affinity: -0.00377  0.15176 (kcal/mol)"
            affinity_pattern = r"Affinity:\s+([-\d.]+)"
            match = re.search(affinity_pattern, log_content)
            if match:
                affinity = float(match.group(1))
        else:
            # Try to find multiple poses first
            affinity_pattern = r"\s+(\d+)\s+([-\d.]+)\s+"
            matches = re.findall(affinity_pattern, log_content)
            
            if matches:
                # Format for regular docking: table with multiple modes
                # Convert to list of tuples (pose, affinity)
                pose_affinities = [(int(pose), float(affinity)) for pose, affinity in matches]
                # Find the best (lowest) affinity and its pose number
                best_pose, best_affinity = min(pose_affinities, key=lambda x: x[1])
                affinity = best_affinity
                pose = best_pose
            else:
                # Try to find single pose format if multiple poses not found
                affinity_pattern = r"Affinity:\s+([-\d.]+)"
                match = re.search(affinity_pattern, log_content)
                if match:
                    affinity = float(match.group(1))
        
        if affinity is not None:
            return affinity, pose
        else:
            error_msg = f"No affinity scores found in log for {ligand_path} with {protein_path}"
            st.session_state.error_log.append(error_msg)
            return None, None
            
    except Exception as e:
        error_msg = f"Exception during docking of {ligand_path} with {protein_path}: {str(e)}"
        st.session_state.error_log.append(error_msg)
        return None, None