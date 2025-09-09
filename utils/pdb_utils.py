import os
import tempfile
import subprocess
from Bio import PDB
from Bio.PDB import is_aa
from Bio.PDB import PDBParser
from collections import defaultdict
from pymol import cmd
import numpy as np
import streamlit as st


def extract_protein_sequence(pdb_path):
    """Extract protein sequence from a PDB file"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)
        
        # Get all amino acid residues
        residues = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in st.session_state.standard_amino_acids:
                        residues.append(residue.get_resname())
        
        return "".join(residues)
    except Exception as e:
        st.session_state.error_log.append(f"Error extracting sequence from {pdb_path}: {str(e)}")
        return None

def extract_ligand_from_complex(complex_path, output_path):
    """Extract ligand from a protein-ligand complex using PyMOL"""
    try:
        # Create a PyMOL script to extract the ligand
        pymol_script = f"""
load {complex_path}
select protein, polymer
remove protein
save {output_path}
quit
"""
        # Write the script to a temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pml', delete=False) as f:
            f.write(pymol_script)
            script_path = f.name
        
        # Run PyMOL in headless mode to execute the script
        result = subprocess.run(['pymol', '-c', script_path], 
                              capture_output=True, text=True)
        
        # Clean up the temporary script
        os.unlink(script_path)
        
        if result.returncode != 0:
            st.session_state.error_log.append(f"PyMOL error extracting ligand from {complex_path}: {result.stderr}")
            return False
        
        # Check if the output file was created and has content
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            return True
        else:
            st.session_state.error_log.append(f"Failed to extract ligand from {complex_path}")
            return False
            
    except Exception as e:
        st.session_state.error_log.append(f"Exception extracting ligand from {complex_path}: {str(e)}")
        return False

def classify_residues(pdb_file):
    """Classify residues in a PDB file as protein or ligand"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    
    protein_chains = defaultdict(set)
    ligand_residues = set()
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.resname.strip()
                is_hetatm = residue.id[0].strip() != ''
                
                if is_aa(residue):
                    protein_chains[chain.id].add(residue.id[1])
                elif resname in {"UNK", "UNL"}:
                    ligand_residues.add((chain.id, resname, residue.id[1]))
                elif is_hetatm and resname not in st.session_state.EXCLUDED_IONS:
                    ligand_residues.add((chain.id, resname, residue.id[1]))
    
    return protein_chains, ligand_residues

def extract_ligand_from_complex_pymol(complex_path, output_path):
    """Extract ligand from a protein-ligand complex using PyMOL"""
    try:
        protein_chains, ligand_res = classify_residues(complex_path)
        
        # Initialize PyMOL
        cmd.reinitialize()
        cmd.load(complex_path, "complex")
        
        # Remove excluded ions and water
        cmd.remove(" or ".join([f"resn {res}" for res in st.session_state.EXCLUDED_IONS]))
        
        # Select and save ligand
        if ligand_res:
            ligand_sel = " or ".join(
                [f"(chain {c} and resn {rn} and resi {rid})" for (c, rn, rid) in ligand_res]
            )
            cmd.select("ligand", f"complex and ({ligand_sel})")
            cmd.save(output_path, "ligand")
            return True
        else:
            st.session_state.error_log.append(f"No ligand found in {complex_path}")
            return False
            
    except Exception as e:
        st.session_state.error_log.append(f"Error extracting ligand from {complex_path}: {str(e)}")
        return False

def parse_residue_string(res_str):
    """Parse a residue string into a list of residue numbers"""
    residues = set()
    parts = res_str.split('+')
    for part in parts:
        part = part.strip()
        if '-' in part:
            start, end = part.split('-')
            try:
                start = int(start)
                end = int(end)
                residues.update(range(start, end+1))
            except ValueError:
                continue
        else:
            try:
                residue = int(part)
                residues.add(residue)
            except ValueError:
                continue
    return sorted(residues)

def residues_to_pymol_format(residues):
    """Convert a list of residues to PyMOL format"""
    return "+".join(map(str, residues))

def calculate_bounding_box(protein_path, residues, padding=2.0):
    """Calculate a bounding box around specified residues"""
    try:
        # Initialize PyMOL
        cmd.reinitialize()
        cmd.load(protein_path, "protein")
        
        # Select residues
        res_str = residues_to_pymol_format(residues)
        cmd.select("pocket", f"resi {res_str}")
        
        # Get center of mass
        center = cmd.centerofmass("pocket")
        
        # Get coordinates of selected atoms
        model = cmd.get_model("pocket")
        coords = np.array([atom.coord for atom in model.atom])
        
        # Calculate bounding box size
        min_coords = np.min(coords, axis=0) - padding
        max_coords = np.max(coords, axis=0) + padding
        size = max_coords - min_coords
        
        # Update center based on bounding box
        center = (min_coords + max_coords) / 2
        
        return center, size
    except Exception as e:
        st.session_state.error_log.append(f"Error calculating bounding box: {str(e)}")
        return None, None

# protein structure
def get_protein_structure(pdb_path):
    """Get protein structure data from PDB file"""
    with open(pdb_path, 'r') as f:
        pdb_data = f.read()
    return pdb_data