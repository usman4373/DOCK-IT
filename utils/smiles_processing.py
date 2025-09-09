from rdkit import Chem
from rdkit.Chem import AllChem
import streamlit as st


def smiles_to_pdb(smiles, output_path, name):
    """Convert SMILES string to PDB format"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        
        with open(output_path, 'w') as f:
            f.write(Chem.MolToPDBBlock(mol))
        return True
    except Exception as e:
        st.session_state.error_log.append(f"Error converting SMILES {name}: {str(e)}")
        return False