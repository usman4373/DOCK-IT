import streamlit as st
import streamlit.components.v1 as components
import py3Dmol


def create_3dmol_view(pdb_data, residues, center, size, padding=2.0):
    """Create a 3Dmol.js view of the protein with highlighted residues and bounding box"""
    view = py3Dmol.view(width=600, height=400)
    view.addModel(pdb_data, 'pdb')
    
    # Style the protein with rainbow colors
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    
    # Highlight selected residues
    if residues:
        res_str = " or ".join([f"{r}" for r in residues])
        view.setStyle({'resi': res_str}, {'stick': {'colorscheme': 'redCarbon'}})
    
    # Add bounding box with addBox method
    if center is not None and size is not None:
        x, y, z = center
        sx, sy, sz = size
        
        # Add a semi-transparent box
        view.addBox({
            'center': {'x': x, 'y': y, 'z': z},
            'dimensions': {'w': sx, 'h': sy, 'd': sz},
            'color': 'red',
            'opacity': 0.8,
            'wireframe': False
        })

        # Create box vertices
        corners = [
            [x-sx/2, y-sy/2, z-sz/2],
            [x+sx/2, y-sy/2, z-sz/2],
            [x+sx/2, y+sy/2, z-sz/2],
            [x-sx/2, y+sy/2, z-sz/2],
            [x-sx/2, y-sy/2, z+sz/2],
            [x+sx/2, y-sy/2, z+sz/2],
            [x+sx/2, y+sy/2, z+sz/2],
            [x-sx/2, y+sy/2, z+sz/2]
        ]
        
        # Draw box edges
        edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),  # bottom
            (4, 5), (5, 6), (6, 7), (7, 4),  # top
            (0, 4), (1, 5), (2, 6), (3, 7)   # sides
        ]
        
        for i, j in edges:
            view.addLine({
                'start': {'x': corners[i][0], 'y': corners[i][1], 'z': corners[i][2]},
                'end': {'x': corners[j][0], 'y': corners[j][1], 'z': corners[j][2]},
                'color': 'black',
                'radius': 0.2
            })
    
    view.zoomTo()
    return view._make_html()
