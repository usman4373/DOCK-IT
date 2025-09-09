# DOCK-IT

A Streamlit application that performs virtual screening with the [GNINA](https://github.com/gnina/gnina) molecular-docking software and conducts protein–ligand interaction analysis using [PLIP](https://github.com/pharmai/plip).


## 📑 Table of contents

1. [📝 Overview](#-overview)
2. [📦 Installation](#-installation)
    - [GNINA Installation](#gnina-installation)
3. [🏃 How to run](#-how-to-run)
4. [📂 Input file formats & expected structure](#-input-file-formats--expected-structure)
5. [🔄 App workflows and usage](#-app-workflows-and-usage)
   - [Blind Docking](#blind-docking)
   - [Site-Specific Docking with Reference Ligand](#site-specific-docking-with-reference-ligand)
   - [Residue Based Docking](#residue-based-docking)
6. [📚 Citation](#-citation)
7. [🤝 Acknowledgments](#-acknowledgments)


## 📝 Overview

<p align="justify"> <b>DOCK-IT</b> is a comprehensive Streamlit-based application that facilitates virtual screening of compounds against protein targets using GNINA, an open-source molecular docking software. The application provides three distinct docking workflows with detailed interaction analysis. Everything is fully automated, eliminating manual steps and significantly reducing the time required for the entire process. </p>

**Key features include:**

- Blind docking across entire protein structures

- Site-specific docking using reference ligands

- Residue-based docking with custom binding site definition

- SMILES to PDB conversion for ligand preparation

- Protein-ligand interaction analysis using PLIP

- Comprehensive results reporting with Excel formatting

- Fully automated workflows


## 📦 Installation

**Create new conda environment (Recommended):**

**From repository root, open terminal and run:**

```
conda create -n dockit python=3.10
```

```
conda activate dockit
```

```
pip install -r requirements.txt
```

**Install PyMOL from `conda-forge` instead of `pip`:**

```
conda install -c conda-forge pymol-open-source
```

- Alternatively, download PyMOL from the [official website](https://www.pymol.org/) and add it to the system's PATH.

### GNINA Installation

Choose one of the following methods to install `GNINA`:

- Install GNINA from source code by following the instructions from the official repository:

    https://github.com/gnina/gnina

- Download GNINA Binary:

    Choose any of the following based on your hardware:

    If you have a GPU:
  
    ```
    wget https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2
    ``` 
    
    If you do not have a GPU:

    ```
    wget https://github.com/gnina/gnina/releases/download/v1.0.1/gnina
    ```

- Once you’ve downloaded the GNINA binary, rename it to `gnina` (removing the version number) and make the binary executable:

```
chmod +x gnina
```


## Add GNINA Binary to System's PATH

- After making `GNINA` binary executable, note the full path to the binary file. For example:
    
    `/home/gray_pc/apps/bin/gnina`

- Open your shell configuration file, `bashrc` or `zshrc` (run in terminal):

```
nano ~/.bashrc
```

- Add the following line at the end (replace `/home/gray_pc/apps/bin/gnina` with your actual path):

```
export PATH="/home/gray_pc/apps/bin/gnina:$PATH"
```
    
- Save the file by:

    `Ctrl+O` → Writes the file.

    Prompt appears asking for the file name.

    `Enter` → Confirms the current file name and proceeds with saving.

    `Ctrl+X` → Closes the editor.

- Reload the shell:

```
source ~/.bashrc   # or source ~/.zshrc
```

### Note: Ensure the GNINA binary is executable (chmod +x) and added to your system’s `PATH`, otherwise the application will not run.


## 🏃 How to run

**Navigate to the directory containing the DOCK-IT `app.py` file and run in terminal:**

```
streamlit run app.py
```

The application will open in your default web browser at `http://localhost:8501`


## 📂 Input file formats & expected structure

- Protein files must be in `.pdb` format, and all protein structures should be stored in a single directory.

- DOCK-IT supports ligand input in two ways:
    
    a. `.pdb` format ligand files stored in a single directory.

    b. A `.csv` file containing SMILES strings along with molecule names or IDs.

- The application will automatically convert SMILES to 3D structures in `.pdb` format

- Reference Ligand Files (for site-specific docking):

    The required format is `.pdb` files of protein–ligand complexes, containing protein structures with their bound reference ligands.

    The application will automatically extract the ligand coordinates from these files to define the binding site.


## 🔄 App workflows and usage

DOCK-IT provides three distinct docking workflows:

### Blind Docking

Screen ligands against the entire protein surface to identify potential binding sites.

**Workflow:**

- Input protein `PDB` files directory

- Input ligands (either as `PDB` files directory or SMILES `CSV` file)

- Set docking parameters (scoring function, exhaustiveness, etc.)

- Run screening

- Review results with detailed interaction reports


### Site-Specific Docking with Reference Ligand

Focus docking on a specific binding site using a known reference ligand.

**Workflow:**

- Input protein `PDB` files directory

- Input reference ligand complexes directory

- Input screening ligands (`PDB` files directory or SMILES `CSV` file)

- The application automatically maps proteins to reference ligands based on sequence similarity

- Set docking parameters including binding box padding around reference ligand

- Run targeted screening

- Review results with detailed interaction reports


### Residue-Based Docking

Define a custom binding site by specifying specific protein residues.

**Workflow:**

- Input protein `PDB` files directory

- Input screening ligands (`PDB` files directory or SMILES `CSV` file)

- Specify residue positions for each protein using format:

  Individual residues: `45+67+23`

  Residue ranges: `79-100` or `20-27+45+63-76`

- Adjust bounding box size with interactive 3D visualization

- Set docking parameters

- Run targeted screening

- Review results with detailed interaction reports


**All workflows include:**

- Automatic pose splitting and complex creation

- Protein-ligand interaction analysis using `PLIP`

- Formatted Excel reports

- Per-protein result organization


# 📚 Citation

If you use DOCK-IT in your research, please cite:

```
DOCK-IT. GitHub: https://github.com/usman4373/DOCK-IT
```

## 🤝 Acknowledgments

- `GNINA` development team for the excellent docking software
- `Streamlit` team for the powerful web application framework
- All open-source libraries that make this project possible
