# mwSuMD - Molecular Dynamics Supervision Tool

## Introduction

Welcome to mwSuMD, a tool for supervising molecular dynamics simulations. This README will guide you through the setup process to ensure a smooth experience with mwSuMD.

## Prerequisites

Before you start, make sure you have the following prerequisites installed on your system:

- Python (>=3.7)
- pip (Python package installer)

## Installation

1. Clone the mwSuMD repository to your local machine:

    ```bash
    git clone https://github.com/pipitoludovico/mwSuMD_multiEngine.git
    ```

2. Navigate to the cloned directory:

    ```bash
    cd mwSuMD
    ```

3. Create a virtual environment to isolate dependencies with venv or conda if you prefer:

    ```bash
    python -m venv mwsumd
    ```

    ```bash
    conda create -n mwsumd 
    ```
    

4. Activate the virtual environment:

    - **On Windows:**
    
        ```bash
        venv\Scripts\activate
        ```
    
    - **On macOS and Linux:**
    
        ```bash
        source venv/bin/activate
        ```
	or with conda in any OS:
        ```bash
        conda activate mwsumd
        ```

5. Install the required dependencies:

    ```bash
    pip install -r requirements.txt
    ```

6. Install mwSuMD as a Python package:

    ```bash
    pip install .
    ```

## Setting up the System

1. Create your working folder, and make a new folder named "system" in. The settings_mwSuMD.inp file must be located outside this folder.

2. Copy your structure inside "system" and run your equilibration there using your preferred molecular dynamics engine among OPENMM, NAMD, ACEMD, or GROMACS. You can use either AMBER or CHARMM force fields, or GROMACS's own OPLS ff.

3. Ensure the "system" folder contains a properly equilibrated system with the necessary binary files generated during equilibration.

## Configuration

1. Modify the `settings_mwSuMD.inp` setting file, selecting one or two metrics for supervision using the MDAnalysis atomic selection syntax. You can supervise distances, contacts, number of hydrogen bonds, number of solvation molecules around your selection.

## Running mwSuMD

1. Add mwSuMD to your system's path for convenience.

2. Launch mwSuMD where `settings_mwSuMD.inp` is located:

    ```bash
    mwSuMD.py
    ```

Enjoy using mwSuMD for monitoring your molecular dynamics simulations! If you have any questions or encounter issues, feel free to reach out to us on the [mwSuMD GitHub repository](https://github.com/your-username/mwSuMD).

