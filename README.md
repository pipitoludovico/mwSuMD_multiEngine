# mwSuMD - Molecular Dynamics Supervision Tool

## Introduction

Welcome to mwSuMD, a tool for supervising molecular dynamics simulations. This README will guide you through the setup process to ensure a smooth experience with mwSuMD. Multiple walker supervised molecular dynamics (mwSuMD) is an adaptive sampling molecular dynamics (MD) technique for studying ligand-receptor binding and unbinding pathways or conformational changes (https://doi.org/10.7554/eLife.96513.2) not accessible to classic MD simulations. It is a scalable multi-engine software meant to run batches of simulations either in parallel or serial. The philosophy behind mwSuMD is to narrow down as much as possible the use of third-party software, limiting the tools only to GNU free licenced software. Furthermore, we decided to open the software to different engines to facilitate its use for those familiar with either NAMD, GROMACS, ACEMD or OPENMM.
The code is available and you can easily customise the files, your inputs and your criteria to conduct your simulations. 


## Prerequisites
Conda is not mandatory but is warmly suggested, to avoid dependency conflicts or versioning issues.

Before you start, make sure you have the following prerequisites installed on your system:

- Python (>=3.7)
- pip (Python package installer)
- VMD
- GPUtil==1.4.0
- MDAnalysis==2.7.0
- numpy==1.24.3
- openmm
- pandas==2.0.3
- setuptools==70.0.1
- openmmplumed (optional)


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
!!! You "./system" folder should include the input file you used to initially equilibrate the system !!!


# mwSuMD looks for file extension to determine the engine you used for the equilibration!
It is important you either keep your input file or even create a tag file for the engine selection.
Each engine compiles the binaries in a slight different way. mwSuMD looks for file extensions that are usually default for that specific engine.

For example: if you used NAMD to produce the binaries, mwSuMD will expect either your .namd input file, or either an empty tag file ending with .namd inside the "./system" folder. The same way, if you used ACEMD you will need a file ending in .inp, an .gro if you used GROMACS, or a .chk file if you used openMM.
The choice of these files was dictated by their simplicity and uniqueness.

1. Create your working folder, and make a new folder named "system" in. The settings_mwSuMD.inp file must be located outside this folder.

2. Copy your structure inside "system" and run your equilibration there using your preferred molecular dynamics engine among OPENMM, NAMD, ACEMD, or GROMACS. You can use either AMBER or CHARMM force fields, or GROMACS's own OPLS ff.

3. Ensure the "system" folder contains a properly equilibrated system with the necessary binary files generated during equilibration.


## Extra Parameters
You can add additional parameter files by either adding your "'.param', '.prm', '.par', '.top', '.rtf', '.str'" files straight into the "./system" folder or, for a tidy working environment, by making an external directory called "parameters" and adding your files there. Remember: your custom parameter files will be added last, hence they will update the builtin files!

## Using PLUMED
It is possible to run MetaDynamics in conjunction with mwSuMD, should you want to. To activate PLUMED, make a folder called "plumed" outside "./system", in your cwd, and write a file called "plumed.inp" in it. Edit your file following your needs and mwSuMD will read that during each step.


## Configuration

1. Modify the `settings_mwSuMD.inp` setting file, selecting one or two metrics for supervision using the MDAnalysis atomic selection syntax. You can supervise distances, contacts, number of hydrogen bonds, number of solvation molecules around your selection.

## Running mwSuMD

1. Add mwSuMD to your system's path for convenience.

2. Launch mwSuMD where `settings_mwSuMD.inp` is located:

    ```bash
    mwSuMD.py
    ```

Enjoy using mwSuMD for monitoring your molecular dynamics simulations! If you have any questions or encounter issues, feel free to reach out to us on the [mwSuMD GitHub repository](https://github.com/your-username/mwSuMD).

## Config File Keyword explanation


Restart [YES/NO]; specify if you wish to continue a previously interrupted mwSuMD simulation. If the “trajectories” folder is present and contains xtc files output from the previous simulation, then is Restart YES by default

WrapEngine = [VMD/MDA]; select how to perform the wrapping of the short simulations before computation of the supervised metric(s) [set MDA to use MDAnalysis, VMD to use VMD]

WrapOn = selection for wrapping

FilterOut = selection for filtering out and reducing the weight of the xtc files stored inside “trajectories”

Output = Output name of the xtc files, without extension

Temperature = Simulation temperature in K

CheckEvery = number of mwSuMD cycles between checks for detecting a stuck simulation, [Default= based on RelaxTime]

RelaxTime = duration of the unsupervised MD performed in case of stuck simulation, in ns; [Default = 5 ns]

Fails = number of times a simulation can be detected as stuck before automatically killing it, [Default = 5]

Tolerance = Data dispersion (multiplied by 100), measured at the checkpoint defined by CheckEvery that is tolerated, otherwise, a simulation is considered stuck and

automatically killed; a low value of 10 or less should have the effect of avoiding that the simulation is considered stuck [Default = 30]

NumberCV = [1/2]; decide if one or two metrics are to be supervised

Metric_1 = [Distance/RMSD/Contacts/HB/Solvation]

Distance: between Sel_1 and Sel_2

RMSD: between Sel_1=selection for superimposition and Sel_2=selection for RMSD Contacts: between Sel_1 and Sel_2

HB: hydrogen bonds between sel_1 and sel_2

Solvation: waters between sel_1 and sel_2

Sel_1 = selection in MDAnalysis selection language

Sel_2 = selection in MDAnalysis selection language

Cutoff_1 = Threshold of Metric_1 to reach (no unit of measure); mwSuMD stops after this vaklue is reached

Transition_1 = [negative/positive]; negative: Metric_1 should decrease; positive: Metric_1 should increase


Metric_2 = [Distance/RMSD/; Contacts/HB/ Solvation]

Distance: between Sel_3 and Sel_4

RMSD: between Sel_3=selection for superimposition and Sel_4=selection for RMSD; Contacts: between Sel_3 and Sel_4

HB: between sel_3 and sel_4

Solvation: waters between sel_3 and sel_4

Sel_3 = selection in MDAnalysis selection language

Sel_4 = selection in MDAnalysis selection language

Cutoff_2 = Threshold of Metric_2 to reach (no unit of measure); mwSuMD stops after it is reached]

Transition_2 = [negative/positive];

negative: Metric_2 should decrease

positive: Metric_2 should increase

Walkers = number of walkers to run; depends on the number of GPUs available for parallel computing

Timewindow = length of each short simulation, in ps; usually 500 – 1000 ps for binding simulations, 20 -100 ps for unbinding simulations, 50 -300 ps for conformational changes

Timestep = integration time step, in fs

Savefreq = time interval between saved frames, in ps; should be 1/3 of Timewindow to obtain 3 frames for each short simulation, 1/3 of Timewindow to obtain 4 frames for each short simulation, etc.


It is possible to run plumed at the same time as the mwSuMD simulations by creating a folder called "plumed" outside "system", where this input file is located, and placing your plumed.inp file in it.

