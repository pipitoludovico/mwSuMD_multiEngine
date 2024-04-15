import os
from subprocess import Popen, DEVNULL
import pkg_resources
import MDAnalysis as Mda
from warnings import filterwarnings

filterwarnings(action='ignore')


class mwInputParser:
    folder = os.getcwd()
    os.makedirs(f'{folder}/trajectories', exist_ok=True)
    inputFile = 'simulation_settings_mwSuMD.inp'
    initialParameters = {'Root': os.getcwd()}
    selection_list = []
    walker_metrics = []
    package_dir = ''
    if 'parameters' not in os.listdir(folder):
        package_dir = pkg_resources.resource_filename('mwSuMD_lib', 'parameters')
        parameterFolderPath = os.path.abspath(package_dir)
    else:
        parameterFolderPath = os.path.abspath('parameters')
    new_value = 0
    max_value = 0
    metric_1 = 0
    metric_2 = 0
    excludedGPUS = []

    def __init__(self):
        self.customInputFileExtension = ('namd', 'inp', 'mdp')
        self.outExtensions = ('cpt', 'cpi', 'coor', 'vel', 'xsc')
        self.fileExtensions = ('.psf', '.pdb', '.mdp', '.gro', '.cpt', 'top', '.prmtop', '.tpr')
        self.initialParametersameter_extensions = ('.param', '.prm', '.par', '.top', '.rtf', '.str')
        self.trajCount = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
        self.allowedMetrics = ("DISTANCE", "RMSD", "CONTACTS", "HB", "SOLVATION")
        if not os.path.isfile(f'{self.folder}/{self.inputFile}'):
            print('Input file for SuMD simulation required')
            quit()
        if not os.path.exists(f'{self.folder}/plumed'):
            self.initialParameters['PLUMED'] = None
        else:
            for file in os.listdir(self.folder + "/plumed"):
                if file.endswith('.inp'):
                    self.initialParameters['PLUMED'] = f'{self.folder}/plumed/{file}'

    def checkEngine(self):
        self.getSystem()
        self.initialParameters['MDEngine'] = 'GROMACS' if any(
            file.endswith('.gro') for file in os.listdir('./system')) else 'NAMD' if any(
            file.endswith('.namd') for file in os.listdir('./system')) else 'ACEMD' if any(file.endswith('.inp') for file in os.listdir("./system")) else "OPENMM"
        if self.initialParameters.get('MDEngine') == 'GROMACS':
            if not any(file.endswith('tpr') for file in os.listdir('./system')):
                try:
                    Popen(f"touch {self.folder}/system/gic.mdp", shell=True, stdout=DEVNULL)
                    Popen(
                        f"gmx grompp -c {self.folder}/system/{self.initialParameters['GRO']} -p {self.folder}/system/{self.initialParameters['TOP']} -f {self.folder}/system/gic.mdp -o {self.folder}/system/{str(self.initialParameters['GRO']).replace('gro', 'tpr')}",
                        shell=True, stdout=DEVNULL)
                except FileNotFoundError:
                    print(
                        "You need a .tpr file with GROMACS for MDAnalysis as GROMACS's .gro or .top lack any meaningful information for MDAnalysis")
                    exit()

    def getSystem(self):
        for ext in self.fileExtensions:
            for file in os.listdir('./system'):
                if file.endswith(ext):
                    self.initialParameters[ext.replace('.', '').upper()] = file

    def getParameters(self):
        if 'Parameters' not in self.initialParameters:
            self.initialParameters['Parameters'] = []
            for params in os.listdir(f'{self.folder}/system'):
                if params.endswith(self.initialParametersameter_extensions):
                    self.initialParameters['Parameters'].append(self.folder + "/system/" + params)
            for dirpath, dirnames, generalParams in os.walk(self.parameterFolderPath):
                for filename in [f for f in generalParams if f.endswith(self.initialParametersameter_extensions)]:
                    self.initialParameters['Parameters'].append(dirpath + "/" + filename)
        else:
            return self.initialParameters['Parameters']

    def getReferencePDB(self):
        def checkRMSDoption():
            if not os.path.isdir(f'{self.folder}/system/reference'):
                print('You need a reference folder with a reference pdb in it if you use the RMSD as a metric.')
                exit()
            if len(os.listdir(f'{self.folder}/system/reference')) > 0:
                files = os.listdir(f'{self.folder}/system/reference')
                reference = next((file for file in files if file.endswith('.pdb')), None)
                if reference:
                    self.initialParameters['REFERENCE'] = reference
                else:
                    print("Put a reference pdb file in the 'reference' folder inside system and rerun.")
                    exit()

        if self.initialParameters.get('Metric_1') == 'RMSD' or self.initialParameters.get('Metric_2') == 'RMSD':
            checkRMSDoption()

    def getForcefields(self):
        self.initialParameters['Forcefield'] = 'CHARMM' \
            if any(fileSys.endswith('.psf') for fileSys in os.listdir(f'{self.folder}/system')) else 'AMBER' \
            if any(fileSys.endswith('.prmtop') for fileSys in os.listdir(f'{self.folder}/system')) else 'GROMOS'

    def getSettingsFromInputFile(self):
        if self.initialParameters.get('MDEngine') != 'GROMACS':
            u = Mda.Universe(f"{self.folder}/system/{self.initialParameters['PDB']}")
        else:
            u = Mda.Universe(f"{self.folder}/system/{self.initialParameters['GRO']}")

        # Default settings:
        self.initialParameters['Metric_1'] = None
        self.initialParameters['Metric_2'] = None
        self.initialParameters['Restart'] = None
        self.initialParameters['CUSTOMFILE'] = None
        self.initialParameters['Timestep'] = 2
        self.initialParameters['Savefreq'] = 20
        self.initialParameters['Wrap'] = 'protein and name CA'
        self.initialParameters['Fails'] = 5
        self.initialParameters['Tolerance'] = 0.3

        for customFile in os.listdir(f"{self.initialParameters['Root']}/system"):
            if customFile.startswith('production') and customFile.endswith(self.customInputFileExtension):
                print("Custom Input File found: ", customFile)
                print(
                    "Warning: Make sure your custom input file is pointing at the binaries in the new restart folder!")
                self.initialParameters['CUSTOMFILE'] = f"{self.folder}/system/{customFile}"
                if self.trajCount == 0 and self.initialParameters['Restart'] == 'NO':
                    for extension in self.outExtensions:
                        Popen(f'cp {self.folder}/system/*.{extension} restarts/previous.{extension}',
                              shell=True)

        with open(self.inputFile, "r") as infile:
            self.initialParameters['NumberCV'] = None
            self.initialParameters['RelaxTime'] = 5
            self.initialParameters['Relax'] = False
            self.initialParameters['CheckEvery'] = None
            for line in infile:
                if line.startswith('#'):
                    continue
                if line.startswith('RelaxTime'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['RelaxTime'] = float(line.split('=')[1].strip())

                if line.startswith('CheckEvery'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['CheckEvery'] = float(line.split('=')[1].strip())

                if line.startswith('Restart'):
                    mwInputParser.initialParameters['Restart'] = line.split('=')[1].strip().upper()

                if line.startswith('Output'):
                    mwInputParser.initialParameters['Output'] = line.split('=')[1].strip()

                if line.startswith('ligand_HB'):
                    mwInputParser.initialParameters['ligand_HB'] = line.split('=')[1].strip()

                if line.startswith('Tolerance'):
                    if line.split('=')[1].strip() != '':
                        if float(line.split("=")[1].strip()) in range(0, 101):
                            self.initialParameters['Tolerance'] = float(line.split('=')[1].strip()) / 100
                        else:
                            raise ValueError("Tolerance range valid: 0 - 100")

                if line.startswith('NumberCV'):
                    if line.split('=')[1].strip() != '':
                        if int(line.split('=')[1].strip()) in range(1, 3):
                            self.initialParameters['NumberCV'] = int(line.split('=')[1].strip())
                        else:
                            raise ValueError("Only NumberCV 1 or 2 allowed")
                    else:
                        raise ValueError("Please input the number of CVs (1 or 2 allowed)")

                if line.startswith('Metric_1'):
                    if line.split('=')[1].strip() != '':
                        if line.split('=')[1].strip().upper() in self.allowedMetrics:
                            self.initialParameters['Metric_1'] = line.split('=')[1].strip().upper()
                        else:
                            raise ValueError("Invalid Metric 1. Only metrics allowed: ", self.allowedMetrics)

                elif line.startswith('Metric_2'):
                    if line.split('=')[1].strip() != '':
                        if line.split('=')[1].strip().upper() in self.allowedMetrics:
                            self.initialParameters['Metric_2'] = line.split('=')[1].strip().upper()
                        else:
                            raise ValueError("Invalid Metric 2. Only metrics allowed: ", self.allowedMetrics)

                if line.startswith('Cutoff_1'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['Cutoff_1'] = float(line.split('=')[1].strip())

                if line.startswith('Cutoff_2'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['Cutoff_2'] = float(line.split('=')[1].strip())

                if line.startswith('Walkers'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['Walkers'] = int(line.split('=')[1].strip())
                    else:
                        raise ValueError("Please set a number of desired walkers in your input file")

                if line.startswith('Timewindow'):
                    if line.split('=')[1].strip() != '' and line.split('=')[1].strip().isdigit():
                        self.initialParameters['Timewindow'] = int(line.split('=')[1].strip())
                    else:
                        raise ValueError("Please set a Timewindow value.")

                if line.startswith('Timestep'):
                    if line.split('=')[1].strip() != '' and line.split('=')[1].strip().isdigit():
                        self.initialParameters['Timestep'] = int(line.split('=')[1].strip())
                    else:
                        raise ValueError("Please set the Timestep in your input file as an integer number.")

                if line.startswith('Savefreq'):
                    if line.split('=')[1].strip() != '' and line.split('=')[1].strip().isdigit():
                        self.initialParameters['Savefreq'] = int(line.split('=')[1].strip())
                    else:
                        raise ValueError("Please set the Savefreq in your input file as an integer number.")

                if line.startswith('Transition_1'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['Transition_1'] = line.split('=')[1].strip().lower()

                if line.startswith('Transition_2'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['Transition_2'] = line.split('=')[1].strip().lower()

                if line.startswith('Wrap'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['Wrap'] = line.split('=')[1].strip()

                if line.startswith('Fails'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['Fails'] = int(line.split('=')[1].strip())

                if line.startswith("Sel_"):
                    if line.split('=')[1].strip() != '':
                        if len(u.select_atoms(f"{line.split('=')[1].strip()}")) != 0:
                            self.selection_list.append(line.split('=')[1].strip())
                        else:
                            raise ValueError(
                                "One of your selection pointed to 0 atoms: please check your selection with your structure file")
        # cleaning Universe as the check is completed
        del u
        if self.initialParameters['NumberCV'] == 2 and (
                not self.initialParameters.get('Metric_1') or not self.initialParameters.get('Metric_2')):
            raise ValueError('Make sure to select both Metric 1 and 2 if NumberCV = 2!')

        if not self.initialParameters.get('Metric_1') and not self.initialParameters.get('Metric_2'):
            raise ValueError(
                "Please make sure if you choose at least one metric to supervise (Distance, Contacts, RMSD, HB)")
        if self.initialParameters.get('NumberCV') == 2 and (
                not self.initialParameters.get('Metric_1') or not self.initialParameters.get('Metric_2')):
            raise ValueError(
                "Please make sure if you use CV2 to specify all the CVs choosing one metric to supervise (Distance, Contacts, RMSD, HB)")
        if not self.initialParameters.get('NumberCV'):
            print("Please set the NumberCV in the input file")
            exit()

    def getRestartOutput(self):
        os.makedirs('restarts', exist_ok=True)
        if self.trajCount != 0 or self.initialParameters['Restart'] == 'YES':
            directory = f'{self.folder}/restarts'
        else:
            directory = f'{self.folder}/system'
        self.initialParameters.update(
            {extension: file for file in os.listdir(directory) for extension in self.outExtensions if
             file.endswith(extension) and 'minimzed' not in file})

    def getSettings(self):
        self.checkEngine(), self.getParameters(), self.getForcefields()
        self.getSettingsFromInputFile(), self.getReferencePDB(), self.getRestartOutput()
        return self.initialParameters, self.selection_list, self.parameterFolderPath
