import os
from subprocess import Popen, DEVNULL
from warnings import filterwarnings

import MDAnalysis as Mda
import pkg_resources

filterwarnings(action='ignore')


class mwInputParser:
    folder = os.getcwd()
    os.makedirs(f'{folder}/trajectories', exist_ok=True)
    inputFile = 'simulation_settings_mwSuMD.inp'
    initialParameters = {'Root': os.getcwd()}
    selection_list = []
    walker_metrics = []
    package_dir = pkg_resources.resource_filename('mwSuMD_lib', 'parameters')
    parameterFolderPath = os.path.abspath(package_dir)
    coordExtensions = ('.pdb', '.gro', '.inpcrd', '.restrt')
    parameterPaths = [parameterFolderPath, ]
    if 'parameters' in os.listdir(folder):
        parameterPaths.append(os.path.join(folder, "parameters"))
    new_value = 0
    max_value = 0
    metric_1 = 0
    metric_2 = 0
    excludedGPUS = []

    def __init__(self):
        self.customInputFileExtension = ('namd', 'inp', 'mdp')
        self.outExtensions = ('cpi', 'coor', 'vel', 'xsc')
        self.fileExtensions = ('.psf', '.pdb', '.mdp', '.gro', '.cpt', 'top', '.prmtop', '.tpr')
        self.initialParametersameter_extensions = ('.param', '.prm', '.par', '.top', '.rtf', '.str')
        self.trajCount = len([traj for traj in os.listdir(f'{self.folder}/trajectories') if traj.endswith('.xtc')])
        self.allowedMetrics = ("DISTANCE", "RMSD", "CONTACTS", "HB", "SOLVATION")
        if not os.path.isfile(f'{self.folder}/{self.inputFile}'):
            print('Input file for SuMD simulation required')
            quit()

    def checkEngine(self):
        self.getSystem()
        self.initialParameters['MDEngine'] = 'GROMACS' if any(
            file.endswith('.gro') for file in os.listdir('./system')) else 'NAMD' if any(
            file.endswith('.namd') for file in os.listdir('./system')) else 'ACEMD' if any(
            file.endswith('.inp') for file in os.listdir("./system")) else "OPENMM" if any(
            file.endswith('.chk') for file in os.listdir("./system")) else None
        if not self.initialParameters['MDEngine']:
            raise FileNotFoundError(
                "No MD Engine detected. Make sure you kept your input setting's file in the system folder")
        if self.initialParameters.get('MDEngine') == 'GROMACS':
            if not any(file.endswith('tpr') for file in os.listdir('./system')):
                try:
                    Popen(f"touch {self.folder}/system/gic.mdp", shell=True, stdout=DEVNULL).wait()
                    Popen(
                        f"gmx grompp -c {self.folder}/system/{self.initialParameters['GRO']} -p {self.folder}/system/{self.initialParameters['TOP']} -f {self.folder}/system/gic.mdp -o {self.folder}/system/{str(self.initialParameters['GRO']).replace('gro', 'tpr')}",
                        shell=True, stdout=DEVNULL).wait()
                except FileNotFoundError:
                    print(
                        "You need a .tpr file with GROMACS for MDAnalysis as GROMACS's .gro or .top lack any meaningful information for MDAnalysis")
                    exit()

    def getSystem(self):
        for ext in self.fileExtensions:
            for file in os.listdir('./system'):
                if file.endswith(ext):
                    self.initialParameters[ext.replace('.', '').upper()] = f"{self.folder}/system/{file}"

    def getParameters(self):
        if not os.path.exists(f'{self.folder}/plumed'):
            self.initialParameters['PLUMED'] = None
        else:
            for file in os.listdir(self.folder + "/plumed"):
                if file.endswith('.inp'):
                    self.initialParameters['PLUMED'] = f'{self.folder}/plumed/{file}'

        if 'Parameters' not in self.initialParameters:
            self.initialParameters['Parameters'] = []
            for main, _, paramFolder in os.walk(f'{self.folder}/system'):
                for params in paramFolder:
                    if params.endswith(self.initialParametersameter_extensions):
                        self.initialParameters['Parameters'].append(f"{main}/{params}")
            for path_ in self.parameterPaths:
                for dirpath, dirnames, generalParams in os.walk(path_):
                    for filename in [f for f in generalParams if f.endswith(self.initialParametersameter_extensions)]:
                        self.initialParameters['Parameters'].append(dirpath + "/" + filename)
        return self.initialParameters['Parameters']

    def getReferencePDB(self) -> None:
        """Check if any folder called ./reference exists inside cwd or system. Else, it will set the reference using the
         first coordinate file found inside ./system"""
        if not os.path.exists(f'{self.folder}/system/reference'):
            systemPath = f"{self.folder}/system/"
        else:
            systemPath = f"{self.folder}/system/reference/"
        pdbFile = [f"{systemPath}" + pdb for pdb in os.listdir(f"{systemPath}") if pdb.endswith(self.coordExtensions)][
            0]
        self.initialParameters['REFERENCE'] = pdbFile

    def getForcefields(self):
        self.initialParameters['Forcefield'] = 'CHARMM' \
            if any(fileSys.endswith('.psf') for fileSys in os.listdir(f'{self.folder}/system')) else 'AMBER' \
            if any(fileSys.endswith('.prmtop') for fileSys in os.listdir(f'{self.folder}/system')) else 'GROMOS'

    def getSettingsFromInputFile(self):
        if self.initialParameters.get('MDEngine') != 'GROMACS':
            u = Mda.Universe(f"{self.initialParameters['PDB']}")
        else:
            u = Mda.Universe(f"{self.initialParameters['GRO']}")

        # Default settings:
        self.initialParameters['NumberCV'] = 1

        self.initialParameters['CheckEvery'] = 100
        self.initialParameters['Temperature'] = 310
        self.initialParameters['WrapEngine'] = "MDA"
        self.initialParameters['WrapOn'] = "protein"
        self.initialParameters['FilterOut'] = "protein"
        self.initialParameters['PROTEIN_RESTRAINTS'] = None
        self.initialParameters['MEMBRANE_RESTRAINTS'] = None
        self.initialParameters['LIGAND_RESNAME'] = None
        self.initialParameters['Metric_1'] = None
        self.initialParameters['Metric_2'] = None
        self.initialParameters['CUSTOMFILE'] = None

        self.initialParameters['Timestep'] = 2
        self.initialParameters['Savefreq'] = 20

        self.initialParameters['Relax'] = False
        self.initialParameters['RelaxTime'] = 10  # duration of the cMD
        self.initialParameters['RelaxTimestep'] = 2  # Timestep in ps
        self.initialParameters['RelaxSavefreq'] = 20
        self.initialParameters['KeepPlumedForRelax'] = False

        self.initialParameters['Wrap'] = 'protein and name CA'
        self.initialParameters['Fails'] = 5
        self.initialParameters['Tolerance'] = 0.3
        self.initialParameters['ACTUAL_DISTANCE'] = None

        for customFile in os.listdir(f"{self.initialParameters['Root']}/system"):
            if customFile.startswith('production') and customFile.endswith(self.customInputFileExtension):
                self.initialParameters['CUSTOMFILE'] = f"{self.folder}/system/{customFile}"
                if self.trajCount == 0:
                    for extension in self.outExtensions:
                        Popen(f'cp {self.folder}/system/*.{extension} restarts/previous.{extension}',
                              shell=True).wait()

        with open(self.inputFile, "r") as infile:

            for line in infile:
                if line.startswith('#'):
                    continue

                if line.startswith("PROTEIN_RESTRAINTS"):
                    if line.split("=")[1].lower().strip() == "no":
                        continue
                    else:
                        tmpline = line.split("=")[1].split(",")
                        cleanLine = [parte.strip() for parte in tmpline]
                        self.initialParameters['PROTEIN_RESTRAINTS'] = cleanLine

                if line.startswith("MEMBRANE_RESTRAINTS"):
                    if line.split("=")[1].lower().strip().startswith('y'):
                        self.initialParameters['MEMBRANE_RESTRAINTS'] = True

                if line.startswith("LIGAND_RESNAMES"):
                    if line.split("=")[1].lower().strip() == "no":
                        continue
                    else:
                        tmpline = line.split("=")[1].split(",")
                        cleanLine = [parte.strip() for parte in tmpline]
                        self.initialParameters['LIGAND_RESNAMES'] = cleanLine

                if line.startswith("WrapEngine"):
                    self.initialParameters['WrapEngine'] = line.split("=")[1].strip()

                if line.startswith('Temperature'):
                    if line.split("=")[1].lower().strip() == "no":
                        continue
                    else:
                        self.initialParameters['Temperature'] = float(line.split('=')[1].strip())

                if line.split("=")[0].strip() == "KeepPlumedForRelax":
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['KeepPlumedForRelax'] = True

                if line.split("=")[0].strip() == "RelaxTime":
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['RelaxTime'] = float(line.split('=')[1].strip())

                if line.startswith('RelaxSavefreq'):
                    if line.split('=')[1].strip() != '' and line.split('=')[1].strip().isdigit():
                        self.initialParameters['RelaxSavefreq'] = int(line.split('=')[1].strip())
                    else:
                        raise ValueError("Please set the Savefreq in your input file as an integer number.")

                if line.split("=")[0].strip() == "RelaxTimestep":
                    if line.split('=')[1].strip() != '' and line.split('=')[1].strip().isdigit():
                        self.initialParameters['RelaxTimestep'] = int(line.split('=')[1].strip())
                    else:
                        raise ValueError("Please set the Timestep in your input file as an integer number.")

                if line.startswith('CheckEvery'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['CheckEvery'] = int(line.split('=')[1].strip())

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

                if line.startswith('WrapOn'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['WrapOn'] = line.split('=')[1].strip()

                if line.startswith('FilterOut'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['FilterOut'] = line.split('=')[1].strip()

                if line.startswith('Fails'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['Fails'] = int(line.split('=')[1].strip())

                if line.startswith("Sel_"):
                    if line.split('=')[1].strip() != '':
                        selection_ = line.split('=')[1].strip()
                        if len(u.select_atoms(selection_)) != 0:
                            self.selection_list.append(line.split('=')[1].strip())
                        else:
                            raise ValueError(f"Your selection {selection_} pointed to 0 atoms: please check your selection with your structure file")
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
        if self.trajCount != 0:
            directory = f'{self.folder}/restarts'
        else:
            directory = f'{self.folder}/system'
        self.initialParameters.update(
            {extension: file for file in os.listdir(directory) for extension in self.outExtensions if
             file.endswith(extension) and 'minimzed' not in file})

    def getSettings(self):
        self.checkEngine(), self.getParameters(), self.getForcefields()
        self.getSettingsFromInputFile(), self.getRestartOutput(), self.getReferencePDB()
        return self.initialParameters, self.selection_list, self.parameterFolderPath
