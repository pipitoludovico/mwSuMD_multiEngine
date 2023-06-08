import argparse
import os
import subprocess
import MDAnalysis as mda


class mwInputParser:
    folder = os.getcwd()
    os.makedirs(f'{folder}/trajectories', exist_ok=True)
    inputFile = 'simulation_settings_mwSuMD.inp'
    initialParameters = {'Root': os.getcwd()}
    selection_list = []
    walker_metrics = []
    parameterFolderPath = os.path.abspath('parameters')
    new_value = 0
    max_value = 0
    metric_1 = 0
    metric_2 = 0
    excludedGPUS = []

    def __init__(self):
        self.customInputFileExtension = ('namd', 'inp', 'mdp')
        self.outExtensions = ('cpt', 'cpi', 'gro', 'tpr', 'coor', 'vel', 'xsc')
        self.fileExtensions = ('.psf', '.pdb', '.mdp', '.gro', '.cpt', 'top', '.prmtop', '.tpr')
        self.initialParametersameter_extensions = ('.param', '.prmtop', '.prm', '.par')
        self.trajCount = len([x for x in os.scandir(f'{self.folder}/trajectories')])
        self.allowedMetrics = ("Distance", "RMSD", "Contacts")
        if not os.path.isfile(f'{self.folder}/{self.inputFile}'):
            print('Input file for SuMD simulation required')
            quit()

    def checkEngine(self):
        self.initialParameters['MDEngine'] = 'GROMACS' if any(
            file.endswith('.gro') for file in os.listdir('./system')) else 'NAMD' if any(
            file.endswith('.namd') for file in os.listdir('./system')) else 'ACEMD'
        self.getSystem()

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
                for reference in os.listdir(f'{self.folder}/system/reference'):
                    if reference.endswith('.pdb') or reference.endswith('.gro'):
                        self.initialParameters['REFERENCE'] = reference
                    else:
                        print("Put a reference pdb file in the 'reference' folder inside system and rerun.")
                        exit()

        if self.initialParameters['NumberCV'] == 1 and self.initialParameters['Metric_1'] == 'RMSD':
            checkRMSDoption()
        if self.initialParameters['NumberCV'] == 2 or self.initialParameters['Metric_1'] == 'RMSD' or \
                self.initialParameters['Metric_2'] == 'RMSD':
            checkRMSDoption()

    def getForcefields(self):
        self.initialParameters['Forcefield'] = 'CHARMM' \
            if any(fileSys.endswith('.psf') for fileSys in os.listdir(f'{self.folder}/system')) else 'AMBER' \
            if any(fileSys.endswith('.prmtop') for fileSys in os.listdir(f'{self.folder}/system')) else 'GROMOS'

    def getSettingsFromInputFile(self):
        u = mda.Universe(f"{self.initialParameters['PDB']}")
        # Default settings:
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
                if self.trajCount == 0 and self.initialParameters[
                    'Restart'] == 'NO':  # first run we sort post-equilibration files
                    for extension in self.outExtensions:
                        subprocess.Popen(f'cp {self.folder}/system/*.{extension} restarts/previous.{extension}',
                                         shell=True)

        with open(self.inputFile, "r") as infile:
            self.initialParameters['RelaxTime'] = 5
            self.initialParameters['Relax'] = False
            for line in infile:
                if line.startswith('#'):
                    continue
                if line.startswith('RelaxTime'):
                    if line.split('=')[1].strip() != '':
                        self.initialParameters['RelaxTime'] = float(line.split('=')[1].strip())

                if line.startswith('Tolerance'):
                    if line.split('=')[1].strip() != '' and float(line.split("=")[1].strip()) in range(0, 101):
                        self.initialParameters['Tolerance'] = float(line.split('=')[1].strip()) / 100
                    else:
                        raise ValueError("Tolerance range valid: 0 - 100")

                if line.startswith('NumberCV'):
                    if line.split('=')[1].strip() != '' and int(line.split('=')[1].strip()) in range(1, 3):
                        self.initialParameters['NumberCV'] = int(line.split('=')[1].strip())
                    else:
                        raise ValueError("Only NumberCV 1 or 2 allowed")

                if line.startswith('Metric_1'):
                    if line.split('=')[1].strip() != '' and line.split('=')[1].strip() in self.allowedMetrics:
                        self.initialParameters['Metric_1'] = line.split('=')[1].strip().upper()
                    else:
                        raise ValueError("Invalid Metric 1. Only metrics allowed: ", self.allowedMetrics)

                if line.startswith('Metric_2'):
                    if line.split('=')[1].strip() != '' and line.split('=')[1].strip() in self.allowedMetrics:
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

                if line.startswith('Sel_'):
                    if line.split('=')[1].strip() != '':
                        if len(u.select_atoms(f"{line.split('=')[1].strip()}")) != 0:
                            self.selection_list.append(line.split('=')[1].strip())
                        else:
                            raise ValueError(
                                "You atom pointed to 0 atoms: please check your selection with your structure file")

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

    def argumentParser(self):

        ap = argparse.ArgumentParser()
        ap.add_argument("-m", '--mode', type=str, default='parallel', required=False,
                        help="specify -m parallel or serial mode [Default = parallel]")
        ap.add_argument('-e', '--exclude', nargs='*', required=False,
                        help=' use -e to exclude a list of GPUs from being used by mwSuMD: e.g. -e 0 3')
        ap.add_argument('-c', '--command', type=str, nargs='?', required=False,
                        help=' use -c to define a specific command you want to use to run your engine.'
                             'Use "" to define the command: "gmx mdrun -deffnm npt -bonded gpu". '
                             'Be aware of the GPU batch division and let mwSuMD sort the GPUs.')
        args = ap.parse_args()

        if 'parallel' in vars(args).values():
            self.initialParameters['Mode'] = 'parallel'
        elif 'serial' in vars(args).values():
            self.initialParameters['Mode'] = 'serial'

        self.initialParameters['COMMAND'] = None
        if args.command is not None:
            self.initialParameters['COMMAND'] = args.command
            print(self.initialParameters['COMMAND'])

        self.initialParameters['EXCLUDED_GPUS'] = None
        if args.exclude is not None and len(args.exclude) != 0:
            self.initialParameters['EXCLUDED_GPUS'] = [int(x) for x in args.exclude]

        if not os.path.isdir(f'{self.folder}/plumed'):
            self.initialParameters['PLUMED'] = None
        else:
            for file in os.listdir(self.folder + "/plumed"):
                if file.endswith('.inp'):
                    self.initialParameters['PLUMED'] = f'{self.folder}/plumed/{file}'

        with open(self.inputFile, "r") as infile:
            for line in infile:
                if line.startswith('#'):
                    continue

                if line.startswith('Restart'):
                    self.initialParameters['Restart'] = line.split('=')[1].strip().upper()

                if line.startswith('Output'):
                    self.initialParameters['Output'] = line.split('=')[1].strip()

                if line.startswith('ligand_HB'):
                    self.initialParameters['ligand_HB'] = line.split('=')[1].strip()

    def getRestartOutput(self):
        os.makedirs('restarts', exist_ok=True)
        if self.trajCount != 0 or self.initialParameters['Restart'] == 'YES':
            directory = f'{self.folder}/restarts'
        else:
            directory = f'{self.folder}/system'
        self.initialParameters.update(
            {extension: file for file in os.listdir(directory) for extension in self.outExtensions if
             file.endswith(extension)})

    def countTraj_logTraj(self, metric):
        """ At what cycle number mwSuMD was stopped? """
        self.trajCount = len([x for x in os.scandir(f'{self.folder}/trajectories')])
        if self.trajCount == 0 and (metric == 0 or metric == 10 ** 6):
            with open('walkerSummary.log', 'w') as logF:
                logF.write('#' * 5 + " Simulation Starts " + '#' * 5 + "\n")
        elif self.trajCount != 0 and (metric == 0 or metric == 10 ** 6):
            with open('walkerSummary.log', 'a') as logF:
                logF.write(str(self.trajCount) + " RESUMED " + str(metric) + "\n")
        elif self.initialParameters['Relax'] is True:
            with open('walkerSummary.log', 'a') as logF:
                logF.write(str(self.trajCount) + " RELAXATION PROTOCOL  " + str(metric) + "\n")
        else:
            with open('walkerSummary.log', 'a') as logF:
                logF.write(str(self.trajCount) + " " + str(metric) + "\n")

    def getSettings(self):
        print("Loading setting parameters...")
        self.checkEngine(), self.getParameters(), self.getForcefields()
        self.getSettingsFromInputFile(), self.getReferencePDB(), self.argumentParser(), self.getRestartOutput()
        return self.initialParameters, self.selection_list, self.parameterFolderPath
