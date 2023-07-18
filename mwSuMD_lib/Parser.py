import argparse
import os
import subprocess
import signal
import MDAnalysis as Mda


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
    from warnings import filterwarnings
    filterwarnings(action='ignore')

    def __init__(self):
        self.customInputFileExtension = ('namd', 'inp', 'mdp')
        self.outExtensions = ('cpt', 'cpi', 'gro', 'tpr', 'coor', 'vel', 'xsc')
        self.fileExtensions = ('.psf', '.pdb', '.mdp', '.gro', '.cpt', 'top', '.prmtop', '.tpr')
        self.initialParametersameter_extensions = ('.param', '.prmtop', '.prm', '.par')
        self.trajCount = len([x for x in os.scandir(f'{self.folder}/trajectories')])
        self.allowedMetrics = ("DISTANCE", "RMSD", "CONTACTS")
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
                files = os.listdir(f'{self.folder}/system/reference')
                reference = next((file for file in files if file.endswith('.pdb')), None)
                if reference:
                    self.initialParameters['REFERENCE'] = reference
                else:
                    print("Put a reference pdb file in the 'reference' folder inside system and rerun.")
                    exit()

        if self.initialParameters['Metric_1'] == 'RMSD' or self.initialParameters['Metric_2'] == 'RMSD':
            checkRMSDoption()

    def getForcefields(self):
        self.initialParameters['Forcefield'] = 'CHARMM' \
            if any(fileSys.endswith('.psf') for fileSys in os.listdir(f'{self.folder}/system')) else 'AMBER' \
            if any(fileSys.endswith('.prmtop') for fileSys in os.listdir(f'{self.folder}/system')) else 'GROMOS'

    def getSettingsFromInputFile(self):
        u = Mda.Universe(f"{self.folder}/system/{self.initialParameters['PDB']}")
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

                if line.startswith(f'Sel_'):
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

        if not self.initialParameters.get('Metric_1') and not self.initialParameters.get('Metric_2'):
            raise ValueError(
                "Please make sure if you choose at least one metric to supervise (Distance, Contacts, RMSD, HB)")
        if self.initialParameters.get('NumberCV') == 2 and (
                not self.initialParameters.get('Metric_1') or not self.initialParameters.get('Metric_2')):
            raise ValueError(
                "Please make sure if you use CV2 to specify all the CVs choosing one metric to supervise (Distance, Contacts, RMSD, HB)")

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
        ap.add_argument("-k", '--kill', required=False, action='store_true', default=False,
                        help="Stop mwSuMD from the current working directory")
        ap.add_argument("-j", '--join', nargs='*', required=False,
                        help="Merge the trajectories from one step to another: e.g. -j 1 10 or -j all to merge every step.")

        args = ap.parse_args()

        if args.kill is True:
            import os

            os.system('val=$(<.mypid ) && kill -9 $val')
            os.kill(os.getpid(), signal.SIGKILL)
        if args.join is not None:
            from mwSuMD_lib import Merger

            merger = Merger.TrajMerger()
            merger.loadTrajectories()
            if args.join == 'all':
                merger.mergeAll()
            if len(args.join) != 0:
                merger.mergeFrom(args.join[0], args.join[1])
            else:
                print(
                    "Error: incorrect arguments for -j. -join needs 2 number to set the starting and ending steps to be merged")
                exit()
            exit()

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
        self.initialParameters.update({extension: file for file in os.listdir(directory) for extension in self.outExtensions if file.endswith(extension) and 'minimzed' not in file})

    def countTraj_logTraj(self, metric):
        """ At what cycle number mwSuMD was stopped? """
        self.trajCount = len([x for x in os.scandir(f'{self.folder}/trajectories')])
        if self.initialParameters.get('Metric_1') is None:
            self.initialParameters['Metric_1'] = 'No Metric 1 chosen'
        if self.initialParameters.get('Metric_2') is None:
            self.initialParameters['Metric_2'] = 'No Metric 2 chosen'
        if self.trajCount == 0 and (metric == 0 or metric == 10 ** 6):
            with open('walkerSummary.log', 'w') as logF:
                logF.write('#' * 5 + " Simulation Starts " + '#' * 5 + "\n")
        elif self.trajCount != 0 and (metric == 0 or metric == 10 ** 6):
            with open('walkerSummary.log', 'a') as logF:
                logF.write(str(self.trajCount) + "\nCheckpoint with Metric 1:" + (
                    self.initialParameters["Metric_1"]) + "\tMetric 2: " + (self.initialParameters[
                                                                                "Metric_2"] + f" and the following selections: \n{self.selection_list}\n"))
        elif self.initialParameters['Relax'] is True:
            with open('walkerSummary.log', 'a') as logF:
                logF.write(str(self.trajCount) + "\nClassic Protocol with Metric 1\t:" + (
                    self.initialParameters["Metric_1"]) + "\tMetric 2:" + (self.initialParameters[
                                                                               "Metric_2"] + f" and the following selections: \n{self.selection_list}\n"))
        else:
            with open('walkerSummary.log', 'a') as logF:
                logF.write(str(self.trajCount) + " " + str(metric) + "\n")

    def getSettings(self):
        print("Loading setting parameters...")
        self.checkEngine(), self.getParameters(), self.getForcefields()
        self.getSettingsFromInputFile(), self.getReferencePDB(), self.argumentParser(), self.getRestartOutput()
        return self.initialParameters, self.selection_list, self.parameterFolderPath
