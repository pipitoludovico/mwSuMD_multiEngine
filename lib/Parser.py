import argparse
import os


class mwInputParser:
    folder = os.getcwd()
    os.makedirs(f'{folder}/trajectories', exist_ok=True)
    inputFile = 'simulation_settings_mwSuMD.inp'
    par = {}
    selection_list = []
    walker_metrics = []
    parPath = os.path.abspath('parameters')
    new_value = 0
    max_value = 0
    metric_1 = 0
    metric_2 = 0
    mode = 'parallel'
    excludedGPUS = []

    def __init__(self):
        self.outExtensions = ('coor', 'vel', 'xsc')
        self.fileExtensions = ('.psf', '.pdb', '.mdp', '.gro', '.cpt', '.itp', 'top', '.prmtop')
        self.charmmParamExt = ('.param', '.prmtop', '.prm')
        self.trajCount = len([x for x in os.scandir(f'{self.folder}/trajectories')])
        if not os.path.isfile(f'{self.folder}/{self.inputFile}'):
            print('Input file for SuMD simulation required')
            quit()

    def checkEngine(self):
        self.par['MDEngine'] = 'GROMACS' if any(file.endswith('.gro') for file in os.listdir('./system')) else 'ACEMD'
        self.getSystem()

    def getSystem(self):
        for ext in self.fileExtensions:
            for file in os.listdir('./system'):
                if file.endswith(ext):
                    self.par[ext.replace('.', '').upper()] = file

    def getParameters(self):
        if 'Parameters' not in self.par:
            self.par['Parameters'] = []
            for params in os.listdir(f'{self.folder}/system'):
                if params.endswith(self.charmmParamExt):
                    self.par['Parameters'].append(params)
            for dirpath, dirnames, generalParams in os.walk(self.parPath):
                for filename in [f for f in generalParams if f.endswith(self.charmmParamExt)]:
                    self.par['Parameters'].append(filename)
        else:
            return self.par['Parameters']

    def getReferencePDB(self):
        if self.par['Metric_1'] or self.par['Metric_2'] == 'RMSD':
            if not os.path.isdir(f'{self.folder}/system/reference'):
                print('You need a reference folder with a reference pdb in it if you use the RMSD as a metric.')
                exit()
            if len(os.listdir(f'{self.folder}/system/reference')) > 0:
                for reference in os.listdir(f'{self.folder}/system/reference'):
                    if reference.endswith('.pdb'):
                        self.par['REFERENCE'] = reference
            else:
                print("Put a reference pdb file in the 'reference' folder inside system and rerun.")
                exit()

    def getForcefields(self):
        self.par['Forcefield'] = 'CHARMM'
        for fileSys in os.listdir(f'{self.folder}/system'):
            if fileSys.endswith('.psf'):
                self.par['Forcefield'] = 'CHARMM'
            if fileSys.endswith('.prmtop'):
                self.par['Forcefield'] = 'AMBER'

    def getMetrics(self):
        # Default settings:
        self.par['Timestep'] = 2
        self.par['Savefreq'] = 20
        self.par['Wrap'] = 'protein and name CA'

        with open(self.inputFile, "r") as infile:
            self.par['RelaxTime'] = 5
            self.par['Relax'] = False
            for line in infile:
                if line.startswith('#'):
                    continue
                if line.startswith('RelaxTime'):
                    if line.split('=')[1].strip() != '':
                        self.par['RelaxTime'] = float(line.split('=')[1].strip())

                if line.startswith('NumberCV'):
                    self.par['NumberCV'] = int(line.split('=')[1].strip())

                if line.startswith('Metric_1'):
                    self.par['Metric_1'] = line.split('=')[1].strip().upper()

                if line.startswith('Metric_2'):
                    self.par['Metric_2'] = line.split('=')[1].strip().upper()

                if line.startswith('Cutoff_1'):
                    self.par['Cutoff_1'] = float(line.split('=')[1].strip())

                if line.startswith('Cutoff_2'):
                    self.par['Cutoff_2'] = float(line.split('=')[1].strip())

                if line.startswith('Walkers'):
                    self.par['Walkers'] = int(line.split('=')[1].strip())

                if line.startswith('Timewindow'):
                    self.par['Timewindow'] = int(line.split('=')[1].strip())

                if line.startswith('Timestep'):
                    self.par['Timestep'] = int(line.split('=')[1].strip())

                if line.startswith('Savefreq'):
                    self.par['Savefreq'] = int(line.split('=')[1].strip())

                if line.startswith('Sel_'):
                    self.selection_list.append(line.split('=')[1].strip())

                if line.startswith('Transition_1'):
                    self.par['Transition_1'] = line.split('=')[1].strip().lower()

                if line.startswith('Transition_2'):
                    self.par['Transition_2'] = line.split('=')[1].strip().lower()

                if line.startswith('Wrap'):
                    self.par['Wrap'] = line.split('=')[1].strip()

    def argumentParser(self):
        ap = argparse.ArgumentParser()
        ap.add_argument("-m", '--mode', type=str, default='parallel', required=False,
                        help="specify -m parallel or serial mode [Default = parallel]")
        ap.add_argument('-e', '--exclude', nargs='*', required=False,
                        help=' use -e to exclude a list of GPUs from being used by mwSuMD: e.g. -e 0 3')
        args = ap.parse_args()
        if 'parallel' in vars(args).values():
            self.mode = 'parallel'
        else:
            self.mode = 'serial'
        if args.exclude is not None and len(args.exclude) != 0:
            self.excludedGPUS = [x for x in args.exclude]

        if not os.path.isdir(f'{self.folder}/plumed'):
            self.par['PLUMED'] = None
        else:
            for file in os.listdir(self.folder + "/plumed"):
                if file.endswith('.inp'):
                    self.par['PLUMED'] = f'{self.folder}/plumed/{file}'

        with open(self.inputFile, "r") as infile:
            for line in infile:
                if line.startswith('#'):
                    continue

                if line.startswith('Restart'):
                    self.par['Restart'] = line.split('=')[1].strip().upper()

                if line.startswith('Output'):
                    self.par['Output'] = line.split('=')[1].strip()

                if line.startswith('ligand_HB'):
                    self.par['ligand_HB'] = line.split('=')[1].strip()

    def getRestartOutput(self):
        os.makedirs('restarts', exist_ok=True)
        if self.trajCount != 0 or self.par['Restart'] == 'YES':
            directory = f'{self.folder}/restarts'
        else:
            directory = f'{self.folder}/system'

        self.par.update({extension: file for file in os.listdir(directory) for extension in self.outExtensions if
                         file.endswith(extension)})

    def countTraj_logTraj(self, metric):
        """ At what cycle number mwSuMD was stopped? """
        self.trajCount = len([x for x in os.scandir(f'{self.folder}/trajectories')])
        if self.par['Restart'] == 'NO':
            if self.trajCount == 0 and (metric == 0 or metric == 10 ** 6):
                with open('walkerSummary.log', 'w') as logF:
                    logF.write('#' * 5 + " Simulation Starts " + '#' * 5 + "\n")
            elif self.trajCount != 0 and (metric == 0 or metric == 10 ** 6):
                with open('walkerSummary.log', 'a') as logF:
                    metric = 'RESUMED'
                    logF.write(str(self.trajCount) + " " + str(metric) + "\n")
            else:
                with open('walkerSummary.log', 'a') as logF:
                    logF.write(str(self.trajCount) + " " + str(metric) + "\n")
        if self.par['Relax'] is True:
            with open('walkerSummary.log', 'a') as logF:
                logF.write(str(self.trajCount) + " RELAXATION PROTOCOL  " + str(metric) + "\n")
        if self.par['Restart'] == 'YES':
            if self.trajCount != 0 and (metric == 0 or metric == 10 ** 6):
                metric = 'RESTART'
                with open('walkerSummary.log', 'a') as logF:
                    logF.write(str(self.trajCount) + " " + str(metric) + "\n")
            else:
                with open('walkerSummary.log', 'a') as logF:
                    logF.write(str(self.trajCount) + " " + str(metric) + "\n")

    def getSettings(self):
        print("Loading setting parameters...")
        self.checkEngine(), self.getParameters(), self.getForcefields()
        self.getMetrics(), self.getReferencePDB(), self.argumentParser(), self.getRestartOutput()
        return self.par, self.selection_list, self.parPath
