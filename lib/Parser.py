import argparse
import os


class mwInputParser:
    folder = os.getcwd()
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
        self.groExtensions = ('.mpd', '.gro', '.cpt', '.itp', 'top')
        self.paramExt = ('.param', '.prmtop', '.prm')
        os.makedirs(f'{self.folder}/trajectories', exist_ok=True)
        self.trajCount = len([x for x in os.scandir(f'{self.folder}/trajectories')])
        if not os.path.isfile(f'{self.folder}/{self.inputFile}'):
            print('Input file for SuMD simulation required')
            quit()

    def checkEngine(self):
        engineCheck = ('psf', 'prmtop')
        self.par = {'MDEngine': 'GROMACS' if file.endswith(engineCheck) else 'ACEMD' for file in
                    os.listdir(f'{self.folder}/system')}
        if self.par['MDEngine'] == 'GROMACS':
            self.getGromacsSystem()

    def getTopology(self):
        self.par['PSF'] = None
        self.par['PDB'] = None
        for topos in os.listdir(f'{self.folder}/system'):
            if topos.endswith('.psf'):
                self.par['PSF'] = topos
            if topos.endswith('.pdb'):
                self.par['PDB'] = topos

    def getGromacsSystem(self):
        for groFile in os.listdir('./system'):
            if groFile.endswith(self.groExtensions):
                self.par[f'GROMACS_{str(self.groExtensions).replace(".", "")}'] = groFile

    def getParameters(self):
        if 'Parameters' not in self.par:
            self.par['Parameters'] = []
            for params in os.listdir(f'{self.folder}/system'):
                if params.endswith(self.paramExt):
                    self.par['Parameters'].append(params)
            for dirpath, dirnames, generalParams in os.walk(self.parPath):
                for filename in [f for f in generalParams if f.endswith(self.paramExt)]:
                    self.par['Parameters'].append(filename)
        else:
            return self.par['Parameters']

    def getReferencePDB(self):
        if self.par['Metric_1'] or self.par['Metric_2'] == 'RMSD':
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
            for line in infile:
                if line.startswith('#'):
                    continue
                if line.startswith('RelaxTime'):
                    self.par['RelaxTime'] = int(line.split('=')[1].strip())

                if line.startswith('NumberCV'):
                    self.par['NumberCV'] = int(line.split('=')[1].strip())

                if line.startswith('Metric_1'):
                    self.par['Metric_1'] = line.split('=')[1].strip().upper()

                if line.startswith('Metric_2'):
                    self.par['Metric_2'] = line.split('=')[1].strip().upper()

                if line.startswith('Cutoff_1'):
                    self.par['Cutoff_1'] = int(line.split('=')[1].strip())

                if line.startswith('Cutoff_2'):
                    self.par['Cutoff_2'] = int(line.split('=')[1].strip())

                if line.startswith('Walkers'):
                    self.par['Walkers'] = int(line.split('=')[1].strip())

                if line.startswith('Timewindow'):
                    self.par['Timewindow'] = line.split('=')[1].strip()

                if line.startswith('Timestep'):
                    self.par['Timestep'] = line.split('=')[1].strip()

                if line.startswith('Savefreq'):
                    self.par['Savefreq'] = line.split('=')[1].strip()

                if line.startswith('Sel_'):
                    self.selection_list.append(line.split('=')[1].strip())

                if line.startswith('Transition_1'):
                    self.par['Transition_1'] = line.split('=')[1].strip().lower()

                if line.startswith('Transition_2'):
                    self.par['Transition_2'] = line.split('=')[1].strip().lower()

                if line.startswith('Wrap'):
                    self.par['Wrap'] = line.split('=')[1].strip()

                if line.startswith('GPU_ID'):  # could go as argv...
                    self.par['GPU_ID'] = line.split('=')[1].strip()

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

        self.par['PLUMED'] = None
        for file in os.listdir(self.folder + "/system"):
            if file.endswith('.inp'):
                self.par['PLUMED'] = f'{self.folder}/system/{file}'

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
        if self.trajCount == 0:
            if self.par['Restart'] == 'YES':
                directory = f'{self.folder}/restarts'
            else:
                directory = f'{self.folder}/system'
            for file in os.listdir(directory):
                for extension in self.outExtensions:
                    if file.endswith(extension):
                        self.par[extension] = file
        if self.trajCount != 0:
            directory = f'{self.folder}/restarts'
            for file in os.listdir(directory):
                for extension in self.outExtensions:
                    if file.endswith(extension):
                        self.par[extension] = file

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
        self.checkEngine(), self.getTopology(), self.getParameters(), self.getForcefields()
        self.getMetrics(), self.getReferencePDB(), self.argumentParser(), self.getRestartOutput()
        return self.par, self.selection_list, self.parPath
