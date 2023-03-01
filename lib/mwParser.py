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
    slope = 0

    def __init__(self):
        self.trajCount = 0
        self.outExtensions = ('coor', 'vel', 'xsc')
        self.groExtensions = ('.mpd', '.gro', '.cpt', '.itp', 'top')
        self.paramExt = ('.param', '.prmtop', '.prm')
        os.makedirs('trajectories', exist_ok=True)
        self.getSettings()
        if not os.path.isfile(self.inputFile):
            print('Input file for SuMD simulation required')
            quit()

    def getRoot(self):
        return self.folder

    def checkEngine(self):
        engineCheck = ('psf', 'prmtop')
        self.par = {'MDEngine': 'GROMACS' if file.endswith(engineCheck) else 'ACEMD' for file in os.listdir('system')}
        if self.par['MDEngine'] == 'GROMACS':
            self.getGromacsSystem()

    def getTopology(self):
        self.par['PSF'] = None
        self.par['PDB'] = None
        for topos in os.listdir('./system'):
            if topos.endswith('.psf'):
                self.par['PSF'] = topos
            if topos.endswith('.pdb'):
                self.par['PDB'] = topos

    def getGromacsSystem(self):
        for groFile in os.listdir('./system'):
            if groFile.endswith(self.groExtensions):
                self.par[f'GROMACS_{str(self.groExtensions).replace(".", "")}'] = groFile

    def getParameters(self):
        self.par['Parameters'] = []
        for params in os.listdir('./system'):
            if params.endswith(self.paramExt):
                self.par['Parameters'].append(params)
        for dirpath, dirnames, generalParams in os.walk(self.parPath):
            for filename in [f for f in generalParams if f.endswith(self.paramExt)]:
                self.par['Parameters'].append(filename)

    def getReferencePDB(self):
        if self.par['Metric_1'] or self.par['Metric_2'] == 'RMSD':
            if len(os.listdir('./system/reference')) > 0:
                for reference in os.listdir('./system/reference'):
                    if reference.endswith('.pdb'):
                        self.par['REFERENCE'] = reference
            else:
                print("Put a reference pdb file in the 'reference' folder inside system and rerun.")
                exit()

    def getForcefields(self):
        self.par['Forcefield'] = 'CHARMM'
        for fileSys in os.listdir('./system'):
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
            for line in infile:
                if line.startswith('#'):
                    continue
                if line.startswith('NumberCV'):
                    self.par['NumberCV'] = int(line.split('=')[1].strip())

                if line.startswith('Metric_1'):
                    self.par['Metric_1'] = line.split('=')[1].strip().upper()

                if line.startswith('Metric_2'):
                    self.par['Metric_2'] = line.split('=')[1].strip()
                    if len(self.par['Metric_2']) == 0:
                        self.par['Metric_2'] = None

                if line.startswith('Cutoff_1'):
                    self.par['Cutoff_1'] = int(line.split('=')[1].strip())

                if line.startswith('Cutoff_2'):
                    self.par['Cutoff_2'] = int(line.split('=')[1].strip())

                if line.startswith('Slope'):
                    self.par['Slope'] = line.split('=')[1].strip().upper()

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
        self.par['PLUMED'] = None
        with open(self.inputFile, "r") as infile:
            for line in infile:
                if line.startswith('#'):
                    continue

                if line.startswith('Restart'):
                    self.par['Restart'] = line.split('=')[1].strip()

                if line.startswith('Output'):
                    self.par['Output'] = line.split('=')[1].strip()

                if line.startswith('PLUMED'):
                    if line.split('=')[1].strip() == '':
                        continue
                    else:
                        self.par['PLUMED'] = line.split("=")[1].strip()
                        self.par['PLUMED'] = os.path.abspath(self.par['PLUMED'])

                if line.startswith('ligand_HB'):
                    self.par['ligand_HB'] = line.split('=')[1].strip()
                    if len(self.par['ligand_HB']) == 0:
                        self.par['ligand_HB'] = None

    def getRestartOutput(self):
        os.makedirs('restarts', exist_ok=True)
        if self.par['Restart'] == 'YES':
            directory = './restarts'
        else:
            directory = './system'
        for file in os.listdir(directory):
            for extension in self.outExtensions:
                if file.endswith(extension):
                    self.par[extension] = file

    def countTraj_logTraj(self, metric=0, slope=0):
        """ At what cycle number mwSuMD was stopped? """
        self.trajCount = len([x for x in os.scandir('trajectories')])
        if self.trajCount == 0:
            with open('walkerSummary.log', 'w') as logF:
                logF.write('#' * 5 + " Simulation Starts " + '#' * 5 + "\n")
        else:
            with open('walkerSummary.log', 'a') as logF:
                logF.write(str(self.trajCount) + " " + str(metric) + " " + str(slope) + "\n")
                logF.close()
        return self.trajCount

    def getSettings(self):
        self.checkEngine(), self.getTopology(), self.getParameters(), self.getForcefields()
        self.getMetrics(), self.getReferencePDB(), self.argumentParser(), self.getRestartOutput()
        self.countTraj_logTraj()
        logF = open("settings.txt", "w")
        for idx, sel in enumerate(self.selection_list):
            logF.write('Metric_%s	%s\n' % (str(idx), sel))
        print(self.par, file=logF)
        logF.close()
        return self.par, self.selection_list, self.parPath
