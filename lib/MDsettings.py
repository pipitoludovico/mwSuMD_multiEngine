import os

from .EngineInputParser import EngineInputs
from .Parser import mwInputParser


class MDsetter(mwInputParser):
    def __init__(self, par):
        super(MDsetter, self).__init__()
        self.par = par
        self.walkers_snaphot = self.par['Walkers']
        self.cycle = len(os.listdir(f'{self.folder}/trajectories'))

    def createInputFile(self):
        if self.par['Relax'] is True:
            self.par['Walkers'] = 1

        for walker in range(1, self.par['Walkers'] + 1):
            os.makedirs('tmp', exist_ok=True)
            os.makedirs('tmp/walker_' + str(walker), exist_ok=True)
            # we check if user-defined production is in the system folder and use it if present
            for file in os.listdir(f'{self.folder}/system'):
                if file.startswith('production'):
                    os.system(f'cp {self.folder}/system/{file} tmp/walker_{walker}')
                else:
                    if self.par['MDEngine'] == 'GROMACS':
                        mw_file = open(f'tmp/walker_{walker}/input_{walker}_{self.cycle}.mdp', 'w')
                    elif self.par['MDEngine'] == 'NAMD':
                        mw_file = open(f'tmp/walker_{walker}/input_{walker}_{self.cycle}.namd', 'w')
                    else:
                        mw_file = open(f'tmp/walker_{walker}/input_{walker}_{self.cycle}.inp', 'w')
                    # we get the inputfile based on the engine and fill the mw_file            
                    txt = EngineInputs(self.par).getInputFile()

                    if self.par['MDEngine'] != 'GROMACS':
                        if self.par['Parameters'] is not None:
                            for e in self.par['Parameters']:
                                mw_file.write(f'parameters		{mwInputParser.parameterFolderPath}/%s\n' % e)
    
                    for line in txt:
                        mw_file.write(line)
    
                    if self.par['PLUMED'] is not None:
                        mw_file.write('plumedFile		%s\n' % self.par['PLUMED'])
                    mw_file.close()
            # if in relaxation mode, we reset the number of walkers back to original
        self.par['Walkers'] = self.walkers_snaphot
        self.cycle += 1
