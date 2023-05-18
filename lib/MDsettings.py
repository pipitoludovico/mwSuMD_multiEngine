import os

from .EngineInputParser import EngineInputs
from .Parser import mwInputParser


class MDsetter(mwInputParser):
    def __init__(self, par):
        super(MDsetter, self).__init__()
        self.setterParameters = par
        self.walkers_snaphot = self.initialParameters['Walkers']
        self.cycle = len(os.listdir(f'{self.folder}/trajectories'))
        self.custom_input = None

    def createInputFile(self):
        if self.initialParameters['Relax'] is True:
            self.setterParameters['Walkers'] = 1

        for walker in range(1, self.setterParameters['Walkers'] + 1):
            os.makedirs('tmp', exist_ok=True)
            os.makedirs('tmp/walker_' + str(walker), exist_ok=True)
            # we check if user-defined production is in the system folder and use it if present
            self.custom_input = self.initialParameters['CUSTOMFILE']

            if self.custom_input is None:
                if self.setterParameters['MDEngine'] == 'GROMACS':
                    mw_file = open(f'tmp/walker_{walker}/input_{walker}_{self.cycle}.mdp', 'w')
                elif self.setterParameters['MDEngine'] == 'NAMD':
                    mw_file = open(f'tmp/walker_{walker}/input_{walker}_{self.cycle}.namd', 'w')
                else:
                    mw_file = open(f'tmp/walker_{walker}/input_{walker}_{self.cycle}.inp', 'w')
                # we get the inputfile based on the engine and fill the mw_file
                txt = EngineInputs(self.setterParameters).getInputFile()

                if self.setterParameters['MDEngine'] != 'GROMACS':
                    # adding parameters lines to input files != GROMACS
                    if self.setterParameters['Parameters'] is not None:
                        for e in self.setterParameters['Parameters']:
                            mw_file.write(f'parameters		%s\n' % e)
                # writing the input file
                for line in txt:
                    mw_file.write(line)
                if self.setterParameters['PLUMED'] is not None:
                    mw_file.write('plumedFile		%s\n' % self.setterParameters['PLUMED'])
                mw_file.close()
                # if in relaxation mode, we reset the number of walkers back to original
            else:
                os.system(f'cp {self.initialParameters["CUSTOMFILE"]} {self.initialParameters["Root"]}/tmp/walker_{walker}')
        self.initialParameters['Walkers'] = self.walkers_snaphot
        self.cycle += 1
