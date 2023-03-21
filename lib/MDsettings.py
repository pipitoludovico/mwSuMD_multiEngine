import os

from .EngineInputParser import EngineInputs
from .Parser import mwInputParser


class MDsetter(mwInputParser):
    def __init__(self, par):
        super(MDsetter, self).__init__()
        self.par = par
        self.walkers_snaphot = self.par['Walkers']
        self.cycle = len(os.listdir(f'{self.folder}/trajectories'))

    def createACEMDinputFile(self):
        if self.par['Relax'] is True:
            self.par['Walkers'] = 1

        for walker in range(1, self.par['Walkers'] + 1):
            os.makedirs('tmp', exist_ok=True)
            os.makedirs('tmp/walker_' + str(walker), exist_ok=True)
            eq_file = open(f'tmp/walker_{walker}/input_{walker}_{self.cycle}.inp', 'w')

            if self.par['Forcefield'] == 'CHARMM':
                txt = EngineInputs(self.par).getInputFile()

                for line in txt:
                    eq_file.write(line)

                if self.par['Parameters'] is not None:
                    for e in self.par['Parameters']:
                        eq_file.write(f'parameters		{mwInputParser.parPath}/%s\n' % e)

                if self.par['PLUMED'] is not None:
                    eq_file.write('plumedFile		%s\n' % self.par['PLUMED'])

                eq_file.close()

            if self.par['Forcefield'] == 'AMBER':
                txt = EngineInputs(self.par).getInputFile()

                for e in txt:
                    eq_file.write(e)

                if self.par['Parameters'] is not None:
                    for e in self.par['Parameters']:
                        eq_file.write(f'parameters		{mwInputParser.parPath}/%s\n' % e)

                if self.par['PLUMED'] is not None:
                    eq_file.write('plumedFile		%s\n' % self.par['PLUMED'])

                eq_file.close()
        self.par['Walkers'] = self.walkers_snaphot
        self.cycle += 1

    def createGROMACSinputFile(self):
        if self.par['Relax'] is True:
            self.par['Walkers'] = 1
        for walker in range(1, self.par['Walkers'] + 1):
            os.makedirs('tmp', exist_ok=True)
            os.makedirs('tmp/walker_' + str(walker), exist_ok=True)
            eq_file = open(f'tmp/walker_{walker}/input_{walker}_{self.cycle}.inp', 'w')
        self.par['Walkers'] = self.walkers_snaphot
