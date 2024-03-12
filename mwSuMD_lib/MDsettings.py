import os

from .EngineInputParser import EngineInputs
from .Parser import mwInputParser
from warnings import filterwarnings

filterwarnings(action='ignore')


class MDsetter(mwInputParser):
    def __init__(self, par):
        super(mwInputParser, self).__init__()
        self.setterParameters = par
        self.walkers_snaphot = self.initialParameters['Walkers']
        self.cycle = len([trajFile for trajFile in os.listdir(f'{self.folder}/trajectories') if trajFile.endswith('xtc')])
        self.custom_input = None

    def createInputFile(self):
        ext = ".mpd" if self.setterParameters['MDEngine'] == 'GROMACS' else ('.namd' if self.setterParameters['MDEngine'] == 'NAMD' else '.inp')
        if self.initialParameters['Relax'] is True:
            self.setterParameters['Walkers'] = 1

        for walker in range(1, self.setterParameters['Walkers'] + 1):
            os.makedirs('tmp', exist_ok=True)
            os.makedirs('tmp/walker_' + str(walker), exist_ok=True)
            # we check if user-defined production is in the system folder and use it if present

            self.custom_input = self.initialParameters['CUSTOMFILE']
            if not self.custom_input:
                txt = EngineInputs(self.setterParameters).getInputFile()
                with open(f'tmp/walker_{walker}/input_{walker}_{self.cycle}{ext}', 'w') as mw_file:
                    # we get the inputfile based on the engine and fill the mw_file
                    if self.setterParameters['MDEngine'] == 'ACEMD' or self.setterParameters['MDEngine'] == 'NAMD':
                        # adding parameters lines to input files != GROMACS
                        if self.setterParameters['Parameters'] is not None:
                            for e in self.setterParameters['Parameters']:
                                if e.endswith('.par') or e.endswith('.prm'):
                                    mw_file.write(f'parameters		%s\n' % e)
                    # writing the input file
                    for line in txt:
                        mw_file.write(line)
                    if self.setterParameters['PLUMED'] is not None:
                        mw_file.write('plumedFile		%s\n' % self.setterParameters['PLUMED'])
            else:
                os.system(f'cp {self.initialParameters["CUSTOMFILE"]} {self.initialParameters["Root"]}/tmp/walker_{walker}')
        self.initialParameters['Walkers'] = self.walkers_snaphot
        self.cycle += 1
