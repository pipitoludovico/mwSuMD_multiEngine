from .FolderOps import FolderOps
from .mwParser import mwInputParser


class MDsetter(mwInputParser):
    def __init__(self):
        super(MDsetter, self).__init__()

    def createInput(self, walker_number, trajCount):
        savefreq = int((int(self.par['Savefreq']) * 1000) / int(self.par['Timestep']))
        eq_file = open(f'tmp/walker_{walker_number}/input_{walker_number}_{trajCount}.inp', 'w')

        if self.par['Forcefield'] == 'CHARMM':
            txt = FolderOps().getAcemdInputFile(self.par, savefreq, trajCount)

            for line in txt:
                eq_file.write(line)

            if self.par['Parameters'] is not None:
                for e in self.par['Parameters']:
                    eq_file.write(f'parameters		{mwInputParser.parPath}/%s\n' % e)

            if self.par['PLUMED'] is not None:
                eq_file.write('plumedFile		%s\n' % self.par['PLUMED'])

            eq_file.close()

        if self.par['Forcefield'] == 'AMBER':
            txt = FolderOps.getAMBERinputFile(self.par, savefreq)

            for e in txt:
                eq_file.write(e)

            if self.par['PLUMED'] is not None:
                eq_file.write('plumedFile		%s\n' % self.par['PLUMED'])

            eq_file.close()
