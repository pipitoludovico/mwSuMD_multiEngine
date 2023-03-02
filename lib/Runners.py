import os

from .Parser import mwInputParser


class Runner(mwInputParser):

    def __init__(self, par):
        self.par = par
        super(Runner, self).__init__()

    def runMD(self):
        if self.par['MDEngine'] == 'ACEMD':
            self.runACEMD(self.trajCount)
        else:
            self.runGROMACS(self.trajCount)

    def runACEMD(self, trajCount):
        for walker in range(1, self.par['Walkers'] + 1):
            os.chdir('tmp/walker_' + str(walker))
            print(f'acemd3 --device {walker - 1} input_{walker}_{trajCount} 1> acemd.log')
            # os.system(f'acemd3 --device{batch[walker-1]} input_{walker}_{trajCount} 1> acemd.log')
            self.wrap()
            print("Wrapped successfully")
            os.chdir(self.folder)

    def runGROMACS(self, trajCount):
        # to do
        for walker in range(1, self.par['Walkers'] + 1):
            os.chdir('tmp/walker_' + str(walker + 1))
            if self.par['PLUMED'] is not None:
                os.system(f'echo test {[walker - 1]} {trajCount}')

    def wrap(self):
        psf = None
        from moleculekit.molecule import Molecule
        if self.par['Forcefield'] == 'CHARMM':
            psf = '../../system/%s' % self.par['PSF']

        elif self.par['Forcefield'] == 'AMBER':
            psf = '../../system/%s' % self.par['PRMTOP']
        xtc = '%s_%s.xtc' % (self.par['Output'], self.trajCount)
        for file in os.listdir(os.getcwd()):
            if file == str(xtc):
                mol = Molecule(psf)
                mol.read(xtc)
                mol.wrap()
                mol.wrap(self.par['Wrap'])
                mol.write('wrapped.xtc')
