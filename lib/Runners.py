import os

from .Parser import mwInputParser


class Runner(mwInputParser):

    def __init__(self, par):
        self.par = par
        super(Runner, self).__init__()

    def runMD(self, GPUbatches):
        if self.par['MDEngine'] == 'ACEMD':
            self.runACEMD(self.trajCount, GPUbatches)
        else:
            self.runGROMACS(self.trajCount, GPUbatches)

    def runACEMD(self, trajCount, GPUbatches):
        for GPUbatch in GPUbatches:
            process = 0
            while process < len(GPUbatch):
                walker_number = process + 1
                os.chdir('tmp/walker_' + str(walker_number))
                command = f'acemd3 --device {GPUbatch[process]} input_{walker_number}_{trajCount}.inp 1> acemd.log'
                process += 1
                print(command)
                try:
                    self.wrap()
                    print("Wrapped successfully")
                except:
                    print(f"Walker {walker_number} wrapping has failed. No results were produced."
                          f"Now moving to the next walker.")
                os.chdir(self.folder)
            process = 0
            print(process)

    def runGROMACS(self, trajCount, GPUbatches):
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
            if file == str(xtc) or file.endswith('.xtc'):
                try:
                    mol = Molecule(psf)
                    mol.read(file)
                    mol.wrap()
                    mol.wrap(self.par['Wrap'])
                    mol.write('wrapped.xtc')
                except:
                    print(f"{file} was used for wrapping")
