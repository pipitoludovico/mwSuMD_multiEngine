import os

import MDAnalysis as Mda
from MDAnalysis import transformations

from .Parser import mwInputParser


class TrajectoryOperator(mwInputParser):
    def __init__(self):
        super(TrajectoryOperator, self).__init__()

    def wrap(self, folder):
        os.chdir(f'{self.folder}/tmp/walker_' + str(folder))
        print('wrapping in ' + os.getcwd())
        ext = ('xtc', 'dcd')
        psf = None

        if self.initialParameters['Forcefield'] == 'CHARMM':
            if self.initialParameters['PSF'] is not None:
                psf = '../../system/%s' % self.initialParameters['PSF']
        elif self.initialParameters['Forcefield'] == 'AMBER':
            psf = '../../system/%s' % self.initialParameters['PRMTOP']
        elif self.initialParameters['Forcefield'] == 'GROMOS':
            for tpr in os.listdir(os.getcwd()):
                if tpr.startswith(self.initialParameters['Output']) and tpr.endswith('.tpr'):
                    psf = tpr

        for trajectory in os.listdir(os.getcwd()):
            if trajectory.startswith(self.initialParameters['Output']) and trajectory.endswith(ext):
                u = Mda.Universe(psf, trajectory)
                selection = u.select_atoms(f"{self.initialParameters['Wrap']}")
                if len(selection.atoms) == 0:
                    print("your wrapping selection selected 0 atoms! using protein and name CA instead...")
                    selection = u.select_atoms('protein and name CA')
                ag = u.atoms
                workflow = (transformations.unwrap(ag),
                            transformations.center_in_box(selection),
                            transformations.wrap(ag, compound='fragments'))
                u.trajectory.add_transformations(*workflow)

                with Mda.Writer('wrapped.xtc', ag) as w:
                    for ts in u.trajectory:
                        w.write(ag)
        os.chdir(self.folder)
