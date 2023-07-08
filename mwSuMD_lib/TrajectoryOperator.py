import os
import signal
import traceback
import MDAnalysis as Mda
from MDAnalysis import transformations

from .Parser import mwInputParser


class TrajectoryOperator(mwInputParser):
    def __init__(self):
        super(mwInputParser, self).__init__()
        self.setting_error = f"File or setting missing. "

    def wrap(self, folder):
        try:
            os.chdir(f'{self.folder}/tmp/walker_' + str(folder))
            ext = ('xtc', 'dcd')
            trajFile = None
            psf = None

            if self.initialParameters['Forcefield'] == 'CHARMM':
                if self.initialParameters['PSF'] is None:
                    print(self.setting_error)
                    os.kill(os.getpid(), signal.SIGKILL)
                    raise FileNotFoundError
                else:
                    psf = '../../system/%s' % self.initialParameters['PSF']
            if self.initialParameters['Forcefield'] == 'AMBER':
                if self.initialParameters['PRMTOP'] is None:
                    print(self.setting_error)
                    os.kill(os.getpid(), signal.SIGKILL)
                    raise FileNotFoundError
                else:
                    psf = '../../system/%s' % self.initialParameters['PRMTOP']
            if self.initialParameters['Forcefield'] == 'GROMOS':
                for new_coords in os.listdir(os.getcwd()):
                    if new_coords.endswith('.tpr'):
                        psf = new_coords

            for trajectory in os.listdir(os.getcwd()):
                if trajectory.endswith(ext):
                    trajFile = trajectory

            try:
                u = Mda.Universe(psf, trajFile)
            except:
                print("No trajectory or psf found. Check your simulation parameters and make sure production went well")
                os.kill(os.getpid(), signal.SIGKILL)

            selection = u.select_atoms(f"{self.initialParameters['Wrap']}")
            if len(selection.atoms) == 0:
                print("your wrapping selection selected 0 atoms! using protein and name CA instead...")
                selection = u.select_atoms('protein and name CA')
            ag = u.atoms
            workflow = (transformations.unwrap(ag),
                        transformations.center_in_box(selection),
                        transformations.wrap(ag, compound="fragments"))
            u.trajectory.add_transformations(*workflow)
            with Mda.Writer('wrapped.xtc', ag) as w:
                for ts in u.trajectory:
                    w.write(ag)
            os.chdir(self.folder)
        except:
            print("Wrapping failed: check your MD results in tmp/walker_x folders")
            os.kill(os.getpid(), signal.SIGKILL)
            raise Exception(traceback.format_exc())
