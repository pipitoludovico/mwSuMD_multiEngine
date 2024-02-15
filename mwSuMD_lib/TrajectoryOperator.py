import os
import MDAnalysis as Mda
from MDAnalysis import transformations

from .Parser import mwInputParser
from .Loggers import Logger
from warnings import filterwarnings

filterwarnings(action='ignore')


class TrajectoryOperator(mwInputParser):
    def __init__(self):
        super(mwInputParser, self).__init__()
        self.setting_error = f"File or setting missing. "
        self.trajCount = len([traj for traj in os.listdir(f'{self.folder}/trajectories') if traj.endswith('.xtc')])

    def wrap(self, folder):
        try:
            os.chdir(f'{self.folder}/tmp/walker_' + str(folder))
            ext = ('xtc', 'dcd')
            trajFile = None
            psf = None

            if self.initialParameters['Forcefield'] == 'CHARMM':
                if self.initialParameters['PSF'] is None:
                    Logger.LogToFile('a', self.trajCount, self.setting_error)
                    raise FileNotFoundError
                else:
                    psf = '../../system/%s' % self.initialParameters['PSF']
            if self.initialParameters['Forcefield'] == 'AMBER':
                if self.initialParameters['PRMTOP'] is None:
                    Logger.LogToFile('a', self.trajCount, self.setting_error)
                    raise FileNotFoundError
                else:
                    psf = '../../system/%s' % self.initialParameters['PRMTOP']
            if self.initialParameters['Forcefield'] == 'GROMOS':
                for new_coords in os.listdir(os.getcwd()):
                    if new_coords.endswith('.tpr'):
                        psf = new_coords
                    else:
                        psf = "../../system/%s" % self.initialParameters.get('TPR')
            for trajectory in os.listdir(os.getcwd()):
                if trajectory.endswith(ext):
                    trajFile = trajectory

            if trajFile.endswith('.dcd'):
                conv = Mda.Universe(psf, trajFile)
                all_atoms = conv.select_atoms('all')
                with Mda.Writer('converted.xtc') as converter:
                    for ts in conv.trajectory:
                        if ts:
                            converter.write(all_atoms)
                trajFile = 'converted.xtc'

            try:
                u = Mda.Universe(psf, trajFile)
            except:
                Logger.LogToFile('a', self.trajCount,
                                 "No trajectory or psf found. Check your simulation parameters and make sure production went well")
                raise FileNotFoundError

            selection = u.select_atoms(f"{self.initialParameters['Wrap']}")
            if len(selection.atoms) == 0:
                Logger.LogToFile('a', self.trajCount,
                                 "your wrapping selection selected 0 atoms! using protein and name CA instead...")
                selection = u.select_atoms('protein and name CA')
            ag = u.atoms
            workflow = (transformations.unwrap(ag),
                        transformations.center_in_box(selection),
                        transformations.wrap(ag, compound="fragments"))
            u.trajectory.add_transformations(*workflow)
            with Mda.Writer('wrapped.xtc', ag) as w:
                for ts in u.trajectory:
                    if ts:
                        w.write(ag)
            os.chdir(self.folder)
        except Exception as e:
            with open('wrappingLog.txt', 'a') as wrapLog:
                wrapLog.write(f"Wrapping Exception in {os.getcwd()}. Exception: {e} ")
            Logger.LogToFile('a', self.trajCount, "Wrapping failed: check your MD results in tmp/walker_x folders")
            raise FileNotFoundError
