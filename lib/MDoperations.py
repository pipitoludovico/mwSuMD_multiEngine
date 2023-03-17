import os

import numpy as np


class MDoperator:
    def __init__(self, par, root, trajCount):
        self.par = par
        self.folder = root
        self.trajCount = trajCount

    def saveStep(self, best_walker):
        _best_walker = best_walker
        if self.par['Relax'] is True:
            best_walker = 1
        else:
            best_walker = _best_walker
        os.makedirs('trajectories', exist_ok=True)
        """Handle the restart files and the xtc storage"""
        for r in range(1, int(self.par['Walkers']) + 1):
            if r == best_walker:
                print("Renaming the restart files and moving the best walker's restart bin files")
                os.chdir('tmp/walker_%s' % str(r))
                # moving the best frame to the trajectory folder
                os.system(f'cp wrapped.xtc {self.folder}/trajectories/step_{self.trajCount}.xtc')
                # moving and renaming the binary files to the restart folder
                if self.par['MDEngine'] == 'ACEMD':
                    os.system(f'cp *.coor {self.folder}/restarts/previous.coor')
                    os.system(f'cp *.xsc {self.folder}/restarts/previous.xsc')
                    os.system(f'cp *.vel {self.folder}/restarts/previous.vel')

                if self.par['MDEngine'] == 'GROMACS':
                    os.system(f'mv *.gro {self.folder}/restarts/previous.coor')
                    os.system(f'mv *.cpt {self.folder}/restarts/previous.cpt')

                if self.par['PLUMED'] is not None:
                    os.system('cp HILLS  %s/restarts/ ' % self.folder)
                    os.system('cp COLVAR  %s/restarts/ ' % self.folder)
                    os.system('cp grid.dat  %s/restarts/ ' % self.folder)
        print("FINISHED SAVING FRAMES")
        os.chdir(self.folder)
        os.system('rm -r tmp')


def checkIfStuck(values, fails):
    if fails > 5:
        with open('walkerSummary.log') as logFile:
            logFile.write('\nSimulation seemed stuck and it has been terminated')
        exit()
    else:
        # we check if the simulation did not change much,
        # and we run the relaxation protocol
        if np.std(values) < 0.6:
            print('\nSimulation might be stuck. Running the relaxation protocol.')
            return True
        else:
            # if it did not, we continue as usual
            return False
