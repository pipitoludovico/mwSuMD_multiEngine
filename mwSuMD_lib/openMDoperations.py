import os

import numpy as np


class MDoperator:
    def __init__(self, par, root):
        self.par = par
        self.folder = root
        self.cycle = len(os.listdir(f'{self.folder}/trajectories'))

    def saveStep(self, best_walker):
        os.makedirs('trajectories', exist_ok=True)
        """Handle the restart files and the xtc storage"""
        for r in range(1, int(self.par['Walkers']) + 1):
            if r == best_walker:
                print("Renaming the restart files and moving the best walker's restart bin files")
                os.chdir('tmp/walker_%s' % str(r))
                os.system(f'cp wrapped.xtc {self.folder}/trajectories/{self.par["Output"]}_step_{self.cycle}.xtc')
                # moving and renaming the binary files to the restart folder
                os.system(f'cp *.chk {self.folder}/restarts/previous.chk')
                os.system(f'cp *.xml {self.folder}/restarts/previous.xml')

                if self.par['PLUMED'] is not None:
                    os.system('cp HILLS  %s/restarts/ ' % self.folder)
                    os.system('cp COLVAR  %s/restarts/ ' % self.folder)
                    os.system('cp grid.dat  %s/restarts/ ' % self.folder)
        print("FINISHED SAVING FRAMES")
        self.cycle += 1
        os.chdir(self.folder)
        os.system('rm -r tmp')
        self.par['Relax'] = False

    def checkIfStuck(self, values, accumulatedFails) -> bool:
        if accumulatedFails > self.par['Fails'] * int(self.par['NumberCV']):
            print('X' * 200)
            print("\nThe simulation is stuck and it has been terminated\n")
            print('X' * 200)
            with open('walkerSummary.log', 'a') as logFile:
                logFile.write('\nSimulation seemed stuck and it has been terminated\n')
                logFile.close()
            exit()
        else:
            failCount = 0
            for vals in values:
                if vals is not None:
                    x = np.array(vals)
                    x_norm = (x - np.min(x)) / (np.max(x) - np.min(x))
                    if np.std(x_norm) < self.par['Tolerance']:
                        failCount += 1
            if failCount == self.par['NumberCV']:
                print('\nSimulation might be stuck. Running the relaxation protocol.')
                print("")
                return True
            else:
                # if it did not, we continue as usual
                return False
