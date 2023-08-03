import os
import subprocess

import numpy as np


class MDoperator:
    def __init__(self, par, root):
        self.par = par
        self.folder = root
        self.initialTrajectoriesInFolder = len(
            [trajFile for trajFile in os.listdir(f'{self.folder}/trajectories') if trajFile.endswith('xtc')])
        self.cycle = None

    def saveStep(self, best_walker):
        self.cycle = self.initialTrajectoriesInFolder
        os.makedirs('trajectories', exist_ok=True)
        """Handle the restart files and the xtc storage"""
        for r in range(1, int(self.par['Walkers']) + 1):
            if r == best_walker:
                os.chdir('tmp/walker_%s' % str(r))
                # moving the best frame to the trajectory folder
                os.system(f'cp wrapped.xtc {self.folder}/trajectories/{self.par["Output"]}_step_{self.cycle}.xtc')
                self.cycle += 1
                # moving and renaming the binary files to the restart folder
                if self.par['MDEngine'] != 'GROMACS':
                    os.system(f'cp *.coor {self.folder}/restarts/previous.coor')
                    os.system(f'cp *.xsc {self.folder}/restarts/previous.xsc')
                    os.system(f'cp *.vel {self.folder}/restarts/previous.vel')

                elif self.par['MDEngine'] == 'GROMACS':
                    os.system(f'cp *.gro {self.folder}/restarts/previous.gro')
                    os.system(f'cp *.cpt {self.folder}/restarts/previous.cpt')
                    os.system(f'cp *.tpr {self.folder}/restarts/previous.tpr')

                elif self.par['PLUMED'] is not None:
                    os.system('cp HILLS  %s/restarts/ ' % self.folder)
                    os.system('cp COLVAR  %s/restarts/ ' % self.folder)
                    os.system('cp grid.dat  %s/restarts/ ' % self.folder)
        print("FINISHED SAVING FRAMES")
        os.chdir(self.folder)
        os.system('rm -r tmp')
        self.par['Relax'] = False

    def prepareTPR(self, walk_count, trajcount, customFile=None):
        gro = f'{self.par["Root"]}/system/{self.par["GRO"]}' if self.initialTrajectoriesInFolder == 0 else f'{self.par["Root"]}/restarts/previous.gro'
        cpt = f'{self.par["Root"]}/system/{self.par["CPT"]}' if self.initialTrajectoriesInFolder == 0 else f'{self.par["Root"]}/restarts/previous.cpt'
        if customFile is None:
            command = f'gmx grompp -f input_{walk_count}_{trajcount}.mdp' \
                      f' -c {gro} -t {cpt} -p {self.folder}/system/{self.par["TOP"]}' \
                      f' -o {self.par["Output"]}_{trajcount}_{walk_count}.tpr -maxwarn 1 &>tpr_log.log'
        else:
            command = f'gmx grompp -f input_{walk_count}_{trajcount}.mdp' \
                      f' -c {gro} -t {cpt} -p {self.folder}/system/{self.par["TOP"]}' \
                      f' -o production.tpr -maxwarn 1 &>tpr_log.log'

        tprPreparation = subprocess.Popen(command, shell=True)
        tprPreparation.wait()

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
