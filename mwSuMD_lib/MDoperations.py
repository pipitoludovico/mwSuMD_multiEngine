import os
import subprocess
import numpy as np
from .Loggers import Logger
from warnings import filterwarnings

filterwarnings(action='ignore')


class MDoperator:
    def __init__(self, par, root):
        self.par = par
        self.folder = root
        self.extensions = ('.coor', '.xsc', '.vel', '.gro', '.cpt', '.tpr')
        self.cycle = len([trajFile for trajFile in os.listdir(f'{self.folder}/trajectories') if trajFile.endswith('xtc')])

    def saveStep(self, best_walker, walker_score, best_metric_result):
        """Handles the restart files and the binary storage for OPENMM"""
        os.chdir('tmp/walker_%s' % str(best_walker))
        check = [binary for binary in os.listdir("./") if binary.endswith(self.extensions)]
        if check:
            # moving the best frame to the trajectory folder
            os.system(f'cp wrapped.xtc {self.folder}/trajectories/{self.par["Output"]}_step_{self.cycle}.xtc')
            if self.par['MDEngine'] != 'GROMACS':
                os.system(f'cp *.coor {self.folder}/restarts/previous.coor')
                os.system(f'cp *.xsc {self.folder}/restarts/previous.xsc')
                os.system(f'cp *.vel {self.folder}/restarts/previous.vel')
            elif self.par['MDEngine'] == 'GROMACS':
                os.system(f'cp *.gro {self.folder}/restarts/previous.gro')
                os.system(f'cp "$(ls -t *.cpt | head -1)" {self.folder}/restarts/previous.cpt')
                os.system(f'cp *.tpr {self.folder}/restarts/previous.tpr')
            elif self.par['PLUMED']:
                os.system('cp HILLS  %s/restarts/ ' % self.folder)
                os.system('cp COLVAR  %s/restarts/ ' % self.folder)
                os.system('cp grid.dat  %s/restarts/ ' % self.folder)
        else:
            Logger.LogToFile('ad', self.cycle, "No binary saved: restarting from last checkpoint.")

        os.chdir(self.folder)
        self.cycle += 1

        if self.par['Relax']:
            with open('walkerSummary.log', 'a') as walkerSummary:
                info_to_write = " RELAXATION PROTOCOL SCORE: " + str(walker_score) + " Metrics: " + str(best_metric_result) + "\n"
                walkerSummary.write(info_to_write)
        else:
            with open('walkerSummary.log', 'a') as walkerSummary:
                info_to_write = str(self.cycle) + " Best Walker: " + str(best_walker) + " Best Metric: " + str(walker_score) + " Last Metric: " + str(best_metric_result) + "\n"
                walkerSummary.write(info_to_write)
        Logger.LogToFile('a', self.cycle, "FINISHED SAVING FRAMES")
        os.system('rm -r tmp')
        self.par['Relax'] = False

    def prepareTPR(self, walk_count, trajcount, customFile=None):
        gro = f'{self.par["Root"]}/system/{self.par["GRO"]}' if self.cycle == 0 else f'{self.par["Root"]}/restarts/previous.gro'
        cpt = f'{self.par["Root"]}/system/{self.par["CPT"]}' if self.cycle == 0 else f'{self.par["Root"]}/restarts/previous.cpt'
        if customFile is None:
            command = f'gmx grompp -f input_{walk_count}_{trajcount}.mdp' \
                      f' -c {gro} -t {cpt} -p {self.folder}/system/{self.par["TOP"]}' \
                      f' -o {self.par["Output"]}_{trajcount}_{walk_count}.tpr -maxwarn 3 &>tpr_log.log'
        else:
            command = f'gmx grompp -f input_{walk_count}_{trajcount}.mdp' \
                      f' -c {gro} -t {cpt} -p {self.folder}/system/{self.par["TOP"]}' \
                      f' -o production.tpr -maxwarn 3 &>tpr_log.log'

        tprPreparation = subprocess.Popen(command, shell=True)
        tprPreparation.wait()

    def checkIfStuck(self, values, accumulatedFails) -> bool:
        if accumulatedFails >= self.par['Fails'] * int(self.par['NumberCV']):
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
