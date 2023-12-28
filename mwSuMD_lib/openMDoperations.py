import os
from .Loggers import Logger

import numpy as np


class MDoperator:
    def __init__(self, par, root):
        self.par = par
        self.folder = root
        self.extensions = ('.chk', '.xml')
        self.cycle = len(
            [trajFile for trajFile in os.listdir(f'{self.folder}/trajectories') if trajFile.endswith('xtc')])

    def saveStep(self, best_walker, walker_score, best_metric_result):
        """Handle the restart files and the xtc storage"""
        os.chdir('tmp/walker_%s' % str(best_walker))
        check = [binary for binary in os.listdir("./") if binary.endswith(self.extensions)]
        if check:
            os.system(f'cp wrapped.xtc {self.folder}/trajectories/{self.par["Output"]}_step_{self.cycle}.xtc')
            # moving and renaming the binary files to the restart folder
            os.system(f'cp *.chk {self.folder}/restarts/previous.chk')
            os.system(f'cp *.xml {self.folder}/restarts/previous.xml')
            Logger.LogToFile('a', self.cycle, "FINISHED SAVING FRAMES")
            if self.par['PLUMED'] is not None:
                os.system('cp HILLS  %s/restarts/ ' % self.folder)
                os.system('cp COLVAR  %s/restarts/ ' % self.folder)
                os.system('cp grid.dat  %s/restarts/ ' % self.folder)
        else:
            Logger.LogToFile('ad', self.cycle, "No binary saved: restarting from last checkpoint.")

        os.chdir(self.folder)
        self.cycle += 1

        if self.par['Relax']:
            with open('walkerSummary.log', 'a') as walkerSummary:
                info_to_write = " RELAXATION PROTOCOL SCORE: " + str(walker_score) + " Metrics: " + str(
                    best_metric_result) + "\n"
                walkerSummary.write(info_to_write)
        else:
            with open('walkerSummary.log', 'a') as walkerSummary:
                info_to_write = "Best Walker: " + str(best_walker) + " Best Metric: " + str(
                    walker_score) + " Last Metric: " + str(best_metric_result) + "\n"
                walkerSummary.write(info_to_write)
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
