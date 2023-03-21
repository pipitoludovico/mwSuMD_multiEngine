import os

import numpy as np


class MDoperator:
    def __init__(self, par, root):
        self.par = par
        self.folder = root
        self.cycle = len(os.listdir(f'{self.folder}/trajectories'))

    def saveStep(self, best_walker):
        _best_walker = best_walker
        if self.par['Relax'] is True:
            best_walker = 1
        os.makedirs('trajectories', exist_ok=True)
        """Handle the restart files and the xtc storage"""
        for r in range(1, int(self.par['Walkers']) + 1):
            if r == best_walker:
                print("Renaming the restart files and moving the best walker's restart bin files")
                os.chdir('tmp/walker_%s' % str(r))
                # moving the best frame to the trajectory folder
                os.system(f'cp wrapped.xtc {self.folder}/trajectories/step_{self.cycle}.xtc')
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
        self.cycle += 1
        os.chdir(self.folder)
        os.system('rm -r tmp')
        self.par['Relax'] = False


def checkIfStuck(values, accumulatedFails) -> bool:
    if accumulatedFails > 2:
        print('X' * 200)
        print("\nThe simulation is stuck and it has been terminated")
        print('X' * 200)
        with open('walkerSummary.log', 'a') as logFile:
            logFile.write('\nSimulation seemed stuck and it has been terminated')
            logFile.close()
        exit()
    else:
        # we check if the simulation did not change much,
        # and we run the relaxation protocol
        if np.std(values) < 0.6:
            print('\nSimulation might be stuck. Running the relaxation protocol.')
            print("")
            return True
        else:
            # if it did not, we continue as usual
            return False


def getSlope(values_metric) -> float:
    """Compute the least square methods on the data
    list provided called by other metrics functions"""
    data = dict(enumerate(values_metric))
    meanTime = np.array(list(data.values())).mean()
    meanDist = np.array(list(data.keys())).mean()
    nume = [(float(value) - meanTime) * (float(key) - meanDist) for key, value in data.items()]
    deNume = [(float(value) - meanTime) ** 2 for value in data.values()]
    try:
        slope = float(np.sum(nume)) / float(np.sum(deNume))
        return slope
    except:
        print("Slope deNumerator was 0.")
        slope = 0
    return slope
