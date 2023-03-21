import os
import subprocess

import lib.MDoperations
from .MDoperations import MDoperator
from .MDsettings import MDsetter
from .Metrics import MetricsParser
from .Parser import mwInputParser
from .Runners import Runner


class Checker(mwInputParser):
    def __init__(self):
        super(Checker).__init__()
        self.bestWalker = None
        self.trajCount = len(os.listdir(f'{self.folder}/trajectories'))

    def checkIfFailed(self, vals, accumulatedFails):
        print('#' * 200)
        print('Checking if trajectory is stuck with values: ' + str(vals) +
              ". Total fails accumulated: " + str(accumulatedFails))
        if lib.MDoperations.checkIfStuck(vals, accumulatedFails) is True:
            self.runRelaxationProtocol()
            accumulatedFails += 1
            print("Number of fails accumulated: " + str(accumulatedFails))
        return accumulatedFails

    def runRelaxationProtocol(self):
        print("")
        print('#' * 200)
        print('Relaxation Protocol begins now:')
        print('#' * 200)
        print("")
        # The relaxation protocol starts here
        self.par['Relax'] = True
        # we create a special input file that has a longer runtime (5ns default or user-defined)
        MDsetter(self.par).createACEMDinputFile()
        # we run this inside walker_1 for convenience
        os.chdir('tmp/walker_1')
        subprocess.Popen(f'acemd3 --device 0 input_1_{self.trajCount}.inp 1> acemd.log', shell=True).wait()
        os.chdir(f'{self.folder}')
        Runner(self.par).wrap(1)
        # we then compute its metrics as a reference
        self.walker_metrics = MetricsParser().getChosenMetrics()
        # then the last coordinate is saved
        MDoperator(self.par, self.folder).saveStep(self.bestWalker)
        self.trajCount += 1
        # we then extract the best metric/score and store it as a reference
        self.bestWalker, self.max_value, self.metric_1, self.metric_2 = MetricsParser().getBestWalker(
            self.walker_metrics)

        if self.par['NumberCV'] == 1:
            MetricsParser().countTraj_logTraj(lib.MDoperations.getSlope(self.walker_metrics))
        else:
            MetricsParser().countTraj_logTraj(lib.MDoperations.getSlope(self.walker_metrics[0]))
        print("\nRelaxation Protocol Ended")
        print('#' * 200)
        print('\n\n')
        # setting our check to False and end the protocol, beginning a new cycle.
        self.par['Relax'] = False