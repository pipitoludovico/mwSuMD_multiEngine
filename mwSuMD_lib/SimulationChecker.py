import os
import subprocess

from .MDoperations import MDoperator
from .MDsettings import MDsetter
from .Metrics import MetricsParser
from .Parser import mwInputParser
from .GPUoperations import ProcessManager
from .TrajectoryOperator import TrajectoryOperator


class Checker(mwInputParser):
    def __init__(self):
        super(mwInputParser, self).__init__()
        self.best_metric_result = None
        self.best_average_metric_2 = None
        self.best_average_metric_1 = None
        self.best_walker_score = None
        self.averages = None
        self.scores = None
        self.best_value = None
        self.bestWalker = None
        self.trajCount = len(os.listdir(f'{self.folder}/trajectories'))
        self.GPUIDs = ProcessManager.getGPUids()

    def checkIfFailed(self, vals1=None, vals2=None, accumulatedFails=0):
        print('#' * 200)
        print('Checking if trajectory is stuck with values: ' + str(vals1) +
              ". Total fails accumulated: " + str(accumulatedFails))
        mdOperator = MDoperator(self.initialParameters, self.folder)
        if mdOperator.checkIfStuck([vals1, vals2], accumulatedFails) is True:
            self.relaxSystem()
            accumulatedFails += 1
        else:
            accumulatedFails += 0
            print("Number of fails accumulated: " + str(accumulatedFails))
        return accumulatedFails

    def relaxSystem(self):
        print('Relaxation Protocol begins now:')
        print('#' * 200)
        # The relaxation protocol starts here
        self.initialParameters['Relax'] = True
        GPU = self.GPUIDs[0]
        # we create a special input file that has a longer runtime (5ns default or user-defined)
        MDsetter(self.initialParameters).createInputFile()
        # we run this inside walker_1 for convenience
        os.chdir('tmp/walker_1')
        for file in os.listdir(os.getcwd()):
            if file.endswith('.inp'):
                subprocess.Popen(f'acemd3 --device {GPU} {file} 1> relax.log', shell=True).wait()
            elif file.endswith('.namd'):
                subprocess.Popen(f'namd3 +p8 +devices {GPU} {file} 1> relax.log', shell=True).wait()
            elif file.endswith('.mdp'):
                subprocess.Popen(
                    f'gmx convert-tpr -s {self.folder}/restarts/previous.tpr -extend {int(self.initialParameters["RelaxTime"] * 1000)} -o {self.initialParameters["Output"]}_{self.trajCount}.tpr &>tpr_log.log',
                    shell=True).wait()
                command = f'gmx mdrun -deffnm {self.initialParameters["Output"]}_{self.trajCount}'
                subprocess.Popen(command, shell=True).wait()
        os.chdir(f'{self.folder}')
        TrajectoryOperator().wrap(1)
        # we then compute its metrics as a reference
        if self.initialParameters['NumberCV'] == 1:
            self.scores, self.averages = MetricsParser().getChosenMetrics()
        else:
            self.walker_metrics, self.averages = MetricsParser().getChosenMetrics()
        # then the last coordinate is saved
        MDoperator(self.initialParameters, self.folder).saveStep(1)
        # we then extract the best metric/score and store it as a reference
        if self.initialParameters['NumberCV'] == 1:
            self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(
                self.scores, self.averages)
        else:
            self.bestWalker, self.best_walker_score, self.best_average_metric_1, self.best_average_metric_2 = MetricsParser().getBestWalker(
                self.walker_metrics[0], self.walker_metrics[1], self.averages[0], self.averages[1])
            self.best_metric_result = [self.best_average_metric_1, self.best_average_metric_2]
        with open('walkerSummary.log', 'a') as walkerSummary:
            info_to_write = str(self.trajCount) + " RELAXATION PROTOCOL SCORE: " + str(self.best_walker_score) + " Metrics: " + str(self.best_metric_result) + "\n"
            walkerSummary.write(info_to_write)
        # self.trajCount += 1
        print("\nRelaxation Protocol Ended")
        print('#' * 200)
        print('\n\n')
        # setting our check to False and end the protocol, beginning a new cycle.
        self.initialParameters['Relax'] = False
