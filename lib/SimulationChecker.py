import os
import subprocess

from .MDoperations import MDoperator
from .MDsettings import MDsetter
from .Metrics import MetricsParser
from .Parser import mwInputParser
from .TrajectoryOperator import TrajectoryOperator


class Checker(mwInputParser):
    def __init__(self):
        super(Checker).__init__()
        self.averages = None
        self.scores = None
        self.best_value = None
        self.bestWalker = None
        self.trajCount = len(os.listdir(f'{self.folder}/trajectories'))

    def checkIfFailed(self, vals, accumulatedFails):
        print('#' * 200)
        print('Checking if trajectory is stuck with values: ' + str(vals) +
              ". Total fails accumulated: " + str(accumulatedFails))
        if MDoperator.checkIfStuck(vals, accumulatedFails) is True:
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
        # we create a special input file that has a longer runtime (5ns default or user-defined)
        MDsetter(self.initialParameters).createInputFile()
        # we run this inside walker_1 for convenience
        os.chdir('tmp/walker_1')
        for file in os.listdir(os.getcwd()):
            if file.endswith('.inp'):
                subprocess.Popen(f'acemd3 --device 0 {file} 1> relax.log', shell=True).wait()
            elif file.endswith('.namd'):
                subprocess.Popen(f'namd3 +p8 +device 0 {file} 1> relax.log', shell=True).wait()
            elif file.endswith('.mdp'):
                subprocess.Popen(f'gmx grompp -f production.mdp'
                                 f' -c {self.folder}/system/{self.initialParameters["GRO"]}'
                                 f' -t {self.folder}/system/{self.initialParameters["CPT"]}'
                                 f' -p {self.folder}/system/{self.initialParameters["TOP"]}'
                                 f' -o {self.initialParameters["Output"]}_{self.trajCount}.tpr').wait()
                command = f'gmx mdrun -deffnm {self.initialParameters["Output"]}_{self.trajCount}'
                subprocess.Popen(command).wait()

        os.chdir(f'{self.folder}')
        TrajectoryOperator().wrap(1)
        # we then compute its metrics as a reference
        if self.initialParameters['NumberCV'] == 1:
            self.scores, self.averages = MetricsParser().getChosenMetrics()
        else:
            self.walker_metrics, self.averages = MetricsParser().getChosenMetrics()
        # then the last coordinate is saved
        MDoperator(self.initialParameters, self.folder).saveStep(1)
        self.trajCount += 1
        # we then extract the best metric/score and store it as a reference
        if self.initialParameters['NumberCV'] == 1:
            self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(
                self.scores, self.averages)
        else:
            self.bestWalker, self.best_walker_score, self.best_average_metric_1, self.best_average_metric_2 = MetricsParser().getBestWalker(
                self.walker_metrics[0], self.walker_metrics[1], self.averages[0], self.averages[1])
            self.best_metric_result = [self.best_average_metric_1, self.best_average_metric_2]

        MetricsParser().countTraj_logTraj(["RELAXATION PROTOCOL SCORE: " + str(self.best_walker_score) + " Metric: "
                                           + str(self.best_metric_result)])
        print("\nRelaxation Protocol Ended")
        print('#' * 200)
        print('\n\n')
        # setting our check to False and end the protocol, beginning a new cycle.
        self.initialParameters['Relax'] = False
