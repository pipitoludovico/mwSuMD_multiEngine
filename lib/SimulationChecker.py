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
        self.bestWalker = None
        self.trajCount = len(os.listdir(f'{self.folder}/trajectories'))

    def checkIfFailed(self, vals, accumulatedFails):
        print('#' * 200)
        print('Checking if trajectory is stuck with values: ' + str(vals) +
              ". Total fails accumulated: " + str(accumulatedFails))
        if MDoperator.checkIfStuck(vals, accumulatedFails) is True:
            self.relaxSystem()
            accumulatedFails += 1
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
                subprocess.Popen(f'namd3 +auto-provide +device 0 {file} 1> relax.log', shell=True).wait()
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
        self.walker_metrics = MetricsParser().getChosenMetrics()
        # then the last coordinate is saved
        MDoperator(self.initialParameters, self.folder).saveStep(1)
        self.trajCount += 1
        # we then extract the best metric/score and store it as a reference
        self.bestWalker, self.max_value, self.metric_1, self.metric_2 = MetricsParser().getBestWalker(
            self.walker_metrics)

        if self.initialParameters['NumberCV'] == 1:
            MetricsParser().countTraj_logTraj(MetricsParser.getSlope(self.walker_metrics))
        else:
            MetricsParser().countTraj_logTraj(MetricsParser.getSlope(self.walker_metrics[0]))
        print("\nRelaxation Protocol Ended")
        print('#' * 200)
        print('\n\n')
        # setting our check to False and end the protocol, beginning a new cycle.
        self.initialParameters['Relax'] = False
