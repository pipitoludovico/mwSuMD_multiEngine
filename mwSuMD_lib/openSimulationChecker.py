import os

from .openMDoperations import MDoperator
from .Metrics import MetricsParser
from .Parser import mwInputParser
from .openRunners import Runner


class Checker(mwInputParser):
    def __init__(self):
        super(mwInputParser, self).__init__()
        self.last_frame_metrics = None
        self.best_metric_result = None
        self.best_average_metric_2 = None
        self.best_average_metric_1 = None
        self.best_walker_score = None
        self.averages = None
        self.scores = None
        self.best_value = None
        self.bestWalker = None
        self.trajCount = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])

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
        self.initialParameters['Relax'] = True
        Runner(self.initialParameters).runAndWrap()
        if self.initialParameters['NumberCV'] == 1:
            self.scores, self.last_frame_metrics = MetricsParser().getChosenMetrics()
        else:
            self.walker_metrics, self.last_frame_metrics = MetricsParser().getChosenMetrics()
        # we then extract the best metric/score and store it as a reference
        if self.initialParameters['NumberCV'] == 1:
            self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(self.scores, self.last_frame_metrics)
        else:
            self.bestWalker, self.best_walker_score, self.best_average_metric_1, self.best_average_metric_2 = MetricsParser().getBestWalker(self.walker_metrics[0], self.walker_metrics[1], self.last_frame_metrics[0], self.last_frame_metrics[1])
            self.best_metric_result = [self.best_average_metric_1, self.best_average_metric_2]
        MDoperator(self.initialParameters, self.folder).saveStep(self.bestWalker)
        with open('walkerSummary.log', 'a') as walkerSummary:
            info_to_write = str(self.trajCount) + " RELAXATION PROTOCOL SCORE: " + str(self.best_walker_score) + " Metrics: " + str(self.best_metric_result) + "\n"
            walkerSummary.write(info_to_write)
        self.trajCount += 1
        print("\nRelaxation Protocol Ended")
        print('#' * 200)
        print('\n\n')
        # setting our check to False and end the protocol, beginning a new cycle.
        self.initialParameters['Relax'] = False
