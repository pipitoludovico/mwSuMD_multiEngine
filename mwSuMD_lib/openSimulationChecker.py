import os

from .openMDoperations import MDoperator
from .Metrics import MetricsParser
from .Parser import mwInputParser
from .openRunners import Runner
from .Loggers import Logger


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
        Logger.LogToFile("w", self.trajCount,
                         "*" * 200 + '\nChecking if trajectory is stuck with values: ' + str(vals1) +
                         ". Total fails accumulated: " + str(accumulatedFails) + "#" * 200)
        mdOperator = MDoperator(self.initialParameters, self.folder)
        if mdOperator.checkIfStuck([vals1, vals2], accumulatedFails) is True:
            Logger.LogToFile("ad", self.trajCount, "\nRUNNING RELAXATION PROTOCOL" + "#" * 200)
            self.relaxSystem()
            accumulatedFails += 1
        else:
            accumulatedFails += 0
            Logger.LogToFile("ad", self.trajCount, "Number of fails accumulated: " + str(accumulatedFails))
        return accumulatedFails

    def relaxSystem(self):
        Logger.LogToFile('a', self.trajCount, 'Relaxation Protocol begins now:\n' + ('#' * 200))
        self.initialParameters['Relax'] = True
        Runner().runAndWrap()
        self.scores = MetricsParser().getChosenMetrics()
        # we then extract the best metric/score and store it as a reference
        self.bestWalker, self.best_walker_score, self.best_metric_result = None, None, None
        self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(self.scores)
        MDoperator(self.initialParameters, self.folder).saveStep(self.bestWalker, self.best_walker_score, self.best_metric_result)

        Logger.LogToFile('ad', self.trajCount, "\nRelaxation Protocol Ended\n" + "#" * 200)
        self.trajCount += 1
        # setting our check to False and end the protocol, beginning a new cycle.
        self.initialParameters['Relax'] = False
