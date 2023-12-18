import time
from os import listdir
from .openMDoperations import MDoperator
from .Metrics import MetricsParser
from .Parser import mwInputParser
from .openRunners import Runner
from .Loggers import Logger


class ProtocolRunner(mwInputParser):
    def __init__(self):
        super(mwInputParser, self).__init__()
        self.best_metric_result = None
        self.best_average_metric_1 = None
        self.best_average_metric_2 = None
        self.best_walker_score = None
        self.scores = None
        self.last_frame_metrics = None
        self.bestWalker = None
        self.trajCount = len([traj for traj in listdir(f'{self.folder}/trajectories') if traj.endswith('.xtc')])

    def runStandardProtocol(self):
        Logger.LogToFile("w", self.trajCount, "*" * 200 + "\nRunning mwSuMD protocol\n" + "#" * 200)
        # create input files per walker
        begin = time.perf_counter()
        # running the simulations
        Runner(self.initialParameters).runAndWrap()
        # compute metrics
        if self.initialParameters['NumberCV'] == 1:
            self.scores, self.last_frame_metrics = MetricsParser().getChosenMetrics()
        else:
            self.walker_metrics, self.last_frame_metrics = MetricsParser().getChosenMetrics()
        # get metrics
        if self.initialParameters['NumberCV'] == 1:
            self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(self.scores, self.last_frame_metrics)
            mwInputParser().countTraj_logTraj(["Best Walker: " + str(self.bestWalker) + " Best Metric: " + str(self.best_walker_score) + " Last Metric: " + str(self.best_metric_result)])
        else:
            self.bestWalker, self.best_walker_score, self.best_average_metric_1, self.best_average_metric_2 = MetricsParser().getBestWalker(self.walker_metrics[0], self.walker_metrics[1], self.last_frame_metrics[0], self.last_frame_metrics[1])
            self.best_metric_result = [self.best_average_metric_1, self.best_average_metric_2]
            mwInputParser().countTraj_logTraj(["Best Walker: " + str(self.bestWalker) + " Score Result: " + str(self.best_walker_score) + " Last Metrics from best: " + str(self.best_metric_result)])
        MDoperator(self.initialParameters, self.folder).saveStep(self.bestWalker)

        end = time.perf_counter()
        final = end - begin
        Logger.LogToFile('ad', self.trajCount, "Cycle completed in:" + str(final))
        return self.bestWalker, self.best_walker_score, self.best_metric_result
