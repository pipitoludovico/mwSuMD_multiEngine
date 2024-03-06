import os.path
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
        if os.path.exists("tmp"):
            os.system("rm -r tmp")
        begin = time.perf_counter()
        # running the simulations
        Runner().runAndWrap()
        # compute metrics
        self.scores = MetricsParser().getChosenMetrics()
        # get metrics
        self.bestWalker, self.best_walker_score, self.best_metric_result = None, None, None
        self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(self.scores)
        MDoperator(self.initialParameters, self.folder).saveStep(self.bestWalker, self.best_walker_score, self.best_metric_result)

        end = time.perf_counter()
        final = end - begin
        Logger.LogToFile('ad', self.trajCount, "Cycle completed in:" + str(final))
        return self.bestWalker, self.best_walker_score, self.best_metric_result
