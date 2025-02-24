import time
import os

from mwSuMD_lib.MetricOperators.MDoperations import MDoperator
from mwSuMD_lib.MDsetters.MDsettings import MDsetter
from mwSuMD_lib.MetricOperators.Metrics import MetricsParser
from mwSuMD_lib.Parsers.InputfileParser import mwInputParser
from .Runners import Runner
from mwSuMD_lib.Utilities.Loggers import Logger
from warnings import filterwarnings

filterwarnings(action='ignore')


class ProtocolRunner(mwInputParser):
    def __init__(self, openMM):
        super(mwInputParser, self).__init__()
        self.openMM = openMM
        self.best_metric_result = None
        self.best_average_metric_2 = None
        self.best_average_metric_1 = None
        self.best_walker_score = None
        self.bestWalker = None
        self.last_frame_metrics = None
        self.scores = None
        self.trajCount = len([traj for traj in os.listdir(f'{self.folder}/trajectories') if traj.endswith('.xtc')])

    def runStandardProtocol(self):
        Logger.LogToFile("w", self.trajCount, "*" * 200 + "\nRunning mwSuMD protocol\n" + "#" * 200)
        # create input files per walker after purging the existing one
        begin = time.perf_counter()
        op = MDoperator(self.initialParameters, self.folder, self.openMM)
        if not self.openMM:
            MDsetter(self.initialParameters).createInputFile()
            Runner(self.initialParameters, self.openMM, op).runMD()
        else:
            Runner(self.initialParameters, self.openMM, op).runAndWrap()

        # compute metrics
        self.scores = MetricsParser().getChosenMetrics()
        # get metrics
        self.bestWalker, self.best_walker_score, self.best_metric_result = None, None, None
        self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(self.scores)
        MDoperator(self.initialParameters, self.folder, self.openMM).saveStep(self.bestWalker, self.best_walker_score, self.best_metric_result)
        os.chdir(self.folder)
        end = time.perf_counter()
        final = end - begin
        Logger.LogToFile('ad', self.trajCount, "Cycle completed in:" + str(final))
        return self.bestWalker, self.best_walker_score, self.best_metric_result
