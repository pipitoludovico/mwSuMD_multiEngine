import time
import os

from .MDoperations import MDoperator
from .MDsettings import MDsetter
from .Metrics import MetricsParser
from .Parser import mwInputParser
from .Runners import Runner
from .Loggers import Logger
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
        if os.path.exists('./tmp'):
            os.system('rm -r ./tmp')
        MDsetter(self.initialParameters).createInputFile()
        if not self.openMM:
            Runner(self.initialParameters, self.openMM).runMD()
        else:
            Runner(self.initialParameters, self.openMM).runAndWrap()

        # compute metrics
        self.scores = MetricsParser().getChosenMetrics()
        # get metrics
        self.bestWalker, self.best_walker_score, self.best_metric_result = None, None, None
        self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(self.scores)
        MDoperator(self.initialParameters, self.folder, self.openMM).saveStep(self.bestWalker, self.best_walker_score, self.best_metric_result)

        end = time.perf_counter()
        final = end - begin
        Logger.LogToFile('ad', self.trajCount, "Cycle completed in:" + str(final))
        return self.bestWalker, self.best_walker_score, self.best_metric_result
