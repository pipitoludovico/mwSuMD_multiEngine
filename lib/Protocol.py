import time
import os

from .MDoperations import MDoperator
from .MDsettings import MDsetter
from .Metrics import MetricsParser
from .Parser import mwInputParser
from .Runners import Runner


class ProtocolRunner(mwInputParser):
    def __init__(self):
        super(ProtocolRunner).__init__()
        self.best_metric_result = None
        self.best_average_metric_2 = None
        self.best_average_metric_1 = None
        self.best_walker_score = None
        self.bestWalker = None
        self.last_frame_metrics = None
        self.scores = None

    def runStandardProtocol(self):
        print("")
        print('#' * 200)
        print("Running mwSuMD protocol")
        print('#' * 200)

        begin = time.perf_counter()
        # purging previous tmp folder:
        if os.path.exists('./tmp'):
            os.system('rm -r ./tmp')
        MDsetter(self.initialParameters).createInputFile()
        Runner(self.initialParameters).runMD()
        # compute metrics
        if self.initialParameters['NumberCV'] == 1:
            self.scores, self.last_frame_metrics = MetricsParser().getChosenMetrics()
        else:
            self.walker_metrics, self.last_frame_metrics = MetricsParser().getChosenMetrics()
        # get metrics
        if self.initialParameters['NumberCV'] == 1:
            self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(
                self.scores, self.last_frame_metrics)
            mwInputParser().countTraj_logTraj(
                ["Best Walker: " + str(self.bestWalker) + " Best Metric: " + str(
                    self.best_walker_score) + " Last Metric: " + str(self.best_metric_result)])
        else:
            self.bestWalker, self.best_walker_score, self.best_average_metric_1, self.best_average_metric_2 = \
                MetricsParser().getBestWalker(self.walker_metrics[0], self.walker_metrics[1],
                                              self.last_frame_metrics[0], self.last_frame_metrics[1])
            self.best_metric_result = [self.best_average_metric_1, self.best_average_metric_2]
            mwInputParser().countTraj_logTraj(["Best Walker: " + str(self.bestWalker) + " Score Result: " + str(
                self.best_walker_score) + " Last Metrics from best: " + str(self.best_metric_result)])

        MDoperator(self.initialParameters, self.folder).saveStep(self.bestWalker)

        end = time.perf_counter()
        final = end - begin
        print("Cycle completed in:" + str(final))
        return self.bestWalker, self.best_walker_score, self.best_metric_result
