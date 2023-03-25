import time

from .MDoperations import MDoperator
from .MDsettings import MDsetter
from .Metrics import MetricsParser
from .Parser import mwInputParser
from .Runners import Runner


class ProtocolRunner(mwInputParser):
    def __init__(self):
        super(ProtocolRunner).__init__()
        self.bestWalker = None

    def runStandardProtocol(self):
        print("")
        print('#' * 200)
        print("Running mwSuMD protocol")
        print('#' * 200)
        # create input files per walker
        begin = time.perf_counter()
        MDsetter(self.initialParameters).createInputFile()
        # if no user-defined mdp file is provided, we will create a default one
        # running the simulations
        Runner(self.initialParameters).runMD()
        # getting the metrics and choose the best one
        self.walker_metrics = MetricsParser().getChosenMetrics()
        self.bestWalker, self.max_value, self.metric_1, self.metric_2 = MetricsParser().getBestWalker(
            self.walker_metrics)
        print("best walker: " + str(self.bestWalker) + " " + str(self.max_value))
        # update values and log them
        MDoperator(self.initialParameters, self.folder).saveStep(self.bestWalker)
        if self.max_value != 0:
            mwInputParser().countTraj_logTraj(self.max_value)
        else:
            mwInputParser().countTraj_logTraj('')
        end = time.perf_counter()
        final = end - begin
        print("Cycle completed in:" + str(final))
        return self.bestWalker, self.max_value, self.metric_1, self.metric_2
