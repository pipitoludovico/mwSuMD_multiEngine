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
        if self.par['MDEngine'] == 'ACEMD':
            MDsetter(self.par).createACEMDinputFile()
        # if no user-defined mdp file is provided, we will create a default one
        if self.par['MDEngine'] == 'GROMACS' and self.par['MDP'] is None:
            MDsetter(self.par).createGROMACSinputFile()
        # running the simulations
        Runner(self.par).runMD()
        # getting the metrics and choose the best one
        self.walker_metrics = MetricsParser().getChosenMetrics()
        self.bestWalker, self.max_value, self.metric_1, self.metric_2 = MetricsParser().getBestWalker(
            self.walker_metrics)
        print("best walker: " + str(self.bestWalker) + " " + str(self.max_value))
        # update values and log them
        MDoperator(self.par, self.folder).saveStep(self.bestWalker)
        mwInputParser().countTraj_logTraj(self.max_value)
        end = time.perf_counter()
        final = end - begin
        print("Cycle completed in:" + str(final))
        return self.bestWalker, self.max_value, self.metric_1, self.metric_2
