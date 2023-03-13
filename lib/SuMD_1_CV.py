from .MDoperations import MDoperator
from .Metrics import MetricsParser
from .MDsettings import MDsetter
from .Runners import Runner
from .Parser import mwInputParser


class suMD1(mwInputParser):

    def __init__(self, par):
        super(suMD1, self).__init__()
        self.par = par
        self.bestWalker = None
        self.metric = 0

    def run_SuMD_1_CV(self):
        suMD1(self.par).countTraj_logTraj(self.max_value)
        if self.par['Transition_1'] == 'positive':
            if self.par['Slope'] == 'NO':
                self.max_value = 0
                while self.max_value < self.par['Cutoff_1']:
                    self.max_value = self.runProtocol()
            if self.par['Slope'] == 'YES':
                self.runSlope()

        if self.par['Transition_1'] == 'negative':
            if self.par['Slope'] == 'NO':
                self.max_value = 10 ** 6
                while self.max_value > self.par['Cutoff_1']:
                    self.max_value = self.runProtocol()
            if self.par['Slope'] == 'YES':
                self.runSlope()

    def runProtocol(self):
        # create input files per walker
        import time
        begin = time.perf_counter()
        MDsetter(self.par).createInput(self.trajCount)
        # running the simulations
        Runner(self.par).runMD()
        # getting the metrics and choose the best one
        self.walker_metrics = MetricsParser(self.par).getChosenMetrics()
        self.bestWalker, self.max_value = MetricsParser(self.par).getBestWalker(self.walker_metrics)
        print("best walker: " + str(self.bestWalker, self.max_value))
        # update values and log them
        MDoperator(self.par, self.folder, self.trajCount).saveStep(self.bestWalker)
        suMD1(self.par).countTraj_logTraj(self.max_value)
        end = time.perf_counter()
        final = end - begin
        print("Cycle completed in:" + str(final))
        return self.max_value

    def runSlope(self):
        max_cycles = 1 / (int(self.par['Timewindow']) / 10 ** 5)  # run for 1 microsecond and then stop
        while self.trajCount < max_cycles:
            self.max_value = self.runProtocol()
