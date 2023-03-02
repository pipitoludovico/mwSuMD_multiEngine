from lib import MDoperations
from lib import Metrics
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
        suMD1(self.par).countTraj_logTraj(self.max_value)
        MDsetter(self.par).createInput(self.trajCount)
        Runner(self.par).runMD()
        self.walker_metrics = Metrics.Metrics(self.par).getChosenMetrics(self.selection_list)
        self.bestWalker, self.max_value = Metrics.MetricsParser(self.par).getBestWalker(self.walker_metrics)
        MDoperations.MDoperations(self.par, self.folder, self.trajCount).saveStep(self.bestWalker)
        self.trajCount += 1
        print("RUN PROT POST COUNT")
        print(self.trajCount)
        return self.max_value

    def runSlope(self):
        max_cycles = 1 / (int(self.par['Timewindow']) / 10 ** 5)  # run for 1 microsecond and then stop
        while self.trajCount < max_cycles:
            self.max_value = self.runProtocol()
