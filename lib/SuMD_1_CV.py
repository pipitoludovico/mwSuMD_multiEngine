from lib import FolderOps
from lib import MDoperations
from lib import metrics
from .MDsettings import MDsetter
from .Runners import Runner
from .mwParser import mwInputParser


class suMD1(mwInputParser):

    def __init__(self, par):
        super(suMD1, self).__init__()
        print('INIT SUMD1')
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
        MDsetter(self.par).createInput(self.trajCount)
        Runner(self.par).runMD()
        self.walker_metrics = metrics.Metrics(self.par).getChosenMetrics(self.selection_list)
        self.bestWalker, self.max_value = metrics.Metrics(self.par).getBestWalker(self.walker_metrics)
        MDoperations.MDoperations(self.par, self.folder, self.trajCount).saveStep(self.bestWalker)
        self.trajCount += 1
        FolderOps.FolderOps(self.par).countTraj_logTraj(self.max_value)
        return self.max_value

    def runSlope(self):
        max_cycles = 1 / (int(self.par['Timewindow']) / 10 ** 5)  # run for 1 microsecond and then stop
        while self.trajCount < max_cycles:
            self.max_value = self.runProtocol()
