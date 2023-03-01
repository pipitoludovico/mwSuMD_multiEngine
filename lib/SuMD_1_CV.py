from .mwParser import mwInputParser
from .MDoperations import MDoperations


class suMD1(mwInputParser):

    def __init__(self):
        super(suMD1, self).__init__()
        # self.mdOpeartions = MDoperations()

    def run_SuMD_1_CV(self):
        if self.par['Slope'] == 'NO':
            self.runSlopeNo()
            print(self.trajCount)
        elif self.par['Slope'] == 'YES':
            self.runSlopeYes()

    def runSlopeNo(self):
        if self.par['Transition_1'] == 'negative':
            self.max_value = 1000000
            while self.max_value > float(self.par['Cutoff_1']):
                self.new_value, self.walker_metrics, self.slope = MDoperations().setAndRun()
                if self.new_value < self.max_value:
                    self.max_value = self.new_value
                    print(self.trajCount, self.new_value, self.slope)
                    mwInputParser().countTraj_logTraj(self.walker_metrics, self.slope)
                    self.trajCount += 1
                    self.walker_metrics.clear()

        if self.par['Transition_1'] == 'positive':
            self.max_value = 0
            while self.max_value < float(self.par['Cutoff_1']):
                MDoperations().setAndRun()
                if self.new_value > self.max_value:
                    self.max_value = self.new_value
                    mwInputParser.countTraj_logTraj(self.new_value, self.slope)
                    self.trajCount += 1
                    self.walker_metrics.clear()

    def runSlopeYes(self):
        max_cycles = 1 / (int(self.par['Timewindow']) / 10 ** 5)  # run for 1 microsecond and then stop
        while self.trajCount < max_cycles:
            MDoperations().setAndRun()
