import os

from .Parser import mwInputParser
from .Protocol import ProtocolRunner
from .SimulationChecker import Checker


class suMD1(mwInputParser):

    def __init__(self, par):
        super(suMD1, self).__init__()
        self.par = par
        self.bestWalker = None
        self.checkVals_1 = []
        self.checkVals_2 = []
        self.fails = 0
        self.cycle = len(os.listdir(f'{self.folder}/trajectories'))

    def run_mwSuMD(self):
        suMD1(self.par).countTraj_logTraj(self.max_value)
        if self.par['NumberCV'] == 1:
            if self.par['Transition_1'] == 'positive':
                self.max_value = 0
                while self.max_value < self.par['Cutoff_1']:
                    self.max_value = self.runProtocol()

            if self.par['Transition_1'] == 'negative':
                self.max_value = 10 ** 6
                while self.max_value > self.par['Cutoff_1']:
                    self.max_value = self.runProtocol()
        else:
            if self.par['Transition_1'] == 'positive' and self.par['Transition_2'] == 'positive':
                self.metric_1, self.metric_2 = 0, 0
                while self.metric_1 < self.par['Cutoff_1'] and self.metric_2 < self.par['Cutoff_2']:
                    self.metric_1, self.metric_2 = self.runProtocol()

            if self.par['Transition_1'] == 'positive' and self.par['Transition_2'] == 'negative':
                self.metric_1, self.metric_2 = 0, 10 ** 6
                while self.metric_1 < self.par['Cutoff_1'] and self.metric_2 > self.par['Cutoff_2']:
                    self.metric_1, self.metric_2 = self.runProtocol()

            if self.par['Transition_1'] == 'negative' and self.par['Transition_2'] == 'negative':
                self.metric_1, self.metric_2 = 10 ** 6, 10 ** 6
                while self.metric_1 > self.par['Cutoff_1'] and self.metric_2 > self.par['Cutoff_2']:
                    self.metric_1, self.metric_2 = self.runProtocol()

            if self.par['Transition_1'] == 'negative' and self.par['Transition_2'] == 'positive':
                self.metric_1, self.metric_2 = 10 ** 6, 0
                while self.metric_1 > self.par['Cutoff_1'] and self.metric_2 < self.par['Cutoff_2']:
                    self.metric_1, self.metric_2 = self.runProtocol()

    def runProtocol(self):
        mwInputParser().getRestartOutput()
        self.bestWalker, self.max_value, self.metric_1, self.metric_2 = ProtocolRunner().runStandardProtocol()
        self.cycle += 1
        self.checkVals_1.append(self.metric_1)
        self.checkVals_2.append(self.metric_2)

        if self.cycle % (self.par['RelaxTime'] * 1000 / self.par['Timewindow']) == 0:
            checker = Checker()
            self.fails = checker.checkIfFailed(self.checkVals_1, self.fails)
            self.fails = checker.checkIfFailed(self.checkVals_2, self.fails)
            self.checkVals_1.clear()
            self.checkVals_2.clear()

        if self.par['NumberCV'] == 1:
            return self.max_value
        else:
            return self.metric_1, self.metric_2
