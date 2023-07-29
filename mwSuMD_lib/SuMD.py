import os

from .Parser import mwInputParser
from .Protocol import ProtocolRunner
from .SimulationChecker import Checker


class suMD1(mwInputParser):

    def __init__(self, par):
        super(mwInputParser, self).__init__()
        self.scores = None
        self.output_to_check = 0
        self.checkVals_1 = []
        self.checkVals_2 = []
        self.parameters = par
        self.bestWalker = None
        self.fails = 0
        self.cycle = len(os.listdir(f'{self.folder}/trajectories'))

    def run_mwSuMD(self):
        suMD1(self.parameters).countTraj_logTraj(self.output_to_check)
        if self.parameters['NumberCV'] == 1:
            if self.parameters['Metric_1']:
                x = 1
            else:
                x = 2
            if self.parameters[f'Transition_{x}'] == 'positive':
                self.output_to_check = 0
                while self.output_to_check < self.parameters[f'Cutoff_{x}']:
                    self.output_to_check = self.runProtocol()

            if self.parameters[f'Transition_{x}'] == 'negative':
                self.output_to_check = 10 ** 6
                while self.output_to_check > self.parameters[f'Cutoff_{x}']:
                    self.output_to_check = self.runProtocol()
        else:
            if self.parameters['Transition_1'] == 'positive' and self.parameters['Transition_2'] == 'positive':
                self.metric_1, self.metric_2 = 0, 0
                while self.metric_1 < self.parameters['Cutoff_1'] and self.metric_2 < self.parameters['Cutoff_2']:
                    self.output_to_check[0], self.output_to_check[1] = self.runProtocol()

            if self.parameters['Transition_1'] == 'positive' and self.parameters['Transition_2'] == 'negative':
                self.metric_1, self.metric_2 = 0, 10 ** 6
                while self.metric_1 < self.parameters['Cutoff_1'] and self.metric_2 > self.parameters['Cutoff_2']:
                    self.output_to_check[0], self.output_to_check[1] = self.runProtocol()

            if self.parameters['Transition_1'] == 'negative' and self.parameters['Transition_2'] == 'negative':
                self.metric_1, self.metric_2 = 10 ** 6, 10 ** 6
                while self.metric_1 > self.parameters['Cutoff_1'] and self.metric_2 > self.parameters['Cutoff_2']:
                    self.output_to_check[0], self.output_to_check[1] = self.runProtocol()

            if self.parameters['Transition_1'] == 'negative' and self.parameters['Transition_2'] == 'positive':
                self.metric_1, self.metric_2 = 10 ** 6, 0
                while self.metric_1 > self.parameters['Cutoff_1'] and self.metric_2 < self.parameters['Cutoff_2']:
                    self.output_to_check[0], self.output_to_check[1] = self.runProtocol()

    def runProtocol(self):
        mwInputParser().getRestartOutput()
        if self.parameters['NumberCV'] == 1:
            self.bestWalker, self.scores, self.output_to_check = ProtocolRunner().runStandardProtocol()
            self.checkVals_1.append(self.output_to_check)
        else:
            self.bestWalker, self.scores, self.output_to_check = ProtocolRunner().runStandardProtocol()
            self.checkVals_1.append(self.output_to_check[0])
            self.checkVals_2.append(self.output_to_check[1])
        self.cycle += 1

        if self.cycle % int(self.parameters['RelaxTime'] * 1000 / self.parameters['Timewindow']) == 0:
            checker = Checker()
            if self.parameters['NumberCV'] == 1:
                self.fails += checker.checkIfFailed(self.checkVals_1, self.fails)
            else:
                self.fails += checker.checkIfFailed(self.checkVals_1, self.checkVals_2, self.fails)
            self.checkVals_1.clear()
            self.checkVals_2.clear()

        return self.output_to_check
