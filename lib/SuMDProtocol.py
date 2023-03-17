import os
import time

import lib.MDoperations
from .MDoperations import MDoperator
from .MDsettings import MDsetter
from .Metrics import MetricsParser
from .Parser import mwInputParser
from .Runners import Runner


class suMD1(mwInputParser):

    def __init__(self, par):
        super(suMD1, self).__init__()
        self.par = par
        self.bestWalker = None
        self.checkVals = []
        self.fails = 0

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
        # create input files per walker
        begin = time.perf_counter()
        MDsetter(self.par).createInput(self.trajCount)
        # running the simulations
        Runner(self.par).runMD()
        # getting the metrics and choose the best one
        self.walker_metrics = MetricsParser(self.par).getChosenMetrics()
        self.bestWalker, self.max_value, self.metric_1, self.metric_2 = MetricsParser(self.par).getBestWalker(
            self.walker_metrics)
        print("best walker: " + str(self.bestWalker) + " " + str(self.max_value))
        # update values and log them
        MDoperator(self.par, self.folder, self.trajCount).saveStep(self.bestWalker)
        self.trajCount += 1
        suMD1(self.par).countTraj_logTraj(self.max_value)
        self.checkVals.append(self.max_value)
        end = time.perf_counter()
        final = end - begin
        print("Cycle completed in:" + str(final))

        if self.trajCount % (self.par['RelaxTime']*1000/self.par['Timewindow']) == 0:
            print('Checking if failed')
            if lib.MDoperations.checkIfStuck(self.checkVals, self.fails) is True:
                self.fails += 1
                self.runRelaxationProtocol()
        self.checkVals.clear()

        if self.par['NumberCV'] == 1:
            return self.max_value
        else:
            return self.metric_1, self.metric_2

    def runRelaxationProtocol(self):
        self.par['Relax'] = True
        # controllare due metriche/score a distanza di RT/TimeWindow
        MDsetter(self.par).createInput(self.trajCount)
        os.chdir('tmp/walker_1')
        os.system(f'acemd3 --device 0 input_1_{self.trajCount}.inp 1> acemd.log')

        self.walker_metrics = MetricsParser(self.par).getChosenMetrics()

        MDoperator(self.par, self.folder, self.trajCount).saveStep(self.bestWalker)
        self.trajCount += 1
        self.bestWalker, self.max_value, self.metric_1, self.metric_2 = MetricsParser(self.par).getBestWalker(
            self.walker_metrics)
        suMD1(self.par).countTraj_logTraj((self.metric_1, self.metric_2))
        print("\nRelaxation Protocol Ended")
        self.par['Relax'] = False
        self.runProtocol()
