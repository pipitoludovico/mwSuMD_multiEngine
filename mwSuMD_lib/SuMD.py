import os
import pandas as pd

from .Parser import mwInputParser
from .Protocol import ProtocolRunner
from .SimulationChecker import Checker


class suMD1(mwInputParser):
    pars = mwInputParser()
    parameters, selection_list, parameterFolderPath = pars.getSettings()
    if not parameters.get('NumberCV'):
        print("Please set the NumberCV in the input file")
    settings_df = pd.DataFrame(sorted(list(parameters.items())), columns=['Setting', 'Parameter'])
    print(settings_df)

    def __init__(self):
        super(mwInputParser, self).__init__()
        self.scores = None
        self.output_to_check = 0
        self.checkVals_1 = []
        self.checkVals_2 = []
        self.bestWalker = None
        self.fails = 0
        self.cycle = len(os.listdir(f'{self.folder}/trajectories'))

    def run_mwSuMD(self):

        suMD1().countTraj_logTraj(self.output_to_check)
        if self.parameters['NumberCV'] == 1:
            if self.parameters['Metric_1']:
                x = 1
            else:
                x = 2
            if self.parameters[f'Transition_{x}'] == 'positive':
                self.output_to_check = 0
                while self.output_to_check < self.parameters[f'Cutoff_{x}']:
                    self.compareAndUpdateSettings()
                    self.output_to_check = self.runProtocol()

            if self.parameters[f'Transition_{x}'] == 'negative':
                self.output_to_check = 10 ** 6
                while self.output_to_check > self.parameters[f'Cutoff_{x}']:
                    self.compareAndUpdateSettings()
                    self.output_to_check = self.runProtocol()
        else:
            if self.parameters['Transition_1'] == 'positive' and self.parameters['Transition_2'] == 'positive':
                self.metric_1, self.metric_2 = 0, 0
                while self.metric_1 < self.parameters['Cutoff_1'] and self.metric_2 < self.parameters['Cutoff_2']:
                    self.compareAndUpdateSettings()
                    self.metric_1, self.metric_2 = self.runProtocol()

            if self.parameters['Transition_1'] == 'positive' and self.parameters['Transition_2'] == 'negative':
                self.metric_1, self.metric_2 = 0, 10 ** 6
                while self.metric_1 < self.parameters['Cutoff_1'] and self.metric_2 > self.parameters['Cutoff_2']:
                    self.compareAndUpdateSettings()
                    self.metric_1, self.metric_2 = self.runProtocol()

            if self.parameters['Transition_1'] == 'negative' and self.parameters['Transition_2'] == 'negative':
                self.metric_1, self.metric_2 = 10 ** 6, 10 ** 6
                while self.metric_1 > self.parameters['Cutoff_1'] and self.metric_2 > self.parameters['Cutoff_2']:
                    self.compareAndUpdateSettings()
                    self.metric_1, self.metric_2 = self.runProtocol()

            if self.parameters['Transition_1'] == 'negative' and self.parameters['Transition_2'] == 'positive':
                self.metric_1, self.metric_2 = 10 ** 6, 0
                while self.metric_1 > self.parameters['Cutoff_1'] and self.metric_2 < self.parameters['Cutoff_2']:
                    self.compareAndUpdateSettings()
                    self.metric_1, self.metric_2 = self.runProtocol()
        print("#" * 200 + "\nTHRESHOLD METRICS REACHED: FINAL RELAXATION PROTOCOL:\n" + "#" * 200)
        Checker().relaxSystem()

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

    def compareAndUpdateSettings(self) -> None:
        parametersSnapshot = self.parameters.copy()
        selectionShapshot = self.selection_list.copy()
        self.selection_list.clear()
        self.parameters, self.selection_list, self.parameterFolderPath = self.pars.getSettings()
        print("snapshot: ", selectionShapshot)
        print("list: ", self.selection_list)
        if parametersSnapshot != self.parameters or selectionShapshot != self.selection_list:
            temp = set(self.selection_list) - set(selectionShapshot)
            selectionShapshot.clear()
            changes = []
            for key, value in self.parameters.items():
                if key in parametersSnapshot and parametersSnapshot[key] != value:
                    changes.append(f"{key}: {value}")
            with open('walkerSummary.log', 'a') as walkCheck:
                walkCheck.write(f"Changed settings during protocol: {changes} {temp}\n")
            changes.clear()
            del changes
            del temp
        del parametersSnapshot, selectionShapshot
