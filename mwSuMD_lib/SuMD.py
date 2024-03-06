import os
import pandas as pd
from .Parser import mwInputParser
from .Protocol import ProtocolRunner
from .SimulationChecker import Checker
from .Loggers import Logger

from warnings import filterwarnings

filterwarnings(action='ignore')


class suMD1(mwInputParser):
    pars = mwInputParser()
    parameters, selection_list, parameterFolderPath = pars.getSettings()
    settings_df = pd.DataFrame(sorted(list(parameters.items())), columns=['Setting', 'Parameter'])

    def __init__(self):
        super(mwInputParser, self).__init__()
        self.scores = None
        self.output_to_check = 0
        self.checkVals_1 = []
        self.checkVals_2 = []
        self.bestWalker = None
        self.fails = 0
        self.cycle = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
        Logger.PrintSettingsToFile("w", self.cycle, str(self.settings_df))
        Logger.countTraj_logTraj(self.initialParameters, self.selection_list)

    def run_mwSuMD(self):
        x = 1
        condition = None
        if self.parameters['NumberCV'] == 1:
            self.output_to_check = self.runProtocol()
            try:
                if self.parameters['Metric_1'] and not self.initialParameters['Metric_2']:
                    x = 1
                if self.initialParameters['Metric_2'] and not self.initialParameters['Metric_1']:
                    x = 2
                if self.parameters[f'Transition_{x}'] == 'positive':
                    condition = self.output_to_check > self.parameters[f'Cutoff_{x}']
                    while not condition and self.parameters['NumberCV'] == 1:
                        self.compareAndUpdateSettings()
                        self.output_to_check = self.runProtocol()
                if self.parameters[f'Transition_{x}'] == 'negative':
                    condition = self.output_to_check < self.parameters[f'Cutoff_{x}']
                    while not condition and self.parameters['NumberCV'] == 1:
                        self.compareAndUpdateSettings()
                        self.output_to_check = self.runProtocol()
            except:
                self.run_mwSuMD()

        if self.parameters['NumberCV'] == 2:
            self.metric_1, self.metric_2 = self.runProtocol()
            try:
                if self.parameters['Transition_1'] == 'positive' and self.parameters['Transition_2'] == 'positive':
                    condition = self.metric_1 > self.parameters['Cutoff_1'] and self.metric_2 > self.parameters[
                        'Cutoff_2']
                    while not condition and self.parameters['NumberCV'] == 2:
                        self.compareAndUpdateSettings()
                        self.metric_1, self.metric_2 = self.runProtocol()

                if self.parameters['Transition_1'] == 'positive' and self.parameters['Transition_2'] == 'negative':
                    condition = self.metric_1 > self.parameters['Cutoff_1'] and self.metric_2 < self.parameters[
                        'Cutoff_2']
                    while not condition and self.parameters['NumberCV'] == 2:
                        self.compareAndUpdateSettings()
                        self.metric_1, self.metric_2 = self.runProtocol()

                if self.parameters['Transition_1'] == 'negative' and self.parameters['Transition_2'] == 'negative':
                    condition = self.metric_1 < self.parameters['Cutoff_1'] and self.metric_2 < self.parameters[
                        'Cutoff_2']
                    while not condition and self.parameters['NumberCV'] == 2:
                        self.compareAndUpdateSettings()
                        self.metric_1, self.metric_2 = self.runProtocol()

                if self.parameters['Transition_1'] == 'negative' and self.parameters['Transition_2'] == 'positive':
                    condition = self.metric_1 < self.parameters['Cutoff_1'] and self.metric_2 > self.parameters[
                        'Cutoff_2']
                    while not condition and self.parameters['NumberCV'] == 2:
                        self.compareAndUpdateSettings()
                        self.metric_1, self.metric_2 = self.runProtocol()
            except:
                self.run_mwSuMD()
        if not condition:
            self.run_mwSuMD()
        Logger.LogToFile('w', self.cycle,
                         "#" * 200 + "\nTHRESHOLD METRICS REACHED: FINAL RELAXATION PROTOCOL:\n" + "#" * 200)
        checker = Checker()
        checker.relaxSystem()

    def runProtocol(self):
        mwInputParser().getRestartOutput()
        self.bestWalker, self.scores, self.output_to_check = ProtocolRunner().runStandardProtocol()

        if self.initialParameters['NumberCV'] == 1:
            metrics = self.output_to_check
            self.checkVals_1.append(self.output_to_check)
        else:
            metrics = list(self.output_to_check.keys())
            self.checkVals_1.append(self.output_to_check[metrics[0]])
            self.checkVals_2.append(self.output_to_check[metrics[1]])
        self.cycle += 1

        if self.parameters['CheckEvery'] is not None:
            if self.cycle % int(self.initialParameters.get('CheckEvery')) == 0:
                checker = Checker()
                if self.parameters['NumberCV'] == 1:
                    self.fails += checker.checkIfFailed(self.checkVals_1, self.fails)
                else:
                    self.fails += checker.checkIfFailed(self.checkVals_1, self.checkVals_2, self.fails)
                self.checkVals_1.clear()
                self.checkVals_2.clear()
        else:
            if self.cycle % int(self.parameters['RelaxTime'] * 1000 / self.parameters['Timewindow']) == 0:
                checker = Checker()
                if self.parameters['NumberCV'] == 1:
                    self.fails += checker.checkIfFailed(self.checkVals_1, self.fails)
                else:
                    self.fails += checker.checkIfFailed(self.checkVals_1, self.checkVals_2, self.fails)
                self.checkVals_1.clear()
                self.checkVals_2.clear()

        if self.parameters['NumberCV'] == 1:
            return float(self.output_to_check)
        else:
            return float(self.output_to_check[metrics[0]]), float(self.output_to_check[metrics[1]])

    def compareAndUpdateSettings(self) -> None:
        parametersSnapshot = self.parameters.copy()
        selectionShapshot = self.selection_list.copy()
        self.selection_list.clear()
        self.parameters, self.selection_list, self.parameterFolderPath = self.pars.getSettings()
        tempParametersSnapshot = self.parameters.copy()
        if self.cycle != 0:
            try:
                del parametersSnapshot['coor'], parametersSnapshot['vel'], parametersSnapshot['xsc']
                del tempParametersSnapshot['coor'], tempParametersSnapshot['vel'], tempParametersSnapshot['xsc']
            except print("\n"):
                del parametersSnapshot['cpt'], parametersSnapshot['gro'], parametersSnapshot['tpr']
                del tempParametersSnapshot['cpt'], tempParametersSnapshot['gro'], tempParametersSnapshot['tpr']

        if parametersSnapshot != tempParametersSnapshot or selectionShapshot != self.selection_list:
            temp = set(self.selection_list) - set(selectionShapshot)
            selectionShapshot.clear()
            changes = []
            for key, value in self.parameters.items():
                if key == 'coor' or key == 'vel' or key == 'xsc':
                    continue
                else:
                    if key in parametersSnapshot and parametersSnapshot[key] != value:
                        changes.append(f"{key}: {value}")
            try:
                with open('walkerSummary.log', 'a') as walkCheck:
                    walkCheck.write(f"Changed settings during protocol: {changes} {temp}\n")
            except Exception as e:
                print(f"Error writing to file: {e}")
                print(f"changes: {changes}")
                print(f"temp: {temp}")

            Logger.PrintSettingsToFile("w", self.cycle, str(self.settings_df))
            changes.clear()
            del changes
            del temp
        del parametersSnapshot, selectionShapshot
