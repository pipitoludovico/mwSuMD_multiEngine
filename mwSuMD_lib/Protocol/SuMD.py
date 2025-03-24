import os

import numpy as np
import pandas as pd
import shutil

from mwSuMD_lib.Parsers.InputfileParser import mwInputParser
from mwSuMD_lib.Protocol import ProtocolSelector
from mwSuMD_lib.MDutils.SimulationChecker import Checker
from mwSuMD_lib.Utilities.Loggers import Logger

from warnings import filterwarnings

filterwarnings(action='ignore')


class suMD1(mwInputParser):
    pars = mwInputParser()
    parameters, selection_list, parameterFolderPath = pars.getSettings()
    settings_df = pd.DataFrame(sorted(list(parameters.items())), columns=['Setting', 'Parameter'])

    def __init__(self, openMM):
        super(mwInputParser, self).__init__()
        self.openMM = openMM
        self.checker = Checker(self.openMM)
        self.ProtocolRunner = ProtocolSelector.ProtocolRunner(self.openMM)
        self.scores = None
        self.output_to_check = 0
        self.checkVals_1 = []
        self.checkVals_2 = []
        self.bestWalker = None
        self.fails = 0
        self.condition = None
        self.cycle = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
        if self.cycle == 0:
            Logger.LogToFile('ad', '',
                             "If you want to use a custom input file for your engine of choice, please, place it in the './system' folder, and call it \n "
                             "production.inp/namd/mdp (according to your engine). MwSuMD will use that instead of the default file. \n"
                             "If you choose to do so, make sure it points to a folder named 'restart' for continuing from the last binaries!\n")
        Logger.PrintSettingsToFile("w", self.cycle, str(self.settings_df))
        Logger.countTraj_logTraj(self.initialParameters, self.selection_list)
        self.run_First()

    def run_First(self):
        # purging old tmp folders to start
        if os.path.exists("./tmp"):
            print("Found existing ./tmp. Purging the folder for the first cycle.")
            shutil.rmtree('./tmp', ignore_errors=True)

        if self.parameters['NumberCV'] == 1:
            self.output_to_check = self.runProtocol()
        if self.parameters['NumberCV'] == 2:
            self.metric_1, self.metric_2 = self.runProtocol()

    def UpdateCondition(self, x=None) -> None:
        if self.parameters['NumberCV'] == 1:
            if self.parameters[f'Transition_{x}'] == 'positive':
                self.condition = self.output_to_check > self.parameters[f'Cutoff_{x}']
            else:
                self.condition = self.output_to_check < self.parameters[f'Cutoff_{x}']
        if self.parameters['NumberCV'] == 2:
            if self.parameters['Transition_1'] == 'positive' and self.parameters['Transition_2'] == 'positive':
                self.condition = self.metric_1 > self.parameters['Cutoff_1'] and self.metric_2 > self.parameters[
                    'Cutoff_2']
            if self.parameters['Transition_1'] == 'positive' and self.parameters['Transition_2'] == 'negative':
                self.condition = self.metric_1 > self.parameters['Cutoff_1'] and self.metric_2 < self.parameters[
                    'Cutoff_2']
            if self.parameters['Transition_1'] == 'negative' and self.parameters['Transition_2'] == 'positive':
                self.condition = self.metric_1 < self.parameters['Cutoff_1'] and self.metric_2 > self.parameters[
                    'Cutoff_2']
            if self.parameters['Transition_1'] == 'negative' and self.parameters['Transition_2'] == 'negative':
                self.condition = self.metric_1 < self.parameters['Cutoff_1'] and self.metric_2 < self.parameters[
                    'Cutoff_2']

    def run_mwSuMD(self):
        x = 1
        if self.parameters['NumberCV'] == 1:
            try:
                if self.parameters['Metric_1'] and not self.initialParameters['Metric_2']:
                    x = 1
                if self.initialParameters['Metric_2'] and not self.initialParameters['Metric_1']:
                    x = 2
                while not self.condition and self.parameters['NumberCV'] == 1:
                    self.compareAndUpdateSettings()
                    self.output_to_check = self.runProtocol()
                    self.UpdateCondition(x)
            except Exception as e:
                print("Conditions changed: rerunning mwSuMD")

                print(repr(e))
                pars = mwInputParser()
                self.parameters, self.selection_list, self.parameterFolderPath = pars.getSettings()
                self.run_mwSuMD()

        if self.parameters['NumberCV'] == 2:
            try:
                while not self.condition and self.parameters['NumberCV'] == 2:
                    self.compareAndUpdateSettings()
                    self.metric_1, self.metric_2 = self.runProtocol()
                    self.UpdateCondition()
            except Exception as e:
                print("Conditions changed: rerunning mwSuMD")
                pars = mwInputParser()
                self.parameters, self.selection_list, self.parameterFolderPath = pars.getSettings()
                print(repr(e))
                self.run_mwSuMD()

        if not self.condition:
            self.run_mwSuMD()
        Logger.LogToFile('w', self.cycle,
                         "#" * 200 + "\nTHRESHOLD METRICS REACHED: FINAL RELAXATION PROTOCOL:\n" + "#" * 200)
        if not self.openMM:
            self.checker.relaxSystemMulti()
        else:
            self.checker.relaxSystem()

    def runProtocol(self):
        mwInputParser().getRestartOutput()
        self.bestWalker, self.scores, self.output_to_check = self.ProtocolRunner.runStandardProtocol()

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
                if self.parameters['NumberCV'] == 1:
                    self.fails += self.checker.checkIfFailed(self.checkVals_1, self.fails)
                else:
                    self.fails += self.checker.checkIfFailed(self.checkVals_1, self.checkVals_2, self.fails)
                self.checkVals_1.clear()
                self.checkVals_2.clear()
        else:
            if self.cycle % int(self.parameters['RelaxTime'] * 1000 / self.parameters['Timewindow']) == 0:
                if self.parameters['NumberCV'] == 1:
                    self.fails += self.checker.checkIfFailed(self.checkVals_1, self.fails)
                else:
                    self.fails += self.checker.checkIfFailed(self.checkVals_1, self.checkVals_2, self.fails)
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
        if not self.openMM:
            if self.cycle != 0:
                try:
                    if self.parameters['MDEngine'] != 'GROMACS':
                        del parametersSnapshot['coor'], parametersSnapshot['vel'], parametersSnapshot['xsc']
                        del self.parameters['coor'], self.parameters['vel'], self.parameters['xsc']
                    else:
                        del parametersSnapshot['CPT'], parametersSnapshot['GRO'], parametersSnapshot['TPR']
                        del self.parameters['CPT'], self.parameters['GRO'], self.parameters['TPR']
                except Exception as e:
                    Logger.LogToFile('w', self.cycle, "#" * 200 + f"\nSetting's dictionary has changed. " + "#" * 200)
                    Logger.LogToFile('w', self.cycle, "#" * 200 + f"\n{repr(e)}" + "#" * 200)

        if parametersSnapshot != self.parameters or selectionShapshot != self.selection_list:
            temp = set(self.selection_list) - set(selectionShapshot)
            changes = []
            for key, value in self.parameters.items():
                if key == "ACTUAL_DISTANCE":
                    self.parameters["ACTUAL_DISTANCE"] = parametersSnapshot['ACTUAL_DISTANCE']
                if key in ('coor', 'vel', ' xsc', 'ACTUAL_DISTANCE'):
                    continue
                else:
                    if parametersSnapshot[key] != value:
                        changes.append(f"{key}: {value}")
            try:
                if len(changes) > 0:
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
