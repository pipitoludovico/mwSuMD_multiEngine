import os

from .InputTemplates import Template
from .Parser import mwInputParser


class EngineInputs(mwInputParser):

    def __init__(self, par):
        super(mwInputParser, self).__init__()
        self.par = par
        self.timewindow = self.initialParameters['Timewindow']
        self.trajCount = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
        mwInputParser().getRestartOutput()

    def getInputFile(self):
        if self.initialParameters['Relax'] is True:
            self.par['Timewindow'] = self.par['RelaxTime'] * 1000
            print("\nTemporary changing the timewindow for relaxation protocol to: " + str(self.par['Timewindow']) + " ps.")

        restartInput = Template().inputFile
        if self.par['Restart'] == 'YES' or self.trajCount != 0:
            restartInput = [line.replace('system', 'restarts')
                            if ('bin' in line or 'extendedSystem' in line) and 'system' in line else line
                            for line in restartInput]

        self.par['Timewindow'] = self.timewindow
        return restartInput
