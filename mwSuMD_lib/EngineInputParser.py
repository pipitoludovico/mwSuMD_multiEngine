import os
from warnings import filterwarnings
from .InputTemplates import Template
from .Parser import mwInputParser
from .Loggers import Logger

filterwarnings(action='ignore')


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
            Logger.LogToFile('ad', self.trajCount, "\nTemporary changing the timewindow for relaxation protocol to: " + str(self.par['Timewindow']) + " ps.")
        restartInput = Template().inputFile
        self.par['Timewindow'] = self.timewindow
        return restartInput
