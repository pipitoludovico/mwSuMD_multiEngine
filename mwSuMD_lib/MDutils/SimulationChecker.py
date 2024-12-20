import os
import subprocess
from subprocess import DEVNULL
from random import choice

from mwSuMD_lib.MetricOperators.MDoperations import MDoperator
from mwSuMD_lib.MDsetters.MDsettings import MDsetter
from mwSuMD_lib.MetricOperators.Metrics import MetricsParser
from mwSuMD_lib.Parsers.InputfileParser import mwInputParser
from mwSuMD_lib.Protocol.Runners import Runner
from mwSuMD_lib.Utilities.ProcessAndGPUutilities import ProcessManager
from mwSuMD_lib.MDutils.TrajectoryWrapper import TrajectoryOperator
from mwSuMD_lib.Utilities.Loggers import Logger
from warnings import filterwarnings

filterwarnings(action='ignore')


class Checker(mwInputParser):
    def __init__(self, openMM):
        super(mwInputParser, self).__init__()
        self.openMM = openMM
        self.best_metric_result = None
        self.best_average_metric_2 = None
        self.best_average_metric_1 = None
        self.best_walker_score = None
        self.averages = None
        self.scores = None
        self.best_value = None
        self.bestWalker = None
        self.trajCount = len([trajFile for trajFile in os.listdir(f'{self.folder}/trajectories') if trajFile.endswith('xtc')])

    def checkIfFailed(self, vals1=None, vals2=None, accumulatedFails=0):
        Logger.LogToFile("w", self.trajCount, "#" * 200 + '\nChecking if trajectory is stuck with values: ' + str(vals1) + ". Total fails accumulated: " + str(accumulatedFails) + "\n" + "#" * 200)
        mdOperator = MDoperator(self.initialParameters, self.folder, self.openMM)
        if mdOperator.checkIfStuck([vals1, vals2], accumulatedFails) is True:
            Logger.LogToFile("ad", self.trajCount, "\nRUNNING RELAXATION PROTOCOL" + "#" * 200)
            if not self.openMM:
                self.relaxSystemMulti()
            else:
                self.relaxSystem()
            accumulatedFails += 1
        else:
            accumulatedFails += 0
            Logger.LogToFile("ad", self.trajCount, "Number of fails accumulated: " + str(accumulatedFails) + "\n")
        return accumulatedFails

    def relaxSystemMulti(self):
        self.trajCount = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
        Logger.LogToFile('ad', self.trajCount, 'Relaxation Protocol begins now:\n' + ('#' * 200))
        # The relaxation protocol starts here
        self.initialParameters['Relax'] = True
        manager = ProcessManager()
        if not self.initialParameters['NOGPU']:
            GPUs = manager.getGPUids()
        else:
            GPUs = self.initialParameters.get('NOGPU')
        # let's exclude the GPU id if we want to keep a GPU for other jobs
        if self.initialParameters['EXCLUDED_GPUS']:
            Logger.LogToFile('ad', self.trajCount, f"EXCLUDED GPU: {self.initialParameters['EXCLUDED_GPUS']}")
            for excluded in self.initialParameters['EXCLUDED_GPUS']:
                GPUs.remove(excluded)
        # strGPU = map(str, GPUs)  # gsd calls an error when using multiple GPUs. See https://github.com/openmm/openmm/issues/2535
        # jointGPUs = ','.join(strGPU)
        jointGPUs = choice(GPUs)
        # we create a special input file that has a longer runtime (5ns default or user-defined)
        MDsetter(self.initialParameters).createInputFile()
        # we run this inside walker_1 for convenience
        os.chdir('tmp/walker_1')
        for file in os.listdir(os.getcwd()):
            if file.endswith('.inp'):
                subprocess.Popen(f'acemd --device {jointGPUs} {file} 1> relax.log', shell=True).wait()
            elif file.endswith('.namd'):
                subprocess.Popen(f'namd3 +p8 +devices {jointGPUs} {file} 1> relax.log', shell=True).wait()
            elif file.endswith('.mdp'):
                subprocess.Popen(
                    f'gmx convert-tpr -s {self.folder}/restarts/previous.tpr -extend {int(self.initialParameters["RelaxTime"] * 1000)} -o {self.initialParameters["Output"]}_{self.trajCount}.tpr > tpr_log.log 2>&1',
                    shell=True, stdout=DEVNULL).wait()
                command = f'gmx mdrun -deffnm {self.initialParameters["Output"]}_{self.trajCount} > relax.log 2>&1'
                subprocess.Popen(command, shell=True, stdout=DEVNULL).wait()
        os.chdir(f'{self.folder}')
        if self.initialParameters['WrapEngine'] == 'MDA':
            TrajectoryOperator().wrap(1)
        else:
            TrajectoryOperator().wrapVMD(1)
        # we then compute its metrics as a reference
        self.scores = MetricsParser().getChosenMetrics()
        # then the last coordinate is saved
        # we then extract the best metric/score and store it as a reference
        self.bestWalker, self.best_walker_score, self.best_metric_result = None, None, None
        self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(self.scores)
        MDoperator(self.initialParameters, self.folder, self.openMM).saveStep(1, self.best_walker_score,
                                                                              self.best_metric_result)

        Logger.LogToFile('ad', self.trajCount, "\nRelaxation Protocol Ended\n" + "#" * 200)
        self.trajCount += 1
        # setting our check to False and end the protocol, beginning a new cycle.
        self.initialParameters['Relax'] = False

    def relaxSystem(self):
        Logger.LogToFile('a', self.trajCount, 'Relaxation Protocol begins now:\n' + ('#' * 200))
        self.initialParameters['Relax'] = True
        op = MDoperator(self.initialParameters, self.folder, self.openMM)
        Runner(self.initialParameters, self.openMM, op).runAndWrap()
        self.scores = MetricsParser().getChosenMetrics()
        # we then extract the best metric/score and store it as a reference
        self.bestWalker, self.best_walker_score, self.best_metric_result = None, None, None
        self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(self.scores)
        MDoperator(self.initialParameters, self.folder, self.openMM).saveStep(self.bestWalker, self.best_walker_score, self.best_metric_result)

        Logger.LogToFile('ad', self.trajCount, "\nRelaxation Protocol Ended\n" + "#" * 200)
        self.trajCount += 1
        # setting our check to False and end the protocol, beginning a new cycle.
        self.initialParameters['Relax'] = False
