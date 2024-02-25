import os
import re
import subprocess
import numpy as np
from random import choice
from .MDsettings import MDsetter
from .Metrics import MetricsParser
from .GPUoperations import ProcessManager
from .TrajectoryOperator import TrajectoryOperator
from .Loggers import Logger
from warnings import filterwarnings

filterwarnings(action='ignore')


class MDoperator:
    def __init__(self, par, root):
        self.par = par
        self.folder = root
        self.extensions = ('.coor', '.xsc', '.vel', '.gro', '.cpt', '.tpr')
        self.cycle = len(
            [trajFile for trajFile in os.listdir(f'{self.folder}/trajectories') if trajFile.endswith('xtc')])
        self.best_metric_result = None
        self.best_average_metric_2 = None
        self.best_average_metric_1 = None
        self.best_walker_score = None
        self.averages = None
        self.scores = None
        self.best_value = None
        self.bestWalker = None
        self.walker_metrics = []

    def saveStep(self, best_walker, walker_score, best_metric_result):
        """Handles the restart files and the binary storage for OPENMM"""
        os.chdir('tmp/walker_%s' % str(best_walker))
        check = [binary for binary in os.listdir("./") if binary.endswith(self.extensions)]
        if check:
            self.cycle += 1
            # moving the best frame to the trajectory folder
            os.system(f'cp wrapped.xtc {self.folder}/trajectories/{self.par["Output"]}_step_{self.cycle}.xtc')
            if self.par['MDEngine'] != 'GROMACS':
                os.system(f'cp *.coor {self.folder}/restarts/previous.coor')
                os.system(f'cp *.xsc {self.folder}/restarts/previous.xsc')
                os.system(f'cp *.vel {self.folder}/restarts/previous.vel')
            if self.par['MDEngine'] == 'GROMACS':
                os.system(f'cp *.gro {self.folder}/restarts/previous.gro')
                os.system(f'cp "$(ls -t *.cpt | head -1)" {self.folder}/restarts/previous.cpt')
                os.system(f'cp *.tpr {self.folder}/restarts/previous.tpr')
            if self.par['PLUMED']:
                with open(self.par['PLUMED'], 'r') as plumedINP:
                    for line in plumedINP.readlines():
                        match = re.search(r'FILE=([^\s]+)', line)
                        if match:
                            outFile = match.group(1)
                            Logger.LogToFile('ad', self.cycle, str(os.getcwd() + outFile))
                            os.system(f'cp {outFile} %s/restarts/' % self.folder)
                for filename in os.listdir("./"):
                    if '.' not in filename:
                        fullname = os.path.join("./", filename)
                        os.system(f'cp {fullname}  %s/restarts/ ' % self.folder)
                    if filename.endswith(".dat"):
                        os.system(f'cp {filename} %s/restarts/' % self.folder)
                os.system('cp plumed.log  %s/restarts/ ' % self.folder)
        else:
            Logger.LogToFile('ad', self.cycle, "No binary saved: restarting from last checkpoint.")
            with open('walkerSummary.log', 'a') as walkerSummary:
                walkerSummary.write(f"No binary produced. Simulation failed at cycle: {self.cycle}")

        os.chdir(self.folder)

        with open('walkerSummary.log', 'a') as walkerSummary:
            if self.par['NumberCV'] == 1:
                if self.par['Relax'] and check:
                    info_to_write = str(self.cycle) + " RELAXATION SCORE: " + str(
                        round(walker_score, 3)) + " Metrics: " + str(
                        best_metric_result) + "\n"
                if not self.par['Relax'] and check:
                    info_to_write = str(self.cycle) + " Best Walker: " + str(best_walker) + " Best Metric: " + str(
                        round(walker_score, 3)) + " Last Metric: " + str(best_metric_result) + "\n"
            if self.par['NumberCV'] == 2:
                if self.par['Relax'] and check:
                    info_to_write = str(self.cycle) + " RELAXATION SCORE: " + str(
                        round(walker_score, 3)) + " Metrics: " + str(best_metric_result) + "\n"
                if not self.par['Relax'] and check:
                    info_to_write = str(self.cycle) + " Best Walker: " + str(best_walker) + " Score Result: " + str(
                        round(walker_score, 3)) + " Last Metrics from best: " + str(best_metric_result) + "\n"
            walkerSummary.write(info_to_write)

        os.system('rm -r tmp')
        self.par['Relax'] = False

    def prepareTPR(self, walk_count, trajcount, customFile=None):
        gro = f'{self.par["Root"]}/system/{self.par["GRO"]}' if self.cycle == 0 else f'{self.par["Root"]}/restarts/previous.gro'
        cpt = f'{self.par["Root"]}/system/{self.par["CPT"]}' if self.cycle == 0 else f'{self.par["Root"]}/restarts/previous.cpt'
        if customFile is None:
            command = f'gmx grompp -f input_{walk_count}_{trajcount}.mdp' \
                      f' -c {gro} -t {cpt} -p {self.folder}/system/{self.par["TOP"]}' \
                      f' -o {self.par["Output"]}_{trajcount}_{walk_count}.tpr -maxwarn 3 &>tpr_log.log'
        else:
            command = f'gmx grompp -f input_{walk_count}_{trajcount}.mdp' \
                      f' -c {gro} -t {cpt} -p {self.folder}/system/{self.par["TOP"]}' \
                      f' -o production.tpr -maxwarn 3 &>tpr_log.log'

        tprPreparation = subprocess.Popen(command, shell=True)
        tprPreparation.wait()

    def checkIfStuck(self, values, accumulatedFails) -> bool:
        if accumulatedFails >= self.par['Fails'] * int(self.par['NumberCV']):
            with open('walkerSummary.log', 'a') as logFile:
                logFile.write(
                    '\nSimulation seemed stuck. It will run the last relaxation protocol and it will be terminated\n')
                logFile.close()
            self.relaxSystem()
            exit()
        else:
            failCount = 0
            for vals in values:
                if vals is not None:
                    x = np.array(vals)
                    x_norm = (x - np.min(x)) / (np.max(x) - np.min(x))
                    if np.std(x_norm) < self.par['Tolerance']:
                        failCount += 1
            if failCount == self.par['NumberCV']:
                print('\nSimulation might be stuck. Running the relaxation protocol.')
                print("")
                return True
            else:
                # if it did not, we continue as usual
                return False

    def relaxSystem(self):
        Logger.LogToFile('ad', self.cycle, 'Relaxation Protocol begins now:\n' + ('#' * 200))
        self.par['Relax'] = True
        manager = ProcessManager()
        GPUs = manager.getGPUids()
        # let's exclude the GPU id if we want to keep a GPU for other jobs
        if self.par['EXCLUDED_GPUS']:
            Logger.LogToFile('ad', self.cycle, f"EXCLUDED GPU: {self.par['EXCLUDED_GPUS']}")
            for excluded in self.par['EXCLUDED_GPUS']:
                GPUs.remove(excluded)
        strGPU = map(str, GPUs)
        # jointGPUs = ','.join(strGPU)
        jointGPUs = str(choice(GPUs))
        # we create a special input file that has a longer runtime (5ns default or user-defined)
        MDsetter(self.par).createInputFile()
        # we run this inside walker_1 for convenience
        os.chdir('tmp/walker_1')
        for file in os.listdir(os.getcwd()):
            if file.endswith('.inp'):
                subprocess.Popen(f'acemd --device {jointGPUs} {file} 1> relax.log', shell=True).wait()
            elif file.endswith('.namd'):
                subprocess.Popen(f'namd3 +p8 +devices {jointGPUs} {file} 1> relax.log', shell=True).wait()
            elif file.endswith('.mdp'):
                subprocess.Popen(
                    f'gmx convert-tpr -s {self.folder}/restarts/previous.tpr -extend {int(self.par["RelaxTime"] * 1000)} -o {self.par["Output"]}_{self.cycle}.tpr &>tpr_log.log',
                    shell=True).wait()
                command = f'gmx mdrun -deffnm {self.par["Output"]}_{self.cycle}'
                subprocess.Popen(command, shell=True).wait()
        os.chdir(f'{self.folder}')
        TrajectoryOperator().wrap(1)
        # we then compute its metrics as a reference
        self.scores = MetricsParser().getChosenMetrics()
        # then the last coordinate is saved
        # we then extract the best metric/score and store it as a reference
        self.bestWalker, self.best_walker_score, self.best_metric_result = None, None, None
        self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(self.scores)

        MDoperator(self.par, self.folder).saveStep(1, self.best_walker_score, self.best_metric_result)

        Logger.LogToFile('ad', self.cycle, "\nRelaxation Protocol Ended\n" + "#" * 200)
        self.cycle += 1
        # setting our check to False and end the protocol, beginning a new cycle.
        self.par['Relax'] = False
