import time
import concurrent.futures
import subprocess
import re
from mwSuMD_lib.Utilities.ProcessAndGPUutilities import ProcessManager
from mwSuMD_lib.MDutils.TrajectoryWrapper import *
from mwSuMD_lib.OpenMMbranch.openMMsetter import *

from signal import signal, SIGPIPE, SIG_DFL
from warnings import filterwarnings

filterwarnings(action='ignore')
signal(SIGPIPE, SIG_DFL)


class Runner(mwInputParser):

    def __init__(self, par, openMM, md_operator):
        self.par = par
        self.openMM = openMM
        self.md_operator = md_operator
        super(mwInputParser, self).__init__()
        self.walk_count = 1
        self.trajCount = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
        self.customProductionFile = None

    def runMD(self):
        if self.initialParameters['CUSTOMFILE'] is not None:
            self.customProductionFile = self.initialParameters['CUSTOMFILE']
        else:
            self.customProductionFile = None
        Logger.LogToFile('a', self.trajCount, "Trajectory count: " + str(self.trajCount))
        self.runSimulation()

        Logger.LogToFile('ad', self.trajCount, "Wrapping results...")
        for i in range(1, self.initialParameters['Walkers'] + 1):
            files = os.listdir(f"./tmp/walker_{i}")
            trajectory = next((file for file in files if file.endswith('.xtc') or file.endswith('dcd')), None)
            if trajectory:
                del files
            else:
                Logger.LogToFile('ad', self.trajCount, "No trajectory found. Check your tmp folder.")
                exit()
        try:
            for walker_wrap in range(1, self.par['Walkers'] + 1):
                if self.par['WrapEngine'] == 'MDA':
                    TrajectoryOperator().wrap(walker_wrap)
                else:
                    TrajectoryOperator().wrapVMD(walker_wrap)
        except Exception as e:
            print("Wrapping failed: ", e)
            exit()

    def runSimulation(self):
        # let's divide the available GPU in batches by the number of walkers
        manager = ProcessManager()
        if not self.initialParameters.get('NOGPU'):
            GPUs = manager.getGPUids()
        else:
            GPUs = self.initialParameters.get('NOGPU')
        # let's exclude the GPU id if we want to keep a GPU for other jobs
        if self.initialParameters['EXCLUDED_GPUS'] is not None:
            Logger.LogToFile("ad", self.trajCount, "\nEXCLUDED GPU:" + str(self.initialParameters['EXCLUDED_GPUS']))
            for excluded in self.initialParameters['EXCLUDED_GPUS']:
                GPUs.remove(excluded)
        GPUbatches, lenIDs = manager.createBatches(walkers=self.par['Walkers'], total_gpu_ids=GPUs)
        Logger.LogToFile('ad', self.trajCount, '#' * 200)
        processes = []
        start_time_serial = time.perf_counter()
        for GPUbatch in GPUbatches:
            for GPU in GPUbatch:
                os.chdir('tmp/walker_' + str(self.walk_count))
                command = self.lauchEngine(self.trajCount, self.walk_count, GPU, self.customProductionFile)
                processes.append(subprocess.Popen(command, shell=True, stdout=DEVNULL))
                self.walk_count += 1
                os.chdir(self.folder)
        for process in processes:
            process.wait()
        end_time_serial = time.perf_counter()
        final_time_serial = end_time_serial - start_time_serial
        Logger.LogToFile("ad", self.trajCount, f"Serial Final Time: {final_time_serial:.2f} seconds")
        self.trajCount += 1

        Logger.LogToFile("ad", self.trajCount, "\nMD Runs completed.\n" + "#" * 200)

    def lauchEngine(self, trajCount, walk_count, GPU, customFile=None):
        command = ''
        core_division = (int(len(os.sched_getaffinity(0)) / self.initialParameters["Walkers"]))
        num_threads = core_division - self.initialParameters["Walkers"]
        offset = num_threads + 1
        taks_master = f'-nt {num_threads} -pin on -pme gpu -nb gpu -bonded gpu -update gpu'
        plumed = f'-plumed {self.par["PLUMED"]}' if self.par['PLUMED'] is not None else ''

        if self.initialParameters.get("PLUMED") and self.trajCount != 0:
            plumedCopy = ''
            with open(self.par['PLUMED'], 'r') as plumedINP:
                for line in plumedINP.readlines():
                    match = re.search(r'FILE=([^\s]+)', line)
                    if match:
                        outFile = match.group(1).split(',')
                        for chunk in outFile:
                            plumedCopy += f"cp ../../restarts/{chunk} .; "

            for filename in os.listdir("../../restarts"):
                if '.' not in filename:
                    fullname = os.path.join("../../restarts", filename)
                    plumedCopy += f'cp {fullname} .;'
            if any(plumedFile.endswith(".dat") for plumedFile in os.listdir("../../restarts")):
                plumedCopy += "cp ../../restarts/*.dat .;"
        else:
            plumedCopy = ''

        if self.par['MDEngine'] == 'GROMACS':
            if not self.initialParameters.get('NOGPU'):
                gpuCall = f'-gpu_id {GPU}'
            else:
                gpuCall = ''
            self.md_operator.prepareTPR(walk_count, trajCount, customFile)
            # standard call: no custom file, no custom command
            if self.initialParameters['COMMAND'] is None and self.customProductionFile is None:
                command = f'gmx mdrun  -v {plumed} -deffnm {self.initialParameters["Output"]}_{trajCount}_{walk_count} {gpuCall} {taks_master} -pinoffset {(offset * GPU)} -nstlist {self.initialParameters["Timewindow"]} > gromacs.log 2>&1'
            # no custom file / custom command
            if self.initialParameters['COMMAND'] is not None and self.customProductionFile is None:
                command = f'{self.initialParameters["COMMAND"]} {gpuCall} -deffnm {self.par["Output"]}_{trajCount}_{walk_count} {plumed} > gromacs.log 2>&1'
            # custom file / no custom command
            if self.initialParameters['COMMAND'] is None and self.customProductionFile is not None:
                command = f'gmx mdrun  -v {plumed} {gpuCall} {taks_master} -pinoffset {(offset * GPU)} -nstlist {self.initialParameters["Timewindow"]} -deffnm production {plumed} > gromacs.log 2>&1'
            # custom file / custom command
            if self.initialParameters['COMMAND'] is not None and self.customProductionFile is not None:
                command = f'{self.initialParameters["COMMAND"]} {gpuCall} -deffnm production {plumed} > gromacs.log 2>&1'

        if self.par['MDEngine'] == 'ACEMD':
            if not self.initialParameters.get('NOGPU'):
                gpuCall = f'--device {GPU}'
            else:
                gpuCall = '--platfrom CPU'
            # standard call: no custom file, no custom command
            if customFile is None and self.initialParameters['COMMAND'] is None:
                command = f'acemd3 {gpuCall} input_{walk_count}_{trajCount}.inp > acemd.log'
            # no custom file / custom command
            if customFile is None and self.initialParameters['COMMAND'] is not None:
                command = f'{self.initialParameters["COMMAND"]} {gpuCall} input_{walk_count}_{trajCount}.inp > acemd.log'
            # custom file / no custom command
            if customFile is not None and self.initialParameters['COMMAND'] is None:
                command = f'acemd3 {gpuCall} production.inp > acemd.log'
            # custom file / custom command
            if customFile is not None and self.initialParameters['COMMAND'] is not None:
                command = f'{self.initialParameters["COMMAND"]} {gpuCall} production.inp > acemd.log'

        if self.par['MDEngine'] == 'NAMD':
            if not self.initialParameters.get('NOGPU'):
                gpuCall = f'+devices {GPU}'
            else:
                gpuCall = ''
            # standard call: no custom file, no custom command
            if customFile is None and self.initialParameters['COMMAND'] is None:
                command = f'namd3 +p8 {gpuCall} input_{walk_count}_{trajCount}.namd > namd.log'
            # no custom file / custom command
            if customFile is None and self.initialParameters['COMMAND'] is not None:
                command = f'{self.initialParameters["COMMAND"]} {gpuCall} input_{walk_count}_{trajCount}.namd > namd.log'
            # custom file / no custom command
            if customFile is not None and self.initialParameters['COMMAND'] is None:
                command = f'namd3 +p8 {gpuCall} {gpuCall} production.namd > namd.log'
            # custom file / custom command
            if customFile is not None and self.initialParameters['COMMAND'] is not None:
                command = f'{self.initialParameters["COMMAND"]} {gpuCall} production.namd > namd.log'
        return plumedCopy + " " + command

    # OpenMM section
    def runAndWrap(self):
        walker_snapshot = self.initialParameters['Walkers']
        if self.initialParameters['Relax'] is True:
            self.initialParameters['Walkers'] = 1
        Logger.LogToFile('a', self.trajCount, "Trajectory count: " + str(self.trajCount))
        self.runSimulationOpen()

        Logger.LogToFile('ad', self.trajCount, "Wrapping results...")

        for i in range(1, self.initialParameters['Walkers'] + 1):
            files = os.listdir(f"./tmp/walker_{i}")
            trajectory = next((file for file in files if file.endswith('.xtc') or file.endswith('dcd')), None)
            if trajectory:
                del files
            else:
                raise Exception("No trajectory found. Check your tmp folder.")

        try:
            processes = []
            if self.par['WrapEngine'] == 'MDA':
                wrapcall = TrajectoryOperator().wrap
            else:
                wrapcall = TrajectoryOperator().wrapVMD
            with concurrent.futures.ProcessPoolExecutor() as executor:
                for walker_wrap in range(1, self.initialParameters['Walkers'] + 1):
                    future = executor.submit(wrapcall, walker_wrap)
                    processes.append(future)
        except Exception as e:
            print("Wrapping failed: ", e)
            exit()
        self.initialParameters['Walkers'] = walker_snapshot

    def runSimulationOpen(self):
        # let's divide the available GPU in batches by the number of walkers
        manager = ProcessManager()
        if not self.initialParameters.get('NOGPU'):
            GPUs = manager.getGPUids()
        else:
            GPUs = self.initialParameters.get('NOGPU')
        # let's exclude the GPU id if we want to keep a GPU for other jobs
        if self.initialParameters['EXCLUDED_GPUS'] is not None:
            Logger.LogToFile('ad', self.trajCount, f"EXCLUDED GPU: {self.initialParameters['EXCLUDED_GPUS']}")
            for excluded in self.initialParameters['EXCLUDED_GPUS']:
                GPUs.remove(excluded)
        GPUbatches, idList = manager.createBatches(walkers=self.initialParameters['Walkers'], total_gpu_ids=GPUs)
        Logger.LogToFile('ad', self.trajCount, '#' * 200)
        start_time_serial = time.perf_counter()
        processes = []
        setter = openMMsetter(self.initialParameters)

        with concurrent.futures.ProcessPoolExecutor() as executor:
            for GPUbatch in GPUbatches:
                for GPU in GPUbatch:
                    folder_path = f'tmp/walker_{self.walk_count}'
                    future = executor.submit(setter.runOPENMM, self.walk_count, GPU, folder_path)
                    processes.append(future)
                    self.walk_count += 1

        end_time_serial = time.perf_counter()
        final_time_serial = end_time_serial - start_time_serial
        Logger.LogToFile("ad", self.trajCount, f"Final Time: {final_time_serial:.2f} seconds")
        self.trajCount += 1

        Logger.LogToFile("ad", self.trajCount, "\nMD Runs completed.\n" + "#" * 200)
