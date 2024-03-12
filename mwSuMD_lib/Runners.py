import time
import concurrent.futures

from .MDoperations import *
from .TrajectoryOperator import *
from mwSuMD_lib.openMMsetter import *


from signal import signal, SIGPIPE, SIG_DFL
from warnings import filterwarnings

filterwarnings(action='ignore')
signal(SIGPIPE, SIG_DFL)


class Runner(mwInputParser):

    def __init__(self, par, openMM):
        self.par = par
        self.openMM = openMM
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
        trajOperator = TrajectoryOperator()

        for i in range(1, self.initialParameters['Walkers'] + 1):
            files = os.listdir(f"./tmp/walker_{i}")
            trajectory = next((file for file in files if file.endswith('.xtc') or file.endswith('dcd')), None)
            if trajectory:
                continue
            else:
                Logger.LogToFile('ad', self.trajCount, "No trajectory found. Check your tmp folder.")
                exit()
        try:
            for walker_wrap in range(1, self.par['Walkers'] + 1):
                trajOperator.wrap(walker_wrap)
        except Exception as e:
            print("Wrapping failed: ", e)
            exit()

    def runSimulation(self):
        # let's divide the available GPU in batches by the number of walkers
        manager = ProcessManager()
        GPUs = manager.getGPUids()
        # let's exclude the GPU id if we want to keep a GPU for other jobs
        if self.initialParameters['EXCLUDED_GPUS'] is not None:
            Logger.LogToFile("ad", self.trajCount, "\nEXCLUDED GPU:" + str(self.initialParameters['EXCLUDED_GPUS']))
            for excluded in self.initialParameters['EXCLUDED_GPUS']:
                GPUs.remove(excluded)
        GPUbatches, lenIDs = manager.createBatches(walkers=self.par['Walkers'], total_gpu_ids=GPUs)
        Logger.LogToFile('ad', self.trajCount, '#' * 200)
        processes = []
        Logger.LogToFile("a", self.trajCount, "Running serial mode")
        start_time_serial = time.perf_counter()
        for GPUbatch in GPUbatches:
            for GPU in GPUbatch:
                os.chdir('tmp/walker_' + str(self.walk_count))
                Logger.LogToFile("ad", self.trajCount, "Running in " + os.getcwd())
                command = self.lauchEngine(self.trajCount, self.walk_count, GPU,
                                           self.customProductionFile)
                processes.append(subprocess.Popen(command, shell=True))
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
                    fullname = os.path.join("./", filename)
                    plumedCopy += f'cp ../../restarts/{fullname} .;'
            if any(plumedFile.endswith(".dat") for plumedFile in os.listdir("../../restarts")):
                plumedCopy += "cp ../../restarts/*.dat .;"
        else:
            plumedCopy = ''

        if self.par['MDEngine'] == 'GROMACS':
            MDoperator(self.initialParameters, self.folder, self.openMM).prepareTPR(walk_count, trajCount, customFile)
            if self.initialParameters['COMMAND'] is None and self.customProductionFile is None:
                command = f'gmx mdrun  -v {plumed} -deffnm {self.initialParameters["Output"]}_{trajCount}_{walk_count} -gpu_id {GPU} {taks_master} -pinoffset {(offset * GPU)} -nstlist {self.initialParameters["Timewindow"]} > gromacs.log'
            elif self.initialParameters['COMMAND'] is not None and self.customProductionFile is None:
                command = f'{self.initialParameters["COMMAND"]} -gpu_id {GPU} -deffnm {self.par["Output"]}_{trajCount}_{walk_count} {plumed} > gromacs.log'
            elif self.initialParameters['COMMAND'] is not None and self.customProductionFile is not None:
                command = f'{self.initialParameters["COMMAND"]} -gpu_id {GPU} -deffnm production {plumed} > gromacs.log'

        if self.par['MDEngine'] == 'ACEMD':
            if customFile is not None and self.initialParameters['COMMAND'] is not None:
                command = f'{self.initialParameters["COMMAND"]} --device {GPU} production.inp > acemd.log'
            if customFile is not None and self.initialParameters['COMMAND'] is None:
                command = f'acemd --device {GPU} production.inp > acemd.log'
            if customFile is None and self.initialParameters['COMMAND'] is None:
                command = f'acemd --device {GPU} input_{walk_count}_{trajCount}.inp > acemd.log'

        if self.par['MDEngine'] == 'NAMD':
            command = f'{self.initialParameters["COMMAND"]} +devices {GPU} input_{walk_count}_{trajCount}.namd > namd.log' \
                if customFile is not None and self.initialParameters['COMMAND'] is not None \
                else f'namd3 +p8 +devices {GPU} input_{walk_count}_{trajCount}.namd > namd.log' \
                if customFile is not None \
                else f'{self.initialParameters["COMMAND"]} +devices {GPU} input_{walk_count}_{trajCount}.namd > namd.log' \
                if self.initialParameters['COMMAND'] is not None \
                else f'namd3 +p8 +devices {GPU} input_{walk_count}_{trajCount}.namd > namd.log'
        return plumedCopy + " " + command

    # OpenMM section
    def runAndWrap(self):
        walker_snapshot = self.initialParameters['Walkers']
        if self.initialParameters['Relax'] is True:
            self.initialParameters['Walkers'] = 1
        Logger.LogToFile('a', self.trajCount, "Trajectory count: " + str(self.trajCount))
        self.runSimulationOpen()

        Logger.LogToFile('ad', self.trajCount, "Wrapping results...")
        trajOperator = TrajectoryOperator()

        for i in range(1, self.initialParameters['Walkers'] + 1):
            files = os.listdir(f"./tmp/walker_{i}")
            trajectory = next((file for file in files if file.endswith('.xtc') or file.endswith('dcd')), None)
            if trajectory:
                continue
            else:
                raise Exception("No trajectory found. Check your tmp folder.")

        try:
            processes = []
            with concurrent.futures.ProcessPoolExecutor() as executor:
                for walker_wrap in range(1, self.initialParameters['Walkers'] + 1):
                    future = executor.submit(trajOperator.wrap, walker_wrap)
                    processes.append(future)
        except Exception as e:
            print("Wrapping failed: ", e)
            exit()
        self.initialParameters['Walkers'] = walker_snapshot

    def runSimulationOpen(self):
        # let's divide the available GPU in batches by the number of walkers
        manager = ProcessManager()
        GPUs = manager.getGPUids()
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
