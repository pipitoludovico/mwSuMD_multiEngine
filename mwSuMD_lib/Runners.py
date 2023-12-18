import multiprocessing as mp
import time

from .GPUoperations import ProcessManager
from .MDoperations import *
from .TrajectoryOperator import *
from .Loggers import Logger

from signal import signal, SIGPIPE, SIG_DFL
from warnings import filterwarnings

filterwarnings(action='ignore')
signal(SIGPIPE, SIG_DFL)


class Runner(mwInputParser):

    def __init__(self, par):
        self.par = par
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
                raise Exception("No trajectory found. Check your tmp folder.")
        try:
            with mp.Pool() as p:
                p.map(trajOperator.wrap, range(1, self.par['Walkers'] + 1))
            p.close()
            p.join()
        except:
            p.close()
            p.join()
            exit()

    def runGPU_batch(self, trajCount, walk_count, GPUbatch, queue):
        processes = []
        for GPU in GPUbatch:
            os.chdir('tmp/walker_' + str(walk_count))
            command = self.lauchEngine(trajCount, walk_count, GPU,
                                       self.customProductionFile)
            process = subprocess.Popen(command, shell=True)
            processes.append(process)
            walk_count += 1
            os.chdir(self.folder)
            Logger.LogToFile('ad', self.trajCount, command)
        # Wait for all subprocesses to finish
        for process in processes:
            process.wait()
        for GPU in GPUbatch:
            queue.put((trajCount, walk_count, GPU))  # Notify completion
        return walk_count  # Return the updated walk_count value

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
        if self.initialParameters['Mode'] == 'parallel':
            Logger.LogToFile("ad", self.trajCount, "\n\n")
            Logger.LogToFile("ad", self.trajCount, '*' * 200)
            Logger.LogToFile("a", self.trajCount, "Running parallel mode")
            manager = mp.Manager()
            q = manager.Queue()
            start_time_parallel = time.perf_counter()
            walk_count = 1  # Initialize the variable
            results = []
            try:
                with mp.Pool(processes=len(lenIDs)) as pool:
                    for GPUbatch in GPUbatches:
                        results.append(
                            pool.apply_async(self.runGPU_batch, args=(self.trajCount, walk_count, GPUbatch, q)))
                        walk_count += len(GPUbatch)
                    for result in results:
                        result.wait()
                        q.get()
                pool.close()
                pool.join()
                self.trajCount += 1
                end_time_parallel = time.perf_counter()
                Logger.LogToFile("ad", self.trajCount,
                                 f"Time taken with multiprocessing: {end_time_parallel - start_time_parallel:.2f} seconds")
            except:
                pool.close()
                pool.join()
                exit()
        else:
            Logger.LogToFile("a", self.trajCount, "Running serial mode")
            start_time_serial = time.perf_counter()
            for GPUbatch in GPUbatches:
                for GPU in GPUbatch:
                    os.chdir('tmp/walker_' + str(self.walk_count))
                    Logger.LogToFile("ad", self.trajCount, "Running in " + os.getcwd())
                    command = self.lauchEngine(self.trajCount, self.walk_count, GPU,
                                               self.customProductionFile)
                    subprocess.Popen(command, shell=True).wait()
                    self.walk_count += 1
                    os.chdir(self.folder)
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
        taks_master = f'-nt {num_threads} -npme 1 -ntmpi {num_threads} -pin on -pme gpu -nb gpu -bonded gpu -update gpu'
        plumed = f'-plumed {self.par["PLUMED"]}' if self.par['PLUMED'] is not None else ''

        if self.par['MDEngine'] == 'GROMACS':
            MDoperator(self.initialParameters, self.folder).prepareTPR(walk_count, trajCount, customFile)
            if self.initialParameters['COMMAND'] is None and self.customProductionFile is None:
                command = f'gmx mdrun  -v {plumed} -deffnm {self.initialParameters["Output"]}_{trajCount}_{walk_count} -gpu_id {GPU} {taks_master} -pinoffset {(offset * GPU)} -nstlist {self.initialParameters["Timewindow"]} &> gromacs.log'
            elif self.initialParameters['COMMAND'] is not None and self.customProductionFile is None:
                command = f'{self.initialParameters["COMMAND"]} -gpu_id {GPU} -deffnm {self.par["Output"]}_{trajCount}_{walk_count} {plumed} &> gromacs.log'
            elif self.initialParameters['COMMAND'] is not None and self.customProductionFile is not None:
                command = f'{self.initialParameters["COMMAND"]} -gpu_id {GPU} -deffnm production {plumed} &> gromacs.log'

        elif self.par['MDEngine'] == 'ACEMD':
            if customFile is not None and self.initialParameters['COMMAND'] is not None:
                command = f'{self.initialParameters["COMMAND"]} --device {GPU} production.inp &> acemd.log'
            if customFile is not None and self.initialParameters['COMMAND'] is None:
                command = f'acemd3 --device {GPU} production.inp &> acemd.log'
            if customFile is None and self.initialParameters['COMMAND'] is None:
                command = f'acemd3 --device {GPU} input_{walk_count}_{trajCount}.inp &> acemd.log'

        elif self.par['MDEngine'] == 'NAMD':
            command = f'{self.initialParameters["COMMAND"]} +devices {GPU} input_{walk_count}_{trajCount}.namd 1> namd.log' \
                if customFile is not None and self.initialParameters['COMMAND'] is not None \
                else f'namd3 +p8 +devices {GPU} input_{walk_count}_{trajCount}.namd 1> namd.log' \
                if customFile is not None \
                else f'{self.initialParameters["COMMAND"]} +devices {GPU} input_{walk_count}_{trajCount}.namd 1> namd.log' \
                if self.initialParameters['COMMAND'] is not None \
                else f'namd3 +p8 +devices {GPU} input_{walk_count}_{trajCount}.namd 1> namd.log'
        return command
