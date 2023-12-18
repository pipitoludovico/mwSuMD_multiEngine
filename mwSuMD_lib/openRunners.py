import multiprocessing as mp
import time

from mwSuMD_lib.TrajectoryOperator import *
from mwSuMD_lib.openMMsetter import *
from mwSuMD_lib.GPUoperations import ProcessManager


class Runner(mwInputParser):

    def __init__(self, par):
        self.par = par
        super(mwInputParser, self).__init__()
        self.walk_count = 1
        self.trajCount = len([x for x in os.scandir(f'{str(self.initialParameters["Root"])}/trajectories')])
        self.customProductionFile = None

    def runAndWrap(self):
        walker_snapshot = self.par['Walkers']
        if self.par['Relax'] is True:
            self.par['Walkers'] = 1
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

        with mp.Pool() as p:
            p.map(trajOperator.wrap, range(1, self.par['Walkers'] + 1))
        p.close()
        p.join()
        self.par['Walkers'] = walker_snapshot

    def runSimulation(self):
        # let's divide the available GPU in batches by the number of walkers
        manager = ProcessManager()
        GPUs = manager.getGPUids()
        # let's exclude the GPU id if we want to keep a GPU for other jobs
        if self.initialParameters.get("EXCLUDED_GPUS"):
            for excluded in self.initialParameters.get("EXCLUDED_GPUS"):
                GPUs.remove(excluded)
        GPUbatches, idList = manager.createBatches(walkers=self.par['Walkers'], total_gpu_ids=GPUs)
        Logger.LogToFile('ad', self.trajCount, '#' * 200)
        if self.initialParameters['Mode'] == 'parallel':
            Logger.LogToFile("ad", self.trajCount, "\n\n")
            Logger.LogToFile("ad", self.trajCount, '*' * 200)
            Logger.LogToFile("a", self.trajCount, "Running parallel mode")
            runner = RunnerOPENMM()
            manager = mp.Manager()
            q = manager.Queue()
            start_time_parallel = time.perf_counter()
            walk_count = 1
            results = []
            with mp.Pool(processes=len(idList)) as pool:
                for GPUbatch in GPUbatches:
                    for GPU in GPUbatch:
                        results.append(pool.apply_async(runner.runOPENMM, args=(walk_count, GPU)))
                        walk_count += 1
                for res in results:
                    res.get()
                while not q.empty():
                    q.get()
            pool.close()
            pool.join()
            self.trajCount += 1
            end_time_parallel = time.perf_counter()
            Logger.LogToFile("ad", self.trajCount,
                             f"Time taken with multiprocessing: {end_time_parallel - start_time_parallel:.2f} seconds")
        else:
            Logger.LogToFile("a", self.trajCount, "Running serial mode")
            start_time_serial = time.perf_counter()
            for GPUbatch in GPUbatches:
                for GPU in GPUbatch:
                    runner = RunnerOPENMM()
                    runner.runOPENMM(self.walk_count, GPU)
                    self.walk_count += 1
            end_time_serial = time.perf_counter()
            final_time_serial = end_time_serial - start_time_serial
            Logger.LogToFile("ad", self.trajCount, f"Serial Final Time: {final_time_serial:.2f} seconds")
            self.trajCount += 1

        Logger.LogToFile("ad", self.trajCount, "\nMD Runs completed.\n" + "#" * 200)
