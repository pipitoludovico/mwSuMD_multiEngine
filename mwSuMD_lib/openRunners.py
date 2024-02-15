import time
import concurrent.futures

from mwSuMD_lib.TrajectoryOperator import *
from mwSuMD_lib.openMMsetter import *
from mwSuMD_lib.GPUoperations import ProcessManager


class Runner(mwInputParser):

    def __init__(self):
        super(mwInputParser, self).__init__()
        self.walk_count = 1
        self.trajCount = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
        self.customProductionFile = None

    def runAndWrap(self):
        walker_snapshot = self.initialParameters['Walkers']
        if self.initialParameters['Relax'] is True:
            self.initialParameters['Walkers'] = 1
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
            for walker_wrap in range(1, self.initialParameters['Walkers'] + 1):
                trajOperator.wrap(walker_wrap)
        except Exception as e:
            print("Wrapping failed: ", e)
            exit()
        self.initialParameters['Walkers'] = walker_snapshot

    def runSimulation(self):
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

        def run_openMM_async(walk_count, GPU_m):
            return setter.runOPENMM(walk_count, GPU_m)

        with concurrent.futures.ThreadPoolExecutor() as executor:
            for GPUbatch in GPUbatches:
                for GPU in GPUbatch:
                    future = executor.submit(run_openMM_async, self.walk_count, GPU)
                    processes.append(future)
                    self.walk_count += 1

        # for GPUbatch in GPUbatches:
        #     for GPU in GPUbatch:
        #         processes.append(setter.runOPENMM(self.walk_count, GPU))
        #         self.walk_count += 1
        end_time_serial = time.perf_counter()
        final_time_serial = end_time_serial - start_time_serial
        Logger.LogToFile("ad", self.trajCount, f"Final Time: {final_time_serial:.2f} seconds")
        self.trajCount += 1

        Logger.LogToFile("ad", self.trajCount, "\nMD Runs completed.\n" + "#" * 200)
