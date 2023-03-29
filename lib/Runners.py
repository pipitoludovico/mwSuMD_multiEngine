import multiprocessing as mp
import time

from .MDoperations import *
from .TrajectoryOperator import *
from .Utilities import ProcessManager


class Runner(mwInputParser):

    def __init__(self, par):
        self.par = par
        super(Runner, self).__init__()
        self.walk_count = 1
        self.trajCount = len([x for x in os.scandir(f'{self.folder}/trajectories')])
        self.customProductionFile = None

    def runMD(self):
        for file in os.listdir(f'{self.folder}/system'):
            if file == 'production.mdp' or file == 'production.namd' or file == 'production.inp':
                self.customProductionFile = file
        print('mwSuMD is working in ' + os.getcwd())
        print("Trajectory count: " + str(self.trajCount))
        self.runSimulation()

        print("Wrapping results...")
        trajOperator = TrajectoryOperator()

        with mp.Pool() as p:
            p.map(trajOperator.wrap, range(1, self.par['Walkers'] + 1))
            p.terminate()
            p.join()

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
            print(command)
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
        for excluded in self.excludedGPUS:
            GPUs.remove(excluded)
        GPUbatches, idList = manager.createBatches(walkers=self.par['Walkers'], total_gpu_ids=GPUs)

        print('#' * 200)
        if self.initialParameters['Mode'] == 'parallel':
            print("\n\n")
            print('*' * 200)
            print("Running parallel mode")
            manager = mp.Manager()
            q = manager.Queue()
            start_time_parallel = time.perf_counter()
            walk_count = 1  # Initialize the variable
            results = []
            with mp.Pool(processes=len(GPUbatches)) as pool:
                for GPUbatch in GPUbatches:
                    results.append(pool.apply_async(self.runGPU_batch, args=(self.trajCount, walk_count, GPUbatch, q)))
                    walk_count += len(GPUbatch)
                for result in results:
                    result.get()
                print(f"Waiting for all processes to finish...")
                while not q.empty():
                    q.get()
                print(f"All batches finished.")
            pool.close()
            pool.terminate()
            self.trajCount += 1
            end_time_parallel = time.perf_counter()
            print(f"Time taken with multiprocessing: {end_time_parallel - start_time_parallel:.2f} seconds")
        else:
            print("Running serial mode")
            start_time_serial = time.perf_counter()
            for GPUbatch in GPUbatches:
                for GPU in GPUbatch:
                    os.chdir('tmp/walker_' + str(self.walk_count))
                    print("we are in " + os.getcwd())
                    command = self.lauchEngine(self.trajCount, self.walk_count, GPU,
                                               self.customProductionFile)
                    subprocess.Popen(command, shell=True).wait()
                    self.walk_count += 1
                    os.chdir(self.folder)
            end_time_serial = time.perf_counter()
            final_time_serial = end_time_serial - start_time_serial
            print("Serial Final Time:")
            print(final_time_serial)
            self.trajCount += 1

        print("\nMD Runs completed.")
        print('#' * 200)

    def lauchEngine(self, trajCount, walk_count, GPU, customFile=None):
        command = ''
        if self.par['MDEngine'] == 'GROMACS':
            if customFile is not None:
                MDoperator(self.initialParameters, self.folder).prepareTPR(walk_count, trajCount, customFile)
                if self.par['PLUMED'] is not None:
                    command = f'%s {self.par["PLUMED"]} -deffnm production' % "gmx mdrun -plumed" \
                        if self.initialParameters['COMMAND'] is None else f'{self.initialParameters["COMMAND"]}'
                else:
                    command = f'%s -deffnm production' % "gmx mdrun " \
                        if self.initialParameters['COMMAND'] is None else f'{self.initialParameters["COMMAND"]}'
            else:
                MDoperator(self.initialParameters, self.folder).prepareTPR(walk_count, trajCount)
                if self.par['PLUMED'] is not None:
                    command = f'%s ' \
                              f'-plumed {self.par["PLUMED"]} ' \
                              f'-deffnm {self.par["Output"]}_{trajCount}_{walk_count}' % "gmx mdrun" \
                        if self.initialParameters['COMMAND'] is None else f'{self.initialParameters["COMMAND"]}'
                else:
                    command = f'%s -gpu_id {GPU} -deffnm {self.par["Output"]}_{trajCount}_{walk_count}' % "gmx mdrun" \
                        if self.initialParameters['COMMAND'] is None else f'{self.initialParameters["COMMAND"]}'

        elif self.par['MDEngine'] == 'ACEMD':
            if customFile is not None:
                command = f'%s --device {GPU} production.inp 1> acemd.log' % "acemd3" \
                    if self.initialParameters['COMMAND'] is None else f'{self.initialParameters["COMMAND"]}'
            else:
                command = f'%s --device {GPU} input_{walk_count}_{trajCount}.inp 1> acemd.log' % "acemd3" \
                    if self.initialParameters['COMMAND'] is None else f'{self.initialParameters["COMMAND"]}'

        elif self.par['MDEngine'] == 'NAMD':
            if customFile is not None:
                command = f'%s +devices {GPU} production.namd 1> namd.log' % "namd3 +auto-provision +setcpuaffinity" \
                    if self.initialParameters['COMMAND'] is None else f'{self.initialParameters["COMMAND"]}'
            else:
                command = f'%s +devices {GPU} input_{walk_count}_{trajCount}.inp 1> namd.log' % \
                          "namd3 +auto-provision +setcpuaffinity" \
                    if self.initialParameters['COMMAND'] is None else f'{self.initialParameters["COMMAND"]}'
        return command
