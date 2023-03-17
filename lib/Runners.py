import multiprocessing as mp
import os
import subprocess
import time

from .Parser import mwInputParser
from .Utilities import ProcessManager


class Runner(mwInputParser):

    def __init__(self, par):
        self.par = par
        super(Runner, self).__init__()
        self.walk_count = 1
        self.trajCount = 0

    def runMD(self):
        # we set the initial trajectory count:
        self.trajCount = len(os.listdir(f'{self.folder}/trajectories'))
        print("trajectories in folder: " + str(self.trajCount))
        if self.par['MDEngine'] == 'ACEMD':
            self.runACEMD()
        else:
            self.runGROMACS()
        with mp.Pool() as p:
            wrappers = p.map(self.wrap, range(1, self.par['Walkers'] + 1))
        wrappers.clear()

    def runGPU_batch(self, trajCount, walk_count, GPUbatch, queue):
        processes = []
        for GPU in GPUbatch:
            os.chdir('tmp/walker_' + str(walk_count))
            command = f'acemd3 --device {GPU} input_{walk_count}_{trajCount}.inp 1> acemd.log'
            print(command)
            process = subprocess.Popen(command, shell=True)
            processes.append(process)

            walk_count += 1
            os.chdir(self.folder)

        # Wait for all subprocesses to finish
        for process in processes:
            process.wait()
        for GPU in GPUbatch:
            queue.put((trajCount, walk_count, GPU))  # Notify completion

        return walk_count  # Return the updated walk_count value

    def runACEMD(self):
        # let's divide the available GPU in batches by the number of walkers
        manager = ProcessManager()
        GPUs = manager.getGPUids()
        # let's exclude the GPU id if we want to keep a GPU for other jobs
        for excluded in self.excludedGPUS:
            GPUs.remove(excluded)
        GPUbatches, idList = manager.createBatches(walkers=self.par['Walkers'], total_gpu_ids=GPUs)

        if self.mode == 'parallel':
            print("\nRunning parallel mode")
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
            # serial version
            start_time_serial = time.perf_counter()
            for GPUbatch in GPUbatches:
                for GPU in GPUbatch:
                    os.chdir('tmp/walker_' + str(self.walk_count))
                    os.system(f'acemd3 --device {GPU} input_{self.walk_count}_{self.trajCount}.inp 1> acemd.log')
                    os.chdir(self.folder)
                    self.walk_count += 1
            end_time_serial = time.perf_counter()
            final_time_serial = end_time_serial - start_time_serial
            print("Serial Final Time:")
            print(final_time_serial)
        print("\nMD Runs completed.")

    def runGROMACS(self):
        # to do
        for walker in range(1, self.par['Walkers'] + 1):
            os.chdir('tmp/walker_' + str(walker + 1))
            if self.par['PLUMED'] is not None:
                os.system(f'echo test {[walker - 1]} {self.trajCount}')

    def wrap(self, folder):
        import MDAnalysis as Mda
        from MDAnalysis import transformations
        os.chdir(f'{self.folder}/tmp/walker_' + str(folder))
        print(os.getcwd())
        ext = ('xtc', 'dcd')
        psf = None

        if self.par['Forcefield'] == 'CHARMM':
            psf = '../../system/%s' % self.par['PSF']
        elif self.par['Forcefield'] == 'AMBER':
            psf = '../../system/%s' % self.par['PRMTOP']

        traj_name = '%s_%s' % (self.par['Output'], self.trajCount)

        for trajectory in os.listdir(os.getcwd()):
            if trajectory.startswith(traj_name) and trajectory.endswith(ext):
                try:
                    u = Mda.Universe(psf, trajectory)
                    prot = u.select_atoms(f"{self.par['Wrap']}")
                    ag = u.atoms
                    workflow = (transformations.unwrap(ag),
                                transformations.center_in_box(prot),
                                transformations.wrap(ag, compound='fragments'))
                    u.trajectory.add_transformations(*workflow)

                    with Mda.Writer('wrapped_MDA.xtc', ag) as w:
                        for ts in u.trajectory:
                            if ts is not None:
                                w.write(ag)
                except:
                    print(f"{trajectory} was used for wrapping")
        os.chdir(self.folder)
