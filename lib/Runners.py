import multiprocessing as mp
import os
import subprocess
import time

import MDAnalysis as Mda
from MDAnalysis import transformations

from .Parser import mwInputParser
from .Utilities import ProcessManager


class Runner(mwInputParser):

    def __init__(self, par):
        self.par = par
        super(Runner, self).__init__()
        self.walk_count = 1
        self.trajCount = len([x for x in os.scandir(f'{self.folder}/trajectories')])

    def runMD(self):
        # we set the initial trajectory count:
        print('runMD working in ' + os.getcwd())
        print("Trajectory count: " + str(self.trajCount))
        if self.par['MDEngine'] == 'ACEMD':
            self.runSimulation()
        if self.par['MDEngine'] == 'GROMACS':
            self.prepareTPR()
            self.runSimulation()
        print("Wrapping results...")
        with mp.Pool() as p:
            p.map(self.wrap, range(1, self.par['Walkers'] + 1))

    def prepareTPR(self):
        if self.par['MDEngine'] == 'GROMACS' and self.par['MDP'] is not None:
            command = (f'gmx grompp -f {self.folder}/system/{(self.par["MDP"])}'
                       f' -c {self.folder}/system/{self.par["GRO"]}'
                       f' -t {self.folder}/system/{self.par["CPT"]}'
                       f' -p {self.folder}/system/{self.par["TOP"]} '
                       f'-o tmp/{self.par["Output"]}_{self.trajCount}.tpr')
            # subprocess.Popen(command)
            print(command)

    def runGPU_batch(self, trajCount, walk_count, GPUbatch, queue):
        processes = []
        for GPU in GPUbatch:
            os.chdir('tmp/walker_' + str(walk_count))
            if self.par['MDEngine'] == 'ACEMD':
                command = f'acemd3 --device {GPU} input_{walk_count}_{trajCount}.inp 1> acemd.log'
                print(command)
                process = subprocess.Popen(command, shell=True)
                processes.append(process)
                walk_count += 1
                os.chdir(self.folder)
            if self.par['MDEngine'] == 'GROMACS':
                # questo si fa in tmp fuori dai walker
                if self.par['MDP'] is not None:
                    command = f'gmx grompp -f {self.folder}/system/{(self.par["MDP"])}' \
                              f' -c {self.par["GRO"]} -t {self.par["CPT"]} -p {self.par["TOP"]}' \
                              f' -o {self.par["Output"]}_{self.trajCount}.tpr'
                    process = subprocess.Popen(command, shell=True)
                    print(process)
                    print(command)
                    process = subprocess.Popen(command, shell=True)
                    processes.append(process)
                    walk_count += 1
                    os.chdir(self.folder)
                # usare l'mdp come templato ma va in NVT e aggiustato il ts, + altri parametri da aggiustare
                if self.par['PLUMED'] is not None:
                    os.system(f'gmx mdrun -plumed {self.par["PLUMED"]} -deffnm {self.par["Output"]}_{self.trajCount}')
                else:
                    os.system(f'gmx mdrun -deffnm ../{self.par["Output"]}_{self.trajCount}')
                    # path out

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

        if self.mode == 'parallel':
            print("")
            print('#' * 200)
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
        print("")
        print('#' * 200)

    def wrap(self, folder):
        os.chdir(f'{self.folder}/tmp/walker_' + str(folder))
        ext = ('xtc', 'dcd')
        psf = None

        if self.par['MDEngine'] == 'ACEMD':
            if self.par['Forcefield'] == 'CHARMM':
                psf = '../../system/%s' % self.par['PSF']

            if self.par['Forcefield'] == 'AMBER':
                psf = '../../system/%s' % self.par['PRMTOP']

            for trajectory in os.listdir(os.getcwd()):
                if trajectory.startswith(self.par['Output']) and trajectory.endswith(ext):
                    u = Mda.Universe(psf, trajectory)
                    prot = u.select_atoms(f"{self.par['Wrap']}")
                    if len(prot.atoms) == 0:
                        print("your wrapping selection selected 0 atoms! using protein and name CA instead...")
                        prot = u.select_atoms('protein and name CA')
                    ag = u.atoms
                    workflow = (transformations.unwrap(ag),
                                transformations.center_in_box(prot),
                                transformations.wrap(ag, compound='fragments'))
                    u.trajectory.add_transformations(*workflow)

                    with Mda.Writer('wrapped.xtc', ag) as w:
                        for ts in u.trajectory:
                            if ts is not None:
                                w.write(ag)
            os.chdir(self.folder)

        if self.par['MDEngine'] == 'GROMACS':
            for trajectory in os.listdir(os.getcwd()):
                if trajectory.startswith(self.par['Output']) and trajectory.endswith(ext):
                    wrapCommand = f"gmx trjconv -s tmp/{self.par['Output']}_{self.trajCount}.tpr" \
                                  f" -f {trajectory} -o wrapped.xtc -pbc atom -center -select {self.par['Wrap']}"
                    subprocess.Popen(wrapCommand, shell=True)
            os.chdir(self.folder)
