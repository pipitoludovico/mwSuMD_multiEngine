import multiprocessing as mp
import os
import subprocess
import sys
import time

from .Parser import mwInputParser
from .Utilities import ProcessManager


class Runner(mwInputParser):

    def __init__(self, par):
        self.par = par
        super(Runner, self).__init__()
        self.walk_count = 1
        self.trajCount = 0  # Initialize the variable

    def runMD(self):
        if self.par['MDEngine'] == 'ACEMD':
            self.runACEMD(self.trajCount)
        else:
            self.runGROMACS(self.trajCount)

    def runGPU_batch(self, trajCount, walk_count, GPUbatch, queue):
        print("runGPU_batch")
        print(GPUbatch)
        processes = []
        for GPU in GPUbatch:
            os.chdir('tmp/walker_' + str(walk_count))
            print(os.getcwd())
            command = f'acemd3 --device {GPU} input_{walk_count}_{trajCount}.inp 1> acemd.log'
            print(command)
            process = subprocess.Popen(command, shell=True)
            processes.append(process)
            try:
                os.chdir(self.folder)
            except:
                print(f"Walker {walk_count} wrapping has failed. No results were produced.")
            walk_count += 1

        # Wait for all subprocesses to finish
        for process in processes:
            process.wait()

        for GPU in GPUbatch:
            queue.put((trajCount, walk_count, GPU))  # Notify completion
        return walk_count  # Return the updated walk_count value

    def runACEMD(self, trajCount):
        manager = ProcessManager()
        GPUs = manager.getGPUids()
        for excluded in self.excludedGPUS:
            GPUs.remove(excluded)
        GPUbatches, idList = manager.createBatches(walkers=self.par['Walkers'], total_gpu_ids=GPUs)
        print(GPUbatches)

        if self.mode == 'parallel':
            manager = mp.Manager()
            q = manager.Queue()
            start_time_parallel = time.perf_counter()
            walk_count = 1  # Initialize the variable
            results = []
            with mp.Pool(processes=len(GPUbatches)) as pool:
                for GPUbatch in GPUbatches:
                    results.append(pool.apply_async(self.runGPU_batch, args=(trajCount, walk_count, GPUbatch, q)))
                    walk_count += len(GPUbatch)

                # Wait for all the processes to finish
                for result in results:
                    walk_count = result.get()

                # Consume queue until all batches have finished
                print(f"Consuming queue...")
                while not q.empty():
                    batch_id = q.get()
                    print(f"Batch {batch_id} finished.")
                    sys.stdout.flush()

                print(f"All batches finished.")
                sys.stdout.flush()

            pool.close()
            pool.join()

            end_time_parallel = time.perf_counter()
            print(f"Time taken with multiprocessing: {end_time_parallel - start_time_parallel:.2f} seconds")

        else:
            # serial version
            start_time_serial = time.perf_counter()
            for GPUbatch in GPUbatches:
                for GPU in GPUbatch:
                    os.chdir('tmp/walker_' + str(self.walk_count))
                    os.system(f'acemd3 --device {GPU} input_{self.walk_count}_{trajCount}.inp 1> acemd.log')
                    self.wrap()
                    os.chdir(self.folder)
                    self.walk_count += 1
            end_time_serial = time.perf_counter()
            final_time_serial = end_time_serial - start_time_serial
            print("Serial Final Time:")
            print(final_time_serial)
        print("\nMD Runs completed.")

    def runGROMACS(self, trajCount):
        # to do
        for walker in range(1, self.par['Walkers'] + 1):
            os.chdir('tmp/walker_' + str(walker + 1))
            if self.par['PLUMED'] is not None:
                os.system(f'echo test {[walker - 1]} {trajCount}')

    def wrap(self):
        ext = ('.xtc', '.dcd')
        psf = None
        from moleculekit.molecule import Molecule
        if self.par['Forcefield'] == 'CHARMM':
            psf = '../../system/%s' % self.par['PSF']

        elif self.par['Forcefield'] == 'AMBER':
            psf = '../../system/%s' % self.par['PRMTOP']
        xtc = '%s_%s.xtc' % (self.par['Output'], self.trajCount)
        for file in os.listdir(os.getcwd()):
            if file == str(xtc) or file.endswith(ext):  # aggiungere
                try:
                    mol = Molecule(psf)
                    mol.read(file)
                    mol.wrap()
                    mol.wrap(self.par['Wrap'])
                    mol.write('wrapped.xtc')
                except:
                    print(f"{file} was used for wrapping")
