import os
import sys
import time
import threading
from multiprocessing import Queue, Process

from .Parser import mwInputParser
from .Utilities import ProcessManager


class Runner(mwInputParser):

    def __init__(self, par):
        self.par = par
        super(Runner, self).__init__()
        self.walk_count = 1

    def runMD(self):
        if self.par['MDEngine'] == 'ACEMD':
            self.runACEMD(self.trajCount)
        else:
            self.runGROMACS(self.trajCount)

    def runACEMD(self, trajCount):
        manager = ProcessManager()
        # GPUs = manager.getGPUids()
        GPUs = [0, 1]  # test
        GPUbatches, idList = manager.createBatches(walkers=3, total_gpu_ids=GPUs)

        if sys.argv[1] == 'parallel':
            # Parallel version:
            start_time = time.perf_counter()
            t = threading.Thread(target=self.runGPU(GPUbatches))
            t.start()
            t.join()
            end_time = time.perf_counter()
            final_time = end_time - start_time
            print("Multiprocess Final Time:")
            print(final_time)

        else:
            # serial version
            start_time_serial = time.perf_counter()
            for GPUbatch in GPUbatches:
                for GPU in GPUbatch:
                    print(GPU, len(GPUbatch))
                    os.chdir('tmp/walker_' + str(self.walk_count))
                    print(os.getcwd())
                    # command = f'acemd3 --device {GPUbatch[process]}
                    # input_{walker_number}_{trajCount}.inp 1> acemd.log'
                    command = f'echo --device {GPU} input_{self.walk_count}_{trajCount}.inp 1> acemd.log'
                    print(command)
                    self.wrap()
                    os.chdir(self.folder)
                    self.walk_count += 1
            end_time_serial = time.perf_counter()
            final_time_serial = end_time_serial - start_time_serial
            print("Serial Final Time:")
            print(final_time_serial)
        print("DONE")

    def runGPU(self, GPUbatches):
        start_time_parallel = time.perf_counter()
        for GPUbatch in GPUbatches:
            for GPU in GPUbatch:
                os.chdir('tmp/walker_' + str(self.walk_count))
                print(os.getcwd())
                # # command = f'acemd3 --device {GPUbatch[process]} input_{walker_number}_{trajCount}.inp 1> acemd.log'
                command = f'echo --device {GPU} input_{self.walk_count}_{self.trajCount}.inp 1> acemd.log'
                print(command)
                try:
                    self.wrap()
                    os.chdir(self.folder)
                except:
                    print(f"Walker {self.walk_count} wrapping has failed. No results were produced.")
                self.walk_count += 1

            end_time_parallel = time.perf_counter()
            final_time_parallel = end_time_parallel - start_time_parallel
            print("Batch processing time: ")
            print(final_time_parallel)

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
