import os
import argparse

import GPUtil


class ProcessManager:
    def __init__(self):
        self.getpid()

    @staticmethod
    def getpid():
        mypid = os.getpid()
        pidFile = open(".mypid", "w")
        pidFile.write(str(mypid))
        pidFile.close()
        os.system('date >> .dates.txt')

    @staticmethod
    def check_kill_process():
        with open(".mypid", "r") as f:
            for line in f:
                pid = line.split()[0].strip()
        os.system('kill %s' % pid)
        quit()

    @staticmethod
    def getGPUids():
        ap = argparse.ArgumentParser()
        ap.add_argument('-e', '--exclude', nargs='*', required=False,
                        help=' use -e to exclude a list of GPUs from being used by mwSuMD: e.g. -e 0 3')
        args = ap.parse_args()
        max_memory = 0.8
        GPUs = GPUtil.getGPUs()
        freeMemory = 0
        gpu_ids = []
        for GPU in GPUs:
            if GPU.memoryUtil > max_memory:
                continue
            if GPU.memoryFree >= freeMemory:
                freeMemory = GPU.memoryFree
                gpu_ids.append(GPU.id)
        for exclusion in args.exclude:
            gpu_ids.remove(int(exclusion))
        if len(gpu_ids) != 0:
            return gpu_ids
        else:
            print("Please leave at least one GPU to run mwSuMD and run again.")
            exit()

    @staticmethod
    def createBatches(walkers, total_gpu_ids):
        quotient, rest = divmod(walkers, len(total_gpu_ids))
        result = quotient * total_gpu_ids + total_gpu_ids[:rest]
        batches = [result[i:i + len(total_gpu_ids)] for i in range(0, len(result), len(total_gpu_ids))]
        return batches, result
