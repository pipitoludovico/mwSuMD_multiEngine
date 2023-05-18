import os

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
        GPUs = GPUtil.getGPUs()
        gpu_ids = []
        for GPU in GPUs:
            gpu_ids.append(GPU.id)
        if len(gpu_ids) != 0:
            print("Available GPUS: ", gpu_ids)
            return gpu_ids
        else:
            return [0]
            # print("Please leave at least one GPU to run mwSuMD and run again.")
            # exit()

    @staticmethod
    def createBatches(walkers, total_gpu_ids):
        quotient, rest = divmod(walkers, len(total_gpu_ids))
        result = quotient * total_gpu_ids + total_gpu_ids[:rest]
        batches = [result[i:i + len(total_gpu_ids)] for i in range(0, len(result), len(total_gpu_ids))]
        return batches, total_gpu_ids
