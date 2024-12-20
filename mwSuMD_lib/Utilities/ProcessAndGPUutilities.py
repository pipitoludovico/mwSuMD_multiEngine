import os
import psutil
import GPUtil

from mwSuMD_lib.Utilities.Loggers import Logger


class ProcessManager:

    @staticmethod
    def checkIfInstanceIsRunning():
        if os.path.exists("./.mypid"):
            oldPID = ProcessManager.get_pid_from_file("./.mypid")
            if psutil.pid_exists(oldPID):
                print("Only one instance of mwSuMD can run in the same folder.")
                exit()
            else:
                ProcessManager.getpid()

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
        trajCount = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
        GPUs = GPUtil.getGPUs()
        gpu_ids = []
        for GPU in GPUs:
            gpu_ids.append(GPU.id)
        if len(gpu_ids) != 0:
            Logger.LogToFile('ad', trajCount, f"Available GPUS: {gpu_ids}")
            return gpu_ids
        else:
            print("Please leave at least one GPU to run mwSuMD and run again.")
            exit()

    @staticmethod
    def createBatches(walkers, total_gpu_ids):
        quotient, rest = divmod(walkers, len(total_gpu_ids))
        result = quotient * total_gpu_ids + total_gpu_ids[:rest]
        batches = [result[i:i + len(total_gpu_ids)] for i in range(0, len(result), len(total_gpu_ids))]
        return batches, total_gpu_ids

    @staticmethod
    def get_pid_from_file(pid_file):
        if os.path.exists(pid_file):
            with open(pid_file, 'r') as f:
                pid = f.read().strip()
                return int(pid)
        return None
