import os
import psutil
import GPUtil
from time import sleep
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
                ProcessManager.writepid()
        else:
            ProcessManager.writepid()

    @staticmethod
    def writepid():
        mypid = os.getpid()
        with open(".mypid", "w") as pidFile:
            pidFile.write(str(mypid))

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

    @staticmethod
    def CheckIfAnyGPUisFree(GPUs: list):
        for device in GPUtil.getGPUs():
            if device.id in GPUs:
                usage = round(device.memoryUsed / device.memoryTotal)
                if usage < 0.7:
                    check = False
                    return check, device.id
        sleep(60)
