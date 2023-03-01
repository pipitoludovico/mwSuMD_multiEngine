import os


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
