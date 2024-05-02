#!/usr/bin/env python3
import os.path

if not os.path.exists('./system'):
    print('\nPlease make your ./system folder with the equilibrated system files and outputs')
    exit()

from mwSuMD_lib import GPUoperations
from mwSuMD_lib import ArgParser
from mwSuMD_lib import Loggers

parser = ArgParser.ArgParser()
openMM = parser.argumentParser()
GPUoperations.ProcessManager()
Loggers.Logger.LogToFile('ad', '',
                         "If you want to use your personal setting for simulating, please, place it the system folder, call it \n "
                         "production.inp/namd/mdp (according to your engine) and mwSuMD will use that instead of the default file. \n"
                         "If you choose to do so, make sure it points to a folder named 'restart' to look for the restart binaries.\n")


def main():
    from mwSuMD_lib import SuMD
    print("\nRunning mwSuMD")
    multiSumd = SuMD.suMD1(openMM)
    multiSumd.run_mwSuMD()
    exit()


if __name__ == '__main__':
    main()
