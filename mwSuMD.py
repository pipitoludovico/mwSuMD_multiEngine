#!/usr/bin/env python3
import warnings
from mwSuMD_lib.Utilities import ProcessAndGPUutilities, Loggers
from mwSuMD_lib.Parsers import CLIparser

warnings.filterwarnings('ignore')

Loggers.Logger.LogToFile('ad', '',
                         "If you want to use your personal setting for simulating, please, place it the system folder, call it \n "
                         "production.inp/namd/mdp (according to your engine) and mwSuMD will use that instead of the default file. \n"
                         "If you choose to do so, make sure it points to a folder named 'restart' to look for the restart binaries.\n")


def main():
    parser = CLIparser.ArgParser()
    openMM = parser.argumentParser()
    ProcessAndGPUutilities.ProcessManager.checkIfInstanceIsRunning()

    from mwSuMD_lib.Protocol import SuMD
    multiSumd = SuMD.suMD1(openMM)
    multiSumd.run_mwSuMD()
    exit()


if __name__ == '__main__':
    main()
