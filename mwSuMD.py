#!/usr/bin/env python3
import warnings
from mwSuMD_lib.Utilities import ProcessAndGPUutilities, Loggers
from mwSuMD_lib.Parsers import CLIparser

warnings.filterwarnings('ignore')

Loggers.Logger.LogToFile('ad', '',
                         "If you want to use a custom input file for your engine of choice, please, place it in the './system' folder, and call it \n "
                         "production.inp/namd/mdp (according to your engine). MwSuMD will use that instead of the default file. \n"
                         "If you choose to do so, make sure it points to a folder named 'restart' for continuing from the last binaries!\n")


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
