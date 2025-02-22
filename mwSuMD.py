#!/usr/bin/env python3
import warnings
from mwSuMD_lib.Utilities import ProcessAndGPUutilities, Loggers
from mwSuMD_lib.Parsers import CLIparser
warnings.filterwarnings('ignore')


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
