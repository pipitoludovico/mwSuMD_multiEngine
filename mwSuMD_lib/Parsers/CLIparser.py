import warnings
import argparse
from os import path
from mwSuMD_lib.Parsers.InputfileParser import mwInputParser
from mwSuMD_lib.Utilities.Loggers import Logger

warnings.filterwarnings(action='ignore')


class ArgParser:
    def __init__(self):
        if not path.exists('./system'):
            print('\nPlease make your ./system folder with the equilibrated system files and outputs')
            exit()

    @staticmethod
    def argumentParser():
        ap = argparse.ArgumentParser()
        ap.add_argument('-e', '--exclude', nargs='*', required=False,
                        help=' use -e to exclude a list of GPUs from being used by mwSuMD: e.g. -e 0 3')
        ap.add_argument('-c', '--command', type=str, nargs='?', required=False,
                        help=' use -c to define a specific command you want to use to run your engine.'
                             'Use "" to define the command: "gmx mdrun -deffnm npt -bonded gpu". '
                             'Be aware of the GPU batch division and let mwSuMD sort the GPUs.')
        ap.add_argument("-k", '--kill', required=False, action='store_true', default=False,
                        help="Stop mwSuMD from the current working directory")
        ap.add_argument("-j", '--join', nargs='*', required=False,
                        help="Merge the trajectories from one step to another: e.g. -j 1 10, -j 10 , or -j all to merge every step.")
        ap.add_argument("-nogpu", '--nogpu', type=int, required=False,
                        help="Add -nogpu and a number (e.g.: -nogpu 4) to not use the gpu and set the batch simulation size")
        ap.add_argument("-openmm", '--openmm', required=False, action='store_true',
                        help="add -openmm to run mwSuMD with OpenMM")
        args = ap.parse_args()

        if args.kill is True:
            import os
            os.system('val=$(<.mypid ) && kill -9 -$val')
            Logger.LogToFile('ad', "", "\nProcess terminated by user.")
            exit()

        if args.join is not None:
            from mwSuMD_lib.MDutils import Merger
            merger = Merger.TrajMerger()
            merger.LoadParameters(args.join)
            merger.Merge()
            exit()

        mwInputParser.initialParameters['COMMAND'] = None
        if args.command is not None:
            mwInputParser.initialParameters['COMMAND'] = args.command
            print(mwInputParser.initialParameters['COMMAND'])

        mwInputParser.initialParameters['EXCLUDED_GPUS'] = None
        if args.exclude is not None and len(args.exclude) != 0:
            mwInputParser.initialParameters['EXCLUDED_GPUS'] = [int(x) for x in args.exclude]

        if args.nogpu:
            mwInputParser.initialParameters['NOGPU'] = [fakeGPU for fakeGPU in range(args.nogpu)]
        else:
            mwInputParser.initialParameters['NOGPU'] = None

        if args.openmm:
            openMM = True
        else:
            openMM = False

        return openMM
