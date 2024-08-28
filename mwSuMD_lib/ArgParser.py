import warnings
import argparse
from .Parser import mwInputParser
from .Loggers import Logger

warnings.filterwarnings(action='ignore')


class ArgParser:
    def __init__(self):
        pass

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
            Logger.LogToFile('ad', "", "\nProcess terminated by user.")
            os.system('val=$(<.mypid ) && kill -15 -$val')

        if args.join is not None:
            from mwSuMD_lib import Merger
            merger = Merger.TrajMerger()
            merger.loadTrajectories()
            if len(args.join) == 1 and 'all' in args.join[0]:
                merger.mergeAll()
            elif len(args.join) == 1 and 'all' not in args.join:
                merger.mergeFromToEnd(args.join[0])
            elif len(args.join) == 2 and 'all' not in args.join:
                merger.mergeFromTo(args.join[0], args.join[1])
            else:
                print(
                    "Error: incorrect arguments for -j. -join needs 2 number to set the starting and ending steps to be merged")
                exit()
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
