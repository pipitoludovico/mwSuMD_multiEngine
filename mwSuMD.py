#!/usr/bin/env python3
import os.path

import pandas as pd

if not os.path.exists('./system'):
    print('\nPlease make your ./system folder with the equilibrated system files and outputs')
    exit()

from mwSuMD_lib import GPUoperations
from mwSuMD_lib import Parser
from mwSuMD_lib import SuMD
from mwSuMD_lib import ArgParser

ArgParser.ArgParser()

GPUoperations.ProcessManager()
pars = Parser.mwInputParser()
settings, selection_list, parameterFolderPath = pars.getSettings()
settings_df = pd.DataFrame(sorted(list(settings.items())), columns=['Setting', 'Parameter'])

print(settings_df)
print("If you want to use your personal setting for simulating, please, place it the system folder, call it \n"
      "production.inp/namd/mdp (according to your engine) and mwSuMD will use that instead of the default file. \n"
      "If you choose to do so, make sure it points to a folder named 'restart' to look for the restart binaries.\n")


def main():
    if settings['NumberCV'] is None:
        print('Check the number of metrics in the setting input file')
        quit()
    else:
        sumd = SuMD.suMD1(settings)
        sumd.run_mwSuMD()
    exit()


if __name__ == '__main__':
    main()
