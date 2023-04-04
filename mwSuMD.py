import pandas as pd

from lib import *

Utilities.ProcessManager()
pars = Parser.mwInputParser()
settings, selection_list, parameterFolderPath = pars.getSettings()
settings_df = pd.DataFrame(sorted(list(settings.items())), columns=['Setting', 'Parameter'])
print(settings_df)

print("If you want to use your personal setting for simulating, please, place it the system folder and call it \n"
      "production.inp/namd/mdp (according to your engine) and mwSuMD will use that instead of the default file.\n")


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
