import pandas as pd

from mwSuMD_lib import *

GPUoperations.ProcessManager()
pars = Parser.mwInputParser()
settings, selection_list, parameterFolderPath = pars.getSettings()
settings_df = pd.DataFrame(sorted(list(settings.items())), columns=['Setting', 'Parameter'])
with open('settings.txt', 'a') as f:
    dfAsString = settings_df.to_string(header=False, index=False)
    f.write(dfAsString)
    f.write("\n")
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
        print("\nRunning final relaxation protocol")
        settings['Relax'] = True
        SimulationChecker.Checker().relaxSystem()
    exit()


if __name__ == '__main__':
    main()
