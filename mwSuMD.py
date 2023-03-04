import pandas as pd

from lib import *

Utilities.ProcessManager()
pars = Parser.mwInputParser()
pars.getSettings()
settings = pd.DataFrame(sorted(list(pars.par.items())), columns=['keys', 'values'])
print(settings)


def main():
    if pars.par['NumberCV'] == 1:
        sumd = SuMD_1_CV.suMD1(pars.par)
        sumd.run_SuMD_1_CV()

    elif pars.par['NumberCV'] == 2:
        print("Coming next...")
    #     SuMD_2_CV(par, selection_list, PARPATH)

    else:
        print('Check the number of metrics in the setting input file')
        quit()


if __name__ == '__main__':
    main()
