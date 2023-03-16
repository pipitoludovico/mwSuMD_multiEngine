import os

import pandas as pd

from lib import *

Utilities.ProcessManager()
pars = Parser.mwInputParser()
pars.getSettings()
settings = pd.DataFrame(sorted(list(pars.par.items())), columns=['keys', 'values'])
print(settings)


def main():
    if pars.par['NumberCV'] is None:
        print('Check the number of metrics in the setting input file')
        quit()
    else:
        sumd = SuMDProtocol.suMD1(pars.par)
        sumd.run_mwSuMD()


if __name__ == '__main__':
    main()
