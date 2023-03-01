import os

from lib import *


def main():
    print("HI")
    parser = mwParser.mwInputParser()
    parser.getSettings()
    print(parser.par)
    mwutils.ProcessManager()
    print("TRJ:")
    print(parser.trajCount)


    if parser.par['NumberCV'] == 1:
        print("")
        SuMD_1_CV.suMD1().run_SuMD_1_CV()

        # print("Creating sorting result's folders...")
        # result_name = parser.par["Output"]
        # os.makedirs(result_name, exist_ok=True)
        # os.system(f'mv restarts trajectories nohup.out settings.txt {result_name}')

    # elif par['NumberCV'] == '2':
    #     SuMD_2_CV(par, selection_list, PARPATH)

    else:
        print('Check the number of metrics in the setting input file')
        quit()


if __name__ == '__main__':
    main()
