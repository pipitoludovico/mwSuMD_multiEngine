from lib import *

pars = mwParser.mwInputParser()
pars.getSettings()
print(pars.par)

def main():
    if pars.par['NumberCV'] == 1:
        sumd = SuMD_1_CV.suMD1(pars.par)
        sumd.run_SuMD_1_CV()
        # SuMD_1_CV.suMD1().run_SuMD_1_CV()

    # elif par['NumberCV'] == '2':
    #     SuMD_2_CV(par, selection_list, PARPATH)

    # else:
    #     print('Check the number of metrics in the setting input file')
    #     quit()


if __name__ == '__main__':
    main()
