import os

from certifi.__main__ import args

from .mwParser import mwInputParser


class Runner(mwInputParser):
    def __init__(self):
        super(mwInputParser).__init__()

    def runACEMD(self, walkers, trajCount, batch):
        print(os.getcwd())
        for walker in range(walkers):
            print(walker + 1)
            os.chdir('tmp/walker_' + str(walker + 1))
            print(self.par)

            # os.system(f'acemd3 --device{batch[walker-1]} input_{walker}_{trajCount} 1> acemd.log')
            # print(f'acemd3 --device {batch[walker-1]} input_{walker+1}_{trajCount} 1> acemd.log')
            os.chdir(self.folder)
