import os
import re

from warnings import filterwarnings

filterwarnings(action='ignore')


class TrajMerger:
    def __init__(self):
        self.trajList = []
        self.sortedTrajs = None
        self.topology = None
        self.coordinates = None
        self.topologyExtensions = ('.psf', '.prmtop', '.tpr')
        self.coordinatesExtensions = ('.pdb', '.innpcr', '.gro')
        self.inputFile = 'simulation_settings_mwSuMD.inp'
        self.outputFileName = 'merged.xtc'
        self.filterSelection = ''
        self.ReadWrap()

    def LoadParameters(self, interval: list):
        for file in os.listdir("system"):
            if file.endswith(self.topologyExtensions):
                self.topology = file
            if file.endswith(self.coordinatesExtensions):
                self.coordinates = file
        for traj in os.listdir('trajectories'):
            if traj.endswith('.xtc'):
                self.trajList.append('trajectories/' + traj)
        natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)]
        self.sortedTrajs = (sorted(self.trajList, key=natsort))

        if len(interval) == 1 and interval[0] != 'all':
            start: int = int(interval[0]) - 1
            self.sortedTrajs = (sorted(self.trajList[start:], key=natsort))
            print(f"Merging from {self.sortedTrajs[0]} to the end.")
            self.outputFileName = f'merged_from_{interval[0]}_to_last_step.xtc'
        elif len(interval) == 1 and interval[0] == 'all':
            print("Merging all steps into one.")
            self.sortedTrajs = (sorted(self.trajList, key=natsort))
            self.outputFileName = f'merged_full_movie.xtc'
        elif len(interval) == 2:
            startIndex = int(interval[0]) - 1
            endIndex = int(interval[1])
            self.sortedTrajs = (sorted(self.trajList[startIndex:endIndex], key=natsort))
            print(f"Merging from {self.sortedTrajs[0]} to {self.sortedTrajs[-1]}.")
            self.outputFileName = f'merged_from_{interval[0]}_to_{interval[1]}.xtc'

    def Merge(self):
        os.system(f"mdconvert {' '.join(self.sortedTrajs)} -f -o {self.outputFileName}")

    def ReadWrap(self) -> None:
        with open(self.inputFile) as infile:
            for line in infile.readlines():
                if line.startswith('#'):
                    continue
                if line.startswith('Wrap'):
                    if line.split('=')[1].strip() != '':
                        self.filterSelection = line.split('=')[1].strip()
