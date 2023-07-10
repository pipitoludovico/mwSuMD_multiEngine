import os
import re
import MDAnalysis as Mda


class TrajMerger:
    def __init__(self):
        self.trajList = []
        self.trajObjects = []
        self.sortedTrajs = None
        self.topology = None

    def loadTrajectories(self):
        for topology in os.listdir('system'):
            if topology.endswith('.psf'):
                self.topology = topology
        for traj in os.listdir('trajectories'):
            if traj.endswith('.xtc'):
                self.trajList.append(traj)

        natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)]
        self.sortedTrajs = (sorted(self.trajList, key=natsort))

    def mergeAll(self):
        Universe = Mda.Universe(f'system/' + self.topology, *self.sortedTrajs)
        atomsel = Universe.select_atoms('all')
        with Mda.Writer('merged_full_movie.xtc') as W:
            for ts in Universe.trajectory:
                W.write(atomsel)

    def mergeFrom(self, start=1, end=-1):
        if start == '' or start is None:
            start = 0
        if end == -1:
            end = 'final_step'
        selection = self.sortedTrajs[start:end]
        Universe = Mda.Universe(f'system/' + self.topology, *selection)
        atomsel = Universe.select_atoms('all')
        with Mda.Writer(f'merged_from_{start}_to_{end}.xtc') as W:
            for ts in Universe.trajectory:
                W.write(atomsel)
