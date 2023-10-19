import os
import re
import MDAnalysis as Mda


class TrajMerger:
    def __init__(self):
        self.trajList = []
        self.trajObjects = []
        self.sortedTrajs = None
        self.topology = None
        self.topologyExtensions = ('.psf', '.prmtop', '.tpr')

    def loadTrajectories(self):
        for topology in os.listdir('system'):
            if topology.endswith(self.topologyExtensions):
                self.topology = topology
        for traj in os.listdir('trajectories'):
            if traj.endswith('.xtc'):
                self.trajList.append('trajectories/' + traj)

        natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)]
        self.sortedTrajs = (sorted(self.trajList, key=natsort))

    def mergeAll(self):
        Universe = Mda.Universe(f'system/' + self.topology, *self.sortedTrajs, refresh_offset=True)
        atomsel = Universe.select_atoms('all')
        with Mda.Writer('merged_full_movie.xtc') as W:
            for ts in Universe.trajectory:
                W.write(atomsel)

    def mergeFromToEnd(self, start=0):
        selection = self.sortedTrajs[int(start):]
        Universe = Mda.Universe(f'system/' + self.topology, *selection, refresh_offset=True)
        atomsel = Universe.select_atoms('all')
        with Mda.Writer(f'merged_from_{start}_to_lastStep.xtc') as W:
            for ts in Universe.trajectory:
                W.write(atomsel)

    def mergeFromTo(self, start=0, end=-1):
        selection = self.sortedTrajs[int(start):int(end)]
        print(self.sortedTrajs[int(start)])
        print(self.sortedTrajs[int(end)])

        Universe = Mda.Universe(f'system/' + self.topology, *selection, refresh_offset=True)
        atomsel = Universe.select_atoms('all')
        with Mda.Writer(f'merged_from_{start}_to_{end}.xtc') as W:
            for ts in Universe.trajectory:
                W.write(atomsel)
