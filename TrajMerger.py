import os
import re

import mdtraj as mdt


class TrajMerger:
    @staticmethod
    def mergeTrajectories():
        trajList = []
        pdbFile = None
        trajObjects = []

        for pdb in os.listdir('system'):
            if pdb.endswith('.pdb'):
                pdbFile = pdb

        for traj in os.listdir('trajectories'):
            trajList.append(traj)

        natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)]
        sortedTrajs = (sorted(trajList, key=natsort))

        for idx, trajectory in enumerate(sortedTrajs):
            print(idx, trajectory)
            print(f'trajectories/' + trajectory)
            _ = mdt.load(f'trajectories/' + trajectory, top= f'system/' + pdbFile)
            trajObjects.append(_)
        merged = mdt.join(trajObjects,check_topology=True, discard_overlapping_frames=True)
        merged.save_xtc('merged.xtc')


TrajMerger.mergeTrajectories()
