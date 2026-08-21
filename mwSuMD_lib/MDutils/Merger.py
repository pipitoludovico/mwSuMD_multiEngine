import os
import re
from warnings import filterwarnings

import numpy as np
from MDAnalysis import Universe
from MDAnalysis.coordinates.XTC import XTCReader, XTCWriter

filterwarnings(action='ignore')


class TrajMerger:
    def __init__(self):
        self.trajList = []
        self.sortedTrajs = None
        self.topology = None
        self.coordinates = None
        # self.topologyExtensions = ('.psf', '.prmtop', '.tpr')
        self.coordinatesExtensions = ('.pdb', '.inpcrd', '.gro')
        self.inputFile = 'simulation_settings_mwSuMD.inp'
        self.outputFileName = 'merged.xtc'
        self.filterSelection = 'all'
        self.ReadWrap()

    def LoadParameters(self, interval: list):
        for traj in os.listdir('trajectories'):
            if traj.endswith('.xtc'):
                self.trajList.append('trajectories/' + traj)
        if not self.trajList:
            print("No .xtc found in ./trajectories: there is nothing to merge.")
            exit(1)
        self.topology = self.FindTopology()
        natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\\d+)', s)]
        self.sortedTrajs = (sorted(self.trajList, key=natsort))
        try:
            if len(interval) == 1 and interval[0] != 'all':
                start: int = int(interval[0]) - 1
                self.sortedTrajs = sorted(self.trajList, key=natsort)[start:]
                print(f"Merging from {self.sortedTrajs[0]} to the end.")
                self.outputFileName = f'merged_from_{interval[0]}_to_last_step.xtc'
            elif len(interval) == 1 and interval[0] == 'all':
                print("Merging all steps into one.")
                self.sortedTrajs = (sorted(self.trajList, key=natsort))
                self.outputFileName = f'merged_full_movie.xtc'
            elif len(interval) == 2:
                startIndex = int(interval[0]) - 1
                endIndex = int(interval[1])
                self.sortedTrajs = sorted(self.trajList, key=natsort)[startIndex:endIndex]
                print(f"Merging from {self.sortedTrajs[0]} to {self.sortedTrajs[-1]}.")
                self.outputFileName = f'merged_from_{interval[0]}_to_{interval[1]}.xtc'
        except Exception as e:
            print(
                'The index of the step returned an error. Make sure the index you used is inside the trajectory folder',
                e)
            exit()

    def FindTopology(self):
        """Pick the topology by atom count rather than by filesystem order. The stored
        trajectories only contain the FilterOut atoms, so the full system topology does not
        match them. This used to keep the last match of os.listdir, which made the choice
        depend on the order the filesystem happened to return."""
        candidates = sorted(found for found in os.listdir("./") if found.endswith(self.coordinatesExtensions))
        if not candidates:
            print("No topology found in the current directory. Expected 'filtered_topology.pdb', "
                  f"or another file ending in {self.coordinatesExtensions}, produced by the wrapping step.")
            exit(1)
        # the file written by the wrapping step is the one meant for this, so try it first
        preferred = [found for found in candidates if found.startswith('filtered_topology')]
        ordered = preferred + [found for found in candidates if found not in preferred]

        expectedAtoms = None
        for trajectory in sorted(self.trajList):
            try:
                expectedAtoms = XTCReader(trajectory).n_atoms
                break
            except Exception:
                # truncated or unreadable file: try the next one before giving up
                print(f"Warning: could not read '{trajectory}' to count its atoms.")
        if expectedAtoms is None:
            print(f"None of the {len(self.trajList)} trajectories in ./trajectories could be read. "
                  "They are empty or truncated: check the runs that produced them.")
            exit(1)
        for candidate in ordered:
            try:
                if Universe(candidate).atoms.n_atoms == expectedAtoms:
                    return candidate
            except Exception:
                continue
        print(f"None of the topologies in the current directory {candidates} has the "
              f"{expectedAtoms} atoms stored in the trajectories. The wrapping step writes a "
              "'filtered_topology.pdb' matching your FilterOut selection: copy it here.")
        exit(1)

    def Merge(self):
        # Load just the topology first
        print(f"File used: {self.topology}")
        u = Universe(self.topology)
        ag = u.select_atoms(self.filterSelection)
        # Reduced coordinates of the previous unwrapped frame. This is kept ACROSS files:
        # MDAnalysis' NoJump resets its own state whenever ts.frame == 0, which is every
        # time a new xtc is loaded, so the jump between the last frame of one step and the
        # first frame of the next was never removed. Steps are contiguous simulations, so
        # that boundary is exactly where a molecule can appear to cross the box.
        previousReduced = None

        with XTCWriter(self.outputFileName, n_atoms=len(ag)) as writer:
            for traj_file in self.sortedTrajs:
                print(f"Processing: {traj_file}")
                try:
                    # Load each trajectory individually to bypass ChainReader offset issues
                    u.load_new(traj_file)
                    for ts in u.trajectory:
                        previousReduced = self.unwrapFrame(ag, ts, previousReduced)
                        writer.write(ag)
                except Exception as e:
                    # Do not carry on: skipping every file in turn used to leave an empty
                    # merged trajectory behind while still exiting successfully.
                    print(f"CRITICAL ERROR in file {traj_file}: {e}")
                    print(f"Topology '{self.topology}' gives {len(ag)} atoms with selection "
                          f"'{self.filterSelection}'. The stored trajectories only hold the FilterOut "
                          "atoms, so a full system topology will not match them.")
                    print(f"'{self.outputFileName}' is incomplete and must not be used.")
                    exit(1)

    @staticmethod
    def unwrapFrame(ag, ts, previousReduced):
        """Remove periodic jumps relative to the previous frame and return the new reduced
        coordinates. Same algorithm as MDAnalysis' NoJump transformation (eq. B6 of
        10.1021/acs.jctc.2c00327), but the reference frame is carried over between files.
        Assumes atoms move by less than half a box length between saved frames."""
        box = ts.triclinic_dimensions
        if box is None:
            raise ValueError(f"No periodic box dimensions at frame {ts.frame}")
        boxInverse = np.linalg.inv(box)
        currentReduced = ag.positions @ boxInverse
        if previousReduced is None or previousReduced.shape != currentReduced.shape:
            unwrappedReduced = currentReduced
        else:
            unwrappedReduced = currentReduced - np.round(currentReduced - previousReduced)
        ag.positions = unwrappedReduced @ box
        return unwrappedReduced

    def ReadWrap(self) -> None:
        with open(self.inputFile) as infile:
            for line in infile.readlines():
                if line.startswith('#'):
                    continue
                if line.startswith('FilterOut'):
                    value = line.split('=')[1].strip()
                    if value:
                        self.filterSelection = value
                        print("USING THIS SELECTION FOR THE TRAJECTORY:", self.filterSelection)
