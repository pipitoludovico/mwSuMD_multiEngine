import os

import mdtraj
import MDAnalysis as Mda

from .Parser import mwInputParser


class Getters(mwInputParser):
    def __init__(self, par):
        super(Getters, self).__init__()
        self.par = par
        self.deNume = None
        self.nume = None
        self.com = None

    def getRMSD(self, sel_1, sel_2):
        import MDAnalysis.analysis.rms

        pdb = f'{self.folder}/system/reference/' + str(self.par['REFERENCE'])
        psf = None
        if self.par['MDEngine'] == 'ACEMD':
            xtc = 'wrapped.xtc'
            if self.par['Forcefield'] == 'CHARMM':
                psf = f'{self.folder}/system/%s' % self.par['PSF']
            elif self.par['Forcefield'] == 'AMBER':
                psf = f'{self.folder}/system/%s.prmtop' % self.par['Topology']
        else:
            xtc = '%s_%s.xtc' % (self.par['Output'], str(self.trajCount))
            psf = '%s.gro' % self.par['Topology']
        u = Mda.Universe(psf, xtc)
        ref = Mda.Universe(pdb)
        R = Mda.analysis.rms.RMSD(u, ref, select="%s" % sel_1, groupselections=["%s" % sel_2])
        R.run()

        rmsd = R.rmsd.T
        data = list(rmsd[3])
        mean_rmsd = sum(data) / len(data)
        last_rmsd = data[-1]
        print(rmsd)
        print(data)
        print(last_rmsd)

        if self.par['NumberCV'] == 1:  # we return the walker score or slope already
            if self.par['Slope'] == 'NO':  # we return the walker score
                distMetric = (mean_rmsd * last_rmsd) ** 0.5
                return distMetric

            elif self.par['Slope'] == 'YES':  # we return the walker slope
                distMetric = self.getSlope(data)
                return distMetric

        elif self.par['NumberCV'] == 2:  # we return all the metric values and the last distance
            return data, last_rmsd

    def getHB_score(self, selection_list=list):
        psf = f"{self.folder}/system/{self.par['PSF']}"
        xtc = 'wrapped.xtc'
        sel_water_oxy = None
        traj = mdtraj.Trajectory(xtc, topology=psf)
        # to get resnam and atom involved in hydrogen bonds

        def label(hbond):
            hbond_label = '%s:%s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
            return hbond_label

        # write HB for each frame
        # another logger?
        def logHB(Frame2hbond_log):
            with open('/hydrogen_bonds.log', 'a') as hbs_log:
                print(Frame2hbond_log, file=hbs_log)

        # debug info about the paqrametrs for HBì_score
        def logScores_detailed(logScoreIndex, logScoreValue, logScoreNumHB_frame, logScoreHBscore):
            import mdtraj
            psf = None
            sel_water_oxy = None

            #logger from Loggers... poi si vede.

            if self.par['Forcefield'] == 'CHARMM':
                psf = '%s.psf' % self.par['Topology']
            elif self.par['Forcefield'] == 'AMBER':
                psf = '%s.prmtop' % self.par['Topology']

            traj = mdtraj.load(xtc, top=psf)

            ligand = traj.topology.select(self.par['ligand_HB'])
            protein = traj.topology.select("protein")
            # figure out how many frames you loaded, this is how many frames we will look for hbonds in
            n_frames = len(traj)
            # print(n_frames)

            """ This set will give us all of the unique hbonds that are made with the ligand, without repeats
            We will want to have this later so we make it not to avoid repeating hbond calculation"""

            all_hbonds_set = set()

            # This list will store all the hbonds made per frame
            hbonds_each_frame = []

            # We want to create a dictionary containing every frame and the ligand hbonds which occur in that frame
            Frame2hbond = {}
            for frame in range(n_frames):
                # The dictionary "words" are the frame number
                Frame2hbond[frame] = []
                # We are doing the hbond analysis frame by frame
                hbonds = mdtraj.baker_hubbard(traj[frame])
                # print(hbonds)
                hbonds_each_frame.append(hbonds)
                # print(hbonds_each_frame)
                # We only care about the hbonds if they involve the ligand
                for hbond_frame in hbonds:
                    if ((hbond_frame[0] in ligand) and (hbond_frame[2] in protein) or (hbond_frame[2] in ligand) and (
                            hbond_frame[0] in protein)):  # ligand is donating or accepting
                        hbond_label_frame = label(hbond_frame)  # get the atom names of the ligand involved
                        # all_hbonds_set.add(tuple(hbond))
                        all_hbonds_set.add(hbond_label_frame)
                        # The dictionary "definitions" are all the hbonds in that frame
                        # Frame2hbond[frame].append(tuple(hbond))
                        Frame2hbond[frame].append(hbond_label_frame)

            print(Frame2hbond)
            logHB(Frame2hbond)
            # let's get just the atom names of the ligand involved without repetitions
            lig_hetatm = []
            for frame, hbs in Frame2hbond.items():
                for hb in hbs:
                    atom = hb.split(':')
                    if 'ZMA' in atom[0]:
                        lig_hetatm.append(atom[0].split('-')[1])
                    elif 'ZMA' in atom[1]:
                        lig_hetatm.append(atom[1].split('-')[1])
            lig_hetatm = list(set(lig_hetatm))
            LIG_HETATM = ' '.join(lig_hetatm)
            print(LIG_HETATM)

            # get number of HB in each frame
            numHB_frame = []
            for frame, hbs in Frame2hbond.items():
                numHB_frame.append(len(hbs))
            print(numHB_frame)

            # let's get contacts between water and ligand's hetatm involved
            from MDAnalysis.analysis import contacts
            u = Mda.Universe(psf, xtc)

            if len(LIG_HETATM) > 0:  # there are hydrogen bonds
                sel_hetatm_lig = "%s and name %s" % (self.par['ligand_HB'], LIG_HETATM)

            else:  # no hydrogen bonds, then use heteroatoms
                sel_hetatm_lig = "%s and (name O* or name N*) " % self.par['ligand_HB']
            if self.par['Forcefield'] == 'CHARMM':
                sel_water_oxy = "resname TIP3 and name O*"
            elif self.par['Forcefield'] == 'AMBER':
                sel_water_oxy = "resname WAT and name O*"
            lig = u.select_atoms(sel_hetatm_lig)
            wat = u.select_atoms(sel_water_oxy)

            timeseries = []
            for ts in u.trajectory:
                # calculate distances between group_a and group_b
                dist = contacts.distance_array(lig.positions, wat.positions)
                # determine which distances <= radius
                n_contacts = contacts.contact_matrix(dist, 5).sum()
                timeseries.append(n_contacts)

            print(timeseries)

            # part re the score
            score_per_frame = []
            for index, value in enumerate(timeseries):
                print("DEBUG:", index, value, numHB_frame[index])
                if numHB_frame[index] > 0:  # there are HB fot this frame
                    HBscore = value ** (1 / numHB_frame[index])
                else:  # no HB fot this frame
                    HBscore = 10
                # some logging for debugging etc
                score_per_frame.append(HBscore)
                logScores_detailed(index, value, numHB_frame[index], HBscore)
            return score_per_frame, score_per_frame[-1]

    def contacts_misc(self, sel_1, sel_2):
        import MDAnalysis
        from MDAnalysis.analysis import contacts

        # Debug: log with crude values for score computation
        psf = None
        xtc = 'wrapped.xtc'

        if self.par['MDEngine'] == 'ACEMD':
            if self.par['Forcefield'] == 'CHARMM':
                psf = f'{self.folder}/system/%s' % self.par['PSF']
            elif self.par['Forcefield'] == 'AMBER':
                psf = f'{self.folder}/system/%s' % self.par['prmtop']

        elif self.par['MDEngine'] == 'GROMACS':
            psf = '%s.gro' % self.par['gro']

        u = MDAnalysis.Universe(psf, xtc)

        sel_1 = u.select_atoms(sel_1)
        sel_2 = u.select_atoms(sel_2)

        # elif self.par['NumberCV'] == '2' and self.par['Metric_1'] == 'Contacts':
        # sel_1 = u.select_atoms(selection_list[0])
        # sel_2 = u.select_atoms(selection_list[1])

        # elif self.par['NumberCV'] == '2' and self.par['Metric_2'] == 'Contacts':
        # sel_1 = u.select_atoms(selection_list[2])
        # sel_2 = u.select_atoms(selection_list[3])

        timeseries = []
        for ts in u.trajectory:
            if ts:
                # calculate distances between sel_1 and sel_2
                dist = contacts.distance_array(sel_1.positions, sel_2.positions)
                # determine which distances <= radius
                n_contacts = contacts.contact_matrix(dist, 3.5).sum()
                # timeseries.append([ts.frame, n_contacts])
                timeseries.append(n_contacts)

        n = len(timeseries)
        mean_contacts = sum(timeseries) / n
        last_contacts = timeseries[-1]

        if self.par['NumberCV'] == 1:  # we return the walker score or slope already
            if self.par['Slope'] == 'NO':  # we return the walker score
                distMetric = (mean_contacts * last_contacts) ** 0.5
                # logContact_misc(timeseries, mean_contacts, last_contacts, distMetric)
                return distMetric

            elif self.par['Slope'] == 'YES':  # we return the walker slope
                distMetric = self.getSlope(timeseries)
                return distMetric

        elif self.par['NumberCV'] == '2':  # we return all the metric values and the last distance
            return timeseries, last_contacts

    def getSlope(self, values_metric):
        import numpy as np
        """Compute the least square methods on the data
        list provided called by other metrics functions"""
        data = dict(enumerate(values_metric))
        # why are you calculating maxes if we don't use them?
        # togliere time e prendere solo maxDist da loggare
        maxDist, time_max = max((dist, value) for dist, value in data.items())
        meanTime = np.array(list(data.values())).mean()
        meanDist = np.array(list(data.keys())).mean()
        self.nume = [(float(value) - meanTime) * (float(key) - meanDist) for key, value in data.items()]
        self.deNume = [(float(value) - meanTime) ** 2 for value in data.values()]
        try:
            slope = float(np.sum(self.nume)) / float(np.sum(self.deNume))
            return slope
        except:
            print("Slope deNumerator was 0.")
            slope = 0
        return slope

    def compute_center_of_mass(self, select=None):
        import MDAnalysis as Mda
        import numpy as np
        if self.par['MDEngine'] == 'ACEMD':
            psf = f"../../system/{self.par['PSF']}" \
                if self.par['Forcefield'] == 'CHARMM' else f"../../system/{self.par['prmtop']}"
        else:
            psf = f"{self.par['gro']}"

        psf = None
        xtc = "wrapped.xtc"
        u = Mda.Universe(psf, xtc)
        sele = u.select_atoms(select)
        arr = np.empty((sele.n_residues, u.trajectory.n_frames, 3))

        for ts in Mda.log.ProgressBar(u.trajectory):
            arr[:, ts.frame] = sele.center_of_mass()
        return arr

    def getDistance(self, sel_1, sel_2):
        import numpy as np
        # Construct path to trajectory and topology files

        c1 = self.compute_center_of_mass(select=sel_1)
        c2 = self.compute_center_of_mass(select=sel_2)

        # Compute distances
        distances = [np.linalg.norm(a - b) * 10 for a, b in zip(c1, c2)]
        mean_distance = np.mean(distances)
        last_distance = distances[-1]

        # Return computed metric(s)
        if self.par['NumberCV'] == 1:
            if self.par['Slope'] == 'NO':
                return (mean_distance * last_distance) ** 0.5
            else:
                return self.getSlope(distances)
        else:
            return distances, last_distance
