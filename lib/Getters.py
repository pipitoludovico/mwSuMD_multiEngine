class Getters:
    def __init__(self):
        self.denume = None
        self.nume = None
        self.com = None

    def getRMSD(self, par, sel_1, sel_2, trajCount):
        import MDAnalysis as Mda
        import MDAnalysis.analysis.rms
        print()
        pdb = f'../../system/reference/' + str(par['REFERENCE'])
        psf = None
        if par['MDEngine'] == 'ACEMD':
            xtc = '%s_%s_wrapped.xtc' % (par['Output'], str(trajCount))
            if par['Forcefield'] == 'CHARMM':
                psf = f'../../system/%s' % par['PSF']
            elif par['Forcefield'] == 'AMBER':
                psf = f'../../system/%s.prmtop' % par['Topology']
        else:
            xtc = '%s_%s.xtc' % (par['Output'], str(trajCount))
            psf = '%s.gro' % par['Topology']

        u = Mda.Universe(psf, xtc)
        ref = Mda.Universe(pdb)
        R = Mda.analysis.rms.RMSD(u, ref, select="%s" % sel_1, groupselections=["%s" % sel_2])
        R.run()
        rmsd = R.rmsd.T

        data = list(rmsd[3])
        mean_rmsd = sum(data) / len(data)
        last_rmsd = rmsd[-1]
        # minimumRmsd = min(data)
        # bestFrame = data.index(minimumRmsd)

        # try:
        #     maxRMSD = max(x for x in data if x < 100)
        #     maxFrame = data.index(maxRMSD)
        # except:
        #     maxRMSD = 0
        #     maxFrame = 0

        if par['NumberCV'] == 1:
            if par['Slope'] == 'NO':  # we return the walker score
                distMetric = (mean_rmsd * last_rmsd) ** 0.5
                return distMetric
                # distMetric = minimumRmsd
                # logRMSD(data, mean_rmsd, last_rmsd, distMetric)  # logger da sistemare a parte
                # print("\nFrame " + str(bestFrame) + " had the lowest RMSD: " + str(distMetric))
                # return bestFrame, distMetric, ""
            if par['Slope'] == 'YES':
                return self.getSlope(data)
        if par['NumberCV'] == 2:
            return data, last_rmsd
                    # distMetric = minimumRmsd
                    # slope = self.getSlope(data)
                    # return bestFrame, distMetric, slope
            # if par['Transition_1'] == 'positive':
            #     if par['Slope'] == 'NO':  # we return the walker score
                    # distMetric = maxRMSD
                    # print("\nFrame " + str(maxFrame) + " had the highest RMSD: " + str(distMetric))
                    # return maxFrame, distMetric, ""
                    # distMetric = (mean_rmsd * last_rmsd) ** 0.5
                    # return distMetric
                # elif par['Slope'] == 'YES':
                    # distMetric = maxRMSD
                    # slope = self.getSlope(data)
                    # return maxFrame, distMetric, slope

    # @staticmethod
    # def getHB_score(par, n, selection_list=list):
    #     print(selection_list)
    #     psf = None
    #     xtc = '%s_%s_wrapped.xtc' % (par['Output'], str(n))
    #     sel_water_oxy = None
    #
    #     # to get resnam and atom involved in hydrogen bonds
    #     def label(hbond):
    #         hbond_label = '%s:%s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
    #         return hbond_label
    #
    #     # write HB for each frame
    #     def logHB(Frame2hbond_log):
    #         with open('/hydrogen_bonds.log', 'a') as hbs_log:
    #             print(Frame2hbond_log, file=hbs_log)
    #
    #     # debug info about the paqrametrs for HBì_score
    #     def logScores_detailed(logScoreIndex, logScoreValue, logScoreNumHB_frame, logScoreHBscore):
    #         psf = None
    #         sel_water_oxy = None
    #         with open('/HB_scores_crude.log', 'a') as scoreF:
    #             scoreF.write(
    #                 'number wat mols: %s, number HB: %s, HB_score: %s\n' % (str(logScoreValue),
    #                                                                         str(logScoreNumHB_frame),
    #                                                                         str(logScoreHBscore)))
    #
    #         if par['Forcefield'] == 'CHARMM':
    #             psf = '%s.psf' % par['Topology']
    #         elif par['Forcefield'] == 'AMBER':
    #             psf = '%s.prmtop' % par['Topology']
    #
    #         traj = mdtraj.load(xtc, top=psf)
    #
    #         ligand = traj.topology.select(par['ligand_HB'])
    #         protein = traj.topology.select("protein")
    #         # figure out how many frames you loaded, this is how many frames we will look for hbonds in
    #         n_frames = len(traj)
    #         # print(n_frames)
    #
    #         """ This set will give us all of the unique hbonds that are made with the ligand, without repeats
    #         We will want to have this later so we make it not to avoid repeating hbond calculation"""
    #
    #         all_hbonds_set = set()
    #
    #         # This list will store all the hbonds made per frame
    #         hbonds_each_frame = []
    #
    #         # We want to create a dictionary containing every frame and the ligand hbonds which occur in that frame
    #         Frame2hbond = {}
    #         for frame in range(n_frames):
    #             # The dictionary "words" are the frame number
    #             Frame2hbond[frame] = []
    #             # We are doing the hbond analysis frame by frame
    #             hbonds = mdtraj.baker_hubbard(traj[frame])
    #             # print(hbonds)
    #             hbonds_each_frame.append(hbonds)
    #             # print(hbonds_each_frame)
    #             # We only care about the hbonds if they involve the ligand
    #             for hbond_frame in hbonds:
    #                 if ((hbond_frame[0] in ligand) and (hbond_frame[2] in protein) or (hbond_frame[2] in ligand) and (
    #                         hbond_frame[0] in protein)):  # ligand is donating or accepting
    #                     hbond_label_frame = label(hbond_frame)  # get the atom names of the ligand involved
    #                     # all_hbonds_set.add(tuple(hbond))
    #                     all_hbonds_set.add(hbond_label_frame)
    #                     # The dictionary "definitions" are all the hbonds in that frame
    #                     # Frame2hbond[frame].append(tuple(hbond))
    #                     Frame2hbond[frame].append(hbond_label_frame)
    #
    #         print(Frame2hbond)
    #         logHB(Frame2hbond)
    #         # let's get just the atom names of the ligand involved without repetitions
    #         lig_hetatm = []
    #         for frame, hbs in Frame2hbond.items():
    #             for hb in hbs:
    #                 atom = hb.split(':')
    #                 if 'ZMA' in atom[0]:
    #                     lig_hetatm.append(atom[0].split('-')[1])
    #                 elif 'ZMA' in atom[1]:
    #                     lig_hetatm.append(atom[1].split('-')[1])
    #         lig_hetatm = list(set(lig_hetatm))
    #         LIG_HETATM = ' '.join(lig_hetatm)
    #         print(LIG_HETATM)
    #
    #         # get number of HB in each frame
    #         numHB_frame = []
    #         for frame, hbs in Frame2hbond.items():
    #             numHB_frame.append(len(hbs))
    #         print(numHB_frame)
    #
    #         # let's get contacts between water and ligand's hetatm involved
    #         from MDAnalysis.analysis import contacts
    #         u = Mda.Universe(psf, xtc)
    #
    #         if len(LIG_HETATM) > 0:  # there are hydrogen bonds
    #             sel_hetatm_lig = "%s and name %s" % (par['ligand_HB'], LIG_HETATM)
    #
    #         else:  # no hydrogen bonds, then use heteroatoms
    #             sel_hetatm_lig = "%s and (name O* or name N*) " % par['ligand_HB']
    #         if par['Forcefield'] == 'CHARMM':
    #             sel_water_oxy = "resname TIP3 and name O*"
    #         elif par['Forcefield'] == 'AMBER':
    #             sel_water_oxy = "resname WAT and name O*"
    #         lig = u.select_atoms(sel_hetatm_lig)
    #         wat = u.select_atoms(sel_water_oxy)
    #
    #         timeseries = []
    #         for ts in u.trajectory:
    #             # calculate distances between group_a and group_b
    #             dist = contacts.distance_array(lig.positions, wat.positions)
    #             # determine which distances <= radius
    #             n_contacts = contacts.contact_matrix(dist, 5).sum()
    #             timeseries.append(n_contacts)
    #
    #         print(timeseries)
    #
    #         # part re the score
    #         score_per_frame = []
    #         for index, value in enumerate(timeseries):
    #             print("DEBUG:", index, value, numHB_frame[index])
    #             if numHB_frame[index] > 0:  # there are HB fot this frame
    #                 HBscore = value ** (1 / numHB_frame[index])
    #             else:  # no HB fot this frame
    #                 HBscore = 10
    #             # some logging for debugging etc
    #             score_per_frame.append(HBscore)
    #             logScores_detailed(index, value, numHB_frame[index], HBscore)
    #         return score_per_frame, score_per_frame[-1]

    def contacts_misc(self, par, sel_1, sel_2, cycle_number):
        import MDAnalysis
        from MDAnalysis.analysis import contacts

        # Debug: log with crude values for score computation
        psf = None
        xtc = None
        if par['MDEngine'] == 'ACEMD':
            xtc = '%s_%s_wrapped.xtc' % (par['Output'], str(cycle_number))

            if par['Forcefield'] == 'CHARMM':
                psf = '%s.psf' % par['Topology']
            elif par['Forcefield'] == 'AMBER':
                psf = '%s.prmtop' % par['Topology']

        elif par['MDEngine'] == 'GROMACS':
            xtc = '%s_%s.xtc' % (par['Output'], str(cycle_number))
            psf = '%s.gro' % par['Topology']

        u = MDAnalysis.Universe(psf, xtc)

        # let's decide what selections to use according to the number of metrics we want to compute
        # if par['NumberCV'] == '1' and par['Metric_1'] == 'Contacts':
        sel_1 = u.select_atoms(sel_1)
        sel_2 = u.select_atoms(sel_2)

        # elif par['NumberCV'] == '2' and par['Metric_1'] == 'Contacts':
        # sel_1 = u.select_atoms(selection_list[0])
        # sel_2 = u.select_atoms(selection_list[1])

        # elif par['NumberCV'] == '2' and par['Metric_2'] == 'Contacts':
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

        if par['NumberCV'] == '1':  # we return the walker score or slope already
            if par['Slope'] == 'NO':  # we return the walker score
                distMetric = (mean_contacts * last_contacts) ** 0.5
                # logContact_misc(timeseries, mean_contacts, last_contacts, distMetric)
                return distMetric

            elif par['Slope'] == 'YES':  # we return the walker slope
                distMetric = self.getSlope(timeseries)
                return distMetric

        elif par['NumberCV'] == '2':  # we return all the metric values and the last distance
            return timeseries, last_contacts

    def getSlope(self, values_metric):
        import numpy as np
        """Compute the least square methods on the data
        list provided called by other metrics functions"""
        data = {k: v for v, k in enumerate(values_metric)}
        # why are you calculating maxes if we don't use them?
        # maxDist, time_max = max((dist, value) for dist, value in data.items())
        meanTime = np.array(list(data.values())).mean()
        meanDist = np.array(list(data.keys())).mean()
        self.nume = [(float(value) - meanTime) * (float(key) - meanDist) for key, value in data.items()]
        self.denume = [(float(value) - meanTime) ** 2 for value in data.values()]
        slope = float(np.sum(self.nume)) / float(np.sum(self.denume))
        return slope

    def compute_center_of_mass_coincise(self, traj, select=None):
        import numpy as np
        # Get masses of atoms based on the selection
        if select is None:
            masses = np.array([a.element.mass for a in traj.top.atoms])
            xyz = traj.xyz
        else:
            atoms_of_interest = traj.topology.select(select)
            masses = np.array([traj.top.atom(i).element.mass for i in atoms_of_interest])
            xyz = traj.xyz[:, atoms_of_interest]

        # Compute center of mass for each frame
        self.com = np.dot(xyz.astype('float64'), masses / masses.sum(axis=0))
        return self.com

    def getDistance_coincise(self, par, sel_1, sel_2, trajCount):
        import mdtraj
        import numpy as np

        # Construct path to trajectory and topology files
        if par['MDEngine'] == 'ACEMD':
            xtc = f"{par['Output']}_{trajCount}_wrapped.xtc"
            topology = f"../../system/{par['Topology']}.psf" \
                if par['Forcefield'] == 'CHARMM' else f"../../system/{par['Topology']}.prmtop"
        else:
            xtc = f"{par['Output']}_{trajCount}.xtc"
            topology = f"{par['Topology']}.gro"

        # Load trajectory and compute center of mass
        traj = mdtraj.load(xtc, top=topology)
        c1 = self.compute_center_of_mass(traj, select=sel_1)
        c2 = self.compute_center_of_mass(traj, select=sel_2)

        # Compute distances
        distances = [np.linalg.norm(a - b) * 10 for a, b in zip(c1, c2)]
        mean_distance = np.mean(distances)
        last_distance = distances[-1]

        # Return computed metric(s)
        if par['NumberCV'] == '1':
            if par['Slope'] == 'NO':
                return (mean_distance * last_distance) ** 0.5
            else:
                return self.getSlope(distances)
        else:
            return distances, last_distance

    def compute_center_of_mass(self, traj, select=None):
        import numpy as np
        """
        Used by distance() function below.
        Compute the center of mass for each frame.
        Parameters
        ----------
        traj : Trajectory
            Trajectory to compute center of mass for
        select : str, optional, default=all
            a mdtraj.Topology selection string that
            defines the set of atoms of which to calculate
            the center of mass, the default is all atoms
        Returns
        -------
        com : np.ndarray, shape=(n_frames, 3)
             Coordinates of the center of mass for each frame
        """
        # create a np array
        self.com = np.empty((traj.n_frames, 3))
        # get all the atoms in the system
        if select is None:
            masses = np.array([a.element.mass for a in traj.top.atoms])
            masses /= masses.sum()
            xyz = traj.xyz
        # get masses of selected atoms
        else:
            atoms_of_interest = traj.topology.select(select)
            masses = np.array([traj.top.atom(i).element.mass for i in atoms_of_interest])
            masses /= masses.sum()
            # compute COM
            xyz = traj.xyz[:, atoms_of_interest]

        for i, x in enumerate(xyz):
            self.com[i, :] = x.astype('float64').T.dot(masses)
        return self.com

    def getDistance(self, par, sel_1, sel_2, trajCount):
        import mdtraj
        import numpy as np
        """Compute distances between COMs"""
        topology = None
        if par['MDEngine'] == 'ACEMD':
            xtc = '%s_%s_wrapped.xtc' % (par['Output'], str(trajCount))
            if par['Forcefield'] == 'CHARMM':
                topology = f'../../system/%s.psf' % par['Topology']
            elif par['Forcefield'] == 'AMBER':
                topology = '../../system/%s.prmtop' % par['Topology']

        else:
            xtc = '%s_%s.xtc' % (par['Output'], str(trajCount))
            topology = '%s.gro' % par['Topology']

        traj = mdtraj.load(xtc, top=topology)
        # let's decide what selections to use according to the number of metrics we want to compute
        # if par['NumberCV'] == '1' and par['Metric_1'] == 'Distance':
        c1 = self.compute_center_of_mass(traj, select=sel_1)
        c2 = self.compute_center_of_mass(traj, select=sel_2)

        # elif par['NumberCV'] == '2' and par['Metric_1'] == 'Distance':
        # c1 = compute_center_of_mass(traj, select=selection_list[0])
        # c2 = compute_center_of_mass(traj, select=selection_list[1])

        # elif par['NumberCV'] == '2' and par['Metric_2'] == 'Distance':
        # c1 = compute_center_of_mass(traj, select=selection_list[2])
        # c2 = compute_center_of_mass(traj, select=selection_list[3])

        distances = []
        # compute distance between elements sam eposiiton in 2 different lists
        for a, b in zip(c1, c2):
            D = np.linalg.norm(a - b) * 10
            distances.append(D)

        n = len(distances)
        mean_distance = sum(distances) / n
        last_distance = distances[-1]

        if par['NumberCV'] == '1':  # we return the walker score or slope already
            if par['Slope'] == 'NO':  # we return the walker score
                distMetric = (mean_distance * last_distance) ** 0.5
                return distMetric

            elif par['Slope'] == 'YES':  # we return the walker slope
                distMetric = self.getSlope(distances, )
                return distMetric

        elif par['NumberCV'] == '2':  # we return all the metric values and the last distance
            return distances, last_distance

