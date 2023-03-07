import MDAnalysis as Mda
import numpy as np
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.lib.distances import distance_array

from .Loggers import Logger
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

        if self.par['NumberCV'] == 1:  # we return the walker score or slope already
            if self.par['Slope'] == 'NO':  # we return the walker score
                distMetric = (mean_rmsd * last_rmsd) ** 0.5
                return distMetric

            elif self.par['Slope'] == 'YES':  # we return the walker slope
                distMetric = self.getSlope(data)
                return distMetric

        elif self.par['NumberCV'] == 2:  # we return all the metric values and the last distance
            return data, last_rmsd

    def contacts_misc(self, sel_1, sel_2):
        import MDAnalysis

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

        timeseries = [(distance_array(sel_1.positions, sel_2.positions, box=u.dimensions) < 3).sum() for ts in
                      u.trajectory if ts is not None]

        mean_contacts = sum(timeseries) / len(timeseries)
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
        """Compute the least square methods on the data
        list provided called by other metrics functions"""
        data = dict(enumerate(values_metric))
        # togliere time e prendere solo maxDist da loggare
        Logger.logSlope(data)
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
        if self.par['MDEngine'] == 'ACEMD':
            psf = f"../../system/{self.par['PSF']}" \
                if self.par['Forcefield'] == 'CHARMM' else f"../../system/{self.par['prmtop']}"
        else:
            psf = f"{self.par['gro']}"

        xtc = "wrapped.xtc"
        u = Mda.Universe(psf, xtc)
        sele = u.select_atoms(select)
        arr = np.empty((sele.n_residues, u.trajectory.n_frames, 3))

        for ts in Mda.log.ProgressBar(u.trajectory):
            arr[:, ts.frame] = sele.center_of_mass()
        return arr

    def getDistance(self, sel_1, sel_2):
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

    def getHB_score(self):
        import time
        start = time.perf_counter()
        xtc = "wrapped.xtc"
        if self.par['MDEngine'] == 'ACEMD':
            psf = f"../../system/{self.par['PSF']}" \
                if self.par['Forcefield'] == 'CHARMM' else f"../../system/{self.par['prmtop']}"
        else:
            psf = f"{self.par['gro']}"

        u = Mda.Universe(psf, xtc)

        lig_sele = u.select_atoms(f"{self.par['ligand_HB']} and name O* or {self.par['ligand_HB']} and name H*")
        water_sele = u.select_atoms('water')

        waterCont = [(distance_array(lig_sele.positions, water_sele.positions, box=u.dimensions) < 3).sum() for ts in
                     u.trajectory if ts is not None]

        hbonds = HydrogenBondAnalysis(universe=u, between=['protein', f'{self.par["ligand_HB"]}'], d_a_cutoff=3,
                                      d_h_a_angle_cutoff=120, update_selections=False)
        hbonds.run(verbose=True)
        mean_contacts = hbonds.count_by_time().mean()
        last_contact = hbonds.count_by_time()[-1]

        distMetric = ((mean_contacts * last_contact) ** 0.5 if self.par['NumberCV'] == 1 and self.par['Slope'] == 'NO'
                      else (((np.mean(waterCont)) * waterCont[-1]) ** 0.5))

        end = time.perf_counter()
        final = end - start
        print(final)
        if distMetric != 0:
            return distMetric
        if self.par['NumberCV'] == 2:
            return waterCont, last_contact
        else:
            # if no H bond is found, we compute the number of water contacts
            # to determine the degree of solvation of the ligand
            distMetric = (((np.mean(waterCont)) * waterCont[-1]) ** 0.5)
            return distMetric
