import os

import MDAnalysis as Mda
import numpy as np
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.lib.distances import distance_array

from .Parser import mwInputParser


class Getters(mwInputParser):
    def __init__(self, par):
        super(Getters, self).__init__()
        self.par = par
        self.deNume = None
        self.nume = None
        self.com = None
        self.trajCount = len(os.listdir(f'{self.folder}/trajectories'))

    def getDistance(self, sel_1, sel_2):
        c1 = self.compute_center_of_mass(select=sel_1)
        c2 = self.compute_center_of_mass(select=sel_2)

        if len(c1) == 0 or len(c2) == 0:
            print(
                f"Your selection {sel_1 + ' ' + sel_2} resulted in 0 atoms."
                f" Please check your selection in the settings and rerun")
            exit()

        # Compute distances
        distances = [np.linalg.norm(a - b) * 10 for a, b in zip(c1, c2)]
        mean_distance = np.mean(distances)
        last_distance = distances[-1]

        return (mean_distance * last_distance) ** 0.5, distances, last_distance

    def getContacts(self, sel_1, sel_2):
        import MDAnalysis

        psf = None
        xtc = 'wrapped.xtc'

        if self.par['MDEngine'] == 'ACEMD':
            if self.par['Forcefield'] == 'CHARMM':
                psf = f'{self.folder}/system/%s' % self.par['PSF']
            elif self.par['Forcefield'] == 'AMBER':
                psf = f'{self.folder}/system/%s' % self.par['PRMTOP']

        elif self.par['MDEngine'] == 'GROMACS':
            psf = '%s.gro' % self.par['gro']

        u = MDAnalysis.Universe(psf, xtc)

        sel_1 = u.select_atoms(sel_1)
        sel_2 = u.select_atoms(sel_2)

        timeseries = [(distance_array(sel_1.positions, sel_2.positions, box=u.dimensions) < 3).sum() for ts in
                      u.trajectory if ts is not None]

        mean_contacts = sum(timeseries) / len(timeseries)
        last_contacts = timeseries[-1]

        distMetric = (mean_contacts * last_contacts) ** 0.5
        return distMetric, timeseries, last_contacts

    def getRMSD(self, sel_1, sel_2):
        import MDAnalysis.analysis.rms

        pdb = f'{self.folder}/system/reference/' + str(self.par['REFERENCE'])
        psf = None
        if self.par['MDEngine'] == 'ACEMD':
            xtc = 'wrapped.xtc'
            if self.par['Forcefield'] == 'CHARMM':
                psf = f'{self.folder}/system/%s' % self.par['PSF']
            elif self.par['Forcefield'] == 'AMBER':
                psf = f'{self.folder}/system/%s' % self.par['PRMTOP']
        else:
            xtc = '%s_%s.xtc' % (self.par['Output'], str(self.trajCount))
            psf = '%s' % self.par['PRMTOP']
        u = Mda.Universe(psf, xtc)
        ref = Mda.Universe(pdb)
        R = Mda.analysis.rms.RMSD(u, ref, select="%s" % sel_1, groupselections=["%s" % sel_2])
        R.run()

        rmsd = R.rmsd.T
        data = list(rmsd[3])
        mean_rmsd = sum(data) / len(data)
        last_rmsd = data[-1]

        distMetric = (mean_rmsd * last_rmsd) ** 0.5
        return distMetric, data, last_rmsd

    def getHB_score(self, sel_1, sel_2):
        print("getting HB in folder in: " + str(os.getcwd()))
        xtc = "wrapped.xtc"
        if self.par['MDEngine'] == 'ACEMD':
            psf = f"../../system/{self.par['PSF']}" \
                if self.par['Forcefield'] == 'CHARMM' else f"../../system/{self.par['PRMTOP']}"
        else:
            psf = f"{self.par['gro']}"

        u = Mda.Universe(psf, xtc)
        if sel_1 or sel_2 is not None:
            lig_sele = u.select_atoms(f"(({str(sel_1)} and type O) or ({str(sel_1)} and type H)) "
                                      f"or (({str(sel_2)} and type O) or ({str(sel_2)} and type H))")
            water_sele = u.select_atoms('(resname SOL and name OW) or (type OH2) or (type H1) or (type H2)')
            if lig_sele.n_atoms == 0:
                print("Your ligand selection produced 0 atoms"
                      "Check if your selection is correct or present in the psf/pdb")
                exit()
            if water_sele.n_atoms == 0:
                print("Warning: no molecule waters were detected."
                      "Make sure your system doesn't have implicit solvent or has not been filtered")

            waterCont = [(distance_array(lig_sele.positions, water_sele.positions, box=u.dimensions) < 3).sum() for ts
                         in
                         u.trajectory if ts is not None]

            hbonds = HydrogenBondAnalysis(universe=u, between=[f'{sel_1}', f'{sel_2}'], d_a_cutoff=3,
                                          d_h_a_angle_cutoff=120, update_selections=False)
            hbonds.run(verbose=False)
            mean_contacts = hbonds.count_by_time().mean()
            last_contact = hbonds.count_by_time()[-1]

            distMetric = ((mean_contacts * last_contact) ** 0.5)  # if self.par['NumberCV'] == 1
            # else (((np.mean(waterCont)) * waterCont[-1]) ** 0.5))

            if distMetric != 0:
                return distMetric, hbonds.count_by_time(), last_contact
            else:
                # if no H bond is found, we compute the number of water contacts
                # to determine the degree of solvation of the ligand
                distMetric = (((np.mean(waterCont)) * waterCont[-1]) ** 0.5)
                return distMetric, waterCont, last_contact

    def compute_center_of_mass(self, select=None):
        if self.par['MDEngine'] == 'ACEMD':
            psf = f"../../system/{self.par['PSF']}" \
                if self.par['Forcefield'] == 'CHARMM' else f"../../system/{self.par['PRMTOP']}"
        else:
            psf = f"../../system{self.par['GRO']}"

        xtc = "wrapped.xtc"
        u = Mda.Universe(psf, xtc)
        sele = u.select_atoms(select)
        arr = np.empty((sele.n_residues, u.trajectory.n_frames, 3))

        # substitute this line if you want to see fancy progress bar
        # for ts in Mda.log.ProgressBar(u.trajectory):
        for ts in u.trajectory:
            arr[:, ts.frame] = sele.center_of_mass()
        return arr
