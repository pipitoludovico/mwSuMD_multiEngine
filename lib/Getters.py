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
        psf = None
        xtc = 'wrapped.xtc'

        if self.initialParameters['Forcefield'] == 'CHARMM':
            if self.initialParameters['PSF'] is not None:

                psf = '../../system/%s' % self.initialParameters['PSF']
        elif self.initialParameters['Forcefield'] == 'AMBER':
            psf = '../../system/%s' % self.initialParameters['PRMTOP']
        elif self.initialParameters['Forcefield'] == 'GROMOS':
            for tpr in os.listdir(os.getcwd()):
                if tpr.startswith(self.initialParameters['Output']) and tpr.endswith('.tpr'):
                    psf = tpr

        u = Mda.Universe(psf, xtc)
        sel1 = u.select_atoms(sel_1)
        sel2 = u.select_atoms(sel_2)

        distances = []
        for ts in u.trajectory:
            print("TRAJECTORY TS COUNTER")

            if ts is not None:
                com1 = sel1.center_of_mass()
                com2 = sel2.center_of_mass()

                distance = Mda.lib.distances.distance_array(com1, com2)[0][0]
                distances.append(distance)
        mean_lin = np.mean(distances)
        distMetric = (mean_lin * distances[-1]) ** 0.5
        print(distMetric)
        return distMetric, distances, distances[-1]

    def getContacts(self, sel_1, sel_2):
        psf = None
        xtc = 'wrapped.xtc'

        if self.initialParameters['Forcefield'] == 'CHARMM':
            if self.initialParameters['PSF'] is not None:
                psf = '../../system/%s' % self.initialParameters['PSF']
        elif self.initialParameters['Forcefield'] == 'AMBER':
            psf = '../../system/%s' % self.initialParameters['PRMTOP']
        elif self.initialParameters['Forcefield'] == 'GROMOS':
            for tpr in os.listdir(os.getcwd()):
                if tpr.startswith(self.initialParameters['Output']) and tpr.endswith('.tpr'):
                    psf = tpr

        u = Mda.Universe(psf, xtc)

        selection_1 = u.select_atoms(sel_1)
        selection_2 = u.select_atoms(sel_2)

        timeseries = [(distance_array(selection_1.positions, selection_2.positions, box=u.dimensions) < 3).sum() for ts
                      in
                      u.trajectory if ts is not None]

        mean_contacts = sum(timeseries) / len(timeseries)
        last_contacts = timeseries[-1]

        distMetric = (mean_contacts * last_contacts) ** 0.5

        return distMetric, timeseries, timeseries[-1]

    def getRMSD(self, sel_1, sel_2):
        import MDAnalysis.analysis.rms

        pdb = f'{self.folder}/system/reference/' + str(self.par['REFERENCE'])
        psf = None
        if self.par['MDEngine'] != 'GROMACS':
            xtc = 'wrapped.xtc'
            if self.par['Forcefield'] == 'CHARMM':
                psf = f'{self.folder}/system/%s' % self.par['PSF']
            if self.par['Forcefield'] == 'AMBER':
                psf = f'{self.folder}/system/%s' % self.par['PRMTOP']
        else:
            xtc = '%s_%s.xtc' % (self.par['Output'], str(self.trajCount))
            psf = '%s' % self.par['PRMTOP']
        u = Mda.Universe(psf, xtc)
        ref = Mda.Universe(pdb)
        R = Mda.analysis.rms.RMSD(u, ref, tol_mass=10, select="%s" % sel_1, groupselections=["%s" % sel_2])
        R.run()
        rmsd = R.rmsd.T
        data = list(rmsd[3])
        mean_rmsd = sum(data) / len(data)
        last_rmsd = data[-1]

        distMetric = (mean_rmsd * last_rmsd) ** 0.5
        return distMetric, data, data[-1]

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
                return distMetric, hbonds.count_by_time(), mean_contacts
            else:
                # if no H bond is found, we compute the number of water contacts
                # to determine the degree of solvation of the ligand
                distMetric = (((np.mean(waterCont)) * waterCont[-1]) ** 0.5)
                return distMetric, waterCont, mean_contacts
