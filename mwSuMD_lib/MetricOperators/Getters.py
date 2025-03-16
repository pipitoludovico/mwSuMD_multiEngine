import os

import MDAnalysis as Mda
import numpy as np
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.lib.distances import distance_array

from mwSuMD_lib.Parsers.InputfileParser import mwInputParser
from mwSuMD_lib.Utilities.Loggers import Logger


class Getters(mwInputParser):
    def __init__(self):
        super(mwInputParser, self).__init__()
        self.deNume = None
        self.nume = None
        self.com = None
        self.trajCount = len([traj for traj in os.listdir(f'{self.folder}/trajectories') if traj.endswith('.xtc')])
        self.selection_error = "One of your selection from setting file was None. Check that your selection matches a part of your system with MDA atomselection language."
        from warnings import filterwarnings
        filterwarnings(action='ignore')

    def GetMetric(self, metric, sel_1, sel_2):
        psf = "filtered.pdb"
        xtc = 'wrapped.xtc'  #
        u = Mda.Universe(psf, xtc)
        if str(metric).startswith('DISTANCE'):
            sel1 = u.select_atoms(sel_1)
            sel2 = u.select_atoms(sel_2)
            distances = []
            for ts in u.trajectory:
                if ts is not None:
                    if len(sel1) == 0 and len(sel2) == 0:
                        Logger.LogToFile('a', self.trajCount, self.selection_error)
                        raise ValueError
                    else:
                        distance = Mda.lib.distances.distance_array(sel1.center_of_geometry(), sel2.center_of_geometry())[0][0]
                        distances.append(distance)
            mean_lin = np.mean(distances)
            distMetric = (mean_lin * distances[-1]) ** 0.5
            return distMetric, distances, distances[-1]

        if str(metric).upper().startswith('CONTACTS'):
            try:
                if len(u.select_atoms(sel_1)) == 0 and len(u.select_atoms(sel_2)) == 0:
                    Logger.LogToFile('a', self.trajCount, self.selection_error)
                    raise ValueError
                else:
                    selection_1 = u.select_atoms(sel_1)
                    selection_2 = u.select_atoms(sel_2)

                    timeseries = [
                        (distance_array(selection_1.positions, selection_2.positions, box=u.dimensions) < 3).sum()
                        for ts in u.trajectory if ts is not None]

                    mean_contacts = sum(timeseries) / len(timeseries)
                    last_contacts = timeseries[-1]

                    distMetric = (mean_contacts * last_contacts) ** 0.5
                    return distMetric, timeseries, timeseries[-1]
            except Exception as e:
                print("Calculation failed:", e)
                Logger.LogToFile('a', self.trajCount, f"Contact calculation in walker {os.getcwd()} failed.")
                raise ArithmeticError

        if str(metric).startswith('RMSD'):
            import MDAnalysis.analysis.rms

            referencePDB = f'{self.folder}/system/reference/' + str(self.initialParameters['REFERENCE'])
            try:
                ref = Mda.Universe(referencePDB)
                if len(u.select_atoms(sel_1)) == 0 or len(u.select_atoms(sel_2)) == 0:
                    Logger.LogToFile('a', self.trajCount, self.selection_error)
                    raise ValueError
                else:
                    R = Mda.analysis.rms.RMSD(u, ref, tol_mass=100, select="%s" % sel_1, groupselections=["%s" % sel_2])
                    R.run()
                    rmsd = R.rmsd.T
                    data = list(rmsd[3])
                    mean_rmsd = sum(data) / len(data)
                    last_rmsd = data[-1]

                    distMetric = (mean_rmsd * last_rmsd) ** 0.5
                    return distMetric, data, data[-1]
            except Exception as e:
                print("RMSD Exception: ", e)
                Logger.LogToFile('a', self.trajCount, f"RMSD calculation in walker {os.getcwd()} failed.")
                raise ArithmeticError

        if str(metric).startswith('HB'):
            try:
                if sel_1 or sel_2 is not None:
                    hbonds = HydrogenBondAnalysis(universe=u, between=[f'{sel_1}', f'{sel_2}'], d_a_cutoff=3,
                                                  d_h_a_angle_cutoff=120, update_selections=False)
                    hbonds.run(verbose=False)
                    mean_contacts = hbonds.count_by_time().mean()
                    last_contact = hbonds.count_by_time()[-1]
                    distMetric = ((mean_contacts * last_contact) ** 0.5)

                    return distMetric, hbonds.count_by_time(), mean_contacts
                else:
                    Logger.LogToFile('a', self.trajCount, self.selection_error)
                    raise ValueError
            except Exception as e:
                print("Contacts Exception: ", e)
                Logger.LogToFile('a', self.trajCount, f"Contacts calculation in walker {os.getcwd()} failed.")
                raise ArithmeticError

        if str(metric).startswith('SOLVATION'):
            try:
                if sel_1 or sel_2 is not None:
                    lig_sele = u.select_atoms(f"({sel_1}) or ({sel_2})")
                    try:
                        water_sele = u.select_atoms('(resname TIP3 and name H*) or (resname TIP3 and name OH2)')
                        if len(water_sele) == 0:
                            water_sele = u.select_atoms(
                                '(resname SOL and name OW) or (type OH2) or (type H1) or (type H2)')
                            if len(water_sele) == 0:
                                water_sele = u.select_atoms('(resname TIP4 and name H*) or (resname TIP4 and name OH2)')
                                if len(water_sele) == 0:
                                    water_sele = u.select_atoms(
                                        '(resname WAT and name H*) or (resname WAT and name OH*)')
                                else:
                                    water_sele = u.select_atoms('not (protein or membrane or ion) ')
                    except:
                        raise SyntaxError('Supported water resnames are TIP3, TIP4, WAT, SOL')
                    if lig_sele.n_atoms == 0:
                        Logger.LogToFile('a', self.trajCount,
                                         "Your ligand selection produced 0 atoms. Check if your selection is correct or present in the psf/pdb")
                        raise SyntaxError
                    if water_sele.n_atoms == 0:
                        Logger.LogToFile('a', self.trajCount,
                                         "Warning: no molecule waters were detected. Make sure your system doesn't have implicit solvent or has not been filtered")
                        raise SyntaxError

                    waterCont = [(distance_array(lig_sele.positions, water_sele.positions, box=u.dimensions) < 3).sum()
                                 for
                                 ts
                                 in
                                 u.trajectory if ts is not None]
                    mean_contacts = np.mean(waterCont)
                    distMetric = (((np.mean(waterCont)) * waterCont[-1]) ** 0.5)
                    return distMetric, waterCont, mean_contacts
                else:
                    Logger.LogToFile('a', self.trajCount, self.selection_error)
                    raise ValueError
            except Exception as e:
                print("Contacts Exception: ", e)
                Logger.LogToFile('a', self.trajCount, f"Contacts calculation in walker {os.getcwd()} failed.")
                raise ArithmeticError
