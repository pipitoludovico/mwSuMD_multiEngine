import MDAnalysis as Mda
import numpy as np
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.lib.distances import distance_array


def getHB_score():
    import time

    start = time.perf_counter()
    par = {'Slope': 'NO', 'ligand_HB': 'resname UNK', 'Forcefield': 'CHARMM', 'NumberCV': 1}
    psf = "ionized.psf"
    xtc = "ionized.xtc"
    u = Mda.Universe(psf, xtc)

    lig_sele = u.select_atoms(f"{par['ligand_HB']} and name O* or {par['ligand_HB']} and name H*")
    water_sele = u.select_atoms('resname TIP3')

    waterCont = [(distance_array(lig_sele.positions, water_sele.positions, box=u.dimensions) < 3).sum() for ts in
                 u.trajectory]

    hbonds = HydrogenBondAnalysis(universe=u, between=['protein', f'{par["ligand_HB"]}'], d_a_cutoff=3,
                                  d_h_a_angle_cutoff=120, update_selections=False)
    hbonds.run(verbose=True)

    mean_contacts = hbonds.count_by_time().mean()
    last_contact = hbonds.count_by_time()[-1]
    if par['NumberCV'] == 1:
        if par['Slope'] == 'NO':
            distMetric = (mean_contacts * last_contact) ** 0.5
            if distMetric != 0:
                print("risultato normale")
                print(distMetric)
            else:
                distMetric = (((np.mean(waterCont)) * 1 / last_contact) ** 0.5)
                print(distMetric)
    else:
        return waterCont, last_contact
    print('if fail:')
    print((np.mean(waterCont)))
    print((np.sum(waterCont)))
    print(waterCont)

    print(((np.mean(waterCont)) * waterCont[-1]) ** 0.5)
    end = time.perf_counter()
    final = end - start
    print("finale:")
    print(final)


getHB_score()
