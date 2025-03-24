# import numpy as np
#
# a, b = np.array([1, 1, 1, ]), np.array([1000, 2000, 3000, ])
# dist = np.linalg.norm(a-b)
# print(dist)
#
#
# temperature = min(max(610 - (dist * 30), 310))
# print(temperature)
import numpy as np

parametersSnapshot = {'Root': '/home/scratch/ludovico/debugs/rmsd_gradient/dehydrox/restrained/610', 'COMMAND': None,
                      'EXCLUDED_GPUS': None, 'NOGPU': None, 'PLUMED': None,
                      'PSF': '/home/scratch/ludovico/debugs/rmsd_gradient/dehydrox/restrained/610/system/ionized.psf',
                      'PDB': '/home/scratch/ludovico/debugs/rmsd_gradient/dehydrox/restrained/610/system/ionized.pdb',
                      'MDEngine': 'OPENMM', 'Parameters': [
        '/home/scratch/ludovico/debugs/rmsd_gradient/dehydrox/restrained/610/system/extraTopPar/dehydrox.par',
        '/home/scratch/ludovico/debugs/rmsd_gradient/dehydrox/restrained/610/system/extraTopPar/dehydrox.top',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/01b_par_all36_na.prm',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/00a_top_all36_prot.rtf',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/07_toppar_water_ions.str',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/05a_top_all36_lipid.rtf',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/10a_top_all36_cgenff.rtf',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/04a_top_all36_carb.rtf',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/06_toppar_all36_carb_glycolipid.str',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/08_toppar_all36_lipid_inositol.str',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/00b_par_all36_prot.prm',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/01a_top_all36_na.rtf',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/04b_par_all36_carb.prm',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/07_toppar_all36_lipid_cholesterol.str',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/05b_par_all36_lipid.prm',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/09_toppar_all36_lipid_sphingo.str',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/10b_par_all36_cgenff.prm'],
                      'Forcefield': 'CHARMM', 'Metric_1': 'DISTANCE', 'Metric_2': 'CONTACTS', 'CUSTOMFILE': None,
                      'Timestep': 4, 'Savefreq': 20, 'Wrap': 'protein and name CA', 'Fails': 10, 'Tolerance': 0.2,
                      'ACTUAL_DISTANCE': [0, 0, 0], 'NumberCV': 2, 'RelaxTime': 5.0, 'Relax': False, 'CheckEvery': 1000,
                      'Temperature': 610.0, 'WrapEngine': 'VMD', 'WrapOn': 'protein',
                      'FilterOut': 'protein or resname UNL', 'PROTEIN_RESTRAINTS': ['protein'],
                      'MEMBRANE_RESTRAINTS': None, 'LIGAND_RESNAME': None, 'LIGAND_RESNAMES': ['NO'],
                      'Output': 'unbinding', 'Cutoff_1': 30.0, 'Transition_1': 'positive', 'Cutoff_2': 1.0,
                      'Transition_2': 'negative', 'Walkers': 1, 'Timewindow': 40,
                      'REFERENCE': '/home/scratch/ludovico/debugs/rmsd_gradient/dehydrox/restrained/610/system/ionized.pdb'}

cane = {'Root': '/home/scratch/ludovico/debugs/rmsd_gradient/dehydrox/restrained/610', 'COMMAND': None,
        'EXCLUDED_GPUS': None, 'NOGPU': None, 'PLUMED': None,
        'PSF': '/home/scratch/ludovico/debugs/rmsd_gradient/dehydrox/restrained/610/system/ionized.psf',
        'PDB': '/home/scratch/ludovico/debugs/rmsd_gradient/dehydrox/restrained/610/system/ionized.pdb',
        'MDEngine': 'OPENMM', 'Parameters': [
        '/home/scratch/ludovico/debugs/rmsd_gradient/dehydrox/restrained/610/system/extraTopPar/dehydrox.par',
        '/home/scratch/ludovico/debugs/rmsd_gradient/dehydrox/restrained/610/system/extraTopPar/dehydrox.top',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/01b_par_all36_na.prm',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/00a_top_all36_prot.rtf',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/07_toppar_water_ions.str',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/05a_top_all36_lipid.rtf',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/10a_top_all36_cgenff.rtf',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/04a_top_all36_carb.rtf',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/06_toppar_all36_carb_glycolipid.str',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/08_toppar_all36_lipid_inositol.str',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/00b_par_all36_prot.prm',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/01a_top_all36_na.rtf',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/04b_par_all36_carb.prm',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/07_toppar_all36_lipid_cholesterol.str',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/05b_par_all36_lipid.prm',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/09_toppar_all36_lipid_sphingo.str',
        '/home/scratch/software/ludovico/mwsumd/multiEngine/mwSuMD_lib/parameters/10b_par_all36_cgenff.prm'],
        'Forcefield': 'CHARMM', 'Metric_1': 'DISTANCE', 'Metric_2': 'CONTACTS', 'CUSTOMFILE': None, 'Timestep': 4,
        'Savefreq': 20, 'Wrap': 'protein and name CA', 'Fails': 10, 'Tolerance': 0.2, 'ACTUAL_DISTANCE': [0, 0, 0],
        'NumberCV': 2, 'RelaxTime': 5.0, 'Relax': False, 'CheckEvery': 1000, 'Temperature': 610.0, 'WrapEngine': 'VMD',
        'WrapOn': 'protein', 'FilterOut': 'protein or resname UNL', 'PROTEIN_RESTRAINTS': ['protein'],
        'MEMBRANE_RESTRAINTS': None, 'LIGAND_RESNAME': None, 'LIGAND_RESNAMES': ['NO'], 'Output': 'unbinding',
        'Cutoff_1': 30.0, 'Transition_1': 'positive', 'Cutoff_2': 1.0, 'Transition_2': 'negative', 'Walkers': 1,
        'Timewindow': 40,
        'REFERENCE': '/home/scratch/ludovico/debugs/rmsd_gradient/dehydrox/restrained/610/system/ionized.pdb'}

print(cane == parametersSnapshot)
