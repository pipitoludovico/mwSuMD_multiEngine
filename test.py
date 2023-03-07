# import MDAnalysis

# bestFrame = 1
#
# universe = MDAnalysis.Universe('NEUTRAL_fis.pdb', 'output_0_wrapped.xtc')
# selectAll = universe.select_atoms('all')
# extensions = ('.pdb', '.coor', '.xtc')  # Creating xtc and coor using our chosen frame:
# for x in extensions:
#     with MDAnalysis.Writer(f'pizza1{x}', selectAll) as bau:
#         for idx, ts in enumerate(universe.trajectory):
#             if idx == bestFrame:
#                 print(bestFrame)
# bau.write(selectAll)
#
# os.system('mv pizza1.pdb focaccia.pdb')
# os.system('mv pizza1.coor focaccia.pdb')
# os.system('mv pizza1.coor focaccia.pdb')
# data = [[0.0, 20.0, 1.7900487442736785, 17.680490294335872]]
# print(data[0][3])

# import os
# import mdtraj as mdt
# os.getcwd()
#
# traj1 = mdt.load('pizza.xtc', top='NEUTRAL_fis.pdb')
# traj2 = mdt.load('pizza1.xtc', top='NEUTRAL_fis.pdb')
# merged = traj1 + traj2
#
#
# merged.save_xtc('spike_hs_merged.xtc')
#
# value = 4.6320443340061805
# value = 10.01262975125534
# test = [(1, 4.6320443340061805), (1, 10.01262975125534)]
# index = [val[1] for val in test].index(value)
#
# print(index)

# monolista = [0]
# print(monolista[0])

# import numpy as np
# vals = [13.84180969, 15.03983472, 10.7068709]
# values_metric = np.asarray(vals)
#
# da qui
# data = {k: v for v, k in enumerate(values_metric)}
# maxDist, time_max = max((dist, value) for dist, value in data.items())
# meanTime = np.array(list(data.values())).mean()
# meanDist = np.array(list(data.keys())).mean()
#
# nume = [(float(value) - meanTime) * (float(key) - meanDist) for key, value in data.items()]
# denume = [(float(value) - meanTime) ** 2 for value in data.values()]
#
#
# slope = float(np.sum(nume)) / float(np.sum(denume))


# {13.84180969: 0, 15.03983472: 1, 10.7068709: 2} <class 'dict'>
# -1.5674693949999998

# import MDAnalysis as mda
# Load the trajectory file and topology file
# u = mda.Universe('NEUTRAL_fis.psf', 'output_0_wrapped.xtc')
#
# dims = u.dimensions
#
# Write the unit cell dimensions to a binary .xsc file
# with open('unitcell.xsc', 'wb') as f:
#     f.write(dims.tobytes())

# par = {'NumberCV': '1'}
#
# if par['NumberCV'] == 1:
#     print("pizza")

# pizza = [22, 135, 9]
#
# resu = max(x for x in pizza if x < 100)
# print(resu)

# batches gpu
# lst = [0, 1]
#
# walkers = 3
# quotient, rest = divmod(walkers, len(lst))
# result = quotient * lst + lst[:rest]
# batches = [result[i:i + len(lst)] for i in range(0, len(result), len(lst))]
# for idx, bat in enumerate(result):
#     print(idx+1, bat)
# print(result)
# print(batches)
# files = len([x for x in os.scandir('trajectories')])
# print(files)

##########################################################
# batch = [[0, 1, 2, 3], [0, 1, 2, 3]]
# count = sum([len(x) for x in batch])
# print(count)
# number_of_walkers = 4

# import multiprocessing as mp
# import os
# batch = [0, 1, 2, 3]
# trajCount = 0
# process = 0
# while process < len(batch):
#     walker_number = process + 1
#     command = f'acemd3 --device {batch[process]} input_{walker_number}_{trajCount}.inp 1> acemd.log'
#     process += 1
#     print(command)
# process = 0
# print(process)
# #############################à
# import time
# import os
# import multiprocessing as mp
#
# print(os.getcwd())
# batch = [0, 1, 2, 3]
#
#
# def runGPU(setting):
#     GPUs = setting[0]
#     walkers = setting[1]
#     for walker in range(1, walkers + 1):
#         os.chdir('tmp/walker_' + str(walker))
#         os.system(f'echo {GPUs[walker - 1]}')
#         os.chdir('../../')
#     command = os.system(f'echo {setting[0]}_{setting[1]} && pwd')
#

# process = 0
# start_time = time.perf_counter()
# with mp.Pool(processes=len(batch)) as pool:
#     pool.map(runGPU, [[batch, len(batch)]])
#     pool.close()
#     pool.join()
# end_time = time.perf_counter()
# final_time = end_time - start_time
# print(final_time)
# command = f'acemd3 --device {batch[process]} input_{walker_number}_{trajCount}.inp 1> acemd.log'
# process += 1
# process = 0
# print("")
#
# par = {}
# ext = ('.psf','.prmtop')
# for file in os.listdir('system'):
#     if file.endswith(ext):
#         par['MDengine'] = 'ACEMD'
#     else:
#         par['MDengine'] = 'GROMACS'


# distMetric = minimumRmsd
# logRMSD(data, mean_rmsd, last_rmsd, distMetric)  # logger da sistemare a parte
# print("\nFrame " + str(bestFrame) + " had the lowest RMSD: " + str(distMetric))
# return bestFrame, distMetric, ""

# divisione dei batch e gpu
# batch = [[0, 1, 2], [0]]
# walker = 1
# for x in batch:
#     gpu = 0
#     while gpu < len(x):
#         print(gpu, walker)
#         gpu += 1
#         walker += 1


# Center of Mass with MDA
# import MDAnalysis as mda
# import numpy as np
#
# selection = 'segid P0 P1 and (resnum 33:239 or resnum 259:323) and name CA'
#
# psf = "NEUTRAL_fis.psf"
# xtc = "wrapped.xtc"
#
# u = mda.Universe(psf, xtc)
# sele = u.select_atoms(selection)
# arr = np.empty((sele.n_residues, u.trajectory.n_frames, 3))
#
# for ts in mda.log.ProgressBar(u.trajectory):
#     arr[:, ts.frame] = sele.center_of_mass()
#
# print(arr)

# formattazione settings
# import pandas as pd
# diz = {'MDEngine': 'ACEMD', 'PSF': 'NEUTRAL_fis.psf', 'PDB': 'NEUTRAL_fis.pdb', 'Parameters':
# ['par_all36_cgenff_empty.prm', 'par_all36_carb_2.prm', 'par_all35_ethers.prm', 'par_all36_na.prm', 'par_all36m_prot.prm',
# 'par_all22_prot.prm', 'par_all36_lipid.prm', 'par_all36_cgenff.prm', 'par_all36_prot.prm', 'par_all36_carb.prm', 'parm14sb_all.prm'],
# 'Forcefield': 'CHARMM', 'Timestep': '4', 'Savefreq': '20', 'Wrap': 'protein and name CA and segid P0 P1',
# 'NumberCV': 1, 'Metric_1': 'CONTACTS', 'Cutoff_1': 3, 'Transition_1': 'positive', 'Slope': 'YES', 'Metric_2': 'RMSD',
# 'Cutoff_2': 3, 'Transition_2': 'negative', 'Walkers': 1, 'Timewindow': '500', 'REFERENCE': 'NEUTRAL_fis.pdb',
# 'PLUMED': None, 'Restart': 'NO', 'Output': 'output', 'ligand_HB': '', 'coor': 'NEUTRAL_fis.coor', 'vel': 'NEUTRAL_fis.vel',
# 'xsc': 'previous.xsc'}
# df = pd.DataFrame(list(diz.items()), columns=['keys', 'values'])

# lista = [22, 23, 29, 26, 27, 27, 31, 35, 23, 21, 29, 26, 32, 29, 28, 22, 30, 31, 26, 36, 29, 26, 22, 30, 25]
#
# dizio = dict((k,v) for k,v in enumerate(lista))
# print(dizio)
# sortDict = sorted(dizio)
# print(sortDict)


# MULTITHREADING
# import threading
# import time
#
# class Cane:
#     def __init__(self):
#         self.val = 0
#
#     def paint(self):
#         print('test' + str(self.val))
#         self.val += 1
#
#     def runT(self):
#         for x in range(5):
#             t = threading.Thread(target=self.paint())
#             t.start()
#             t.join()
#
# a = Cane()
# start_time_serial = time.perf_counter()
# a.runT()
# end_time_serial = time.perf_counter()
# final_time_serial = end_time_serial - start_time_serial
# print(final_time_serial)
# print('b')
# b = Cane()
# start = time.perf_counter()
# b.paint()
# b.paint()
# b.paint()
# b.paint()
# b.paint()
# end = time.perf_counter()
# final_time_serial = end - start
# print(final_time_serial)
