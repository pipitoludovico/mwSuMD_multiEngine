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

# from mwSuMD import pars
# print(pars)

# def getSlope(self, values_metric):
#     """Compute the least square methods on the data
#     list provided called by other metrics functions"""
#     data = dict(enumerate(values_metric))
#     # togliere time e prendere solo maxDist da loggare
#     Logger.logSlope(data)
#     meanTime = np.array(list(data.values())).mean()
#     meanDist = np.array(list(data.keys())).mean()
#     self.nume = [(float(value) - meanTime) * (float(key) - meanDist) for key, value in data.items()]
#     self.deNume = [(float(value) - meanTime) ** 2 for value in data.values()]
#     try:
#         slope = float(np.sum(self.nume)) / float(np.sum(self.deNume))
#         return slope
#     except:
#         print("Slope deNumerator was 0.")
#         slope = 0
#     return slope
#
# cv2 = ([13, 13, 13], [2.2666170606145175, 2.2666170606145175, 2.2666170606145175],
#        [[13, 17, 15, 12, 16, 16, 16, 18, 13, 14, 14, 14, 16, 12, 15, 11, 14, 15, 14, 15, 14, 13, 13, 16, 13],
#         [13, 17, 15, 12, 16, 16, 16, 18, 13, 14, 14, 14, 16, 12, 15, 11, 14, 15, 14, 15, 14, 13, 13, 16, 13],
#         [13, 17, 15, 12, 16, 16, 16, 18, 13, 14, 14, 14, 16, 12, 15, 11, 14, 15, 14, 15, 14, 13, 13, 16, 13]], [
#            [2.4609648162393216, 2.3908314660477332, 2.466340781125959, 2.5656839415595396, 2.4910222696124085,
#             2.440992108558348, 2.5034504632229435, 2.296762625931633, 2.341785948102454, 2.4385512518704,
#             2.2956150559258255, 2.2627887514507483, 2.3116477575872176, 2.3093611324461505, 2.309551536019565,
#             2.367078671134977, 2.227416655191998, 2.309515315268195, 2.2924107175271904, 2.389840092025217,
#             2.2853145037172635, 2.2543465294538674, 2.2914293407246924, 2.3252992338595564, 2.2666170606145175],
#            [2.4609648162393216, 2.3908314660477332, 2.466340781125959, 2.5656839415595396, 2.4910222696124085,
#             2.440992108558348, 2.5034504632229435, 2.296762625931633, 2.341785948102454, 2.4385512518704,
#             2.2956150559258255, 2.2627887514507483, 2.3116477575872176, 2.3093611324461505, 2.309551536019565,
#             2.367078671134977, 2.227416655191998, 2.309515315268195, 2.2924107175271904, 2.389840092025217,
#             2.2853145037172635, 2.2543465294538674, 2.2914293407246924, 2.3252992338595564, 2.2666170606145175],
#            [2.4609648162393216, 2.3908314660477332, 2.466340781125959, 2.5656839415595396, 2.4910222696124085,
#             2.440992108558348, 2.5034504632229435, 2.296762625931633, 2.341785948102454, 2.4385512518704,
#             2.2956150559258255, 2.2627887514507483, 2.3116477575872176, 2.3093611324461505, 2.309551536019565,
#             2.367078671134977, 2.227416655191998, 2.309515315268195, 2.2924107175271904, 2.389840092025217,
#             2.2853145037172635, 2.2543465294538674, 2.2914293407246924, 2.3252992338595564, 2.2666170606145175]])

# print(*cv2[2])

# result = []
# for sublist in cv2[2]:
#     result += [*sublist]
# print(result)

# avg = np.average(result)
# print(avg)

# a = [[-9.470752089136488], [-9.470752089136488], [-9.470752089136488], [22]]
# b = [[-3.785051307914562], [-3.785051307914562], [-3.785051307914562], [200]]

# score_sum = [(x[0] + y[0]) for x, y in zip(a, b)]
# max_score = max(score_sum)
# list_conv = np.array(score_sum).tolist()
# test = list(score_sum)
# max_index = test.index(max_score) + 1
# print(max_score, max_index)
#
# import argparse
#
# import GPUtil
#
#
# def getGpus():
#     ap = argparse.ArgumentParser()
#     ap.add_argument('-e', '--exclude', nargs='*', required=False,
#                     help=' use -e to exclude a list of GPUs from being used by mwSuMD: e.g. -e 0 3')
#     args = ap.parse_args()
#
#     max_memory = 0.8
#     GPUs = GPUtil.getGPUs()
#     freeMemory = 0
#     gpu_ids = []
#     for GPU in GPUs:
#         if GPU.memoryUtil > max_memory:
#             continue
#         if GPU.memoryFree >= freeMemory:
#             freeMemory = GPU.memoryFree
#             gpu_ids.append(GPU.id)
#     gpu_ids = [0, 1, 2]
#     if args.exclude is not None:
#         for x in args.exclude:
#             print(x)
#             gpu_ids.remove(int(x))
#     print('after')
#     print(gpu_ids)
#
# getGpus()

# data = [[13, 13, 13], [2.2666170606145175, 2.2666170606145175, 2.2666170606145175]]
# print(data[0])

# import argparse
# def argumentParser():
#     excludeGPUS = []
#     ap = argparse.ArgumentParser()
#     ap.add_argument("-m", '--mode', type=str, default='parallel', required=False,
#                     help="specify -m parallel or serial mode [Default = parallel]")
#     ap.add_argument('-e', '--exclude', nargs='*', required=False,
#                     help=' use -e to exclude a list of GPUs from being used by mwSuMD: e.g. -e 0 3')
#     args = ap.parse_args()
#     if 'parallel' in vars(args).values():
#         mode = 'parallel'
#     else:
#         mode = 'serial'
#     if args.exclude is not None and len(args.exclude) != 0:
#         excludeGPUS = [x for x in args.exclude]
#
#
# argumentParser()

# import multiprocessing as mp
# from multiprocessing import Manager, Queue
# import os
#
#
# class Test:
#     def __init__(self):
#         pass
#
#     def wrap(self, fold_num, q):
#         os.chdir('tmp/walker_' + str(fold_num))
#         print("I'm into " + str(os.getcwd()))
#         for file in os.listdir(os.getcwd()):
#             if file.endswith('.xtc'):
#                 print(file)
#         q.put(fold_num)
#         os.chdir('../../')
#
#     def run(self):
#         print("pooling")
#         manager = Manager()
#         q = manager.Queue()
#         with mp.Pool() as pool:
#             results = pool.starmap(self.wrap, [(x, q) for x in range(1, 8)])
#             for x in results:
#                 print(q.get())

# processes = [mp.Process(target=self.wrap, args=(x, q)) for x in range(1, 8)]
# for proc in processes:
#     proc.start()
# for proc in processes:
#     proc.join()
# while not q.empty():
#     print(q.get())


# test = Test()
# test.run()
#
# import numpy as np
#
# k2 = ([13, 13, 13], [13, 13, 13],
#       [[13, 17, 15, 12, 16, 16, 16, 18, 13, 14, 14, 14, 16, 12, 15, 11, 14, 15, 14, 15, 14, 13, 13, 16, 13],
#        [13, 17, 15, 12, 16, 16, 16, 18, 13, 14, 14, 14, 16, 12, 15, 11, 14, 15, 14, 15, 14, 13, 13, 16, 13],
#        [13, 17, 15, 12, 16, 16, 16, 18, 13, 14, 14, 14, 16, 12, 15, 11, 14, 15, 14, 15, 14, 13, 13, 16, 13]],
#       [[13, 17, 15, 12, 16, 16, 16, 18, 13, 14, 14, 14, 16, 12, 15, 11, 14, 15, 14, 15, 14, 13, 13, 16, 13],
#        [13, 17, 15, 12, 16, 16, 16, 18, 13, 14, 14, 14, 16, 12, 15, 11, 14, 15, 14, 15, 14, 13, 13, 16, 13],
#        [13, 17, 15, 12, 16, 16, 16, 18, 13, 14, 14, 14, 16, 12, 15, 11, 14, 15, 14, 15, 14, 13, 13, 16, 13]])
#
#
# def getBestWalker(walkers_metrics):
#     par = {'Walkers': 3, 'Transition_1': 'positive', 'Transition_2': 'negative'}
#     # we create the "allMetrics_1 and 2", walking metrics list 1 and 2
#     print("METRICS FROM GETBEST")
#     print(walkers_metrics)
#     # we create the "allMetrics_1 and 2", walking metrics list 1 and 2
#     allMetricLists = (walkers_metrics[i] for i in [2, 3])
#     # we calculate the averages for each element in the sublist
#     metric_averages = [np.average([item for sublist in upperList for item in sublist]) for upperList in
#                        allMetricLists]
#
#     scores_wm = [
#         ([(i - metric_averages[0]) * (100 / metric_averages[0])]) if par[
#                                                                          'Transition_1'] == 'positive' else (
#             [-(i - metric_averages[0]) * (100 / metric_averages[0])]) for i in walkers_metrics[0]]
#     scores_wm2 = [
#         ([(i - metric_averages[1]) * (100 / metric_averages[1])]) if par[
#                                                                          'Transition_2'] == 'positive' else (
#             [-(i - metric_averages[1]) * (100 / metric_averages[1])]) for i in walkers_metrics[1]]
#
#     score_sum = [(x[0] + y[0]) for x, y in zip(scores_wm, scores_wm2)]
#     list(score_sum)
#     print(score_sum)
#     max_score = max(score_sum)
#     max_index = score_sum.index(max_score) + 1
#     print("BEST RESULTS FROM INSIDE")
#     print(max_index, max_score)
#     return max_index, max_score
#
#
# getBestWalker(k2)

# import numpy as np
#
# values = [np.random.uniform(low=5.5, high=10, size=(100,))]
#
#
# def check(values):
#         print(np.std(values))
#         if np.std(values) < 0.6:
#             print('done')
#             return True
#
#
# check(values)

# x = 90.94356089922285
# y = 113.44331582425268
# s1 = 3
# s2 = 200
#
# k = 0
# while x > s1 and y < s2:
#     print('OK')
#     k += 1
#     if k == 10:
#         break

# u = Mda.Universe('system/NEUTRAL_fis.psf', 'system/output_0_wrapped.xtc')
# lig_sele = u.select_atoms('segid P2')
# print(lig_sele.n_atoms)

# tested with xtc and dcd and it works.
# def wrapMDA():
#     os.chdir('system')
#     import MDAnalysis as Mda
#     from MDAnalysis import transformations
#     ext = ('dcd', 'xtc')
#     traj_name = 'NEUTRAL_fis'
#     for trajectory in os.listdir(os.getcwd()):
#         if trajectory.startswith(traj_name) and trajectory.endswith(ext):
#             u = Mda.Universe('NEUTRAL_fis.psf', trajectory)
#             prot = u.select_atoms("protein and name CA")
#             ag = u.atoms
#             workflow = (transformations.unwrap(ag),
#                         transformations.center_in_box(prot),
#                         transformations.wrap(ag, compound='fragments'))
#             u.trajectory.add_transformations(*workflow)
#
#             with Mda.Writer('wrapped_MDA_TEST_DCD.xtc', ag) as w:
#                 for ts in u.trajectory:
#                     if ts is not None:
#                         w.write(ag)
#
#
# wrapMDA()


# par = {'Relax': True, 'TW': 400}
#
# _tw = par['TW']
#
# for x in range(1, 11):
#     if par['Relax'] is True:
#         par['TW'] = 999
#     if x == 5:
#         par['Relax'] = False
#         par['TW'] = _tw
#         # par['TW'] = 999999999
#
# print(par['TW'])

# Autoext getter
# par = {}
# groExtensions = ('.psf', '.pdb', '.mdp', '.gro', '.cpt', '.itp', 'top')
#
# def getThings():
#     for ext in groExtensions:
#         for file in os.listdir('./system'):
#             if file.endswith(ext):
#                 par[ext.replace('.', '').upper()] = file
#
# getThings()
# print(par['MDP'])

# wrapping AMBER
# def wrapMDA():
#     import MDAnalysis as Mda
#     from MDAnalysis import transformations
#     ext = ('xtc', 'dcd')
#     traj_name = 'neutral'
#     for trajectory in os.listdir(os.getcwd()):
#         if trajectory.startswith(traj_name) and trajectory.endswith(ext):
#             print(trajectory)
#             u = Mda.Universe('neutral.prmtop', trajectory)
#             prot = u.select_atoms("segid P0")
#             if len(prot.atoms) == 0:
#                 print("your wrapping selection selected 0 atoms! using protein and name CA instead...")
#                 prot = u.select_atoms('protein and name CA')
#             ag = u.atoms
#             workflow = (transformations.unwrap(ag),
#                         transformations.center_in_box(prot),
#                         transformations.wrap(ag, compound='fragments'))
#             u.trajectory.add_transformations(*workflow)
#
#             with Mda.Writer('wrapped_MDA.xtc', ag) as w:
#                 for ts in u.trajectory:
#                     if ts is not None:
#                         w.write(ag)
#
# wrapMDA()

#
# class Test:
#     def __init__(self):
#         self.par = {}
#
#     def getDistance(self, sel_1, sel_2):
#         print("getting c1")
#         c1 = self.compute_center_of_mass(select=sel_1)
#         print("getting c2")
#         c2 = self.compute_center_of_mass(select=sel_2)
#         print("sel 1")
#         print(sel_1)
#         print(len(c1))
#         print("sel 2")
#         print(sel_2)
#         print(len(c2))
#
#         if len(c1) == 0 or len(c2) == 0:
#             print(
#                 f"Your selection {sel_1 + ' ' + sel_2} resulted in 0 atoms."
#                 f" Please check your selection in the settings and rerun")
#             exit()
#
#         # Compute distances
#         distances = [np.linalg.norm(a - b) * 10 for a, b in zip(c1, c2)]
#         mean_distance = np.mean(distances)
#         last_distance = distances[-1]
#         print('results')
#         print((mean_distance * last_distance) ** 0.5, distances, last_distance)
#         return (mean_distance * last_distance) ** 0.5, distances, last_distance
#
#     def compute_center_of_mass(self, select=None):
#         xtc = "wrapped.xtc"
#         u = Mda.Universe('neutral.prmtop', xtc)
#         sele = u.select_atoms(select)
#         arr = np.empty((sele.n_residues, u.trajectory.n_frames, 3))
#
#         for ts in u.trajectory:
#             arr[:, ts.frame] = sele.center_of_mass()
#         return arr
#
#
# a = Test()
# a.getDistance('resname MOL', 'protein and resnum 90')

# x = 0
#
# if x % 3 == 0:
#     print('ok')

# vals = [6.4415219420149965, 1.4936400458098136, 4.8526007560352316]
# fails = 0
#
# from lib.MDoperations import checkIfStuck
#
# print(checkIfStuck(vals, fails))

# extrapolating namd input file to list
# namdInputFile = []
#
# with open('./production.inp', 'r') as file:
#     for line in file.readlines():
#         namdInputFile.append(line)
#
# print(namdInputFile)
# setting NAMD input file for PRODUCTION
# folder = os.getcwd()
# par = {'MDEngine': 'NAMD', 'PDB': 'neutral.pdb', 'PSF': 'neutral.psf', 'TOP': 'neutral.prmtop',
#        'PRMTOP': 'neutral.prmtop',
#        'Parameters': ['neutral.prmtop', 'par_all36_cgenff_empty.prm', 'par_all36_carb_2.prm', 'par_all35_ethers.prm',
#                       'par_all36_na.prm', 'par_all36m_prot.prm', 'par_all22_prot.prm', 'par_all36_lipid.prm',
#                       'par_all36_cgenff.prm', 'par_all36_prot.prm', 'par_all36_carb.prm', 'parm14sb_all.prm'],
#        'Forcefield': 'CHARMM', 'Timestep': 4, 'Savefreq': 20, 'Wrap': 'protein and name CA and segid P0 P1',
#        'RelaxTime': 5, 'Relax': False, 'NumberCV': 2, 'Metric_1': 'DISTANCE', 'Cutoff_1': 3.0,
#        'Transition_1': 'negative', 'Metric_2': 'RMSD', 'Cutoff_2': 20.0, 'Transition_2': 'positive', 'Walkers': 3,
#        'Timewindow': 40, 'REFERENCE': 'neutral.pdb', 'PLUMED': None, 'Restart': 'NO', 'Output': 'output',
#        'ligand_HB': '', 'coor': 'previous.coor', 'vel': 'previous.vel', 'xsc': 'previous.xsc'}
#
# print((len(os.listdir(f"{folder}/trajectories"))))
# saveFrequency = 40
#
# inputFile = ['structure               %s\n' % par['PSF'],
#              'coordinates             %s\n' % par['PDB'],
#              'set xscfile [open %s]\n' % f"../../system/{par['xsc']}" if par['Restart'] == 'NO'
#              else f"../../restart/{par['xsc']}" if (len(os.listdir(f'{folder}/trajectories'))) != 0
#              else f"../../restart/{par['xsc']}",
#              'proc get_first_ts { xscfile } {\n',
#              '  set fd [open $xscfile r]\n',
#              '  gets $fd\n',
#              '  gets $fd\n',
#              '  gets $fd line\n',
#              '  set ts [lindex $line 0]\n',
#              '  close $fd\n',
#              '  return $ts\n',
#              '}\n',
#              'set firsttime [get_first_ts %s]\n' % f"{folder}/system/{par['xsc'].replace('.xsc','')}.restart.xsc",
#              'firsttimestep\t   $firsttime\n',
#              'set temp                310;\n',
#              'outputName              %s\n' % par['Output'],
#              f'binCoordinates          %s{par["coor"]} \n' % '../../system/',
#              'binVelocities           ../../system/%s;     # velocities from last run (binary)\n' % par['vel'],
#              'extendedSystem          ../../system/%s;     # cell dimensions from last run (binary)\n' % par['xsc'],
#              'restartfreq             %s;               \n' % str(saveFrequency*100),
#              'dcdfreq                 50000;\n',
#              'dcdUnitCell             yes;                # the file will contain unit cell info\n',
#              'xstFreq\t    %s;               \n' % str(saveFrequency*100),
#              'outputEnergies          %s;               \n' % str(saveFrequency*100),
#              'outputTiming\t    %s;               \n' % str(saveFrequency*100),
#              '# Force-Field Parameters\n',
#              "%s;                 \n" % 'paraTypeCharmm          on' if par['Forcefield'] == 'CHARMM' else 'amber on',
#              '# FF PARAMETERS - These are specified by CHARMM\n',
#              'exclude                 scaled1-4           # non-bonded exclusion policy"\n',
#              '1-4scaling              1.0\n', 'switching               on\n',
#              'vdwForceSwitching       on;                 # New option for force-based switching of vdW\n',
#              'cutoff                  12.0;               # may use smaller, maybe 10., with PME\n',
#              'switchdist              10.0;               # cutoff - 2.\n',
#              'pairlistdist            16.0;               # stores the all the pairs with in the distance\n',
#              'stepspercycle           5;                  # SET TO 5 (or lower than 20 if HMR)\n',
#              'pairlistsPerCycle       2;                  # 2 is the default\n',
#              'timestep                2.0;                # fs/step SET 4 is you use HMR\n',
#              'rigidBonds              all;                # Bound constraint\n',
#              'nonbondedFreq           1;                  # nonbonded forces every step\n',
#              'fullElectFrequency      1;                  # PME every step\n', '\n',
#              'wrapWater               on;                 # wrap water to central cell\n',
#              'wrapAll                 on;                 # wrap other molecules too\n',
#              '# PME (for full-system periodic electrostatics)\n',
#              'PME                     yes;\n',
#              'PMEInterpOrder          6;                  # interpolation order (spline order 6 in charmm)\n',
#              'PMEGridSpacing          1.0;                # maximum PME grid space / used to calculate grid size\n',
#              '\n', '# Constant Pressure Control (variable volume)\n',
#              'useGroupPressure        yes;                # use a hydrogen-group based pseudo-molecular viral\n',
#              'useFlexibleCell         yes;                 # yes for anisotropic system like membrane\n',
#              'useConstantRatio        yes;                 # keeps the ratio of the unit cell.\n',
#              '\n', '# Constant Temperature Control\n',
#              'langevin                on;                 # langevin dynamics\n',
#              'langevinDamping         1.0;                # damping coefficient of 1/ps (keep low) 0.5 if HMR\n',
#              'langevinTemp            $temp;              # random noise at this level\n',
#              "langevinHydrogen        off;                # don't couple bath to hydrogens\n", '\n',
#              '# Constant pressure\n',
#              'langevinPiston          off;                # Nose-Hoover Langevin piston pressure control\n',
#              'langevinPistonTarget    1.01325;            # target pressure in bar 1atm = 1.01325bar\n',
#              'langevinPistonPeriod    50.0;               # oscillation period in fs\n',
#              'langevinPistonDecay     25.0;               # oscillation decay time.\n',
#              'langevinPistonTemp      $temp;              # coupled to heat bath\n', '\n',
#              'margin 1\n',
#              '\n', 'numsteps                %s;             \n' % int(par['Timewindow'] / (par['Timestep'] / 1000000)),
#              'run                     %s;             \n' % int(par['Timewindow'] / (par['Timestep'] / 1000000))]

# gro = ['title                   = OPLS Lysozyme NPT equilibration \n',
#        '; Run parameters\n',
#        'integrator              = md        ; leap-frog integrator\n',
#        'nsteps                  = 500000    ; 2 * 500000 = 1000 ps (1 ns)\n',
#        'dt                      = 0.002     ; 2 fs\n',
#        '; Output control\n',
#        'nstxout                 = 0         ; suppress bulky .trr file by specifying \n',
#        'nstvout                 = 0         ; 0 for output frequency of nstxout,\n',
#        'nstfout                 = 0         ; nstvout, and nstfout\n',
#        'nstenergy               = 5000      ; save energies every 10.0 ps\n',
#        'nstlog                  = 5000      ; update log file every 10.0 ps\n',
#        'nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps\n',
#        'compressed-x-grps       = System    ; save the whole system\n', '; Bond parameters\n',
#        'continuation            = yes       ; Restarting after NPT \n',
#        'constraint_algorithm    = lincs     ; holonomic constraints \n',
#        'constraints             = h-bonds   ; bonds involving H are constrained\n',
#        'lincs_iter              = 1         ; accuracy of LINCS\n',
#        'lincs_order             = 4         ; also related to accuracy\n',
#        '; Neighborsearching\n',
#        'cutoff-scheme           = Verlet    ; Buffered neighbor searching\n',
#        'ns_type                 = grid      ; search neighboring grid cells\n',
#        'nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme\n',
#        'rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)\n',
#        'rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)\n',
#        '; Electrostatics\n',
#        'coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics\n',
#        'pme_order               = 4         ; cubic interpolation\n',
#        'fourierspacing          = 0.16      ; grid spacing for FFT\n',
#        '; Temperature coupling is on\n',
#        'tcoupl                  = V-rescale             ; modified Berendsen thermostat\n',
#        'tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate\n',
#        'tau_t                   = 0.1     0.1           ; time constant, in ps\n',
#        'ref_t                   = 300     300           ; reference temperature, one for each group, in K\n',
#        '; Pressure coupling is on\n',
#        'pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT\n',
#        'pcoupltype              = isotropic             ; uniform scaling of box vectors\n',
#        'tau_p                   = 2.0                   ; time constant, in ps\n',
#        'ref_p                   = 1.0                   ; reference pressure, in bar\n',
#        'compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1\n',
#        '; Periodic boundary conditions\n',
#        'pbc                     = xyz       ; 3-D PBC\n',
#        '; Dispersion correction\n',
#        'DispCorr                = EnerPres  ; account for cut-off vdW scheme\n',
#        '; Velocity generation\n',
#        'gen_vel                 = no        ; Velocity generation is off \n']
#
# with open('prod.mdp', 'r') as test:
#     for line in test.readlines():
#         gro.append(line)
#
# print(gro)

# wrapping GROMACS

import os
def wrapMDA():
    import MDAnalysis as Mda
    from MDAnalysis import transformations
    ext = ('xtc', 'dcd')
    traj_name = '1'
    for trajectory in os.listdir(os.getcwd()):
        if trajectory.startswith(traj_name) and trajectory.endswith(ext):
            print(trajectory)
            u = Mda.Universe('1.tpr', trajectory)
            prot = u.select_atoms("protein")
            if len(prot.atoms) == 0:
                print("your wrapping selection selected 0 atoms! using protein and name CA instead...")
                prot = u.select_atoms('protein and name CA')
            ag = u.atoms
            workflow = (transformations.unwrap(ag),
                        transformations.center_in_box(prot),
                        transformations.wrap(ag, compound='fragments'))
            u.trajectory.add_transformations(*workflow)

            with Mda.Writer('wrapped_MDA.xtc', ag) as w:
                for ts in u.trajectory:
                    if ts is not None:
                        w.write(ag)

wrapMDA()