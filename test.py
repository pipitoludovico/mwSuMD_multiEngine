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
#
# import os
# def wrapMDA():
#     import MDAnalysis as Mda
#     from MDAnalysis import transformations
#     ext = ('xtc', 'dcd')
#     traj_name = '1'
#     for trajectory in os.listdir(os.getcwd()):
#         if trajectory.startswith(traj_name) and trajectory.endswith(ext):
#             print(trajectory)
#             u = Mda.Universe('1.tpr', trajectory)
#             prot = u.select_atoms("protein")
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
# import numpy as np
#
# arr1 = np.full((9, 2, 3), 16.62242288)
# arr1[:, :, 1] = 66.61325028
# arr1[:, :, 2] = 132.32931789
#
# arr2 = np.full((9, 2, 3), 14.36602568)
# arr2[:, :, 1] = 62.69833784
# arr2[:, :, 2] = 131.89819077
#
# c1 = arr1
# c2 = arr2
#
# distances = [np.linalg.norm(a - b) * 10 for a, b in zip(c1, c2)]
# mean_distance = np.mean(distances)
# last_distance = distances[-1]
#
# print(distances, mean_distance, last_distance)


# def compute_center_of_mass(select=None):
#     import MDAnalysis as Mda
#     psf = 'NEUTRAL_fis.psf'
#     xtc = "wrapped.xtc"
#     u = Mda.Universe(psf, xtc)
#     sele = u.select_atoms(select)
#     arr = np.empty((sele.n_residues, u.trajectory.n_frames, 3))
#
#     # substitute this line if you want to see fancy progress bar
#     # for ts in Mda.log.ProgressBar(u.trajectory):
#     for ts in u.trajectory:
#         arr[:, ts.frame] = sele.center_of_mass()
#     print(arr)
#     return arr
#
#
# def getDistance(sel_1='segid P0', sel_2='segid P2'):
#     c1 = compute_center_of_mass(select=sel_1)
#     c2 = compute_center_of_mass(select=sel_2)
#
#     if len(c1) == 0 or len(c2) == 0:
#         print(
#             f"Your selection {sel_1 + ' ' + sel_2} resulted in 0 atoms."
#             f" Please check your selection in the settings and rerun")
#         exit()
#
#     # Compute distances
#     distances_a = [np.linalg.norm(a - b) * 10 for a, b in zip(c1, c2)]
#     mean_distance_a = np.mean(distances)
#     last_distance_a = distances[-1]
#
#     print((mean_distance_a * last_distance_a) ** 0.5, distances_a, last_distance_a)


# def comp2():
#     import MDAnalysis as Mda
#     import numpy as np
#
#     # Load your trajectory and topology
#     u = Mda.Universe("NEUTRAL_fis.psf", "wrapped.xtc")
#
#     # Define your two selections
#     sel1 = u.select_atoms("segid P0")
#     sel2 = u.select_atoms("segid P2")
#
#     distances = []
#     distanceS_np = []
#     for ts in u.trajectory:
#         # Compute the center of mass of each selection
#         com1 = sel1.center_of_mass()
#         com2 = sel2.center_of_mass()
#
#         distance = Mda.lib.distances.distance_array(com1, com2)[0][0]
#
#         distances_np = [np.linalg.norm(a - b) * 10 for a, b in zip(com1, com2)]
#         distanceS_np.append(distances_np)
#
#         distances.append(distance)
#     # distances = np.array(distances)
#     mean_distance = np.mean(distances)
#     last_distance = distances[-1]
#     print(distances)
#     print(last_distance)
#     print(distanceS_np)
#     print(mean_distance)

# mean_distance = np.mean(distances)
# sel1_pos = [sel1.positions for ts in u.trajectory if ts is not None]
# sel2_pos = [sel2.positions for ts in u.trajectory if ts is not None]

# Compute the distance between the centers of mass

# distance = mda.lib.distances.distance_array(com1, com2)
# print(distance)
# print(sel1_pos)
# print(sel2_pos)
# print(distances, last_distance, mean_distance, (mean_distance * last_distance) ** 0.5)

# print(f"The distance between the centers of mass of the two selections is {distance:.3f} Å.")


# print("COMP2 ")
# comp2()

# TEST STDV for checker
# vals = [1042.1788427529493, 1188.3683748604185, 1231.2038868066752, 1257.244082669362, 1249.1861324877304,
#         1272.7317349357868, 1225.9001828756932, 1227.6679315096676, 1344.876566152387, 671.7112380960203,
#         1595.517741891683, 1621.003152386816, 1568.7403368719906, 1579.8772863631048, 1576.051984570963,
#         1578.0473450206298, 1545.6718431715312, 1295.4830651306302, 1424.487573827727, 1425.1332176774808,
#         1436.2839044548994, 1403.1385924102328, 1397.7928848689737, 1459.6992362493745, 1422.2882951690428]
#
# import numpy as np
#
# std = np.std(vals)
# avg = np.average(vals)
# print(std, avg)
# import numpy as np
#
# initialParameters = {'Transition_1': 'positive', 'Transition_2': 'positive'}
#
# test = ([666.6166123583164, 728.0935421312167, 726.8064722293553], [0, 0, 0],
#         [[[86.48604948230368, 56.14707281802907, 666.6166123583164]],
#          [[64.74193515110267, 60.36174557705699, 728.0935421312167]],
#          [[90.40869781133853, 77.98607617612205, 726.8064722293553]]], [[[0, 0]], [[0, 0]], [[0, 0]]])
#
# allMetricLists = (test[i] for i in [2, 3])
# # we calculate the averages for each element in the sublist
# metric_averages = [np.average([item for sublist in upperList for item in sublist]) for upperList in
#                    allMetricLists]
# # print(metric_averages)
# scores_wm = [
#     ([(i - metric_averages[0]) * (100 / metric_averages[0])]) if initialParameters[
#                                                                      'Transition_1'] == 'positive' else (
#         [-(i - metric_averages[0]) * (100 / metric_averages[0])]) for i in test[0]]
# scores_wm2 = [
#     ([(i - metric_averages[1]) * (100 / metric_averages[1])]) if initialParameters[
#                                                                      'Transition_2'] == 'positive' else (
#         [-(i - metric_averages[1]) * (100 / metric_averages[1])]) for i in test[1]]
#
# import numpy as np
#
# if any(np.isnan(val).any for val in [scores_wm, scores_wm2]):
#     print("nan")


# score_sum = [(x[0] + y[0]) for x, y in zip(scores_wm, scores_wm2)]
# list(score_sum)
# max_score = max(score_sum)
# max_index = score_sum.index(max_score) + 1
# pizza = [[73.68924751717799, 74.52631161362292], [70.89893277755907, 72.13195844938991], [73.54088838232087, 73.86075976198593]]
#
# minore = (min(min(x for x in pizza)))
#
# for idx, lista in enumerate(pizza):
#     if minore in lista:
#         print(idx)
#
# idx_val = [idx +1 for idx, lista in enumerate(pizza) if minore in lista]
#
# print(idx_val)
# pizza2 = ['Score of the best 1.5247694223395782', ' Metrics 1:  0', ' Metrics 2: 0', ' Best Metrics: 67.09209972567722']
#
# print("i valori sono: " + " ".join(pizza2))

# import numpy as np
#
# walkers_metrics = [[68.40681676365901, 66.98947894167819, 67.0462169266036, 66.09016507158805],  # DISTANCE
#                    [64.03433644363041, 62.20138715941137, 61.4184966301327, 59.4528484778242],
#                    [67.20794704993656, 67.94094116284039, 69.59623880159859, 68.68284323180802]], \
#     [[1.7798778082582722, 1.5911268346305518, 1.701341448376421, 1.6158511790391388],  # RMSD
#      [1.5198264916843358, 1.5283420344411225, 1.6476582638556438, 1.5569977220128732],
#      [1.6316067378354169, 1.5273519297784623, 1.6640788946351763, 1.5680215720067414]]
#
# allMetricLists = [walkers_metrics[i] for i in [0, 1]]
#
# averages_metric_1 = [np.average(results) for results in walkers_metrics[0]]
# averages_metric_2 = [np.average(results) for results in walkers_metrics[1]]
#
# scores_wm = [([(i - min(averages_metric_1)) * (100 / min(averages_metric_1))]) for i in walkers_metrics[0]]
# scores_wm2 = [([(i - np.average(averages_metric_1)) * (100 / np.average(averages_metric_1)) for i in walkers_metrics[0]])]
#
#
# print(averages_metric_1, averages_metric_2)
# for x in scores_wm:
#     print(x)
#
# print("")
# for y in scores_wm2:
#     print(y)

# walkers_metrics = [[67.55080172781814, 68.57198471483285, 65.68296837125814, 66.72526699638205],
#                    [69.81932494267221, 71.11762158331567, 72.16996068106717, 72.917310659952],
#                    [62.864904136240746, 60.77628648401, 59.26766231092518, 57.204254997916934]], [
#     [1.4574887950163649, 1.5563416332636282, 1.6905481817541355, 1.6517788841510965],
#     [1.6005589431898808, 1.6115067918659267, 1.5959415322050126, 1.7800017850111511],
#     [1.5601941768899734, 1.6392860180628905, 1.643573651531415, 1.629790253289279]]

# walkers_metrics = ([666.6166123583164, 728.0935421312167, 726.8064722293553], [0, 0, 0],
#         [[[86.48604948230368, 56.14707281802907, 666.6166123583164]],
#          [[64.74193515110267, 60.36174557705699, 728.0935421312167]],
#          [[90.40869781133853, 77.98607617612205, 726.8064722293553]]], [[[0, 0]], [[0, 0]], [[0, 0]]])

# print(walkers_metrics[0])

# allMetricLists = (walkers_metrics[i] for i in [0, 1])
# # we calculate the averages for each element in the sublist
# metric_averages = [np.average([item for sublist in upperList for item in sublist]) for upperList in
#                    allMetricLists]
#
# scores_wm = [([(i - metric_averages[0]) * (100 / metric_averages[0])]) for i in walkers_metrics[0]]
# scores_wm2 = [([(i - metric_averages[1]) * (100 / metric_averages[1])]) for i in walkers_metrics[1]]
#
# score_sum = [(x[0] + y[0]) for x, y in zip(scores_wm, scores_wm2)]
#
# listaConvertita = np.array(score_sum).tolist()
# print(listaConvertita)
#
# massimoScore = (max(sottolista) for sottolista in listaConvertita)
# print(*massimoScore)
# list(score_sum)
# max_score = max(score_sum)
# max_index = score_sum.index(max_score) + 1


# averages per frame walker 1 - Distance:
# [68.23039075764112, 65.43151479943582, 65.59349050450218]
# averages per frame walker 2 - Distance:
# [74.9461452237978, 73.83383057460237, 74.44472274755637]
# averages per frame walker 3 - Distance:
# [73.0118791985253, 68.12678017859682, 75.07388666968403]
# averages per frame walker 1 - RMSD:
# [1.518071303385798, 1.5012931293695824, 1.4881288621530513]
# averages per frame walker 2 - RMSD:
# [1.7496403751351113, 1.8159470142191956, 1.735967463609014]
# averages per frame walker 3 - RMSD:
# [1.5757336715415255, 1.7539980741201482, 1.779617380929741]
# Average Metrics 1 of each walker
# [66.4184653538597, 74.40823284865219, 72.07084868226872]
# Average Metrics 2 of each walker
# [1.5024977649694773, 1.7671849509877735, 1.7031163755304715]
# average poll:
# 36.31172432937805
# SCORE SUMS
# -12.950322868920805 9.781879562271655 3.1684433066491935
#
# initialParameters = {"Transition_1": 'negative', 'Transition_2': 'positive'}
# all_metrics = (
#     [[61.182863915064786, 62.61750030774978, 65.09235961562598],
#      [65.65477800598221, 61.87346259963871, 58.129623807388384],
#      [66.06689021141422, 72.9892678497545, 72.88424573411918]],
#     [[1.4401303820477205, 1.379706569886924, 1.32476885713243],
#      [1.580342428681843, 1.358962083714633, 1.5385896250442517],
#      [1.792912866549573, 1.525100782898183, 1.6315333686322762]])
#
# all_avgs = ([62.96424127948018, 61.88595480433643, 70.64680126509597],
#             [1.381535269689025, 1.4926313791469095, 1.6498490060266775])


# def prova(*args):
#     average_poll = np.average(args[0] + args[1])
#     # average_poll = np.average(all_metrics_per_frame_1 + all_metrics_per_frame_2)
#     print("average poll:")
#     print(average_poll)
#     metric_1_scores = []
#     metric_2_scores = []
#     # for average in avg_metrics_1:
#     for average in args[2]:
#         score_1 = ((average / average_poll) - 1)
#         if initialParameters['Transition_1'] == 'negative':
#             abs(score_1)
#             metric_1_scores.append(abs(score_1))
#         else:
#             metric_1_scores.append(score_1)
#         print("SCORE 1:")
#         print(score_1)
#     # for average in avg_metrics_2:
#     for average in args[3]:
#         score_2 = ((average / average_poll) - 1)
#         if initialParameters['Transition_2'] == 'negative':
#             abs(score_2)
#             metric_2_scores.append(abs(score_2))
#         else:
#             metric_2_scores.append(score_2)
#         print("SCORE 2:")
#         print(score_2)
#     score_sum = [(x + y) for x, y in zip(metric_1_scores, metric_2_scores)]
#     list(score_sum)
#     print("SCORE SUM")
#     print(*score_sum)
#     max_score = max(score_sum)
#     print('MAX')
#     print(max_score)
#     max_index = score_sum.index(max_score) + 1
#     print('\nCV2 BEST METRIC RESULTS')
#     print(max_index, max_score, args[2][max_index - 1], args[3][max_index - 1])
# print("RISULTATI DEL CV1 BEST SCORE")


# ordered multiprocessing outputs
# import multiprocessing as mp
#
#
# def quadrato(q, x):
#     print((x * x))
#     q.put((x * x))
#
# pizza = []
# manager = mp.Manager()
# q = manager.Queue()
# results = []
# with mp.Pool() as pool:
#     for x in range(1, 11):
#         results.append(pool.apply_async(quadrato, args=(q, x)))
#         aaa = q.get()
#         print(x, aaa)
#         pizza.append(aaa)
# pool.terminate()
# print(pizza)


import numpy as np

values = [np.random.uniform(low=4.826, high=5.118,  size=(20,))]




def check(values):
    x = np.array(values)
    x_norm = (x-np.min(x))/(np.max(x)-np.min(x))
    print(values)
    # print(*x_norm)
    print((np.std(x_norm)))
    if np.std(x_norm) < 0.3:
        print('stuck')
        return True


check(values)
