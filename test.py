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
#
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

import multiprocessing as mp
import os


class Test:
    def __init__(self):
        pass

    def wrap(self, fold_num):
        for x in fold_num:
            os.chdir('tmp/walker_' + str(x))
            print("I'm into " + str(os.getcwd()))
            for file in os.listdir(os.getcwd()):
                if file.endswith('.xtc'):
                    print(file)
            os.chdir('../../')

    def run(self):
        print("pooling")
        with mp.Pool(7) as p:
            p.map(self.wrap, [range(1, 8)])
            p.close()
            p.join()

test = Test()
test.run()