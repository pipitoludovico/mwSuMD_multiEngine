# import os
# import re
#
# import numpy as np
# from .Metrics import MetricsParser
# from .openRunners import Runner
# from .Loggers import Logger
#
#
# class MDoperator:
#     def __init__(self, par, root):
#         self.par = par
#         self.folder = root
#         self.extensions = ('.chk', '.xml')
#         self.cycle = len(
#             [trajFile for trajFile in os.listdir(f'{self.folder}/trajectories') if trajFile.endswith('xtc')])
#         self.best_metric_result = None
#         self.best_average_metric_2 = None
#         self.best_average_metric_1 = None
#         self.best_walker_score = None
#         self.averages = None
#         self.scores = None
#         self.best_value = None
#         self.bestWalker = None
#         self.walker_metrics = []
#
#     def saveStep(self, best_walker, walker_score, best_metric_result):
#         """Handle the restart files and the xtc storage"""
#         os.chdir('tmp/walker_%s' % str(best_walker))
#         check = [binary for binary in os.listdir("./") if binary.endswith(self.extensions)]
#         if check:
#             self.cycle += 1
#             os.system(f'cp wrapped.xtc {self.folder}/trajectories/{self.par["Output"]}_step_{self.cycle}.xtc')
#             # moving and renaming the binary files to the restart folder
#             os.system(f'cp *.chk {self.folder}/restarts/previous.chk')
#             os.system(f'cp *.xml {self.folder}/restarts/previous.xml')
#             Logger.LogToFile('a', self.cycle, "FINISHED SAVING FRAMES")
#             if self.par['PLUMED']:
#                 with open(self.par['PLUMED'], 'r') as plumedINP:
#                     for line in plumedINP.readlines():
#                         match = re.search(r'FILE=([^\s]+)', line)
#                         if match:
#                             outFile = match.group(1)
#                             Logger.LogToFile('ad', self.cycle, str(os.getcwd() + outFile))
#                             os.system(f'cp {outFile} %s/restarts/' % self.folder)
#                 for filename in os.listdir("./"):
#                     if '.' not in filename:
#                         fullname = os.path.join("./", filename)
#                         os.system(f'cp {fullname}  %s/restarts/ ' % self.folder)
#                     if filename.endswith(".dat"):
#                         os.system(f'cp {filename} %s/restarts/' % self.folder)
#                 os.system('cp plumed.log  %s/restarts/ ' % self.folder)
#         else:
#             Logger.LogToFile('ad', self.cycle, "No binary saved: restarting from last checkpoint.")
#             with open('walkerSummary.log', 'a') as walkerSummary:
#                 walkerSummary.write(f"No binary produced. Simulation failed at cycle: {self.cycle}")
#
#         os.chdir(self.folder)
#
#         with open('walkerSummary.log', 'a') as walkerSummary:
#             if self.par['NumberCV'] == 1:
#                 if self.par['Relax'] and check:
#                     info_to_write = str(self.cycle) + " RELAXATION SCORE: " + str(walker_score) + " Metrics: " + str(
#                         best_metric_result) + "\n"
#                 if not self.par['Relax'] and check:
#                     info_to_write = str(self.cycle) + " Best Walker: " + str(best_walker) + " Best Metric: " + str(
#                         walker_score) + " Last Metric: " + str(best_metric_result) + "\n"
#             if self.par['NumberCV'] == 2:
#                 if self.par['Relax'] and check:
#                     info_to_write = str(self.cycle) + " RELAXATION SCORE: " + str(walker_score) + " Metrics: " + str(
#                         best_metric_result) + "\n"
#                 if not self.par['Relax'] and check:
#                     info_to_write = str(self.cycle) + " Best Walker: " + str(best_walker) + " Score Result: " + str(
#                         walker_score) + " Last Metrics from best: " + str(best_metric_result) + "\n"
#             walkerSummary.write(info_to_write)
#
#         os.system('rm -r tmp')
#         self.par['Relax'] = False
#
#     def checkIfStuck(self, values, accumulatedFails) -> bool:
#         if accumulatedFails >= self.par['Fails'] * int(self.par['NumberCV']):
#             with open('walkerSummary.log', 'a') as logFile:
#                 logFile.write(
#                     '\nSimulation seemed stuck. It will run the last relaxation protocol and it will be terminated\n')
#                 logFile.close()
#             self.relaxSystem()
#             exit()
#         else:
#             failCount = 0
#             for vals in values:
#                 if vals is not None:
#                     x = np.array(vals)
#                     x_norm = (x - np.min(x)) / (np.max(x) - np.min(x))
#                     if np.std(x_norm) < self.par['Tolerance']:
#                         failCount += 1
#             if failCount == self.par['NumberCV']:
#                 print('\nSimulation might be stuck. Running the relaxation protocol.')
#                 print("")
#                 return True
#             else:
#                 # if it did not, we continue as usual
#                 return False
#
#     def relaxSystem(self):
#         Logger.LogToFile('a', self.cycle, 'Relaxation Protocol begins now:\n' + ('#' * 200))
#         self.par['Relax'] = True
#         Runner().runAndWrap()
#         self.scores = MetricsParser().getChosenMetrics()
#         # we then extract the best metric/score and store it as a reference
#         self.bestWalker, self.best_walker_score, self.best_metric_result = None, None, None
#         self.bestWalker, self.best_walker_score, self.best_metric_result = MetricsParser().getBestWalker(self.scores)
#         MDoperator(self.par, self.folder).saveStep(self.bestWalker, self.best_walker_score, self.best_metric_result)
#
#         Logger.LogToFile('ad', self.cycle, "\nRelaxation Protocol Ended\n" + "#" * 200)
#         self.cycle += 1
#         # setting our check to False and end the protocol, beginning a new cycle.
#         self.par['Relax'] = False
