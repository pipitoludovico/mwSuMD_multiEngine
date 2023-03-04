from multiprocessing import Process, Queue
import sys
from threading import Thread

from .Getters import *
from .Parser import *


class MetricsParser(mwInputParser):
    def __init__(self, par):
        self.par = par
        super(MetricsParser, self).__init__()
        import warnings
        warnings.filterwarnings(action='ignore')
        self.walkers_metrics = []

    def getChosenMetrics(self):
        if sys.argv[1] == 'serial':
            self.getChosenMetricsSer()
        else:
            self.getChosenMetricsMP()
        return self.walkers_metrics

    def getChosenMetricsSer(self):
        """Compute single metric and return a list of values per frame"""
        threads = []
        for walkFolderNumber in range(1, int(self.par['Walkers'] + 1)):
            threads.append(Thread(target=self.calculateMetrics, args=(walkFolderNumber, self.selection_list)))
        for process in threads:
            process.start()
            process.join()
        print('Serial metrics:')
        print(self.walkers_metrics)

    def getChosenMetricsMP(self):
        q = Queue()
        processes = []
        returns = []
        for walker in range(1, self.par['Walkers'] + 1):
            p = Process(target=self.calculateMetricsMP, args=(q, walker, self.selection_list))
            processes.append(p)
            p.start()
        for process in processes:
            return_ = q.get()
            returns.append(return_)
        for process in processes:
            process.join()
        print('Parallel metrics:')
        for ret in returns:
            self.walkers_metrics.append(ret[0])
        returns.clear()

    def calculateMetricsMP(self, queue, walker, selection_list):
        walkers_metrics = []
        import glob
        os.chdir(f'tmp/walker_' + str(walker))
        print("CALCULATED METRICS IN: " + str(os.getcwd()))
        print(f'Walking into walker: ' + str(walker))
        if (glob.glob('*.xtc') and glob.glob('*.coor')) or (glob.glob('*.xtc') and glob.glob('*.gro')):
            if self.par['Metric_1'] == 'DISTANCE':
                # AGGIUNGERE LOGGERS PER METRIC
                distMetric = Getters(self.par).getDistance(selection_list[0], selection_list[1])
                walkers_metrics.append(distMetric)

            elif self.par['Metric_1'] == 'CONTACTS':
                contactMetric = Getters(self.par).contacts_misc(selection_list[0], selection_list[1])
                walkers_metrics.append(contactMetric)

            elif self.par['Metric_1'] == 'RMSD':
                distMetric = Getters(self.par).getRMSD(selection_list[0], selection_list[1])
                walkers_metrics.append(distMetric)
            else:
                print("NO METRICS FOUND")
                exit()
        os.chdir(self.folder)
        queue.put(walkers_metrics)

    def calculateMetrics(self, walkFolderNumber, selection_list):
        import glob
        os.chdir(f'tmp/walker_' + str(walkFolderNumber))
        print("CALCULATE METRICS IN")
        print(os.getcwd())
        print(f'Walking into walker: ' + str(walkFolderNumber))
        if (glob.glob('*.xtc') and glob.glob('*.coor')) or (glob.glob('*.xtc') and glob.glob('*.gro')):
            if self.par['Metric_1'] == 'DISTANCE':
                # AGGIUNGERE LOGGERS PER METRIC
                distMetric = Getters(self.par).getDistance(selection_list[0], selection_list[1])
                self.walkers_metrics.append(distMetric)

            elif self.par['Metric_1'] == 'CONTACTS':
                contactMetric = Getters(self.par).contacts_misc(selection_list[0], selection_list[1])
                self.walkers_metrics.append(contactMetric)

            elif self.par['Metric_1'] == 'RMSD':
                distMetric = Getters(self.par).getRMSD(selection_list[0], selection_list[1])
                self.walkers_metrics.append(distMetric)
            else:
                print("NO METRICS FOUND")
                exit()
        os.chdir(self.folder)

    def getBestWalker(self, walkers_metrics):
        index = 0
        value = 0

        if self.par['Transition_1'] == 'positive':  # we want the metric to increase
            value = max([i for i in walkers_metrics if i is not None])
            # max_value = max(walkers_metrics) # so we take the maximum value
            index = walkers_metrics.index(value)

        elif self.par['Transition_1'] == 'negative':  # we want the metric to dencrease
            value = min([i for i in walkers_metrics if i is not None])
            # max_value = min(walkers_metrics) # so we take the minimum value
            index = walkers_metrics.index(value)

        return index + 1, value
