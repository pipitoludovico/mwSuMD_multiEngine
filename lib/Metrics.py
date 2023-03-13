import sys
from multiprocessing import Process, Queue
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
            threads.append(Thread(target=self.calculateMetrics, args=(None, walkFolderNumber, self.selection_list)))
        for process in threads:
            process.start()
            process.join()
        print('Serial metrics:')
        print(self.walkers_metrics)

    def calculateMetrics(self, queue, walkFolderNumber, selection_list):
        walkers_metrics = []
        import glob
        os.chdir(f'tmp/walker_' + str(walkFolderNumber))
        print("CALCULATE METRICS IN")
        print(f'Walking into walker: ' + str(walkFolderNumber))
        if (glob.glob('*.xtc') and glob.glob('*.coor')) or (glob.glob('*.xtc') and glob.glob('*.gro')):
            if self.par['Metric_1'] == 'DISTANCE':
                # AGGIUNGERE LOGGERS PER METRIC
                distMetric = Getters(self.par).getDistance(selection_list[0], selection_list[1])
                self.walkers_metrics.append(distMetric)
                walkers_metrics.append(distMetric)

            elif self.par['Metric_1'] == 'CONTACTS':
                contactMetric = Getters(self.par).contacts_misc(selection_list[0], selection_list[1])
                self.walkers_metrics.append(contactMetric)
                walkers_metrics.append(contactMetric)

            elif self.par['Metric_1'] == 'RMSD':
                distMetric = Getters(self.par).getRMSD(selection_list[0], selection_list[1])
                self.walkers_metrics.append(distMetric)
                walkers_metrics.append(distMetric)
            else:
                print("NO METRICS FOUND")
                exit()
        os.chdir(self.folder)
        if queue is not None:
            queue.put(walkers_metrics)

    def getChosenMetricsMP(self):
        q = Queue()
        processes = []
        returns = []
        for walker in range(1, self.par['Walkers'] + 1):
            p = Process(target=self.calculateMetrics, args=(q, walker, self.selection_list))
            processes.append(p)
            p.start()
        for process in processes:
            return_ = q.get()
            returns.append(return_)
        for process in processes:
            process.join()
        print('Parallel metrics:')
        print(returns)
        for ret in returns:
            self.walkers_metrics.append(ret[0])
        returns.clear()

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
        #
        # print('INTO THE FUNCTION: orignal 1 and 2 + all metrics 1 and 2', walkers_metrics_1, walkers_metrics_2,
        #       allMetric_1, allMetric_2)
        # scores_walkers_metric_1 = []
        # scores_walkers_metric_2 = []
        #
        # metric_1_avg = sum(allMetric_1) / len(allMetric_1)
        # metric_2_avg = sum(allMetric_2) / len(allMetric_2)
        # print('AVERAGE M1', metric_1_avg)
        # print('AVERAGE M2', metric_2_avg)
        #
        # # final score computation
        # for i in walkers_metrics_1:
        #     if i != None:
        #         score_1 = (i - metric_1_avg) * (100 / metric_1_avg)
        #
        #         if par['Transition_1'] == 'negative':
        #             score_1 = score_1 * -1
        #         print("original value - score metric 1", i, score_1)
        #         scores_walkers_metric_1.append(score_1)
        #     elif i == None:
        #         score_1 = None
        # for i in walkers_metrics_2:
        #     if i != None:
        #         score_2 = (i - metric_2_avg) * (100 / metric_2_avg)
        #
        #         if par['Transition_2'] == 'negative':
        #             score_2 = score_2 * -1
        #         print("original value - score metric 2", i, score_2)
        #         scores_walkers_metric_2.append(score_2)
        #     elif i == None:
        #         score_2 = None
        # # if par['Transition_2'] == None:
        # # score_2 = score_2 * -1
        # # print("original value - score metric 2",i, score_2)
        # # scores_walkers_metric_2.append(score_2)
        #
        # # put the final scores for the two metrics together for each walker and sum them up
        # print("TWO scores LISTS", scores_walkers_metric_1, scores_walkers_metric_2)
        # scores_list = []
        # for (score1, score2) in zip(scores_walkers_metric_1, scores_walkers_metric_2):
        #     if score1 != None and score2 != None:
        #         scores_list.append(score1 + score2)
        #     elif score1 == None or score2 == None:
        #         scores_list.append(None)
        # # get best 	walker accoriding to the final scores sum
        # max_value = max([i for i in scores_list if i is not None])
        # # max_value = max(scores_list)
        # max_index = scores_list.index(max_value) + 1
        # print("SCORES FINAL", scores_list, max_index)
        # return scores_list, max_index
