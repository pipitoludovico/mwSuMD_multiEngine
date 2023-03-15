from multiprocessing import Process, Queue

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
        self.createMetricList()
        return self.walkers_metrics

    def createMetricList(self):
        walker_metrics_1 = []
        walker_metrics_2 = []
        all_m_1 = []
        all_m_2 = []
        q = Queue()
        processes = [Process(target=self.calculateMetricsMP, args=(q, walker, self.selection_list))
                     for walker in range(1, self.par['Walkers'] + 1)]
        for process in processes:
            process.start()
        for _ in processes:
            ret = q.get()
            if self.par['NumberCV'] == 1:
                self.walkers_metrics.append(ret[0])
            else:
                walker_metrics_1.append(ret[0][0])
                walker_metrics_2.append(ret[1][0])
                all_m_1.append(ret[2])
                all_m_2.append(ret[3])
        for process in processes:
            process.join()
        if self.par['NumberCV'] == 2:
            self.walkers_metrics = walker_metrics_1, walker_metrics_2, all_m_1, all_m_2

    def calculateMetricsMP(self, queue, walker, selection_list):
        import glob
        walker_metric = []
        walkers_metrics_1 = []
        walkers_metrics_2 = []
        allMetric_1 = []
        allMetric_2 = []
        HB_1 = []
        HB_2 = []
        print(f'\nWalking into walker: ' + str(walker))
        os.chdir(f'tmp/walker_' + str(walker))
        print("CALCULATING METRICS IN: " + str(os.getcwd()))
        if (glob.glob('*.xtc') and glob.glob('*.coor')) or (glob.glob('*.xtc') and glob.glob('*.gro')):
            if self.par['Metric_1'] == 'DISTANCE':
                # AGGIUNGERE LOGGERS Pwalkers_metrics_2ER METRIC
                distMetric, distances, lastDist = Getters(self.par).getDistance(selection_list[0], selection_list[1])
                walker_metric.append(distMetric)
                walkers_metrics_1.append(lastDist)
                allMetric_1 = allMetric_1 + distances
            elif self.par['Metric_2'] == 'DISTANCE':
                distMetric, distances, lastDist = Getters(self.par).getDistance(selection_list[2], selection_list[3])
                walkers_metrics_2.append(lastDist)
                allMetric_2 = allMetric_2 + distances

            if self.par['Metric_1'] == 'CONTACTS':
                contactMetric, timeseries, lastContacts = Getters(self.par).getContacts(selection_list[0],
                                                                                        selection_list[1])
                walker_metric.append(contactMetric)
                walkers_metrics_1.append(lastContacts)
                allMetric_1 = allMetric_1 + timeseries
            elif self.par['Metric_2'] == 'CONTACTS':
                contactMetric, timeseries, lastContacts = Getters(self.par).getContacts(selection_list[0],
                                                                                        selection_list[1])
                walkers_metrics_2.append(lastContacts)
                allMetric_2 = allMetric_2 + timeseries

            if self.par['Metric_1'] == 'RMSD':
                distMetric, data, lastRMSD = Getters(self.par).getRMSD(selection_list[0], selection_list[1])
                walker_metric.append(distMetric)
                walkers_metrics_1.append(lastRMSD)
                allMetric_1 = allMetric_1 + data
            elif self.par['Metric_2'] == 'RMSD':
                distMetric, data, lastRMSD = Getters(self.par).getRMSD(selection_list[2], selection_list[3])
                walkers_metrics_2.append(lastRMSD)
                allMetric_2 = allMetric_2 + data

            elif self.par['ligand_HB'] is not None:
                HB_score, last_HB = Getters(self.par).getHB_score()
                HB_1.append(HB_score)
                HB_2.append(last_HB)
                allMetric_2 = allMetric_2 + HB_score
            else:
                print("NO METRICS FOUND, please check your setting file")
                exit()
        os.chdir(self.folder)
        if self.par['NumberCV'] == 1:
            queue.put(walker_metric)
        else:
            queue.put([walkers_metrics_1, walkers_metrics_2, allMetric_1, allMetric_2])

    def getBestWalker(self, walkers_metrics):
        if self.par['NumberCV'] == 1:
            metric_values = [m for m in walkers_metrics if m is not None]
            print(metric_values)
            if not metric_values:
                return None, None
            if self.par['Transition_1'] == 'positive':
                best_value = max(metric_values)
            else:
                best_value = min(metric_values)
            best_walker = walkers_metrics.index(best_value) + 1
            return best_walker, best_value
        else:
            # we create the "allMetrics_1 and 2", walking metrics list 1 and 2
            allMetricLists = (walkers_metrics[i] for i in [2, 3])
            # we calculate the averages for each element in the sublist
            metric_averages = [np.average([item for sublist in upperList for item in sublist]) for upperList in
                               allMetricLists]

            scores_wm = [
                ([(i - metric_averages[0]) * (100 / metric_averages[0])]) if self.par[
                                                                                 'Transition_1'] == 'positive' else (
                    [-(i - metric_averages[0]) * (100 / metric_averages[0])]) for i in walkers_metrics[0]]
            scores_wm2 = [
                ([(i - metric_averages[1]) * (100 / metric_averages[1])]) if self.par[
                                                                                 'Transition_2'] == 'positive' else (
                    [-(i - metric_averages[1]) * (100 / metric_averages[1])]) for i in walkers_metrics[1]]

            score_sum = [(x[0] + y[0]) for x, y in zip(scores_wm, scores_wm2)]
            list(score_sum)
            max_score = max(score_sum)
            max_index = score_sum.index(max_score) + 1
            return max_index, max_score
