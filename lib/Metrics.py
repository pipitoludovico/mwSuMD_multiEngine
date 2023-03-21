import multiprocessing as mp
from multiprocessing import Manager

from .Getters import *
from .Parser import *


class MetricsParser(mwInputParser):
    def __init__(self):
        super(MetricsParser, self).__init__()
        import warnings
        warnings.filterwarnings(action='ignore')
        self.walkers_metrics = []
        self.walkers_number_snapshot = self.par['Walkers']

    def getChosenMetrics(self):
        self.createMetricList()
        print(self.walkers_metrics)
        print("METRIC RETURNIN IN FOLDER" + os.getcwd())
        return self.walkers_metrics

    def createMetricList(self):
        if self.par['Relax'] is True:
            self.par['Walkers'] = 1
        walker_metrics_1 = []
        walker_metrics_2 = []
        all_m_1 = []
        all_m_2 = []
        manager = Manager()
        q = manager.Queue()
        print("Calculating Metrics:")
        with mp.Pool() as pool:
            results = pool.starmap(self.calculateMetricsMP,
                                   [(q, walker, self.selection_list) for walker in range(1, self.par['Walkers'] + 1)])
            for _ in results:
                ret = q.get()
                if self.par['NumberCV'] == 1:
                    self.walkers_metrics.append(ret[0])
                else:
                    print("PRINTING RETS")
                    print(ret)
                    walker_metrics_1.append(ret[0][0])
                    walker_metrics_2.append(ret[1][0])
                    all_m_1.append(ret[2])
                    all_m_2.append(ret[3])
        pool.join()
        if self.par['NumberCV'] == 2:
            self.walkers_metrics = walker_metrics_1, walker_metrics_2, all_m_1, all_m_2
        self.par['Walkers'] = self.walkers_number_snapshot

    def calculateMetricsMP(self, queue, walker, selection_list):
        import glob
        walker_metric = []
        walkers_metrics_1 = []
        walkers_metrics_2 = []
        allMetric_1 = []
        allMetric_2 = []
        os.chdir(f'tmp/walker_' + str(walker))
        if (glob.glob('*.xtc') and glob.glob('*.coor')) or (glob.glob('*.xtc') and glob.glob('*.gro')):
            if self.par['Metric_1'] == 'DISTANCE':
                distMetric, distances, lastDist = Getters(self.par).getDistance(selection_list[0], selection_list[1])
                walker_metric.append(distMetric)
                walkers_metrics_1.append(lastDist)
                allMetric_1.append(distances)
            if self.par['Metric_2'] == 'DISTANCE':
                distMetric, distances, lastDist = Getters(self.par).getDistance(selection_list[2], selection_list[3])
                walkers_metrics_2.append(lastDist)
                allMetric_2.append(distances)

            if self.par['Metric_1'] == 'CONTACTS':
                contactMetric, timeseries, lastContacts = Getters(self.par).getContacts(selection_list[0],
                                                                                        selection_list[1])
                walker_metric.append(contactMetric)
                walkers_metrics_1.append(lastContacts)
                allMetric_1.append(timeseries)
            if self.par['Metric_2'] == 'CONTACTS':
                contactMetric, timeseries, lastContacts = Getters(self.par).getContacts(selection_list[0],
                                                                                        selection_list[1])
                walkers_metrics_2.append(lastContacts)
                allMetric_2.append(timeseries)

            if self.par['Metric_1'] == 'RMSD':
                distMetric, data, lastRMSD = Getters(self.par).getRMSD(selection_list[0], selection_list[1])
                walker_metric.append(distMetric)
                walkers_metrics_1.append(lastRMSD)
                allMetric_1.append(data)
            if self.par['Metric_2'] == 'RMSD':
                distMetric, data, lastRMSD = Getters(self.par).getRMSD(selection_list[2], selection_list[3])
                walkers_metrics_2.append(lastRMSD)
                allMetric_2.append(data)

            if self.par['Metric_1'] == 'HB':
                HB_score, hbData, last_HB = Getters(self.par).getHB_score(selection_list[0], selection_list[1])
                walker_metric.append(HB_score)
                walkers_metrics_1.append(HB_score)
                allMetric_1.append(hbData)
            if self.par['Metric_2'] == 'HB':
                HB_score, hbData, last_HB = Getters(self.par).getHB_score(selection_list[2], selection_list[3])
                walkers_metrics_2.append(last_HB)
                allMetric_2.append(hbData)
            # wrapping up all the computed metrics into shared memory values
            if self.par['NumberCV'] == 1:
                queue.put(walker_metric)
            if self.par['NumberCV'] == 2:
                queue.put([walkers_metrics_1, walkers_metrics_2, allMetric_1, allMetric_2])
        os.chdir(self.folder)

    def getBestWalker(self, walkers_metrics):
        if self.par['NumberCV'] == 1:
            metric_values = [m for m in walkers_metrics if m is not None]
            if self.par['Transition_1'] == 'positive':
                best_value = max(metric_values)
            else:
                best_value = min(metric_values)
            best_walker = walkers_metrics.index(best_value) + 1
            return best_walker, best_value, walkers_metrics[-1], walkers_metrics[-1]
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
            if max_score is 'nan':
                raise ZeroDivisionError("The reference value you chose did not produce results (averages = 0)."
                                        "Please try a different selection and start again")
            max_index = score_sum.index(max_score) + 1
            return max_index, max_score, walkers_metrics[0][-1], walkers_metrics[1][-1]
