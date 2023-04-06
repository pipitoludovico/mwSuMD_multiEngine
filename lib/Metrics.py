import multiprocessing as mp
from multiprocessing import Manager

from .Getters import *
from .Loggers import Logger
from .Parser import *


class MetricsParser(mwInputParser):
    def __init__(self):
        super(MetricsParser, self).__init__()
        import warnings
        warnings.filterwarnings(action='ignore')
        self.walkers_metrics = []
        self.metric_in_last_frames = []
        self.score_metrics = []
        self.walkers_number_snapshot = self.initialParameters['Walkers']

    def getChosenMetrics(self):
        self.createMetricList()
        if self.initialParameters['NumberCV'] == 1:
            return self.score_metrics, self.metric_in_last_frames
        else:
            return self.walkers_metrics, self.metric_in_last_frames

    def createMetricList(self):
        if self.initialParameters['Relax'] is True:
            self.initialParameters['Walkers'] = 1
        last_frame_metric_1, last_frame_metric_2, all_m_1, all_m_2 = [], [], [], []
        manager = Manager()
        q = manager.Queue()
        print("Calculating Metrics:")
        results = []
        with mp.Pool() as pool:
            for x in range(1, self.initialParameters['Walkers'] + 1):
                results.append(pool.apply(self.calculateMetricsMP, args=(q, x, self.selection_list)))
                ret = q.get()
                if self.initialParameters['NumberCV'] == 1:
                    self.score_metrics.append(*ret[0])
                    self.metric_in_last_frames.append(*ret[1])
                else:
                    last_frame_metric_1.append(*ret[0][0])
                    last_frame_metric_2.append(*ret[1][0])
                    all_m_1.append(*ret[0][1])
                    all_m_2.append(*ret[1][1])
        if self.initialParameters['NumberCV'] == 2:
            self.metric_in_last_frames = last_frame_metric_1, last_frame_metric_2,
            self.walkers_metrics = all_m_1, all_m_2
        self.initialParameters['Walkers'] = self.walkers_number_snapshot

    def calculateMetricsMP(self, queue, walker, selection_list):
        score_metric_1, score_metric_2, last_metrics_1, last_metrics_2, allMetric_1, allMetric_2 = \
            [], [], [], [], [], []
        os.chdir(f'tmp/walker_' + str(walker))
        print(os.getcwd())
        for productions in os.listdir(os.getcwd()):
            if productions == 'wrapped.xtc':
                if self.initialParameters['Metric_1'] == 'DISTANCE':
                    distMetric, distances, lastDist = Getters(self.initialParameters).getDistance(selection_list[0],
                                                                                                  selection_list[1])
                    score_metric_1.append(distMetric)
                    last_metrics_1.append(lastDist)
                    allMetric_1.append(distances)
                if self.initialParameters['Metric_2'] == 'DISTANCE':
                    distMetric, distances, lastDist = Getters(self.initialParameters).getDistance(selection_list[2],
                                                                                                  selection_list[3])
                    last_metrics_2.append(lastDist)
                    allMetric_2.append(distances)

                if self.initialParameters['Metric_1'] == 'CONTACTS':
                    contactMetric, timeseries, lastContacts = Getters(self.initialParameters).getContacts(
                        selection_list[0], selection_list[1])
                    score_metric_1.append(contactMetric)
                    last_metrics_1.append(lastContacts)
                    allMetric_1.append(timeseries)
                if self.initialParameters['Metric_2'] == 'CONTACTS':
                    contactMetric, timeseries, lastContacts = Getters(self.initialParameters).getContacts(
                        selection_list[0], selection_list[1])
                    last_metrics_2.append(lastContacts)
                    allMetric_2.append(timeseries)

                if self.initialParameters['Metric_1'] == 'RMSD':
                    rmsdMetric, data, lastRMSD = Getters(self.initialParameters).getRMSD(selection_list[0],
                                                                                         selection_list[1])
                    score_metric_1.append(rmsdMetric)
                    last_metrics_1.append(lastRMSD)
                    allMetric_1.append(data)
                if self.initialParameters['Metric_2'] == 'RMSD':
                    rmsdMetric, data, lastRMSD = Getters(self.initialParameters).getRMSD(selection_list[2],
                                                                                         selection_list[3])
                    last_metrics_2.append(lastRMSD)
                    allMetric_2.append(data)

                if self.initialParameters['Metric_1'] == 'HB':
                    HB_score, hbData, lastHB = Getters(self.initialParameters).getHB_score(selection_list[0],
                                                                                           selection_list[1])
                    score_metric_1.append(HB_score)
                    last_metrics_1.append(lastHB)
                    allMetric_1.append(hbData)
                if self.initialParameters['Metric_2'] == 'HB':
                    HB_score, hbData, lastHB = Getters(self.initialParameters).getHB_score(selection_list[2],
                                                                                           selection_list[3])
                    last_metrics_2.append(lastHB)
                    allMetric_2.append(hbData)

                if self.initialParameters['NumberCV'] == 1:
                    # Queue
                    queue.put((score_metric_1, last_metrics_1, allMetric_1))
                    Logger(self.initialParameters['Root']).logData(
                        1, walker, self.initialParameters['Metric_1'], allMetric_1, np.mean(allMetric_1),
                        allMetric_1[0][-1], score_metric_1)
                if self.initialParameters['NumberCV'] == 2:
                    score_metric_1.clear()
                    score_metric_2.clear()
                    # Queue
                    queue.put([[last_metrics_1, allMetric_1], [last_metrics_2, allMetric_2]])

                    Logger(self.initialParameters['Root']).logData(1, walker, self.initialParameters['Metric_1'],
                                                                   allMetric_1, np.mean(allMetric_1), allMetric_1[-1],
                                                                   last_metrics_1)
                    Logger(self.initialParameters['Root']).logData(2, walker, self.initialParameters['Metric_2'],
                                                                   allMetric_2, np.mean(allMetric_2), allMetric_2[-1],
                                                                   last_metrics_2)
                os.chdir(self.folder)

    def getBestWalker(self, *args):
        """Returns the walker with the best metric"""
        if self.initialParameters['NumberCV'] == 1:
            metric_scores = args[0]
            if self.initialParameters['Transition_1'] == 'positive':
                best_metric = max(metric_scores)
            else:
                best_metric = min(metric_scores)
            best_walker_index = metric_scores.index(best_metric) + 1
            print('\nCV1 Best Index \t Best Metric \t Associated Last Metric')
            print(best_walker_index, best_metric, args[1][best_walker_index - 1])
            return best_walker_index, best_metric, args[1][best_walker_index - 1]

        else:
            average_poll_1 = np.average(args[0])
            average_poll_2 = np.average(args[1])
            print(average_poll_1, average_poll_2)

            metric_1_scores = []
            metric_2_scores = []
            for walker in args[0]:
                score_1 = ((walker[-1] / average_poll_1) - 1)
                if self.initialParameters['Transition_1'] == 'negative':
                    score_1 = score_1 * -1
                    metric_1_scores.append(score_1)
                else:
                    metric_1_scores.append(score_1)
            for walker in args[1]:
                score_2 = ((walker[-1] / average_poll_2) - 1)
                if self.initialParameters['Transition_2'] == 'negative':
                    score_2 = score_2 * -1
                    metric_2_scores.append(score_2)
                else:
                    metric_2_scores.append(score_2)
            score_sum = [(x + y) for x, y in zip(metric_1_scores, metric_2_scores)]
            list(score_sum)
            max_score = max(score_sum)
            max_index = score_sum.index(max_score) + 1
            return max_index, max_score, args[2][max_index - 1], args[3][max_index - 1]
    #
    # @staticmethod
    # def getSlope(values_metric) -> float:
    #     """Compute the least square methods on the data
    #     list provided called by other metrics functions"""
    #     data = dict(enumerate(values_metric))
    #     meanTime = np.array(list(data.values())).mean()
    #     meanDist = np.array(list(data.keys())).mean()
    #     nume = [(float(value) - meanTime) * (float(key) - meanDist) for key, value in data.items()]
    #     deNume = [(float(value) - meanTime) ** 2 for value in data.values()]
    #     try:
    #         slope = float(np.sum(nume)) / float(np.sum(deNume))
    #         return slope
    #     except:
    #         print("Slope deNumerator was 0.")
    #         slope = 0
    #     return slope
