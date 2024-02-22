import traceback

from .Getters import *
from .Loggers import Logger
from .Parser import *

from signal import signal, SIGPIPE, SIG_DFL

filterwarnings(action='ignore')
signal(SIGPIPE, SIG_DFL)


class MetricsParser(mwInputParser):

    def __init__(self):
        super(mwInputParser, self).__init__()
        self.scores = {}
        self.walkers_number_snapshot = self.initialParameters['Walkers']
        self.trajCount = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])

    def getChosenMetrics(self):
        self.createMetricList()
        return self.scores

    def createMetricList(self):
        if self.initialParameters['Relax'] is True:
            self.initialParameters['Walkers'] = 1
            Logger.LogToFile('a', self.trajCount,
                             f"Relax mode: calculating metrics: {self.initialParameters['Relax']}, Walker set to relax: {self.initialParameters['Walkers']}")
        subprocesses = []
        try:
            for walker in range(1, self.initialParameters.get('Walkers') + 1):
                Logger.LogToFile('a', self.trajCount, f"Calculating metrics in walker #: {walker}")
                subprocesses.append(self.calculateMetrics(walker, self.selection_list))
            self.initialParameters['Walkers'] = self.walkers_number_snapshot
        except FileNotFoundError:
            Logger.LogToFile('a', self.trajCount,
                             "\nMetric calculation failed. Check if all the simulations ended well.")
            exit()

    def calculateMetrics(self, walker, selection_list):
        os.chdir(f'tmp/walker_{walker}')
        try:
            for productions in os.listdir(os.getcwd()):
                if productions == 'wrapped.xtc':
                    metric_1 = self.initialParameters.get('Metric_1')
                    metric_2 = self.initialParameters.get('Metric_2')
                    number_cv = self.initialParameters.get('NumberCV')
                    if metric_1 and (not metric_2 or number_cv == 1):
                        self.calculateMetric(metric_1, selection_list[:2], walker)

                    if not metric_1 and metric_2 and number_cv == 1:
                        self.calculateMetric(metric_2, selection_list[:2], walker)

                    if metric_1 and metric_2 and number_cv == 2:
                        if metric_1 == metric_2:
                            metric_2 += "_2"
                        print("SELECTION LIST 1:", self.selection_list[:2])
                        print("SELECTION LIST 2:", self.selection_list[2:])
                        self.calculateMetric(metric_1, selection_list[:2], walker)
                        self.calculateMetric(metric_2, selection_list[2:], walker)

                    os.chdir(self.folder)
                    for idx, report_metric in enumerate([metric_1, metric_2]):
                        if report_metric:
                            Logger(self.initialParameters['Root']).logData(idx, walker, report_metric,
                                                                           self.scores[walker][report_metric][
                                                                               'allMetricValues'],
                                                                           np.mean(self.scores[walker][report_metric][
                                                                                       'allMetricValues']),
                                                                           self.scores[walker][report_metric][
                                                                               'lastValue'],
                                                                           self.scores[walker][report_metric][
                                                                               'scoreMetric'])
        except Exception:
            Logger.LogToFile('a', self.trajCount,
                             "Metric collection and calculation failed. Make sure the trajectories produced make sense and were successfully completed")
            Logger.LogToFile('a', self.trajCount, traceback.format_exc())
            raise IOError

    def calculateMetric(self, metric, selection_list, walker):
        getter = Getters()
        scoreMetric, allMetricValues, lastValue = getter.GetMetric(metric, selection_list[0], selection_list[1])
        if walker not in self.scores:
            self.scores[walker] = {}
        if metric not in self.scores[walker]:
            self.scores[walker][metric] = {}

        self.scores[walker][metric]['scoreMetric'] = scoreMetric
        self.scores[walker][metric]['lastValue'] = lastValue
        self.scores[walker][metric]['allMetricValues'] = allMetricValues

    def getBestWalker(self, scores: dict):
        """Returns the walker with the best metric"""
        print("SCORES ", scores)
        try:
            if self.initialParameters['NumberCV'] == 1:

                if self.initialParameters['Transition_1'] == 'positive':
                    best_metric_key = max(scores, key=lambda k: scores[k][next(iter(scores[k]))]['allMetricValues'])
                    best_metric_value = scores[best_metric_key][next(iter(scores[best_metric_key]))]['scoreMetric']
                    metric_used = next(iter(scores[best_metric_key]))
                else:
                    best_metric_key = min(scores, key=lambda k: scores[k][next(iter(scores[k]))]['allMetricValues'])
                    best_metric_value = scores[best_metric_key][next(iter(scores[best_metric_key]))]['scoreMetric']
                    metric_used = next(iter(scores[best_metric_key]))
                return best_metric_key, best_metric_value, scores[best_metric_key][metric_used]['allMetricValues'][-1]
            else:
                last_values_walker_1 = []
                last_values_walker_2 = []
                all_metric_new_1 = []
                all_metric_new_2 = []
                scores_walkers_metric_1 = []
                scores_walkers_metric_2 = []

                last_values = {}
                for walker, result in scores.items():
                    last_values[walker] = {}
                    for metric, data in result.items():
                        last_values[walker][metric] = data['lastValue']

                for walker, results in scores.items():
                    for idx, key in enumerate(scores[walker].keys()):
                        if idx == 0:
                            last_values_walker_1.append(scores[walker][key]['lastValue'])
                            all_metric_new_1.append(scores[walker][key]['allMetricValues'])
                        if idx == 1:
                            last_values_walker_2.append(scores[walker][key]['lastValue'])
                            all_metric_new_2.append(scores[walker][key]['allMetricValues'])

                average_last_values_metric1_and_2_walker1 = (np.average(all_metric_new_1))
                average_last_values_metric1_and_2_walker2 = (np.average(all_metric_new_2))

                # print("INTO original 1", last_values_walker_1)
                # print("INTO original 2", last_values_walker_2)
                # print("INTO all metrics 1", allMetric_1)
                # print("INTO all metrics 2", allMetric_2)

                for i in last_values_walker_1:  # wuesto ha last RMSD e last DISTANCE
                    # score_1 = (i - metric_1_avg) * (100 / metric_1_avg)
                    score_1 = (i - average_last_values_metric1_and_2_walker1) * (
                            100 / average_last_values_metric1_and_2_walker1)
                    if self.initialParameters['Transition_1'] == 'negative':
                        score_1 = score_1 * -1
                    # print("SCORE 1, ", score_1, " AVERAGE", average_last_values_metric1_and_2_walker1)
                    scores_walkers_metric_1.append(score_1)
                for i in last_values_walker_2:  # wuesto ha last RMSD e last DISTANCE
                    score_2 = (i - average_last_values_metric1_and_2_walker2) * (
                            100 / average_last_values_metric1_and_2_walker2)
                    if self.initialParameters['Transition_2'] == 'negative':
                        score_2 = score_2 * -1
                    # print("SCORE 2, ", score_2, " AVERAGE", average_last_values_metric1_and_2_walker2)
                    scores_walkers_metric_2.append(score_2)

                # print("AVERAGE M1", average_last_values_metric1_and_2_walker1)
                # print("AVERAGE M2", average_last_values_metric1_and_2_walker2)

                scores_list = []
                for (score1, score2) in zip(scores_walkers_metric_1, scores_walkers_metric_2):
                    scores_list.append(score1 + score2)
                max_value = max([i for i in scores_list if i is not None])
                max_index = scores_list.index(max_value) + 1
                # print(max_index, max_value, (last_values_walker_1, last_values_walker_2))
                return max_index, max_value, last_values[max_index]
        except FileExistsError:
            Logger.LogToFile('a', self.trajCount,
                             "Not all MDs produced the same amount of frames. Please check your MD results and restart.")
            Logger.LogToFile('ad', self.trajCount, traceback.format_exc())
