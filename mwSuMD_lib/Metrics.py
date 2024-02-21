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
                        self.calculateMetric(metric_1, selection_list[:2], walker)
                        self.calculateMetric(metric_2, selection_list[:2], walker)

                    os.chdir(self.folder)
                    for idx, report_metric in enumerate([metric_1, metric_2]):
                        if report_metric:
                            Logger(self.initialParameters['Root']).logData(idx, walker, report_metric,
                                                                           self.scores[walker][report_metric][
                                                                               'allMetricValues'],
                                                                           np.mean(self.scores[walker][report_metric][
                                                                                       'allMetricValues']),
                                                                           self.scores[walker][report_metric]['lastValue'],
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
                averages = {}
                metric_scores = []
                last_values = {}

                for walker, result in scores.items():
                    averages[walker] = {}
                    last_values[walker] = {}
                    for metric, data in result.items():
                        all_metric_values = data['allMetricValues']
                        average = np.average(all_metric_values)
                        averages[walker][metric] = average
                        if self.initialParameters['Transition_1'] == 'negative':
                            last_values[walker][metric] = min(all_metric_values)
                        else:
                            last_values[walker][metric] = max(all_metric_values)

                    metric_keys = list(averages[walker].keys())
                    numerator = averages[walker][metric_keys[0]]
                    denominator = averages[walker][metric_keys[1]]
                    score = (numerator / denominator) - 1

                    if self.initialParameters['Transition_1'] == 'negative':
                        score = score * -1
                        metric_scores.append(score)
                    else:
                        metric_scores.append(score)
                max_index = metric_scores.index(max(metric_scores))
                max_score = max(metric_scores)
                return max_index+1, max_score, last_values[max_index+1]
        except FileExistsError:
            Logger.LogToFile('a', self.trajCount,
                             "Not all MDs produced the same amount of frames. Please check your MD results and restart.")
            Logger.LogToFile('ad', self.trajCount, traceback.format_exc())
