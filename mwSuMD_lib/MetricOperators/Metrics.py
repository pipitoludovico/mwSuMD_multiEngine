import traceback

from mwSuMD_lib.MetricOperators.Getters import *
from mwSuMD_lib.Utilities.Loggers import Logger
from mwSuMD_lib.Parsers.InputfileParser import *

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
        for walker in range(1, self.initialParameters.get('Walkers') + 1):
            try:
                Logger.LogToFile('a', self.trajCount, f"Calculating metrics in walker #: {walker}")
                subprocesses.append(self.calculateMetrics(walker, self.selection_list))
            except FileNotFoundError:
                Logger.LogToFile('a', self.trajCount, "\nMetric calculation failed. Check if all the simulations ended well.")
                os.chdir(self.folder)
                Logger.LogToFile('a', self.trajCount, f"\nAn error occurred in calculating the metrics. Check the simulations inside walker {walker}")
                exit()
            self.initialParameters['Walkers'] = self.walkers_number_snapshot

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
                        self.calculateMetric(metric_2, selection_list[2:], walker)

                    os.chdir(self.folder)
                    for idx, report_metric in enumerate([metric_1, metric_2]):
                        if report_metric in self.scores[walker]:
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
        try:
            if self.initialParameters['NumberCV'] == 1:

                if self.initialParameters['Transition_1'] == 'positive':
                    best_walker = max(scores, key=lambda k: max(scores[k].values(), key=lambda x: x['scoreMetric'])[
                        'scoreMetric'])
                    best_metric_value = max(scores[best_walker].values(), key=lambda x: x['scoreMetric'])['scoreMetric']
                    metric_used = list(scores[best_walker].keys())[0]
                else:
                    best_walker = min(scores, key=lambda k: max(scores[k].values(), key=lambda x: x['scoreMetric'])[
                        'scoreMetric'])
                    best_metric_value = min(scores[best_walker].values(), key=lambda x: x['scoreMetric'])['scoreMetric']
                    metric_used = list(scores[best_walker].keys())[0]
                return best_walker, best_metric_value, scores[best_walker][metric_used]['allMetricValues'][-1]
            else:
                last_values = {}
                for walker, result in scores.items():
                    last_values[walker] = {}
                    for metric, data in result.items():
                        last_values[walker][metric] = data['lastValue']

                all_metric_new_1 = []
                all_metric_new_2 = []
                last_metric_dict = {}

                for walker, results in scores.items():
                    last_metric_dict[walker] = {}
                    for metric_idx, (metric_name, calculatedMetrics_) in enumerate(scores[walker].items()):
                        if metric_idx % 2 == 0:
                            all_metric_new_2.append(calculatedMetrics_['allMetricValues'])
                        else:
                            all_metric_new_1.append(calculatedMetrics_['allMetricValues'])
                        last_metric_dict[walker][metric_name] = calculatedMetrics_['lastValue']
                average_last_values_metric1 = (np.average(all_metric_new_1))
                average_last_values_metric2 = (np.average(all_metric_new_2))

                scores_ = {}
                for walker, last_metrics_names in last_values.items():
                    scores_[walker] = {}
                    for metric_idx, metric_name in enumerate(last_metrics_names.items()):
                        avg = average_last_values_metric1 if metric_idx % 2 == 0 else average_last_values_metric2
                        constant = 1 if self.initialParameters[f'Transition_{metric_idx + 1}'] == 'positive' else -1
                        scores_[walker][metric_name[0]] = (constant * ((metric_name[1] - avg) * (100 / avg)))

                final_result = {}

                if len(last_values) == len(scores_):
                    for walker, metric_names_ in scores_.items():
                        final_result[walker] = 0
                        for metric_k, scores__ in enumerate(metric_names_.items()):
                            final_result[walker] += scores__[1]

                best_walker = max(final_result, key=lambda k: final_result[k])

                return best_walker, round(final_result[best_walker], 3), last_values[best_walker]
        except FileExistsError:
            Logger.LogToFile('a', self.trajCount,
                             "Not all MDs produced the same amount of frames. Please check your MD results and restart.")
            Logger.LogToFile('ad', self.trajCount, traceback.format_exc())
