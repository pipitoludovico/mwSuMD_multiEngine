import multiprocessing as mp
import os
import traceback
from multiprocessing import Manager
import time

from .Getters import *
from .Loggers import Logger
from .Parser import *


class MetricsParser(mwInputParser):
    def __init__(self):
        super(mwInputParser, self).__init__()
        from warnings import filterwarnings
        filterwarnings(action='ignore')
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
        manager = None
        q = None
        if self.initialParameters['Relax'] is True:
            self.initialParameters['Walkers'] = 1
            print("Relax mode: calculating metrics: ", self.initialParameters['Relax'], "Walker set to relax: ",
                  self.initialParameters['Walkers'])
        last_frame_metric_1, last_frame_metric_2, all_m_1, all_m_2 = [], [], [], []

        manager = Manager()
        q = manager.Queue()
        results = []
        try:
            with mp.Pool() as pool:
                for x in range(1, self.initialParameters['Walkers'] + 1):
                    results.append(pool.apply(self.calculateMetricsMP, args=(q, x, self.selection_list)))
                    ret = q.get()
                    if self.initialParameters['NumberCV'] == 1:
                        self.score_metrics.append(*ret[0])
                        self.metric_in_last_frames.append(*ret[1])
                    else:
                        last_frame_metric_1.append(ret[0][0])
                        last_frame_metric_2.append(ret[0][1])
                        all_m_1.append(ret[1][0])
                        all_m_2.append(ret[1][1])
            if self.initialParameters['NumberCV'] == 2:
                self.metric_in_last_frames = last_frame_metric_1, last_frame_metric_2,
                self.walkers_metrics = all_m_1, all_m_2
            self.initialParameters['Walkers'] = self.walkers_number_snapshot
            q.shutdown()
        except:
            print("\nMetric calculation failed. Check if all the simulations ended well.")
            raise Exception(traceback.format_exc())

    def calculateMetricsMP(self, queue, walker, selection_list):
        score_metrics, last_metrics, all_metrics = [], [], []
        os.chdir(f'tmp/walker_{walker}')
        try:
            for productions in os.listdir(os.getcwd()):
                if productions == 'wrapped.xtc':
                    metric_1 = self.initialParameters.get('Metric_1')
                    metric_2 = self.initialParameters.get('Metric_2')
                    number_cv = self.initialParameters.get('NumberCV')

                    if metric_1 and (not metric_2 or number_cv == 1):
                        self.calculateMetric(metric_1, selection_list[:2], score_metrics, last_metrics, all_metrics)

                    if not metric_1 and metric_2 and number_cv == 1:
                        self.calculateMetric(metric_2, selection_list[2:], score_metrics, last_metrics, all_metrics)

                    if metric_1 and metric_2 and number_cv == 2:
                        self.calculateMetric(metric_1, selection_list[:2], score_metrics, last_metrics, all_metrics)
                        self.calculateMetric(metric_2, selection_list[2:], score_metrics, last_metrics, all_metrics)

                    if number_cv == 1:
                        queue.put((score_metrics, last_metrics, all_metrics))
                        Logger(self.initialParameters['Root']).logData(
                            1, walker, metric_1, all_metrics, np.mean(all_metrics), all_metrics[0][-1], score_metrics)

                    if number_cv == 2:
                        queue.put([last_metrics, all_metrics])
                        Logger(self.initialParameters['Root']).logData(1, walker, metric_1, all_metrics,
                                                                       np.mean(all_metrics), last_metrics,
                                                                       score_metrics)
                        Logger(self.initialParameters['Root']).logData(2, walker, metric_2, all_metrics,
                                                                       np.mean(all_metrics), last_metrics,
                                                                       score_metrics)
        except Exception:
            print(
                "Metric collection and calculation failed. Make sure the trajectories produced make sense and were successfully completed")
            print(traceback.format_exc())
            raise
        finally:
            os.chdir(self.folder)

    def calculateMetric(self, metric, selection_list, score_metrics, last_metrics, all_metrics):
        if metric == 'DISTANCE':
            distMetric, distances, lastDist = Getters(self.initialParameters).getDistance(selection_list[0],
                                                                                          selection_list[1])
            score_metrics.append(distMetric)
            last_metrics.append(lastDist)
            all_metrics.append(distances)
        elif metric == 'CONTACTS':
            contactMetric, timeseries, lastContacts = Getters(self.initialParameters).getContacts(selection_list[0],
                                                                                                  selection_list[1])
            score_metrics.append(contactMetric)
            last_metrics.append(lastContacts)
            all_metrics.append(timeseries)
        elif metric == 'RMSD':
            rmsdMetric, data, lastRMSD = Getters(self.initialParameters).getRMSD(selection_list[0], selection_list[1])
            score_metrics.append(rmsdMetric)
            last_metrics.append(lastRMSD)
            all_metrics.append(data)
        elif metric == 'HB':
            HB_score, hbData, lastHB = Getters(self.initialParameters).getHB_score(selection_list[0], selection_list[1])
            score_metrics.append(HB_score)
            last_metrics.append(lastHB)
            all_metrics.append(hbData)

    def getBestWalker(self, *args):
        """Returns the walker with the best metric"""
        try:
            if self.initialParameters['NumberCV'] == 1:
                metric_scores = args[0]
                if self.initialParameters['Transition_1'] == 'positive':
                    best_metric = max(metric_scores)
                else:
                    best_metric = min(metric_scores)
                best_walker_index = metric_scores.index(best_metric) + 1
                return best_walker_index, best_metric, args[1][best_walker_index - 1]
            else:
                average_poll_1 = np.average(args[0])
                average_poll_2 = np.average(args[1])

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
                print("")
                return max_index, max_score, args[2][max_index - 1], args[3][max_index - 1]
        except:
            print("Not all MDs produced the same amount of frames. Please check your MD results and restart.")
            raise Exception(traceback.format_exc())
