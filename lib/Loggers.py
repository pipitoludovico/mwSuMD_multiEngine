import csv
import os


class Logger:
    def __init__(self, root):
        self.root = root
        os.makedirs('logs', exist_ok=True)

    def logData(self, CV, walker, metricName, data, mean_of_data, last_data, scoreMetric):
        os.makedirs(f'{self.root}/reports', exist_ok=True)
        cycle = len(os.listdir(f'{self.root}/trajectories'))
        with open(f'{self.root}/reports/datalogger_METRIC_{CV}.log', 'a') as logFile:
            logFile.write(f'Cycle number: {cycle} '
                          f'Walker: {walker}, Metric: {metricName} '
                          f'All data per frame: {data} '
                          f'Mean of data: {mean_of_data} '
                          f'Last Data: {last_data} '
                          f'Score Metric: {scoreMetric}\n')
            logFile.close()

    @staticmethod
    def logStep(par, trajCount, best_walker, walkers_metrics, walkers_metrics_1=None, walkers_metrics_2=None):
        """Write log file with minimum info"""
        if walkers_metrics_1 is None:
            walkers_metrics_1 = []
        walkers_metrics_1.append(best_walker)

        with open('walkerSummary.log', 'a') as logF:
            logF.write('#####     Step number %s     ######\n' % str(trajCount))

            if par['NumberCV'] == '1':
                write = csv.writer(logF)
                write.writerow(walkers_metrics)

            elif par['NumberCV'] == '2':
                write = csv.writer(logF)
                write.writerow(walkers_metrics_1)
                write.writerow(walkers_metrics_2)
                write.writerow(walkers_metrics)

        with open('walkerSummary.log', 'a') as logF:
            logF.write('#####     Step number %s     ######\n' % str(trajCount))

            if par['NumberCV'] == '1':
                write = csv.writer(logF)
                write.writerow(walkers_metrics)

            elif par['NumberCV'] == '2':
                write = csv.writer(logF)
                write.writerow(walkers_metrics_1)
                write.writerow(walkers_metrics_2)
                write.writerow(walkers_metrics)

    #
    # @staticmethod
    # def logContact_misc(timeseries, mean_contacts, last_contacts, distMetric):
    #     with open('/contacts_crude.log', 'a') as distF:
    #         distF.write('contacts: %s, mean: %s, last: %s, score: %s\n' % (
    #             " ".join(map(str, timeseries)), str(mean_contacts), str(last_contacts), str(distMetric)))
    #
    # @staticmethod
    # def logDistance(logDistances, log_mean_distance, log_last_distance, log_distMetric):
    #     with open('/distance_crude.log', 'a') as distF:
    #         distF.write('distances: %s, mean: %s, last: %s, score: %s\n' % (
    #             " ".join(map(str, logDistances)), str(log_mean_distance),
    #             str(log_last_distance), str(log_distMetric)))
    #
    # @staticmethod
    # def HbLogger(nContacts, frame, HBscore):
    #     with open('/HB_scores_crude.log', 'a') as scoreF:
    #         scoreF.write(
    #             'number wat mols: %s, number HB: %s, HB_score: %s\n' % (str(nContacts),
    #                                                                     str(frame),
    #                                                                     str(HBscore)))

    @staticmethod
    def logSlope(data):
        print("Loggers data:")
        print(data)
        # maxDist = max((dist, value) for dist, value in data.items())
        mD = max(dist for dist in data.values())
        with open('slope_logs', 'w') as slopeLog:
            for frames in range(1, len(data)):
                slopeLog.write(str(frames) + "\t" + str(mD))
