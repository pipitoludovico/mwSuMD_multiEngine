import csv


class Logger:
    def __init__(self):
        pass

    @staticmethod
    def logRMSD(RMSD_data, RMSD_mean_rmsd, RMSD_last_rmsd, RMSD_distMetric):
        with open('rmsd_crude.log', 'a') as distF:
            distF.write('RMSD: %s, mean: %s, last: %s, score: %s\n' % (
                " ".join(map(str, RMSD_data)), str(RMSD_mean_rmsd), str(RMSD_last_rmsd), str(RMSD_distMetric)))

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

    @staticmethod
    def logContact_misc(timeseries, mean_contacts, last_contacts, distMetric):
        with open('/contacts_crude.log', 'a') as distF:
            distF.write('contacts: %s, mean: %s, last: %s, score: %s\n' % (
                " ".join(map(str, timeseries)), str(mean_contacts), str(last_contacts), str(distMetric)))

    @staticmethod
    def logDistance(logDistances, log_mean_distance, log_last_distance, log_distMetric):
        with open('/distance_crude.log', 'a') as distF:
            distF.write('distances: %s, mean: %s, last: %s, score: %s\n' % (
                " ".join(map(str, logDistances)), str(log_mean_distance),
                str(log_last_distance), str(log_distMetric)))
