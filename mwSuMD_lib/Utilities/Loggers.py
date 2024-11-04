import csv
import os


class Logger:
    def __init__(self, root):
        self.root = root
        os.makedirs(f'{self.root}/reports', exist_ok=True)

    def logData(self, CV, walker, metricName, data, mean_of_data, last_data, scoreMetric):
        cycle = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
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

    @staticmethod
    def logSlope(data):
        # maxDist = max((dist, value) for dist, value in data.items())
        mD = max(dist for dist in data.values())
        with open('slope_logs', 'w') as slopeLog:
            for frames in range(1, len(data)):
                slopeLog.write(str(frames) + "\t" + str(mD))

    @staticmethod
    def PrintSettingsToFile(mode, cycle, settings_df):
        if mode == "w":
            with open('settings.txt', 'w') as settingsFile:
                settingsFile.write(f"Settings at cycle: {cycle}.\n")
                settingsFile.write(settings_df)
                settingsFile.write('\n')
        if mode == "a":
            with open('settings.txt', 'a') as settingsFile:
                settingsFile.write(f"\nSettings at cycle: {cycle}.\n")
                settingsFile.write(settings_df)
                settingsFile.write('\n')

    @staticmethod
    def LogToFile(mode, cycle, message):
        if mode == "w":
            with open('LogFile.txt', 'w') as settingsFile:
                settingsFile.write(f"This cycle started at step #{cycle}.\n")
                settingsFile.write(message + "\n")
        if mode == "a":
            with open('LogFile.txt', 'a') as settingsFile:
                settingsFile.write(f"\nSettings at cycle: {cycle}.\n")
                settingsFile.write(message + "\n")
        if mode == "ad":
            with open('LogFile.txt', 'a') as settingsFile:
                settingsFile.write(message + "\n")

    @staticmethod
    def countTraj_logTraj(initialParameters, selection_list):
        """ At what cycle number mwSuMD was stopped? """
        trajCount = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
        if trajCount == 0:
            with open('walkerSummary.log', 'w') as logF:
                logF.write('#' * 5 + " Simulation Starts " + '#' * 5 + "\n")
        else:
            with open('walkerSummary.log', 'a') as logF:
                logF.write(str(trajCount) + " |Checkpoint| Metric 1: " + (str(initialParameters["Metric_1"])) + " Metric 2: " + (str(initialParameters["Metric_2"]) + f" |Selections: {selection_list} |\n"))
