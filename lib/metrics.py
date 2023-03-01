from .Getters import *
from .mwParser import *


class Metrics(mwInputParser):
    def __init__(self):
        super(Metrics).__init__()
        import warnings
        warnings.filterwarnings(action='ignore')
        self.slopeData = {}
        self.walkers_metrics = []

    def getBestWalker(self, walkers_metrics):
        index = 0
        value = 0
        slope = None

        if self.par['Transition_1'] == 'positive':  # we want the metric to increase
            value = max([i[1] for i in walkers_metrics if i is not None])
            index = [val[1] for val in walkers_metrics].index(value)
            slope = walkers_metrics[index]

        elif self.par['Transition_1'] == 'negative':  # we want the metric to decrease
            value = min([i[1] for i in walkers_metrics if i is not None])
            index = [val[1] for val in walkers_metrics].index(value)
            slope = walkers_metrics[index]

        best_walker = index + 1
        print(f"\nWalker {best_walker} scored: {value}, with a slope of: {slope[2]}")
        return best_walker, value, slope[2]

    def getChosenMetrics(self, selection_list, trajCount):
        import glob
        import os
        """Compute single metric and return a list of values per frame"""
        os.chdir(self.folder)
        os.chdir('tmp')
        for walkFolderNumber in range(1, int(self.par['Walkers']) + 1):
            os.chdir(f'walker_' + str(walkFolderNumber))
            if (glob.glob('*.xtc') and glob.glob('*.coor')) or (glob.glob('*.xtc') and glob.glob('*.gro')):
                if self.par['Slope'] == 'NO':
                    if self.par['Metric_1'] == 'Distance':
                        distMetric = Getters().getDistance(self.par, selection_list[0], selection_list[1], trajCount)
                        self.walkers_metrics.append(distMetric)

                if self.par['Metric_1'] == 'Contacts':
                    contactMetric = Getters().contacts_misc(self.par, selection_list[0], selection_list[1], trajCount)
                    self.walkers_metrics.append(contactMetric)

                if self.par['Metric_1'] == 'RMSD':
                    distMetric = Getters().getRMSD(self.par, selection_list[0], selection_list[1], trajCount)
                    self.walkers_metrics.append(distMetric)
                    # bestFrame, rmsdMetric, slope = Getters().getRMSD(self.par, selection_list[0], selection_list[1],
                    #                                                  trajCount)
                    # self.walkers_metrics.append((bestFrame, rmsdMetric, slope))
            else:
                print("NO METRICS FOUND")
                exit()
            os.chdir('../')
        os.chdir(self.folder)
        return self.walkers_metrics
