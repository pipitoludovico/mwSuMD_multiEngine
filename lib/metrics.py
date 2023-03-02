from .Getters import *
from .mwParser import *


class Metrics(mwInputParser):
    def __init__(self, par):
        self.par = par
        super(Metrics, self).__init__()
        import warnings
        warnings.filterwarnings(action='ignore')
        self.walkers_metrics = []

    def getChosenMetrics(self, selection_list):
        import glob
        import os
        """Compute single metric and return a list of values per frame"""
        for walkFolderNumber in range(1, int(self.par['Walkers'] + 1)):
            os.chdir(f'tmp/walker_' + str(walkFolderNumber))
            print(f'walker_' + str(walkFolderNumber))
            if (glob.glob('*.xtc') and glob.glob('*.coor')) or (glob.glob('*.xtc') and glob.glob('*.gro')):
                if self.par['Slope'] == 'NO':
                    if self.par['Metric_1'] == 'DISTANCE':
                        distMetric = Getters(self.par).getDistance(selection_list[0], selection_list[1])
                        self.walkers_metrics.append(distMetric)

                    elif self.par['Metric_1'] == 'CONTACTS':
                        contactMetric = Getters(self.par).contacts_misc(selection_list[0], selection_list[1])
                        self.walkers_metrics.append(contactMetric)

                    elif self.par['Metric_1'] == 'RMSD':
                        distMetric = Getters(self.par).getRMSD(selection_list[0], selection_list[1])
                        print(distMetric)
                        self.walkers_metrics.append(distMetric)
                    else:
                        print("NO METRICS FOUND")
                        exit()
                    os.chdir('../')
                os.chdir(self.folder)
        return self.walkers_metrics


    def getBestWalker(self, walkers_metrics):
        index = 0
        value = 0

        if self.par['Transition_1'] == 'positive':  # we want the metric to increase
            value = max([i for i in walkers_metrics if i is not None])
            # max_value = max(walkers_metrics) # so we take the maximum value
            index = walkers_metrics.index(value)

        elif self.par['Transition_1'] == 'negative':  # we want the metric to dencrease
            value = min([i for i in walkers_metrics if i is not None])
            # max_value = min(walkers_metrics) # so we take the minimum value
            index = walkers_metrics.index(value)

        return index + 1, value
