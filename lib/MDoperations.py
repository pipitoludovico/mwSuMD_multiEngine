import os
import MDAnalysis

from .mwParser import mwInputParser
from .MDsettings import MDsetter
from .FolderOps import FolderOps
from .Runners import Runner
from . metrics import Metrics


class MDoperations(mwInputParser):
    def __init__(self):
        super(MDoperations, self).__init__()
        self.folderOps = FolderOps()
        # self.batches = self.folderOps.createBatches(self.par['Walkers'], FolderOps.getGPUs())

    def setAndRun(self):
        # run acemd3
        self.runMD(self.trajCount, self.par['Walkers'])
        self.new_value, self.walker_metrics, self.slope = self.saveSteps_and_Logs_and_value()
        return self.new_value, self.walker_metrics, self.slope

    def runMD(self, trajCount, number_of_walkers):
        os.chdir(self.folder)
        self.folderOps.checkCycle(trajCount)
        for walk_number in range(1, int(number_of_walkers) + 1):
            os.makedirs('tmp/walker_' + str(walk_number), exist_ok=True)
            MDsetter().createInput(walk_number, trajCount)
            if self.par['PLUMED'] is not None:
                try:
                    os.system('cp restarts/HILLS .')
                    os.system('cp restarts/grid.dat .')
                    os.system('cp restarts/COLVAR .')
                except:
                    continue
            # lancia contemporaneamente il batch di folder 1 2 3 4,5 etc etc, finché
            # non esaurisci il NUMERO DI WALKERS usando i batches creati nel misc.
            # os.chdir('tmp/walker_' + str(walk_number))

            # qui andrà multiprocessing
        gpus = self.folderOps.getGPUs()
        for batch in FolderOps().createBatches(self.par['Walkers'], gpus):
            print(batch)
            if self.par['MDEngine'] == 'ACEMD':
                Runner().runACEMD(self.par['Walkers'], trajCount, batch)
            if self.par['MDEngine'] == 'GROMACS':
                print('test')
                # Runner.runGROMACS(trajCount, batch)

    def saveSteps_and_Logs_and_value(self):
        # walkers_metrics_lists = Metrics().getChosenMetrics(self.selection_list, self.trajCount)
        # list best frame, list rmsd per walker folder...
        # best_walker, new_value, self.slope = Metrics().getBestWalker(walkers_metrics_lists)  # idx+1 value
        # chosenFrame = [walkers_metrics_lists[best_walker - 1][0]]
        # MDoperations().saveStep(best_walker, self.trajCount, self.folder, chosenFrame)
        self.new_value, self.walkers_metrics_lists, self.slope = 2, [6.5], -55
        return self.new_value, self.walkers_metrics_lists, self.slope

    def saveStep(self, best_walker, trajCount, mothfolder, bestFrame):
        """Handle the restart files and the xtc storage"""
        for r in range(1, int(self.par['Walkers']) + 1):
            if r == best_walker:
                print("Renaming the restart files and moving the best walker's restart bin files")
                os.chdir('tmp/walker_%s' % str(r))

                os.system('rm *.coor')  # removing last coordinate
                universe = MDAnalysis.Universe(f'{mothfolder}/system/{self.par["PDB"]}',
                                               f'{self.par["Output"]}_{str(trajCount)}_wrapped.xtc')
                selectAll = universe.select_atoms('all')
                extensions = ('.coor', '.xtc')  # Creating xtc and coor using our chosen frame:
                for extension in extensions:
                    with MDAnalysis.Writer(f'best_frame_{str(trajCount)}{extension}', selectAll) as WRITER:
                        for idx, ts in enumerate(universe.trajectory):
                            if idx == bestFrame[0]:
                                WRITER.write(selectAll)

                # moving the best frame to the trajectory folder
                os.system(f'cp best_frame_{str(trajCount)}.xtc {mothfolder}/trajectories')
                # moving and renaming the binary files to the restart folder
                if self.par['MDEngine'] == 'ACEMD':
                    os.system(f'mv *.coor {mothfolder}/restarts/previous.coor')
                    os.system(f'mv *.xsc {mothfolder}/restarts/previous.xsc')
                    os.system(f'mv *.vel {mothfolder}/restarts/previous.vel')

                if self.par['MDEngine'] == 'GROMACS':
                    os.system(f'mv *.gro {mothfolder}/restarts/previous.coor')
                    os.system(f'mv *.cpt {mothfolder}/restarts/previous.cpt')

                if self.par['PLUMED'] is not None:
                    os.system('cp HILLS  %s/restarts/ ' % mothfolder)
                    os.system('cp COLVAR  %s/restarts/ ' % mothfolder)
                    os.system('cp grid.dat  %s/restarts/ ' % mothfolder)
                print(os.getcwd())
                print("FINISHED SAVING FRAMES")
        os.chdir(self.folder)
        os.system('rm -r tmp')
