import os


class MDoperator:
    def __init__(self, par, root, trajCount):
        self.par = par
        self.folder = root
        self.trajCount = trajCount
        self.walkers_metrics_lists = []

    def saveStep(self, best_walker):
        """Handle the restart files and the xtc storage"""
        for r in range(1, int(self.par['Walkers']) + 1):
            if r == best_walker:
                print("Renaming the restart files and moving the best walker's restart bin files")
                os.chdir('tmp/walker_%s' % str(r))
                # moving the best frame to the trajectory folder
                os.system(f'cp wrapped.xtc {self.folder}/trajectories/step_{self.trajCount}.xtc')
                # moving and renaming the binary files to the restart folder
                if self.par['MDEngine'] == 'ACEMD':
                    os.system(f'cp *.coor {self.folder}/restarts/previous.coor')
                    os.system(f'cp *.xsc {self.folder}/restarts/previous.xsc')
                    os.system(f'cp *.vel {self.folder}/restarts/previous.vel')

                if self.par['MDEngine'] == 'GROMACS':
                    os.system(f'mv *.gro {self.folder}/restarts/previous.coor')
                    os.system(f'mv *.cpt {self.folder}/restarts/previous.cpt')

                if self.par['PLUMED'] is not None:
                    os.system('cp HILLS  %s/restarts/ ' % self.folder)
                    os.system('cp COLVAR  %s/restarts/ ' % self.folder)
                    os.system('cp grid.dat  %s/restarts/ ' % self.folder)
            print("FINISHED SAVING FRAMES")
        os.chdir(self.folder)
        os.system('rm -r tmp')
