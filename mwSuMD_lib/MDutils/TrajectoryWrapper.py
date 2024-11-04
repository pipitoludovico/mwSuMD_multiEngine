import os

import MDAnalysis as Mda
from MDAnalysis import transformations
import traceback

from mwSuMD_lib.Parsers.InputfileParser import mwInputParser
from mwSuMD_lib.Utilities.Loggers import Logger
from warnings import filterwarnings

filterwarnings(action='ignore')


class TrajectoryOperator(mwInputParser):
    def __init__(self):
        super(mwInputParser, self).__init__()
        self.setting_error = f"File or setting missing. "
        self.trajCount = len([traj for traj in os.listdir(f'{self.folder}/trajectories') if traj.endswith('.xtc')])
        if self.initialParameters['Forcefield'] == 'CHARMM':
            if self.initialParameters['PSF'] is None:
                Logger.LogToFile('a', self.trajCount, self.setting_error)
                raise FileNotFoundError
            else:
                self.topology = '../../system/%s' % self.initialParameters['PSF']
                self.coordinates = '../../system/%s' % self.initialParameters['PDB']
        if self.initialParameters['Forcefield'] == 'AMBER':
            if self.initialParameters['PRMTOP'] is None:
                Logger.LogToFile('a', self.trajCount, self.setting_error)
                raise FileNotFoundError
            else:
                self.topology = '../../system/%s' % self.initialParameters['PRMTOP']
                self.coordinates = '../../system/%s' % self.initialParameters['INPCRD']

        if self.initialParameters['Forcefield'] == 'GROMOS':
            for new_coords in os.listdir(os.getcwd()):
                if new_coords.endswith('.tpr'):
                    self.topology = new_coords
                else:
                    self.topology = "../../system/%s" % self.initialParameters.get('TPR')
            self.coordinates = '../../system/%s' % self.initialParameters['GRO']

    def wrap(self, folder):
        try:
            os.chdir(f'{self.folder}/tmp/walker_' + str(folder))
            ext = ('xtc', 'dcd')
            trajFile = None
            for trajectory in os.listdir(os.getcwd()):
                if trajectory.endswith(ext):
                    trajFile = trajectory

            if trajFile.endswith('.dcd'):
                conv = Mda.Universe(self.topology, trajFile)
                all_atoms = conv.select_atoms('all')
                with Mda.Writer('converted.xtc') as converter:
                    for ts in conv.trajectory:
                        if ts:
                            converter.write(all_atoms)
                trajFile = 'converted.xtc'

            try:
                u = Mda.Universe(self.topology, trajFile)
            except:
                Logger.LogToFile('a', self.trajCount, "No trajectory or psf found. Check your simulation parameters and make sure production went well")
                raise FileNotFoundError

            selection = u.select_atoms(f"{self.initialParameters['Wrap']}")
            if len(selection.atoms) == 0:
                Logger.LogToFile('a', self.trajCount,
                                 "your wrapping selection selected 0 atoms! using protein and name CA instead...")
                selection = u.select_atoms('protein and name CA')
            ag = u.atoms
            workflow = [transformations.unwrap(selection), transformations.center_in_box(selection), transformations.wrap(ag, compound="fragments")]
            u.trajectory.add_transformations(*workflow)
            with Mda.Writer('wrapped.xtc', ag) as w:
                for ts in u.trajectory:
                    if ts:
                        w.write(ag)
            os.chdir(self.folder)
        except Exception as e:
            with open('wrappingLog.txt', 'a') as wrapLog:
                wrapLog.write(f"Wrapping Exception in {os.getcwd()}. Exception: {e} ")
            Logger.LogToFile('a', self.trajCount, "Wrapping failed: check your MD results in tmp/walker_x folders")
            raise FileNotFoundError

    def wrapVMD(self, folder):
        try:
            ext = ('xtc', 'dcd')
            trajFile = None
            os.chdir(f'{self.folder}/tmp/walker_' + str(folder))
            for trajectory in os.listdir(os.getcwd()):
                if trajectory.endswith(ext):
                    trajFile = trajectory

            unwrapSel = f"{self.initialParameters['Wrap']}"
            topExt = 'psf' if self.topology.endswith('psf') else 'parm7' if self.topology.endswith('prmtop') else "top"
            pdbExt = "pdb" if self.coordinates.endswith('pdb') else 'rst7' if self.coordinates.endswith('inpcrd') else "gro"
            extension = 'xtc' if trajFile.endswith('xtc') else 'dcd'
            txt = ['package require psfgen\n',
                   'package require pbctools',
                   'resetpsf',
                   f'mol new {self.topology} type {topExt}',
                   f'mol addfile {self.coordinates} type {pdbExt}',
                   f'mol addfile {trajFile} type {extension} first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all',
                   f'set sel [atomselect top "{unwrapSel}"]',
                   "animate goto 0",
                   f'pbc join residue -all -sel "{unwrapSel}"',
                   f'pbc wrap -center com -centersel "{unwrapSel}" -compound residue -all',

                   "proc align { rmolid smolid2 seltext } {",
                   "set ref [atomselect $rmolid $seltext frame 0]",
                   "set sel [atomselect $smolid2 $seltext]",
                   "set all [atomselect $smolid2 all]",
                   "set n [molinfo $smolid2 get numframes]",
                   "",
                   "for { set i 1 } { $i < $n } { incr i } {",
                   "$sel frame $i",
                   "$all frame $i",
                   "$all move [measure fit $sel $ref]",
                   "}",
                   "$ref delete",
                   "$all delete",
                   "$sel delete",
                   "return",
                   "}"
                   f'align 0 0 {unwrapSel}'
                   f'pbc unwrap -sel "{unwrapSel}"',
                   "$sel writepsf filtered.psf",
                   "$sel writepsf ../../filtered_topology.psf" if not os.path.exists("../../filtered_topology.psf") else "",
                   "$sel writepdb ../../filtered_topology.pdb" if not os.path.exists("../../filtered_topology.pdb") else "",
                   'animate write dcd wrapped.dcd beg 1 end -1 waitfor all sel $sel',
                   'quit']
            with open("filterTrj.vmd", "w") as vmdscr:
                for line in txt:
                    vmdscr.write(f'{line}\n')
            os.system('vmd -dispdev text -e filterTrj.vmd > vmd.log 2>&1')
            os.system('mdconvert -o wrapped.xtc -f wrapped.dcd > mdconvert.log')
            os.chdir(self.folder)
        except RuntimeError:
            print(traceback.format_exc())