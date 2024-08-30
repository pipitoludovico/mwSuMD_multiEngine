import os

from .Parser import mwInputParser
from warnings import filterwarnings

filterwarnings(action='ignore')


class Template(mwInputParser):
    def __init__(self):
        super(mwInputParser, self).__init__()
        self.trajCount = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
        self.root = self.initialParameters['Root']
        if self.trajCount == 0:
            self.sys_folder = "../../system"
        else:
            self.sys_folder = "../../restarts"

        if self.initialParameters['MDEngine'] == 'ACEMD':
            self.inputFile = [
                'restart\toff\n',
                'minimize        0\n',
                'run            %sps\n' % self.initialParameters.get('Timewindow'),
                'timeStep        %s\n' % self.initialParameters.get('Timestep'),
                '%s' % f'parmfile ../../system/{self.initialParameters.get("PRMTOP")}\n' if self.initialParameters[
                                                                                                'Forcefield'] == 'AMBER' else '\n',
                f'structure\t../../system/{self.initialParameters.get("PSF")}\n' if self.initialParameters.get(
                    'PSF') is not None else '\n',
                "plumedFile\t%s\n" % (self.initialParameters.get('PLUMED')) if self.initialParameters.get(
                    'PLUMED') is not None else "\n",
                'coordinates             ../../system/%s\n' % self.initialParameters.get('PDB'),
                'temperature     %s\n' % self.initialParameters.get('Temperature'),
                'PME             on\n',
                'cutoff          9.0\n',
                'switchDistance  7.5\n',
                'thermostat      on\n',
                'thermostatDamping       0.1\n',
                'thermostatTemperature   %s\n' % self.initialParameters.get('Temperature'),
                'barostat                off\n',
                'trajectoryFile          %s_%s.xtc\n' % (self.initialParameters.get('Output'), str(self.trajCount)),
                'trajectoryPeriod               %s\n' % int(
                    (self.initialParameters.get("Savefreq") / (self.initialParameters.get('Timestep')) * 10 ** 3)),
                f'binCoordinates          %s/%s\n' % (self.sys_folder, self.initialParameters.get('coor')),
                f'extendedSystem          %s/%s\n' % (self.sys_folder, self.initialParameters.get('xsc')),
                'binVelocities           %s/%s\n' % (self.sys_folder, self.initialParameters.get('vel'))]

        if self.initialParameters['MDEngine'] == 'NAMD':
            self.inputFile = ['structure               ../../system/%s\n' % self.initialParameters.get('PSF'),
                              'coordinates             ../../system/%s\n' % self.initialParameters.get('PDB'),
                              'outputname              %s\n' % self.initialParameters.get('Output'),
                              'binCoordinates\t%s/%s\n' % (self.sys_folder, self.initialParameters.get("coor")),
                              f'binVelocities\t%s/%s\n' % (self.sys_folder, self.initialParameters.get("vel")),
                              f'extendedSystem\t%s/%s\n' % (self.sys_folder, self.initialParameters.get("xsc")),
                              f'set xscfile [open %s/%s]\n' % (self.sys_folder, self.initialParameters.get("xsc")),
                              "plumed on\n" if self.initialParameters.get("PLUMED") is not None else "\n",
                              "plumedfile\t%s\n" % self.initialParameters.get("PLUMED") if self.initialParameters.get(
                                  'PLUMED') is not None else "\n",
                              'proc get_first_ts { xscfile } {\n',
                              '  set fd [open $xscfile r]\n',
                              '  gets $fd\n',
                              '  gets $fd\n',
                              '  gets $fd line\n',
                              '  set ts [lindex $line 0]\n',
                              '  close $fd\n',
                              '  return $ts\n',
                              '}\n',
                              f'set firsttime [get_first_ts %s{self.initialParameters.get("xsc")}]\n' % self.sys_folder,
                              'firsttimestep\t$firsttime\n',
                              'set temp\t%s;' % self.initialParameters.get('Temperature'),
                              'outputEnergies %s\n' % int(self.initialParameters.get("Savefreq") / (
                                          self.initialParameters.get('Timestep') / 10 ** 3)),
                              'dcdfreq\t%s;\t\n' % int(self.initialParameters.get("Savefreq") / (
                                          self.initialParameters.get('Timestep') / 10 ** 3)),
                              'dcdUnitCell\tyes;\n',
                              'xstFreq\t%s;\t\n' % int(
                                  self.initialParameters["Savefreq"] / (self.initialParameters['Timestep'] / 10 ** 3)),
                              '# Force-Field Parameters\n',
                              "%s;\n" % 'paraTypeCharmm\ton' if self.initialParameters[
                                                                    'Forcefield'] == 'CHARMM' else 'amber\ton\n',
                              '%s;\n' % f'parmfile ../../system/{self.initialParameters["PRMTOP"]}' if
                              self.initialParameters["Forcefield"] == "AMBER" else "\n",
                              'exclude                 scaled1-4\n',
                              '1-4scaling              1.0\n',
                              'switching               on\n',
                              'vdwForceSwitching       on;                 # \n',
                              'cutoff                  10.0;               # may use smaller, maybe 10., with PME\n',
                              'switchdist              8.0;               # cutoff - 2.\n',
                              'pairlistdist            16.0;               # \n',
                              'stepspercycle           5;                  # SET TO 5 (or lower than 20 if HMR)\n',
                              'pairlistsPerCycle       2;                  # 2 is the default\n',
                              'timestep\t%s;                # fs/step SET 4 is you use HMR\n' % self.initialParameters[
                                  'Timestep'],
                              'rigidBonds              all;                # Bound constraint\n',
                              '%s' % 'useSettle\ton' if self.initialParameters['Forcefield'] == 'AMBER' else '',
                              '%s' % 'rigidTolerance 1.0e-8' if self.initialParameters['Forcefield'] == 'AMBER' else '',
                              '%s' % 'zeromomentum on' if self.initialParameters['Forcefield'] == 'AMBER' else '',
                              '%s' % 'ljcorrection on' if self.initialParameters['Forcefield'] == 'AMBER' else '',
                              'nonbondedFreq           1;                  # nonbonded forces every step\n',
                              'fullElectFrequency      1;                  # PME every step\n', '\n',
                              'wrapWater               on;                 # wrap water to central cell\n',
                              'wrapAll                 on;                 # wrap other molecules too\n',
                              '# PME (for full-system periodic electrostatics)\n',
                              'PME                     yes;\n',
                              'PMEInterpOrder          6;\n',
                              'PMEGridSpacing          1.0;\n',
                              '\n',
                              'useGroupPressure        yes;                \n',
                              'useFlexibleCell         yes;                \n',
                              'useConstantRatio        yes;                 # keeps the ratio of the unit cell.\n',
                              '\n', '# Constant Temperature Control\n',
                              'langevin                on;                 \n',
                              'langevinDamping         0.5;                \n',
                              'langevinTemp            $temp;              # random noise at this level\n',
                              "langevinHydrogen        off;                # don't couple bath to hydrogens\n", '\n',
                              '# Constant pressure\n',
                              'langevinPiston          off;                \n',
                              'langevinPistonTarget    1.01325;            \n',
                              'langevinPistonPeriod    50.0;               \n',
                              'langevinPistonDecay     25.0;               \n',
                              'langevinPistonTemp      $temp;              \n',
                              '\n',
                              '\n',
                              'numsteps\t%s;\n' % int(self.initialParameters['Timewindow'] / (self.initialParameters['Timestep'] / 10 ** 3)),
                              'run\t%s;\n' % int(self.initialParameters['Timewindow'] / (self.initialParameters['Timestep'] / 10 ** 3))]

        if self.initialParameters['MDEngine'] == 'GROMACS':
            self.inputFile = ['title                   = %s\n' % self.initialParameters['Output'],
                              '; Run parameters\n',
                              'integrator              = md        ; leap-frog integrator\n',
                              'nsteps                  = %s    ; ts (ps) * ns = Timewindow (ps)\n' %
                              int(self.initialParameters['Timewindow'] / (self.initialParameters['Timestep'] / 1000)),
                              'dt                      = %s     ; Timestep/1000 \n' % float(
                                  self.initialParameters['Timestep'] / 10 ** 3),
                              '; Output control\n',
                              'nstxout                 = 0         ; suppress bulky .trr file by specifying \n',
                              'nstvout                 = 0         ; 0 for output frequency of nstxout,\n',
                              'nstfout                 = 0         ; nstvout, and nstfout\n',
                              'nstenergy               = %s      ; savefrequency\n' % int(
                                  self.initialParameters["Savefreq"] / (self.initialParameters['Timestep'] / 10 ** 3)),
                              'nstlog                  = %s      ; update log file every 10.0 ps\n'
                              % int(
                                  self.initialParameters["Savefreq"] / (self.initialParameters['Timestep'] / 10 ** 3)),
                              'nstxout-compressed      = %s      ; save compressed coordinates every 10.0 ps\n'
                              % int(
                                  self.initialParameters["Savefreq"] / (self.initialParameters['Timestep'] / 10 ** 3)),
                              'compressed-x-grps       = System    ; save the whole system\n', '; Bond parameters\n',
                              'continuation            = no       ; Restarting after NPT \n',
                              'constraint_algorithm    = SHAKE     ; holonomic constraints \n',
                              'constraints             = none   ; bonds involving H are constrained\n',
                              'lincs_iter              = 1         ; accuracy of LINCS\n',
                              'lincs_order             = 4         ; also related to accuracy\n',
                              '; Neighborsearching\n',
                              'cutoff-scheme           = Verlet    ; Buffered neighbor searching\n',
                              'ns_type                 = grid      ; search neighboring grid cells\n',
                              'nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme\n',
                              'rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)\n',
                              'rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)\n',
                              '; Electrostatics\n',
                              'coulombtype             = PME       ; Particle Mesh Ewald for long-range\n',
                              'pme_order               = 4         ; cubic interpolation\n',
                              'fourierspacing          = 0.16      ; grid spacing for FFT\n',
                              '; Temperature coupling is on\n',
                              'tcoupl                  = V-rescale             ; modified Berendsen thermostat\n',
                              'tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate\n',
                              'tau_t                   = 0.1     0.1           ; time constant, in ps\n',
                              'ref_t\t= %s     ; reference temperature, one for each group, in K%s\n' % (self.initialParameters.get('Temperature'), self.initialParameters.get('Temperature')),
                              '; Pressure coupling is on\n',
                              'pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT\n',
                              'pcoupltype              = isotropic             ; uniform scaling of box vectors\n',
                              'tau_p                   = 2.0                   ; time constant, in ps\n',
                              'ref_p                   = 1.0                   ; reference pressure, in bar\n',
                              'compressibility\t= 4.5e-5\t; isothermal compressibility of water, bar^-1\n',
                              '; Periodic boundary conditions\n',
                              'pbc                     = xyz       ; 3-D PBC\n',
                              '; Dispersion correction\n',
                              'DispCorr                = EnerPres  ; account for cut-off vdW scheme\n',
                              '; Velocity generation\n',
                              'gen_vel                 = no        ; Velocity generation is off \n']
