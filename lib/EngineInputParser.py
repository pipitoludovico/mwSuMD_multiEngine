import os

from .Parser import mwInputParser


class EngineInputs(mwInputParser):

    def __init__(self, par):
        super(EngineInputs, self).__init__()
        self.par = par
        self.new_value = None
        self.max_value = None
        self.savefreq = int((int(self.par['Savefreq']) * 1000) / int(self.par['Timestep']))
        self.timewindow = self.par['Timewindow']
        self.trajCount = len(os.listdir(f'{self.folder}/trajectories'))
        mwInputParser().getRestartOutput()
        if 'coor' not in self.par:
            self.par['coor'] = ''
        if 'vel' not in self.par:
            self.par['vel'] = ''
        if 'xsc' not in self.par:
            self.par['xsc'] = ''
        if self.par['Forcefield'] == 'CHARMM':
            self.charmmAcemdInput = ['restart off\n',
                                     'minimize        0\n',
                                     'run            %sps\n' % self.par['Timewindow'],
                                     'timeStep        %s\n' % self.par['Timestep'],
                                     'structure               ../../system/%s\n' % self.par['PSF'],
                                     'coordinates             ../../system/%s\n' % self.par['PDB'],
                                     'temperature     310\n',
                                     'PME             on\n',
                                     'cutoff          9.0\n',
                                     'switchDistance  7.5\n',
                                     'thermostat      on\n',
                                     'thermostatDamping       0.1\n',
                                     'thermostatTemperature   310\n',
                                     'barostat                off\n',
                                     'trajectoryFile          %s_%s.xtc\n' % (self.par['Output'], str(self.trajCount)),
                                     'trajectoryPeriod               %s\n' % str(self.savefreq),
                                     f'binCoordinates          ../../system/%s\n' % self.par['coor'],
                                     f'extendedSystem          ../../system/%s\n' % self.par['xsc'],
                                     'binVelocities           ../../system/%s\n' % self.par['vel']]
        if self.par['Forcefield'] == 'AMBER':
            self.amberAcemdInput = ['restart off\n',
                                    'minimize        0\n',
                                    'run            %sps\n' % self.par['Timewindow'],
                                    'timeStep        %s\n' % self.par['Timestep'],
                                    'parmfile               ../../system/%s\n' % self.par['PRMTOP'],
                                    'coordinates             ../../system/%s\n' % self.par['PDB'],
                                    'temperature     310\n',
                                    'PME             on\n',
                                    'cutoff          9.0\n',
                                    'switchDistance  7.5\n',
                                    'thermostat      on\n',
                                    'thermostatDamping       0.1\n',
                                    'thermostatTemperature   310\n',
                                    'barostat                off\n',
                                    'trajectoryFile          %s_%s.xtc\n' % (self.par['Output'], str(self.trajCount)),
                                    'trajectoryPeriod               %s\n' % str(self.savefreq),
                                    'binCoordinates          ../../system/%s\n' % self.par['coor'],
                                    'extendedSystem          ../../system/%s\n' % self.par['xsc'],
                                    'binVelocities           ../../system/%s\n' % self.par['vel']]

    def getInputFile(self):
        if self.par['Relax'] is True:
            self.par['Timewindow'] = self.par['RelaxTime'] * 1000
            print("\nTemporary changing the timewindow for relaxation protocol to: " + str(self.par['Timewindow']))

        restartInput = self.charmmAcemdInput if self.par['Forcefield'] == 'CHARMM' else self.amberAcemdInput
        if self.par['Restart'] == 'YES' or self.trajCount != 0:
            restartInput = [line.replace('system', 'restarts')
                            if ('bin' in line or 'extendedSystem' in line) and 'system' in line else line
                            for line in restartInput]

        self.par['Timewindow'] = self.timewindow
        return restartInput
