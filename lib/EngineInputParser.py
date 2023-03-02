from .Parser import mwInputParser


class EngineInputs(mwInputParser):

    def __init__(self, par):
        super(EngineInputs, self).__init__()
        self.par = par
        self.new_value = None
        self.max_value = None

    def getAcemdInputFile(self, savefreq):
        acemdInput = ['restart off\n',
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
                      'trajectoryPeriod               %s\n' % str(savefreq),
                      f'binCoordinates          ../../system/%s\n' % self.par['coor'],
                      f'extendedSystem          ../../system/%s\n' % self.par['xsc'],
                      'binVelocities           ../../system/%s\n' % self.par['vel']]

        # If we decide to restart, we modify the input file to look straight into the restart folder
        if self.par['Restart'] == 'YES' or self.trajCount != 0:
            print(self.trajCount)
            acemdRestartInput = []
            for line in acemdInput:
                if line.startswith('bin') or line.startswith('extendedSystem'):
                    if 'system' in line:
                        acemdRestartInput.append(line.replace('system', 'restarts'))
                else:
                    acemdRestartInput.append(line)
            return acemdRestartInput
        else:
            return acemdInput

    @staticmethod
    def getAMBERinputFile(par=None, savefreq=None, trajCount=None):
        amberInput = ['restart off\n',
                      'minimize        0\n',
                      'run            %sps\n' % par['Timewindow'],
                      'timeStep        %s\n' % par['Timestep'],
                      'parmfile               %s.prmtop\n' % par['Topology'],
                      'coordinates             %s.pdb\n' % par['Topology'],
                      'temperature     310\n',
                      'PME             on\n',
                      'cutoff          9.0\n',
                      'switchDistance  7.5\n',
                      'thermostat      on\n',
                      'thermostatDamping       0.1\n',
                      'thermostatTemperature   310\n',
                      'barostat                off\n',
                      'trajectoryFile          %s_%s.xtc\n' % (par['Output'], str(trajCount)),
                      'trajectoryPeriod               %s\n' % str(savefreq),
                      'binCoordinates          ../../restarts/%s\n' % par['coor'],
                      'extendedSystem          ../../restarts/%s.xsc\n' % par['xsc'],
                      'binVelocities           ../../restarts/%s.vel\n' % par['vel']]
        return amberInput
