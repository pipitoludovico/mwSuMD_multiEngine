import os

import GPUtil

from .mwParser import mwInputParser


class FolderOps(mwInputParser):

    def __init__(self, par):
        super(FolderOps, self).__init__()
        print('INIT FOLDEROPS')
        self.par = par
        self.new_value = None
        self.max_value = None

    def checkCycle(self, trajCount):
        outExtensions = ('coor', 'vel', 'xsc')
        if trajCount != 0:
            for file in os.listdir('restarts'):
                for extension in outExtensions:
                    if file.endswith(extension):
                        self.par[extension] = file

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
        if self.par['Restart'] == 'YES':
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

    @staticmethod
    def getGPUs():
        max_memory = 0.8
        GPUs = GPUtil.getGPUs()
        freeMemory = 0
        gpu_ids = []
        for GPU in GPUs:
            if GPU.memoryUtil > max_memory:
                continue
            if GPU.memoryFree >= freeMemory:
                freeMemory = GPU.memoryFree
                gpu_ids.append(GPU.id)
        return gpu_ids

    @staticmethod
    def createBatches(walkers, total_gpu_ids):
        quotient, rest = divmod(walkers, len(total_gpu_ids))
        result = quotient * total_gpu_ids + total_gpu_ids[:rest]
        batches = [result[i:i + len(total_gpu_ids)] for i in range(0, len(result), len(total_gpu_ids))]

        return batches
