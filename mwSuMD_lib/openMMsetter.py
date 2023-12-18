import os
import numpy as np
import openmm
import openmm.app as app
from openmm.app import *
from openmm.unit import *
from openmm import *
import pkg_resources

from .Parser import mwInputParser
from .GPUoperations import ProcessManager
from .Loggers import Logger


class RunnerOPENMM(mwInputParser):
    def __init__(self):
        super(mwInputParser, self).__init__()
        self.parameterFolderPath = None
        self.trajCount = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
        self.timeWindow = None
        if 'parameters' not in os.listdir(os.getcwd()):
            package_dir = pkg_resources.resource_filename('mwSuMD_lib', 'parameters')
            self.parameterFolderPath = os.path.abspath(package_dir)
        else:
            self.parameterFolderPath = os.path.abspath('parameters')

    def runOPENMM(self, walker_folder, gpu):
        if self.initialParameters['Relax'] is True:
            manager = ProcessManager()
            GPUs = manager.getGPUids()
            if self.initialParameters.get("EXCLUDED_GPUS"):
                for excluded in self.initialParameters.get("EXCLUDED_GPUS"):
                    GPUs.remove(excluded)
            GPUbatches, idList = manager.createBatches(walkers=self.initialParameters['Walkers'], total_gpu_ids=GPUs)
            strGPU = map(str, idList)
            gpu = ','.join(strGPU)
        Logger.LogToFile("ad", self.trajCount, "Running in " + os.getcwd())

        self.timeWindow = self.initialParameters['Timewindow']
        if self.initialParameters.get('Relax') is True:
            self.initialParameters['Timewindow'] = int(self.initialParameters['RelaxTime'] * 1000)
        number_of_steps = (int(self.initialParameters['Timewindow'] / (self.initialParameters['Timestep'] / 1000)))
        saveFreq = (int(self.initialParameters['Savefreq'] / (self.initialParameters['Timestep'] / 1000)))
        ts = float(self.initialParameters.get('Timestep') / 1000)

        if self.trajCount == 0 and self.initialParameters['Restart'] == "NO":  # if we start from 0 we check in system
            with open('./system/equilibration_checkpnt.xml') as inputChk:
                system = XmlSerializer.deserialize(inputChk.read())
        else:
            with open('./restarts/previous.xml') as inputChk:  # else we look in the restarts
                system = XmlSerializer.deserialize(inputChk.read())

        os.makedirs(f"tmp/walker_{walker_folder}", exist_ok=True)
        if self.initialParameters['Relax'] is True:
            number_of_steps = (
                int((self.initialParameters['RelaxTime'] * 1000) / (self.initialParameters['Timestep'] / 1000)))

        if self.initialParameters['Forcefield'] == 'GROMOS':
            gro = GromacsGroFile(f"{self.initialParameters['Root']}/system/{self.initialParameters['GRO']}")
            top = GromacsTopFile(f"{self.initialParameters['Root']}/system/{self.initialParameters['TOP']}",
                                 periodicBoxVectors=gro.getPeriodicBoxVectors(),
                                 includeDir='/usr/local/gromacs/share/gromacs/top')
            system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9 * nanometer,
                                      switchDistance=0.75 * nanometer,
                                      constraints=HBonds)
            platform = openmm.Platform.getPlatformByName('CUDA')
            properties = {'DeviceIndex': f'{str(gpu)}', 'Precision': 'mixed'}
            integrator = openmm.LangevinIntegrator(310 * kelvin, 1 / picosecond, ts * picoseconds)
            sim = app.Simulation(prmtop.topology, system, integrator, platform, properties)
            sim.context.setPositions(gro.positions)

        if self.initialParameters['Forcefield'] == 'CHARMM':
            pdb = app.PDBFile(f"{self.initialParameters['Root']}/system/{self.initialParameters['PDB']}")
            psf = app.CharmmPsfFile(f"{self.initialParameters['Root']}/system/{self.initialParameters['PSF']}")
            pos = pdb.positions.value_in_unit(unit.nanometers)
            boxLength = np.max(pos, axis=0) - np.min(pos, axis=0)
            x, y, z = boxLength[0], boxLength[1], boxLength[2]
            psf.setBox(x * unit.nanometers, y * unit.nanometers, z * unit.nanometers)
            integrator = openmm.LangevinIntegrator(310 * unit.kelvin, 1 / unit.picosecond, ts * unit.picoseconds)
            platform = openmm.Platform.getPlatformByName('CUDA')
            properties = {'DeviceIndex': str(gpu)}
            sim = app.Simulation(psf.topology, system, integrator, platform, properties)

        else:  # if FF = AMBER
            params = f"{self.initialParameters['Root']}/system/{self.initialParameters['PRMTOP']}"
            prmtop = AmberPrmtopFile(params)
            inpcrd_ = [inpcrd for inpcrd in os.listdir(f"{self.initialParameters['Root']}/system/") if
                       inpcrd.endswith('.inpcrd')][0]
            inpcrdPATH = f"{self.initialParameters['Root']}/system/" + inpcrd_
            inpcrd = AmberInpcrdFile(inpcrdPATH)
            platform = openmm.Platform.getPlatformByName('CUDA')
            properties = {'DeviceIndex': f'{str(gpu)}', 'Precision': 'mixed'}
            integrator = openmm.LangevinIntegrator(310 * kelvin, 1 / picosecond, ts * picoseconds)
            sim = app.Simulation(prmtop.topology, system, integrator, platform, properties)
            if inpcrd.boxVectors is not None:
                sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        if self.trajCount == 0 and self.initialParameters['Restart'] == "NO":  # if we start from 0 we check in system
            sim.loadCheckpoint(f"{self.initialParameters['Root']}/system/equilibration_checkpnt.chk")
        else:
            sim.loadCheckpoint("./restarts/previous.chk")
        sim.reporters.append(app.StateDataReporter(
            f"{self.initialParameters['Root']}/tmp/walker_{walker_folder}/openMM_{walker_folder}.log", saveFreq,
            step=True, potentialEnergy=True, temperature=True))
        sim.reporters.append(app.DCDReporter(
            f"{self.initialParameters['Root']}/tmp/walker_{walker_folder}/{self.initialParameters['Output']}_{walker_folder}.dcd",
            saveFreq, enforcePeriodicBox=True))
        sim.reporters.append(app.CheckpointReporter(
            f"{self.initialParameters['Root']}/tmp/walker_{walker_folder}/{self.initialParameters['Output']}_{walker_folder}.chk",
            saveFreq))
        sim.step(number_of_steps)

        # writing checkpoint and state
        sim.reporters.append(app.CheckpointReporter('restart.chk', saveFreq))
        with open(f'tmp/walker_{walker_folder}/restart.xml', 'w') as output:
            output.write(XmlSerializer.serialize(system))
        sim.reporters.clear()
        self.initialParameters['Timewindow'] = self.timeWindow
