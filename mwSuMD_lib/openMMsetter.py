import os

import numpy as np
from openmm import XmlSerializer, LangevinIntegrator, Platform
import openmm.app as app
from openmm.unit import *
import pkg_resources


class openMMsetter:
    def __init__(self, par):
        self.initialParameters = par
        self.parameterFolderPath = None
        self.trajCount = len([traj for traj in os.listdir('./trajectories') if traj.endswith('.xtc')])
        self.timeWindow = None
        if 'parameters' not in os.listdir(os.getcwd()):
            package_dir = pkg_resources.resource_filename('mwSuMD_lib', 'parameters')
            self.parameterFolderPath = os.path.abspath(package_dir)
        else:
            self.parameterFolderPath = os.path.abspath('parameters')

    def runOPENMM(self, walker_folder, gpu):
        sim = None
        self.timeWindow = self.initialParameters['Timewindow']
        if self.initialParameters.get('Relax') is True:
            self.initialParameters['Timewindow'] = int(self.initialParameters['RelaxTime'] * 1000)
        number_of_steps = (int(self.initialParameters['Timewindow'] / (self.initialParameters['Timestep'] / 1000)))
        saveFreq = (int(self.initialParameters['Savefreq'] / (self.initialParameters['Timestep'] / 1000)))
        ts = float(self.initialParameters.get('Timestep') / 1000)

        if self.trajCount == 0:  # if we start from 0 we check in system
            with open('./system/equilibration_checkpnt.xml') as inputChk:
                system = XmlSerializer.deserialize(inputChk.read())
        else:
            with open('./restarts/previous.xml') as inputChk:  # else we look in the restarts
                system = XmlSerializer.deserialize(inputChk.read())

        os.makedirs(f"tmp/walker_{walker_folder}", exist_ok=True)
        if self.initialParameters['Relax'] is True:
            number_of_steps = int(
                (self.initialParameters['RelaxTime'] * 1000) / (self.initialParameters['Timestep'] / 1000))

        integrator = LangevinIntegrator(310 * kelvin, 1 / picosecond, ts * picoseconds)
        platform = Platform.getPlatformByName('CUDA')
        properties = {'DeviceIndex': str(gpu), 'Precision': 'mixed'}

        if self.initialParameters['Forcefield'] == 'GROMOS':
            gro = app.GromacsGroFile(f"{self.initialParameters['Root']}/system/{self.initialParameters['GRO']}")
            top = app.GromacsTopFile(f"{self.initialParameters['Root']}/system/{self.initialParameters['TOP']}",
                                     periodicBoxVectors=gro.getPeriodicBoxVectors(),
                                     includeDir='/usr/local/gromacs/share/gromacs/top')
            sim = app.Simulation(top.topology, system, integrator, platform, properties)
            sim.context.setPositions(gro.positions)

        if self.initialParameters['Forcefield'] == 'CHARMM':
            pdb = app.PDBFile(f"{self.initialParameters['Root']}/system/{self.initialParameters['PDB']}")
            psf = app.CharmmPsfFile(f"{self.initialParameters['Root']}/system/{self.initialParameters['PSF']}")
            pos = pdb.positions.value_in_unit(nanometers)
            boxLength = np.max(pos, axis=0) - np.min(pos, axis=0)
            x, y, z = boxLength[0], boxLength[1], boxLength[2]
            psf.setBox(x * nanometers, y * nanometers, z * nanometers)
            sim = app.Simulation(psf.topology, system, integrator, platform, properties)

        if self.initialParameters['Forcefield'] == 'AMBER':
            params = f"{self.initialParameters['Root']}/system/{self.initialParameters['PRMTOP']}"
            prmtop = app.AmberPrmtopFile(params)
            inpcrd_ = [inpcrd for inpcrd in os.listdir(f"{self.initialParameters['Root']}/system/") if
                       inpcrd.endswith('.inpcrd')][0]
            inpcrdPATH = f"{self.initialParameters['Root']}/system/" + inpcrd_
            inpcrd = app.AmberInpcrdFile(inpcrdPATH)
            sim = app.Simulation(prmtop.topology, system, integrator, platform, properties)
            if inpcrd.boxVectors is not None:
                sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

        if self.trajCount == 0 and self.initialParameters['Restart'] == "NO":  # if we start from 0 we check in system
            sim.loadCheckpoint(f"{self.initialParameters['Root']}/system/equilibration_checkpnt.chk")
        else:
            sim.loadCheckpoint("./restarts/previous.chk")
        total_steps = int(number_of_steps)

        sim.reporters.append(app.StateDataReporter(
            f"{self.initialParameters['Root']}/tmp/walker_{walker_folder}/openMM_{walker_folder}.log", saveFreq,
            step=True, totalSteps=total_steps, remainingTime=True, potentialEnergy=True, temperature=True))
        sim.reporters.append(app.DCDReporter(
            f"{self.initialParameters['Root']}/tmp/walker_{walker_folder}/{self.initialParameters['Output']}_{walker_folder}.dcd",
            saveFreq, enforcePeriodicBox=True))
        sim.reporters.append(app.CheckpointReporter(
            f"{self.initialParameters['Root']}/tmp/walker_{walker_folder}/{self.initialParameters['Output']}_{walker_folder}.chk",
            saveFreq))
        sim.step(number_of_steps)

        # writing checkpoint and state
        with open(f'tmp/walker_{walker_folder}/restart.xml', 'w') as output:
            output.write(XmlSerializer.serialize(system))
        sim.reporters.clear()
        self.initialParameters['Timewindow'] = self.timeWindow
