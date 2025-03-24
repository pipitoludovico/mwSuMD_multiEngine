import os.path

import numpy as np
import pkg_resources
from subprocess import Popen, DEVNULL
from openmm import *
from openmm import XmlSerializer, Platform, LangevinMiddleIntegrator, CustomExternalForce
import openmm.app as app
from openmm.unit import *
import MDAnalysis as Mda

from mwSuMD_lib.Utilities.Loggers import Logger


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
            self.parameterFolderPath = os.path.abspath('../parameters')

    def runOPENMM(self, walker_folder, gpu, folder_path):
        coord = None
        pdbREF: str = self.initialParameters['REFERENCE']
        psfREF: str = self.initialParameters['PSF']
        refUni = Mda.Universe(psfREF, pdbREF)
        ligand = refUni.select_atoms('resname UNL LIG UNK')
        initialDistance = ligand.center_of_geometry()

        def ApplyRestraints(system_, coords, restraintIndexesLocal):
            """ APPLIES RESTRAINTS TO INDEXES. MUST BE RUN BEFORE BUILDING THE SYSTEM"""
            restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
            # unit.kilocalories_per_mole/unit.nanometers**2
            restraint.addGlobalParameter('k', 5.0 * kilocalories_per_mole / nanometers ** 2)
            restraint.addPerParticleParameter('x0')
            restraint.addPerParticleParameter('y0')
            restraint.addPerParticleParameter('z0')
            for atom in restraintIndexesLocal:
                restraint.addParticle(atom, coords.positions[atom].value_in_unit(nanometers))
            system_.addForce(restraint)

        def GetRestraintIndex(protRestraints, membRestraints, ligandNames, FF, **kwargs) -> list:
            PDBPATH_ = ""
            if FF == "CHARMM":
                PDBPATH_ = self.initialParameters['PDB']
            if FF == "GROMOS":
                PDBPATH_ = self.initialParameters['GRO']
            if FF == "AMBER":
                try:
                    tentativePath = self.initialParameters['INPCRD']
                    if os.path.exists(tentativePath):
                        PDBPATH_ = tentativePath
                    else:
                        PDBPATH_ = self.initialParameters['PDB']
                        if os.path.exists(PDBPATH_):
                            PDBPATH_ = tentativePath
                except:
                    raise FileNotFoundError(
                        'Could not find any INPCRD or PDB FILE FOR THE COORDINATES WHILE DETECTING AMBER FF')

            proteinRes = []
            u = Mda.Universe(PDBPATH_)
            if protRestraints:
                for protRest in protRestraints:
                    selection = u.select_atoms(protRest)
                    proteinRes += list(selection.indices)

            membraneResnames = ["POPC", "POPE", "POPI", "CHL1", "SOPE", "POPS", "SSM"]

            membraneRes = []
            if membRestraints:
                for membSel in membraneResnames:
                    selection = u.select_atoms(f"resname {membSel}")
                    membraneRes += list(selection.indices)

            ligandRes = []
            if ligandNames:
                for ligname in ligandNames:
                    selection = u.select_atoms(f"resname {ligname}")
                    ligandRes += list(selection.indices)
            if kwargs['walker'] == 1 and kwargs['cycle'] == 0:
                if len(ligandRes) > 0:
                    print("Number of restrained ligand atoms: ", len(ligandRes))
                if len(membraneRes) > 0:
                    print("Number of restrained membrane atoms: ", len(proteinRes))
                if len(proteinRes) > 0:
                    print("Number of restrained protein atoms: ", len(proteinRes))

            restraintIndexes = list(
                set([atomIndex for component in (proteinRes, membraneRes, ligandRes) if component is not None for
                     atomIndex in component]))

            return restraintIndexes

        try:
            self.timeWindow = self.initialParameters['Timewindow']
            number_of_steps = int(self.initialParameters['Timewindow'] / (self.initialParameters['Timestep'] / 1000))
            saveFreq = (int(self.initialParameters['Savefreq'] / (self.initialParameters['Timestep'] / 1000)))
            ts = float(self.initialParameters.get('Timestep') / 1000)
            # cutoff di 30
            if self.initialParameters['Temperature'] is None:
                if self.initialParameters['ACTUAL_DISTANCE'] is None:
                    self.initialParameters['ACTUAL_DISTANCE'] = initialDistance

                a, b = initialDistance, np.array([self.initialParameters['ACTUAL_DISTANCE']])
                dist = np.linalg.norm(a - b)
                n = 610 - (dist * 30)

                temperature = max(min(n, n), 310)
                if walker_folder ==1:
                    print("Temperature: ", temperature, "at a distance of  ", dist, "A")
            else:
                temperature = self.initialParameters.get('Temperature', 310)

            integrator = LangevinMiddleIntegrator(temperature * kelvin, 1 / picosecond, ts * picoseconds)
            platform = Platform.getPlatformByName('CUDA')
            if not self.initialParameters['NOGPU']:
                properties = {'DeviceIndex': str(gpu), 'Precision': 'mixed'}
            else:
                properties = {'Precision': 'mixed'}

            if self.initialParameters.get('Relax') is True:
                self.initialParameters['Timewindow'] = int(self.initialParameters['RelaxTime'] * 1000)
                number_of_steps = int(
                    (self.initialParameters['RelaxTime'] * 1000) / (self.initialParameters['Timestep'] / 1000))
                saveFreq = 5000

            os.makedirs(folder_path, exist_ok=True)
            os.chdir(folder_path)
            p_top, charmm, params = None, None, None

            if self.initialParameters['Forcefield'] == 'GROMOS':
                coord = app.GromacsGroFile(f"{self.initialParameters['GRO']}")
                p_top = app.GromacsTopFile(f"{self.initialParameters['TOP']}",
                                           periodicBoxVectors=gro.getPeriodicBoxVectors(),
                                           includeDir='/usr/local/gromacs/share/gromacs/top')
            if self.initialParameters['Forcefield'] == 'CHARMM':
                charmm = True
                PDBPATH = self.initialParameters['PDB']
                coord = app.PDBFile(PDBPATH)
                PSFPATH = self.initialParameters['PSF']
                p_top = app.CharmmPsfFile(PSFPATH)
                pos = coord.positions.value_in_unit(nanometers)
                boxLength = np.max(pos, axis=0) - np.min(pos, axis=0)
                x, y, z = boxLength[0], boxLength[1], boxLength[2]
                p_top.setBox(x * nanometers, y * nanometers, z * nanometers)
                defaultParams = sorted(list(self.initialParameters['Parameters']))
                try:
                    params = app.CharmmParameterSet(*defaultParams)
                except Exception as e:
                    print(repr(e))
                    Logger.LogToFile('a', self.trajCount, repr(e))

            if self.initialParameters['Forcefield'] == 'AMBER':
                p_top = app.AmberPrmtopFile(self.initialParameters['PRMTOP'])

            restraintIdxs = GetRestraintIndex(self.initialParameters['PROTEIN_RESTRAINTS'],
                                              self.initialParameters['MEMBRANE_RESTRAINTS'],
                                              self.initialParameters['LIGAND_RESNAME'],
                                              self.initialParameters['Forcefield'], walker=walker_folder, cycle=self.trajCount)
            if charmm:
                system = p_top.createSystem(params, nonbondedMethod=app.PME, nonbondedCutoff=0.9 * nanometer,
                                            switchDistance=0.75 * nanometer, constraints=app.HBonds, rigidWater=True,
                                            hydrogenMass=4 * amu)
            else:
                system = p_top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=0.9 * nanometer,
                                            switchDistance=0.75 * nanometer, constraints=app.HBonds, rigidWater=True,
                                            hydrogenMass=4 * amu)
            ApplyRestraints(system, coords=coord, restraintIndexesLocal=restraintIdxs)

            def CheckPlumed():
                plumedCopy = ""
                if self.initialParameters.get("PLUMED"):
                    for filename in os.listdir("../../restarts"):
                        if '.' not in filename:
                            fullname = os.path.join("../../restarts", filename)
                            plumedCopy += f'cp {fullname} .;'
                    if any(plumedFile.endswith(".dat") for plumedFile in os.listdir("../../restarts")):
                        plumedCopy += "cp ../../restarts/*.dat .;"
                    Popen(plumedCopy, shell=True, stdout=DEVNULL)
                    try:
                        from openmmplumed import PlumedForce
                        plumedPath = self.initialParameters.get("PLUMED")
                        script = ""
                        with open(plumedPath, 'r') as pd:
                            for line in pd.readlines():
                                script += line
                        system.addForce(PlumedForce(script))
                    except Exception as excep:
                        Logger.LogToFile('a', self.trajCount, repr(excep))
                        raise ModuleNotFoundError

            if self.initialParameters.get("PLUMED") and self.trajCount != 0:
                for file in os.listdir('../../restarts'):
                    if "." not in file:
                        os.system(f"cp ../../restarts/{file} .")

            CheckPlumed()
            sim = app.Simulation(p_top.topology, system, integrator, platform, properties)
            _PATH = 'system' if self.trajCount == 0 else 'restarts'
            xmlPATH = [os.path.abspath(os.path.join(self.initialParameters['Root'], _PATH, file)) for file in
                       os.listdir(os.path.join(self.initialParameters['Root'], _PATH)) if file.endswith(".xml")][0]
            chkPATH = [os.path.abspath(os.path.join(self.initialParameters['Root'], _PATH, file)) for file in
                       os.listdir(os.path.join(self.initialParameters['Root'], _PATH)) if file.endswith(".chk")][0]
            try:
                if xmlPATH:
                    try:
                        sim.loadState(xmlPATH)
                    except:
                        if walker_folder == 1 and self.trajCount == 0:
                            print("A new NVT system was built with default parameters.")
            except Exception as e:
                print(repr(e))
                Logger.LogToFile('a', self.trajCount, repr(e))
            sim.loadCheckpoint(chkPATH)
            sim.context.reinitialize(True)
            sim.context.setStepCount(0)
            sim.reporters.append(
                app.StateDataReporter(f"{self.initialParameters['Root']}/{folder_path}/openMM_{walker_folder}.log",
                                      saveFreq, step=True, totalSteps=number_of_steps, remainingTime=True,
                                      potentialEnergy=True, speed=True, temperature=True))
            try:
                sim.reporters.append(app.XTCReporter(
                    f"{self.initialParameters['Root']}/{folder_path}/{self.initialParameters['Output']}_{walker_folder}.xtc",
                    saveFreq, enforcePeriodicBox=True))
            except Exception as e:
                if self.trajCount == 0 and walker_folder == 1:
                    print(repr(e))
                    print(
                        'Using DCD reported instead of the new XTC reporter. XTC reporter is available with OpenMM>=8.0.0')
                sim.reporters.append(app.DCDReporter(
                    f"{self.initialParameters['Root']}/{folder_path}/{self.initialParameters['Output']}_{walker_folder}.dcd",
                    saveFreq, enforcePeriodicBox=True))
            sim.reporters.append(app.CheckpointReporter(
                f"{self.initialParameters['Root']}/{folder_path}/{self.initialParameters['Output']}_{walker_folder}.chk",
                saveFreq))
            sim.step(number_of_steps)

            final_state = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True,
                                               getParameters=True)
            with open(f'final_state.xml', 'w') as output:
                output.write(XmlSerializer.serialize(final_state))
            sim.reporters.clear()
            self.initialParameters['Timewindow'] = self.timeWindow
            os.chdir(self.initialParameters.get('Root'))
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)
            print(repr(e))
