import numpy as np
import pkg_resources

from openmm import *
import openmm as mm
from openmm import XmlSerializer, LangevinIntegrator, Platform
import openmm.app as app
from openmm.unit import *

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
        try:
            self.timeWindow = self.initialParameters['Timewindow']
            number_of_steps = int(self.initialParameters['Timewindow'] / (self.initialParameters['Timestep'] / 1000))
            saveFreq = (int(self.initialParameters['Savefreq'] / (self.initialParameters['Timestep'] / 1000)))
            ts = float(self.initialParameters.get('Timestep') / 1000)
            temperature = self.initialParameters.get('Temperature', 310)
            integrator = LangevinIntegrator(temperature * kelvin, 1 / picosecond, ts * picoseconds)
            platform = Platform.getPlatformByName('CUDA')
            if not self.initialParameters['NOGPU']:
                properties = {'DeviceIndex': str(gpu), 'Precision': 'mixed'}
            else:
                properties = {'Precision': 'mixed'}

            if self.initialParameters.get('Relax') is True:
                self.initialParameters['Timewindow'] = int(self.initialParameters['RelaxTime'] * 1000)
                number_of_steps = int((self.initialParameters['RelaxTime'] * 1000) / (self.initialParameters['Timestep'] / 1000))
                saveFreq = 50000

            os.makedirs(folder_path, exist_ok=True)
            os.chdir(folder_path)
            p_top, charmm, params = None, None, None

            if self.initialParameters['Forcefield'] == 'GROMOS':
                gro = app.GromacsGroFile(f"{self.initialParameters['Root']}/system/{self.initialParameters['GRO']}")
                p_top = app.GromacsTopFile(f"{self.initialParameters['Root']}/system/{self.initialParameters['TOP']}",
                                           periodicBoxVectors=gro.getPeriodicBoxVectors(),
                                           includeDir='/usr/local/gromacs/share/gromacs/top')
            if self.initialParameters['Forcefield'] == 'CHARMM':
                charmm = True
                PDBPATH = self.initialParameters['Root'] + f"/system/{self.initialParameters['PDB']}"
                pdb = app.PDBFile(PDBPATH)
                PSFPATH = self.initialParameters['Root'] + f"/system/{self.initialParameters['PSF']}"
                p_top = app.CharmmPsfFile(PSFPATH)
                pos = pdb.positions.value_in_unit(nanometers)
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
                p_top = app.AmberPrmtopFile(
                    f"{self.initialParameters['Root']}/system/{self.initialParameters['PRMTOP']}")
            if charmm:
                system = p_top.createSystem(params, nonbondedMethod=app.PME, nonbondedCutoff=0.9 * nanometer, switchDistance=0.75 * nanometer, constraints=app.HBonds, rigidWater=True, hydrogenMass=4 * amu)
            else:
                system = p_top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=0.9 * nanometer, switchDistance=0.75 * nanometer, constraints=app.HBonds, rigidWater=True, hydrogenMass=4 * amu)

            def CheckPlumed():
                if self.initialParameters.get("PLUMED"):
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
            xmlPATH = [os.path.abspath(os.path.join(self.initialParameters['Root'], _PATH, file)) for file in os.listdir(os.path.join(self.initialParameters['Root'], _PATH)) if file.endswith(".xml")][0]
            chkPATH = [os.path.abspath(os.path.join(self.initialParameters['Root'], _PATH, file)) for file in os.listdir(os.path.join(self.initialParameters['Root'], _PATH)) if file.endswith(".chk")][0]
            try:
                sim.loadState(xmlPATH)
                sim.loadCheckpoint(chkPATH)
            except Exception as e:
                print(repr(e))
                Logger.LogToFile('a', self.trajCount, repr(e))
            # checking parameters
            p = sim.context.getParameters()

            for f in system.getForces():
                if isinstance(f, mm.MonteCarloBarostat):
                    f.setFrequency(0)
            for k in p:
                if k == "k":
                    sim.context.setParameter('k', 0)
                    sim.context.reinitialize(True)
            sim.context.setStepCount(0)
            sim.reporters.append(app.StateDataReporter(f"{self.initialParameters['Root']}/{folder_path}/openMM_{walker_folder}.log", saveFreq, step=True, totalSteps=number_of_steps, remainingTime=True, potentialEnergy=True, speed=True, temperature=True))
            try:
                sim.reporters.append(app.XTCReporter(f"{self.initialParameters['Root']}/{folder_path}/{self.initialParameters['Output']}_{walker_folder}.xtc", saveFreq, enforcePeriodicBox=True))
            except Exception as e:
                if self.trajCount == 0 and walker_folder == 1:
                    print(repr(e))
                    print('Using DCD reported instead of the new XTC reporter. XTC reporter is available with OpenMM>=8.0.0')
                sim.reporters.append(app.DCDReporter(f"{self.initialParameters['Root']}/{folder_path}/{self.initialParameters['Output']}_{walker_folder}.dcd", saveFreq, enforcePeriodicBox=True))
            sim.reporters.append(app.CheckpointReporter(f"{self.initialParameters['Root']}/{folder_path}/{self.initialParameters['Output']}_{walker_folder}.chk", saveFreq))
            sim.step(number_of_steps)

            final_state = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True, getParameters=True)
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
