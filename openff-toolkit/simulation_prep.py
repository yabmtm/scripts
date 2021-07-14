from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from simtk import openmm, unit
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator
from openbabel import openbabel
import parmed
import os, re, sys, subprocess

class Simulation:
    def __init__(self, project_directory):
        self.solvent_padding = 10.0 * unit.angstrom
        self.box_size = openmm.vec3.Vec3(3.4,3.4,3.4)*unit.nanometers
        self.ionic_strength = 100 * unit.millimolar # 100
        self.pressure = 1.0 * unit.atmospheres
        self.collision_rate = 91.0 / unit.picoseconds
        self.temperature = 300.0 * unit.kelvin
        self.timestep = 4.0 * unit.femtoseconds
        self.nsteps_equil = 5000 # test
        self.protein_forcefield = 'amber14/protein.ff14SB.xml'
        self.small_molecule_forcefield = 'openff-1.3.0' # 'gaff-2.11'
        self.solvation_forcefield = 'amber14/tip3p.xml'
        self.write_gromacs_files = True
        self.write_openmm_files = True
        self.hydrogen_mass = 4.0*unit.amu
        self.removeCMMotion = False
        self.rigidWater = False
        self.project_directory = project_directory

    def prepare_system(self,receptor_pdb=None, ligand_mol2=None):

        if receptor_pdb:
            receptor_structure = parmed.load_file(receptor_pdb)

        if ligand_mol2:
            for output in ['pdb','sdf']:
                obConversion = openbabel.OBConversion()
                obConversion.SetInAndOutFormats("mol2", output)
                mol = openbabel.OBMol()
                obConversion.ReadFile(mol, ligand_mol2)
                outMDL = obConversion.WriteFile(mol,f'{self.project_directory}/LIG_temp.{output}')
            ligand_structure = parmed.load_file(f'{self.project_directory}/LIG_temp.pdb')
            ligand = Molecule.from_file(f'{self.project_directory}/LIG_temp.sdf')
            os.remove(f'{self.project_directory}/LIG_temp.sdf')
            os.remove(f'{self.project_directory}/LIG_temp.pdb')
        if receptor_pdb and ligand_mol2:
            structure = ligand_structure + receptor_structure
        elif receptor_pdb:
            structure = receptor_structure
        elif ligand_mol2:
            structure = ligand_structure
        else:
            raise(f'prepare_system() requires receptor_pdb and/or ligand_sdf')

        barostat = openmm.MonteCarloBarostat(self.pressure, self.temperature)    
        parmed_forcefield_kwargs = {'removeCMMotion': self.removeCMMotion, 'ewaldErrorTolerance': 5e-04,
            'nonbondedMethod': app.PME, 'constraints': False, 'rigidWater': self.rigidWater, 'hydrogenMass': self.hydrogen_mass}
        openmm_forcefield_kwargs = {'removeCMMotion': self.removeCMMotion, 'ewaldErrorTolerance': 5e-04,
            'nonbondedMethod': app.PME, 'constraints': True, 'rigidWater': self.rigidWater, 'hydrogenMass': self.hydrogen_mass}

        if ligand_structure:
            parmed_system_generator = SystemGenerator(forcefields=[self.protein_forcefield, self.solvation_forcefield],
                barostat=barostat, periodic_forcefield_kwargs=parmed_forcefield_kwargs, molecules=[ligand],
                small_molecule_forcefield=self.small_molecule_forcefield, cache=f'{self.project_directory}/LIG.json')
            openmm_system_generator = SystemGenerator(forcefields=[self.protein_forcefield, self.solvation_forcefield],
                barostat=barostat, periodic_forcefield_kwargs=openmm_forcefield_kwargs, molecules=[ligand],
                small_molecule_forcefield=self.small_molecule_forcefield, cache=f'{self.project_directory}/LIG.json')
        else:
            parmed_system_generator = SystemGenerator(forcefields=[self.protein_forcefield, self.solvation_forcefield],
                barostat=barostat, periodic_forcefield_kwargs=parmed_forcefield_kwargs)
            openmm_system_generator = SystemGenerator(forcefields=[self.protein_forcefield, self.solvation_forcefield],
                barostat=barostat, periodic_forcefield_kwargs=openmm_forcefield_kwargs)

        modeller = app.Modeller(structure.topology, structure.positions)
        modeller.addSolvent(openmm_system_generator.forcefield, model=re.sub('.*/','',
            re.sub('\..*','',self.solvation_forcefield)),
            padding=self.solvent_padding, ionicStrength=self.ionic_strength)
            # padding=solvent_padding for R/RL, boxSize=box_size for L

        parmed_system = parmed_system_generator.create_system(modeller.topology)
        openmm_system = openmm_system_generator.create_system(modeller.topology)

        integrator = openmm.LangevinIntegrator(self.temperature, self.collision_rate, self.timestep)

        # minimize and equilibrate
        print('Minimizing...')
        platform = openmm.Platform.getPlatformByName('CUDA')
        platform.setPropertyDefaultValue('Precision', 'mixed')
        platform.setPropertyDefaultValue('CudaDeviceIndex', '0')
        try:
            context = openmm.Context(openmm_system, integrator, platform)
        except Exception as e:
            platform = openmm.Platform.getPlatformByName('OpenCL')
            context = openmm.Context(openmm_system, integrator, platform)
        context.setPositions(modeller.positions)

        state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
        with open(f'{self.project_directory}/pre-min-state.xml','w') as f:
            f.write(openmm.XmlSerializer.serialize(state))
        openmm.LocalEnergyMinimizer.minimize(context)

        print(f'Equilibrating using {self.nsteps_equil} timesteps...')
        integrator.step(self.nsteps_equil)

        state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
        parmed_system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
        parmed_system = parmed.openmm.load_topology(modeller.topology,
            parmed_system, xyz=state.getPositions(asNumpy=True))

        if self.write_gromacs_files:
            print('Saving Gromacs files...')
            parmed_system.save(f'{self.project_directory}/conf.gro', overwrite=True)
            parmed_system.save(f'{self.project_directory}/topol.top', overwrite=True)

        if self.write_openmm_files:
            with open(f'{self.project_directory}/integrator.xml', 'w') as f:
                f.write(openmm.XmlSerializer.serialize(integrator))
            with open(f'{self.project_directory}/post-equil-state.xml','w') as f:
                f.write(openmm.XmlSerializer.serialize(state))
            with open(f'{self.project_directory}/system.xml','w') as f:
                f.write(openmm.XmlSerializer.serialize(parmed_system))


