#!/usr/bin/env python

from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
from simtk import openmm, unit
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator
import pandas as pd
import parmed
import glob, os, sys, subprocess
from openbabel import openbabel

water_model = 'tip3p'
solvent_padding = 10.0 * unit.angstrom
box_size = openmm.vec3.Vec3(3.7,3.7,3.7)*unit.nanometers
ionic_strength = 100 * unit.millimolar # 100
pressure = 1.0 * unit.atmospheres
collision_rate = 91.0 / unit.picoseconds
temperature = 300.0 * unit.kelvin
timestep = 4.0 * unit.femtoseconds
nsteps_equil = 25000 # test

protein_forcefield = 'amber14/protein.ff14SB.xml'
small_molecule_forcefield = 'openff-1.1.0'
#small_molecule_forcefield = 'gaff-2.11' # only if you really like atomtypes
solvation_forcefield = 'amber14/tip3p.xml'

openmm_write_cutoff = 0 # only write XML files for first n ligands (bc they big)

def prepare_RL_system():
    RL_complex_structure = ligand_structure + receptor_structure

    barostat = openmm.MonteCarloBarostat(pressure, temperature)
    parmed_forcefield_kwargs = {'removeCMMotion': False, 'ewaldErrorTolerance': 5e-04,
        'nonbondedMethod': app.PME, 'constraints': False, 'rigidWater': False, 'hydrogenMass': 3.0*unit.amu}
    parmed_system_generator = SystemGenerator(forcefields=[protein_forcefield,solvation_forcefield],
        barostat=barostat, forcefield_kwargs=parmed_forcefield_kwargs, molecules=[ligand],
        small_molecule_forcefield=small_molecule_forcefield, cache=f'{RL_output_prefix}/LIG{ligand_ndx}.json')
    openmm_forcefield_kwargs = {'removeCMMotion': False, 'ewaldErrorTolerance': 5e-04,
        'nonbondedMethod': app.PME, 'constraints': True, 'rigidWater': True, 'hydrogenMass': 3.0*unit.amu}
    openmm_system_generator = SystemGenerator(forcefields=[protein_forcefield,solvation_forcefield],
        barostat=barostat, forcefield_kwargs=openmm_forcefield_kwargs, molecules=[ligand],
        small_molecule_forcefield=small_molecule_forcefield, cache=f'{RL_output_prefix}/LIG{ligand_ndx}.json')

    modeller = app.Modeller(RL_complex_structure.topology, RL_complex_structure.positions)
    modeller.addSolvent(openmm_system_generator.forcefield, model='tip3p',
        padding=solvent_padding, ionicStrength=ionic_strength) # padding=solvent_padding for RL, boxSize=box_size for L

    parmed_system = parmed_system_generator.create_system(modeller.topology)
    openmm_system = openmm_system_generator.create_system(modeller.topology)

    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)

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
    openmm.LocalEnergyMinimizer.minimize(context)

    print('Equilibrating...')
    integrator.step(nsteps_equil)

    print('Saving RL GMX files...')
    state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
    parmed_system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    parmed_system = parmed.openmm.load_topology(modeller.topology,
        parmed_system, xyz=state.getPositions(asNumpy=True))
    parmed_system.save(f'{RL_output_prefix}/conf.gro', overwrite=True)
    parmed_system.save(f'{RL_output_prefix}/topol.top', overwrite=True)

    if 1: # write serialized xml files too
        with open(f'{RL_output_prefix}/integrator.xml', 'w') as f:
            f.write(openmm.XmlSerializer.serialize(integrator))
        with open(f'{RL_output_prefix}/state.xml','w') as f:
            f.write(openmm.XmlSerializer.serialize(state))
        with open(f'{RL_output_prefix}/system.xml','w') as f:
            f.write(openmm.XmlSerializer.serialize(parmed_system))

for ligand_ndx in range(1,2): # just process the example directory

    print(f'Processing RUN{ligand_ndx}...')

    if 1:
        ligand_mol2 = glob.glob(f'LIG{ligand_ndx}/ligand.mol2')[0]
        receptor_file = glob.glob(f'LIG{ligand_ndx}/receptor.pdb')[0]
        RL_output_prefix = f'LIG{ligand_ndx}'
        for output in ['pdb','sdf']:
            obConversion = openbabel.OBConversion()
            obConversion.SetInAndOutFormats("mol2", output)
            mol = openbabel.OBMol()
            obConversion.ReadFile(mol, ligand_mol2)
            outMDL = obConversion.WriteFile(mol,f'{RL_output_prefix}/LIG_temp.{output}')
        ligand_structure = parmed.load_file(f'{RL_output_prefix}/LIG_temp.pdb')
        ligand = Molecule.from_file(f'{RL_output_prefix}/LIG_temp.sdf')
        ligand_file = glob.glob(f'{RL_output_prefix}/LIG*.sdf')[0]
        ligand_file_pdb = glob.glob(f'{RL_output_prefix}/LIG*pdb')[0]
        ligand = Molecule.from_file(ligand_file)
        receptor = app.PDBFile(receptor_file)
        receptor_structure = parmed.load_file(receptor_file)
        ligand_structure = parmed.load_file(ligand_file_pdb)

    # prepare RL system
        print(f'Processing RL System for RUN{ligand_ndx}')
        prepare_RL_system()
