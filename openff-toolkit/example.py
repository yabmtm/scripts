from simulation_prep import Simulation

project_directory = './LIG1'
s = Simulation(project_directory)
s.prepare_system(receptor_pdb=f'{project_directory}/receptor.pdb', ligand_mol2=f'{project_directory}/ligand.mol2')
#s.prepare_system(ligand_mol2=f'{project_directory}/ligand.mol2')
