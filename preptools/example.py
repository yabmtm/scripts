from simulation_prep import Simulation

project_directory = './LIG1'
s = Simulation(project_directory)
s.prepare_system(receptor_pdb=f'{project_directory}/nsp16.pdb', ligand_mol2=f'{project_directory}/LIG_h.mol2')
#s.prepare_system(ligand_mol2=f'{project_directory}/LIG_h.mol2')
