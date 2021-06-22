import glob,os
from simulation_prep import Simulation

for system_type in ['RL','L']:
  for i in range(24):
    try:
      project_directory = f'../../ROS_{system_type}/RUN{i}'
      if os.path.exists(f'{project_directory}/conf.gro'):
        continue
      print(f'Processing {project_directory}...')
      s = Simulation(project_directory)
      if system_type == 'RL':
        s.prepare_system(receptor_pdb=glob.glob(f'{project_directory}/*.pdb')[0], ligand_mol2=glob.glob(f'{project_directory}/*.mol2')[0])
      if system_type == 'L':
        s.prepare_system(ligand_mol2=glob.glob(f'{project_directory}/*.mol2')[0])
    except Exception as e:
      print(e)
      continue
