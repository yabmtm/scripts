#!/usr/bin/python

from __future__ import print_function
import glob, os, sys, subprocess
import matplotlib
import mdtraj as md
import msmbuilder.utils
import numpy as np
from itertools import combinations
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from msmbuilder.cluster import KCenters, KMedoids
from msmbuilder.decomposition import tICA
from msmbuilder.featurizer import AtomPairsFeaturizer
from msmbuilder.msm import ContinuousTimeMSM, implied_timescales, MarkovStateModel
from sklearn.externals import joblib
from sklearn.pipeline import Pipeline
# import colormaps as cmaps

''' This script will build a markov state model and output a tica landscape and
    a plot of the timescalse vs. lagtimes. It calculates pairwise distances, tica
    components, and IPTs, based off of a sorted list of input trajectories (traj*.xtc)
    and a structure file (structure.gro). '''

# TODO: add interactivity

tica_lagtime = 1 # determine from implied timescales
tica_components = 4
n_clusters = 10
n_timescales = 8
time_step = 0.1
lagtimes = np.array([1,2,3,4,8,10])
tica_metric = 'dihedrals'
cluster_method = 'kcenters'

tica_1, tica_2, tica_3, tica_4 = [], [], [], []
distance_data = []

def calculate_distances():
    print("Calculating distances...")
    traj_files = sorted(glob.glob("traj*xtc"))
    traj = [ md.load(filename, top='structure.gro') for filename in traj_files ]
    indices = [ a.index for a in traj[0].topology.atoms if a.element.symbol != 'NA' and a.element.symbol != 'CL' ]
#        indices = traj[i].topology.select('name==CA')   
    pairs = list(combinations(indices, 2))
    features = AtomPairsFeaturizer(pairs)
    transformed_data = features.fit_transform(traj)
    for i in range(len(transformed_data)):
        np.save('out_' + str(i) + '.npy', transformed_data[i])

def calculate_dihedrals():
    print("Calculating dihedrals...")
    traj_files = sorted(glob.glob("traj*xtc"))
    traj = [ md.load(filename, top='structure.gro') for filename in traj_files ]
    indices = list(traj[0].topology.select('backbone'))
    dihedral_quartets = np.array([indices[i:i+4] for i in range(len(indices)-4)])
    for i in range(len(traj)):
        thetas = md.compute_dihedrals(traj[i], dihedral_quartets)
        cos_sin_dihedrals = np.hstack([np.cos(thetas), np.sin(thetas)])
        np.save('out_' + str(i) + '.npy', thetas)
    
def calculate_tica_components():
    print("Calculating tICA components...")
    in_files = glob.glob("out*npy")
    loaded_files = [ np.load(filename) for filename in in_files ]
    tica = tICA(lag_time=tica_lagtime,
        n_components=int(tica_components)).fit_transform(loaded_files)
    np.save('lag_%d_comp_%d.npy' %(tica_lagtime, tica_components), tica)
    tica_data = 'data_lag_%d_comp_%d' %(tica_lagtime, tica_components)
    joblib.dump(tica, tica_data)
    data = np.load('lag_%d_comp_%d.npy' %(tica_lagtime, tica_components))

    for i in range(len(glob.glob('out*npy'))): # extract the four tICA components
        for j in range(len(data[i])):
            tica_1.append(data[i][j][0])
            tica_2.append(data[i][j][1])
            tica_3.append(data[i][j][2])
            tica_4.append(data[i][j][3])

# Clustering via KCenters
    if cluster_method == 'kcenters':
        print("Clustering via KCenters...")
        clusters = KCenters(n_clusters)
    elif cluster_method == 'kmeans':
        print("Clustering via KMeans...")
        clusters = KMeans(n_clusters)
    else:
        sys.exit("Invalid cluster_method. Use kmeans or kcenters.")
    sequences = clusters.fit_transform(tica)
    np.save('lag_%d_clusters_%d_sequences.npy' %(tica_lagtime, n_clusters), sequences)
    np.save('lag_%d_clusters_%d_center.npy' %(tica_lagtime, n_clusters),
        clusters.cluster_centers_)
    cluster_data = 'lag_%d_clusters_%d.pkl' %(tica_lagtime, n_clusters)
    joblib.dump(sequences, cluster_data)

 # Determining cluster populations
    print("Determining cluster populations...")
    counts = np.array([len(np.where(np.concatenate(sequences)==i)[0]) for i in range(n_clusters)]) # how many frames are in each cluster
    normalized_counts =  counts/float(counts.sum())
    percentages = [ i*100 for i in normalized_counts ]

# Plotting the tICA components
    print("Plotting tICA components with cluster centers...")
    plt.figure(0) # plotting tica_1, tica_2
    plt.hexbin(tica_1, tica_2, bins='log') #, cmap=cmaps.viridis
    x_centers = [clusters.cluster_centers_[i][0] for i in range(len(clusters.cluster_centers_))]
    y_centers = [clusters.cluster_centers_[i][1] for i in range(len(clusters.cluster_centers_))]
    plt.plot(x_centers, y_centers, 'wo')
    for label, x, y in zip(["%.4f"%i for i in percentages], x_centers, y_centers): # adds percentage contribution for each cluster
        plt.annotate(
          label,
          xy = (x, y), xytext = (-20, 20),
          textcoords = 'offset points', ha = 'right', va = 'bottom',
          bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
          arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    plt.savefig('tica_1_2.png')
    plt.figure(1) # plotting tica_1, tica_3
    plt.hexbin(tica_1, tica_3, bins='log')
    x_centers = [clusters.cluster_centers_[i][0] for i in range(len(clusters.cluster_centers_))]
    y_centers = [clusters.cluster_centers_[i][2] for i in range(len(clusters.cluster_centers_))]
    plt.plot(x_centers, y_centers, 'wo')
    for label, x, y in zip([ "%.4f"%i for i in percentages], x_centers, y_centers):
        plt.annotate(
          label,
          xy = (x, y), xytext = (-20, 20),
          textcoords = 'offset points', ha = 'right', va = 'bottom',
          bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
          arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    plt.savefig('tica_1_3.png')
    plt.figure(2) # plotting tica_2, tica_3
    plt.hexbin(tica_2, tica_3, bins='log')
    x_centers = [clusters.cluster_centers_[j][1] for j in range(len(clusters.cluster_centers_))]
    y_centers = [clusters.cluster_centers_[j][2] for j in range(len(clusters.cluster_centers_))]
    plt.plot(x_centers, y_centers, 'wo')
    for label, x, y in zip(["%.4f"%i for i in percentages], x_centers, y_centers):
        plt.annotate(
          label,
          xy = (x, y), xytext = (-20, 20),
          textcoords = 'offset points', ha = 'right', va = 'bottom',
          bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
          arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    plt.savefig('tica_2_3.png')


   # Determining cluster entropy ( this yields errors for me )
    # print("Determining cluster entropy")
    # cluster_entropy = (-1.0*normalized_counts*np.log(normalized_counts)).sum()
    # np.savetxt('cluster_entropy.dat', cluster_entropy)

  
 # Determining the cluster populations and writing out PDBs for cluster centers
    print("Determining cluster populations...")
    counts = np.array([len(np.where(np.concatenate(sequences)==i)[0]) for i in range(n_clusters)]) # how many frames are in each cluster
    normalized_counts =  counts/float(counts.sum())
    np.savetxt('populations.dat', normalized_counts)
    print("Performing cluster analytics and saving center PDBs...\n")
    for i in range(len(glob.glob("traj*xtc"))):
        n_snapshots = len(clusters.distances_[i])
        cluster_indices = np.arange(n_snapshots)[ (clusters.distances_[i] < 1e-6) ] # frames that have centers
        cluster_labels = sequences[i][cluster_indices] # number of cluster
	if cluster_indices.size != 0: # print only the trajectories that have cluster centers
            for j in range(len(cluster_labels)): # for each cluster center found in this trajectory
                print('Cluster center', cluster_labels[j], 'was found in trajectory', str(i) + '.')
                print('It is found on frame', cluster_indices[j], 'and has a relative population of',
                  "%.4f"%percentages[cluster_labels[j]], '%.')

        xtcfile = sorted(glob.glob("traj*xtc"))[i]
        for j in range(len(cluster_indices)): # actually saving the snapshots
            cluster_traj = md.load_frame(xtcfile, cluster_indices[j], top='structure.gro')
            cluster_traj.save_pdb('state_%d.pdb' %cluster_labels[j]+1)


   # Calculating IPTs
    print("\nCalculating Implied Timescales...")
    timescales = implied_timescales(sequences, lagtimes, n_timescales=n_timescales,
        msm=MarkovStateModel(verbose=False))
    
    implied_timescale_data = 'ipt_lag_%d_clusters_%d.pkl' %(tica_lagtime, n_clusters)
    joblib.dump(timescales, implied_timescale_data)
    numpy_timescale_data = 'lag_%d_clusters_%d_timescales.npy' %(tica_lagtime, n_clusters)
    np.savetxt('lagtimes.txt', lagtimes)
    np.save(numpy_timescale_data, timescales)
   
# Plotting IPTs (lagtimes and timescales)
    print("Plotting Implied Timescales...")
    for i in range(n_timescales):
	plt.figure(42)
	plt.plot(lagtimes * time_step, timescales[:, i] * time_step, 'o-')
	plt.yscale('log')
	plt.xlabel('lagtime (ns)')
	plt.ylabel('Implied timescales (ns)')
	plt.savefig('lag_%d_clusters_%d_.png' %(tica_lagtime, n_clusters))

	
if __name__ == '__main__':
#    calculate_distances()
#    calculate_dihedrals()
    calculate_tica_components()
