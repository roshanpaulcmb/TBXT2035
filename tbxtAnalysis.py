#### IMPORTS ####

import argparse as arg
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import colormaps
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import mdtraj as mdt
from scipy.ndimage import gaussian_filter
import os
import time
import umap
import hdbscan
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist

def makeU(files):
    print("\nMaking universe..")
    data = { }
    data["u"] = mda.Universe(files["pdb"], files["dcd"])
    data["protein"] = data["u"].select_atoms("protein")
    
    data["t"] = mdt.load(files["dcd"], top = files["pdb"])
    data["t"].superpose(data["t"], 0) # align with 0th frame
    return data

def align(data, verbose = True):
    aligner = align.AlignTraj(data["protein"], 
                              data["protein"],
                              select='backbone',
                              filename='alignedTraj.dcd',
                              in_memory=False).run(verbose = verbose)
    
    data["protein"] = mda.Universe(data["pdb"], 'alignedTraj.dcd')
    
    return data

def calcRmsd(data, verbose = True):
    rFirst = rms.RMSD(data["protein"],
                 data["protein"],
                 select = "name CA",
                 ref_frame = 0).run(verbose = verbose)
    rFirstM = rFirst.results.rmsd.T
    
    rLast = rms.RMSD(data["protein"],
                 data["protein"],
                 select = "name CA",
                 ref_frame = -1).run(verbose = verbose)
    rLastM = rLast.results.rmsd.T
    
    # Visualization
    fig, ax = plt.subplots(figsize=(8,6))

    inferno_cmap = colormaps["inferno"]

    ax.plot(rFirstM[1], rFirstM[2], color=inferno_cmap(0.3), label="First Frame Reference")
    ax.plot(rLastM[1],  rLastM[2],  color=inferno_cmap(0.7), label="Last Frame Reference")

    ax.set_xlabel("Frame")
    ax.set_ylabel("RMSD ($\AA$)")
    ax.set_title("RMSD Over Trajectory")
    ax.legend()
    plt.savefig("rmsd.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    
    return data

def calcRmsf(data, ref_frame = 0, select = "name CA", verbose = True):
    ## Calculate the RMSF for the selection
    data["proteinCA"] = data["protein"].select_atoms("name CA")
    
    r = rms.RMSF(data["proteinCA"]).run(verbose = verbose)
    
    data["rmsfCA"] = r.results.rmsf
    np.save('rmsfCA.npy', data["rmsfCA"]) #alpha carbon only rmsf

    plt.plot(data["proteinCA"].resids, data["rmsf"])
    plt.xlabel('Residue number')
    plt.ylabel('RMSF ($\AA$)')
    plt.savefig("rmsf.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    
    r = rms.RMSF(data["protein"]).run(verbose = verbose)
    data["rmsfProtein"] = r.results.rmsf
    data["protein"].atoms.tempfactors = data["rmsfProtein"]
    np.save('rmsfProtein.npy', data["rmsfProtein"])
    data["protein"].atoms.write('rmsfProtein.pdb') # pdb with tempfactors set for every atom
    
    return data

def runPca(data, verbose = True):
    print("Running PCA..")
    atoms = data["t"].topology.select("name CA")
    coords = data["t"].xyz[:, atoms, :].reshape(len(data["t"]), -1) 
    # flatten to 2d matrix for PCA
    # flattened rows = frames: [x1 y1 z1 x2 y2 z2 ... xn yn zn]
    # columns = features
    
    nFeatures = coords.shape[1]
    coordsCentered = coords - coords.mean(axis=0) # centered matrix for PCA
    cov = (coordsCentered.T @ coordsCentered) / (nFeatures - 1)
    
    eigvals, eigvecs = np.linalg.eigh(cov)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    
    data["pc"] = coordsCentered@eigvecs # [frames x features] @ [features @ features]
    return data

def calcSfe(data, bins = 200, sigma = 2, temperature = 310):
    print("\nCalculating surface free energy..")
    kB = 1.380649E-23 # J/K
    T = 273.15 + temperature # K
    
    states, xedges, yedges = np.histogram2d(data["pc"][:,0],
                                       data["pc"][:,1],
                                       bins=bins)
    prob = states / np.sum(states) # probability of a state (bin)
    F = -kB*T*np.log(prob, where=(prob>0)) # absolute free energy
    deltaF = F - np.nanmin(F) # relative free energy, sets min to 0

    data["deltaFsmooth"] = gaussian_filter(deltaF, sigma = sigma)
    
    # Visualization
    plt.contourf(data["deltaFsmooth"].T, cmap='inferno')
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.colorbar(label="ΔF (kcal/mol)")
    plt.title("Surface Free Energy Plot")
    plt.savefig("sfePCA.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    return data

def runUmap(data, n_neighbors=90, min_dist=0.6, n_components=3):
    print("\nRunning UMAP..")
    start = time.perf_counter()
    
    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        metric='euclidean',
        random_state=42
    )
    
    data["umap"] = reducer.fit_transform(data["pc"][:, :10]) # using first 10 pcs
    
    end = time.perf_counter()
    print(f"Finished finding contacts in {end - start} seconds")
    
    # Visualization
    plt.figure(figsize=(8,6))
    plt.scatter(data["umap"][:,0], data["umap"][:,1], 
                c = np.arange(len(data["umap"])), 
                cmap="inferno", s=3)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.colorbar(label="Frame")
    plt.title("UMAP of PCA Projection")
    plt.savefig(f"umapN{n_neighbors}D{min_dist}.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    
    return data

def clusterUmap(data, min_cluster_size=500):
    print("\nClustering UMAP with no noise..")

    # HDBSCAN clustering
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=1,
        cluster_selection_epsilon=0.01
    )
    labels = clusterer.fit_predict(data["umap"])

    # Assign noise points (-1) to nearest cluster
    noise_idx = np.where(labels == -1)[0]
    if len(noise_idx) > 0:
        unique_labels = np.unique(labels[labels != -1])
        centroids = np.array([data["umap"][labels == lbl].mean(axis=0) for lbl in unique_labels])

        for i in noise_idx:
            nearest = np.argmin(cdist([data["umap"][i]], centroids))
            labels[i] = unique_labels[nearest]

    data["cluster"] = labels

    # Visualization
    plt.figure(figsize=(8,6))
    unique_clusters = np.unique(labels)
    cmap = plt.get_cmap("tab10")

    for idx, cluster_id in enumerate(unique_clusters):
        cluster_points = data["umap"][labels == cluster_id]
        plt.scatter(
            cluster_points[:,0], 
            cluster_points[:,1], 
            color=cmap(idx % 10), 
            s=4, 
            label=f"Cluster {cluster_id}"
        )

    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title("Clustered UMAP")
    plt.legend(markerscale=3, fontsize=10)
    plt.savefig("umapClustered.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()

    return data

### RUN ####

def runAll(pdb = "system.pdb", 
           dcd = "trajectory.dcd"):
    
    files = {
        "pdb": pdb,
        "dcd": dcd
    }
    
    data = makeU(files)
    data = align(data)
    data = calcRmsd(data)
    data = calcRmsf(data)
    data = runPca(data)
    data = calcSfe(data)
    data = runUmap(data)
    data = clusterUmap(data)

    return None

if __name__ == "__main__":
    runAll()