import argparse

import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDFile
import MDAnalysis.transformations as transformations
from MDAnalysis.analysis import align

import mdtraj as mdt
import os

#### ARGPARSE ####

parser = argparse.ArgumentParser(description = "Tools to edit DCD")
parser.add_argument("--pdb", required = True, help = "PDB file")
parser.add_argument("--dcd", required = True, help = "DCD file")
parser.add_argument("--outputName", required = True,
                    type = str, default = "system",
                    help="Base string name for ALL files. Do not append .pdb or .dcd")
parser.add_argument("--dcdStitch", help = "Addable DCD file for stitching")

args = parser.parse_args()

#### FUNCTIONS ####


def setup(pdb = args.pdb, dcd = args.dcd, loadingMdt = False):
    print("\nSetting up..")
    data = { }
    data["pdb"] = pdb
    data["dcd0"] = dcd
    
    if loadingMdt:
        data["t0"] = mdt.load(data["dcd0"], top = data["pdb"])
    return data

# Do this first if you need to stitch dcds
# Future: make function be able to take multiple dcds
def stitch(data, dcd1 = args.dcdStitch, 
           outputName = args.outputName + "Stitched.dcd",
           save = False):
    print("\nStitching .dcd files..")
    if dcd1 is None:
        raise FileNotFoundError("Use --dcdStitch to stitch dcd to initial dcd")
    if ".dcd" not in dcd1:
        raise ValueError("Input must contain '.dcd'")
    
    data["dcd1"] = dcd1
    data["t1"] = mdt.load(data["dcd1"], top = data["pdb"])
    data["t0"] = mdt.join(data["t0"], data["t1"], 
                          discard_overlapping_frames = True)
    
    if save:
        data["t0"].save(outputName)
        data["dcd0"] = outputName
    return data

def downsample(data,
               start = 0, end = -1, step = 10,
               outputName = args.outputName + "Downsampled.dcd"):
    print("\nCreating downsized .dcd file..")
    data["t0"][start:end:step].save(outputName)
    return data


# For aligning mdtraj trajectories
def alignTraj(data, key = "t0", refFrame = 0):
    print("\nAligning mdtraj trajectory..")
    data[key].superpose(data[key], refFrame)
    return data

def stripWater(data, 
               outputName = args.outputName + "Stripped",
               select = "backbone",
               verbose = True,
               save = False):
    print("\nStripping water from .pdb and .dcd files..")
    data["u"] = mda.Universe(data["pdb"], data["dcd0"])
    data["protein"] = data["u"].select_atoms("protein")
    
    print("\nAligning mdanalysis universe..")
    aligner = align.AlignTraj(data["protein"], 
                              data["protein"],
                              select = select,
                              outputName = args.outputName + 'Aligned.dcd',
                              in_memory=False).run(verbose = verbose)
    
    # Write over aligned dcd
    data["dcd0"] = args.outputName + 'Aligned.dcd'
    data["protein"] = mda.Universe(data["pdb"], data["dcd0"])
    
    if save:
        data["protein"].write(f'{outputName}.pdb')
        data["protein"].write(f'{outputName}.dcd', frames = "all")
    
    return data

def runAll():
    # use whatever functions you need to here
    data = setup(loadingMdt = True)
    # stitch(data)
    # downsample(data)
    # align(data)
    stripWater(data, save = True)
    
    return None


if __name__ == "__main__":
    runAll()