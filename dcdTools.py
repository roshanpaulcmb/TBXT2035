import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDFile
import MDAnalysis.transformations as transformations

import mdtraj as mdt
import os

def setup(pdb = "system.pdb", dcd = 'trajectory.dcd'):
    data = { }
    data["pdb"] = pdb
    data["dcd0"] = dcd
    
    data["t0"] = mdt.load(data["dcd0"], top = data["pdb"])
    return data

def downsample(data, key = "t0",
               start = 0, end = -1, step = 10, 
               outputDcd = "downsampled.dcd"):
    data[key][start:end:step].save(outputDcd)
    return data


# make function be able to take multiple dcds
def stitch(data, dcd1, 
           outputDcd = "stitched.dcd"):
    data["dcd1"] = dcd1
    data["t1"] = mdt.load(data["dcd1"], top = data["pdb"])
    data["stitched"] = outputDcd
    data["tS"] = mdt.join(data["t0"], data["t1"], 
                          discard_overlapping_frames = True)
    
    data["tS"].save(outputDcd)
    return data
    
def align(data, key = "t0"):
    data[key].superpose(data[key], 0)
    return data

def runAll():
    # use whatever functions you need to here
    data = setup("<YOUR PDB>", "<YOUR DCD>")
    align()
    downsample()    
    return None


if __name__ == "__main__":
    runAll()

# def stitch_dcds(data, dcd2, outputDcd = "stitched.dcd"):
#     """
#     Stitches two DCD files together into a single output DCD file.

#     Args:
#         dcd_file1 (str): Path to the first DCD file.
#         dcd_file2 (str): Path to the second DCD file.
#         pdb_file (str): Path to a PDB file containing the topology for both DCDs.
#         output_dcd (str): Path for the output concatenated DCD file.
#     """
#     # 1. Load the first trajectory
#     # The PDB file is required to provide the topology (atom count, bonds, etc.)
#     u1 = mda.Universe(pdb_file, dcd_file1)
    
#     # 2. Load the second trajectory
#     u2 = mda.Universe(pdb_file, dcd_file2)
    
#     # Check if atom counts match, which is necessary for concatenation
#     if u1.atoms.n_atoms != u2.atoms.n_atoms:
#         print("Error: The two DCD files have different numbers of atoms.")
#         return

#     # 3. Open the output file for writing using the MDAnalysis DCDFile writer
#     # The 'w' mode handles writing the header correctly.
#     with DCDFile(output_dcd, 'w', u1.atoms.n_atoms) as w:
#         # 4. Iterate through all frames in the first universe and write them
#         for ts in u1.trajectory:
#             # Write the current timestep (coordinates) to the new DCD file
#             w.write(ts.positions)
        
#         # 5. Iterate through all frames in the second universe and append them
#         for ts in u2.trajectory:
#             w.write(ts.positions)
            

#     u = mda.Universe("system.pdb", "trajectory.dcd")

#     # align to frame 0 of the trajectory
#     align.AlignTraj(u, u, select="protein and backbone", in_memory=True).run()

#     with mda.Writer("aligned.dcd", u.atoms.n_atoms) as W:
#         for ts in u.trajectory:
#             W.write(u.atoms)

#     print(f"Successfully concatenated {dcd_file1} and {dcd_file2} into {output_dcd}")

# # --- Example Usage ---
# # Replace with your actual file names
# dcd_file_1 = 'traj_part1.dcd'
# dcd_file_2 = 'traj_part2.dcd'
# pdb_topology_file = 'system_topology.pdb' # A PDB, PSF, or other topology file
# output_concatenated = 'merged_trajectory.dcd'

# # Create dummy files for testing if they don't exist
# if not os.path.exists(dcd_file_1) or not os.path.exists(dcd_file_2):
#     print("Example files not found. Please replace the variables with your actual file paths.")
# else:
#     stitch_dcds(dcd_file_1, dcd_file_2, pdb_topology_file, output_concatenated)