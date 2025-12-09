from MDAnalysis.analysis import align

u = mda.Universe("system.pdb", "trajectory.dcd")

# align to frame 0 of the trajectory
align.AlignTraj(u, u, select="protein and backbone", in_memory=True).run()

with mda.Writer("aligned.dcd", u.atoms.n_atoms) as W:
    for ts in u.trajectory:
        W.write(u.atoms)
