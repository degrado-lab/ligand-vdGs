'''
Generate density map of FG center of masses for viewing in pymol.

After generating the .dx file, view in pymol by: 
load $YOURFILE.dx
isosurface hotspots, $YOURFILE, 0.05 # adjust threshold as needed
'''

import os
import prody as pr
from scipy.stats import gaussian_kde
import numpy as np

CG = 'N=C(N)N'
pdbname = '6br5'

com_coords = []
for pdbfile in os.listdir('output'):
    if pdbfile.startswith(CG):
        par = pr.parsePDB(os.path.join('output', pdbfile))
        hetatms = par.hetatm
        assert len(set(hetatms.getResindices())) == 1
        com = pr.calcCenter(hetatms)
        com_coords.append(com)
com_coords = np.array(com_coords)

kde = gaussian_kde(com_coords.T, bw_method=0.2) # try different bandwidths

# Create 3d grid and evaluate kde
# Define grid bounds
xmin, ymin, zmin = com_coords.min(axis=0) - 5
xmax, ymax, zmax = com_coords.max(axis=0) + 5

# Create grid
x, y, z = np.mgrid[xmin:xmax:50j, ymin:ymax:50j, zmin:zmax:50j]
positions = np.vstack([x.ravel(), y.ravel(), z.ravel()])

# Evaluate KDE
density = kde(positions).reshape(x.shape)

def write_dx(filename, grid, values, origin):
    with open(filename, "w") as f:
        nx, ny, nz = grid.shape
        f.write("object 1 class gridpositions counts {} {} {}\n".format(nx, ny, nz))
        f.write("origin {} {} {}\n".format(*origin))
        f.write("delta {} 0 0\n".format((xmax - xmin) / (nx - 1)))
        f.write("delta 0 {} 0\n".format((ymax - ymin) / (ny - 1)))
        f.write("delta 0 0 {}\n".format((zmax - zmin) / (nz - 1)))
        f.write("object 2 class gridconnections counts {} {} {}\n".format(nx, ny, nz))
        f.write("object 3 class array type double rank 0 items {} data follows\n".format(nx * ny * nz))

        flat = values.ravel(order='F')  # PyMOL uses Fortran-style order
        for i in range(0, len(flat), 3):
            f.write(" ".join(f"{v:.6e}" for v in flat[i:i+3]) + "\n")

        f.write("object 'density' class field\n")

write_dx(f"{pdbname}_{CG}.dx", density, density, [xmin, ymin, zmin])
