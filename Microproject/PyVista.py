import pyvista as pv
import numpy as np

cave_pts = np.load('cave_pts.npy')
points = pv.wrap(cave_pts)
print(cave_pts.size)
surf = points.reconstruct_surface(nbr_sz=3600, sample_spacing=0.8, progress_bar=True)

pl = pv.Plotter(shape=(1, 2))
pl.add_mesh(points)
pl.add_title('Point Cloud of 3D Surface')
pl.subplot(0, 1)
pl.add_mesh(surf, color=True, show_edges=True)
pl.add_title('Reconstructed Surface')
pl.show()
