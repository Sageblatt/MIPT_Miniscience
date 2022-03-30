import numpy as np
import vtk


cave_pts = np.load('cave_pts.npy')
barrel_pts = np.load('barrel_pts.npy')

unstructuredGrid = vtk.vtkUnstructuredGrid()
points = vtk.vtkPoints()

depth = vtk.vtkDoubleArray()
depth.SetName("Depth")


for i in cave_pts:
    points.InsertNextPoint(i[2], i[0], i[1])
    depth.InsertNextValue(i[1])

ply = vtk.vtkPolyData()
ply.SetPoints(points)
ply.GetPointData().AddArray(depth)


delny = vtk.vtkDelaunay3D()
delny.SetInputData(ply)
delny.BoundingTriangulationOn()
# delny.SetTolerance(0.01)
delny.SetAlpha(12.1)
delny.BoundingTriangulationOff()
# delny.Update()


unstructuredGrid.SetPoints(points)
unstructuredGrid.GetPointData().AddArray(depth)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetInputConnection(delny.GetOutputPort())
writer.SetFileName("Delaunay.vtu")
writer.Write()

# Delaunay using mpl and scipy
# import matplotlib.pyplot as plt
# from scipy.spatial import Delaunay
# cave_pts = np.array([*cave_pts, *barrel_pts])
# fig = plt.figure(figsize=(4, 4))
# ax = fig.add_subplot(111, projection='3d')
# tri = Delaunay(cave_pts, furthest_site=True, incremental=True)
# ax.plot_trisurf(cave_pts[:, 2], cave_pts[:, 0], cave_pts[:, 1],
#                 triangles=tri.simplices, cmap=plt.cm.Spectral)
# plt.show()
