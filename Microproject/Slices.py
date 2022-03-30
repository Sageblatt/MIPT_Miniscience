import numpy as np
import vtk

cave_pts = np.load('cave_pts.npy')
bar_for_cave_pts = np.load('bar_for_cave_pts.npy')
barrel_pts = np.load('barrel_pts.npy')

unstructuredGrid = vtk.vtkUnstructuredGrid()
points = vtk.vtkPoints()

depth = vtk.vtkDoubleArray()
depth.SetName("Depth")


for i in cave_pts:
    points.InsertNextPoint(i[2], i[0], i[1])
    depth.InsertNextValue(i[1])
for i in barrel_pts:
    points.InsertNextPoint(i[2], i[0], i[1])
    depth.InsertNextValue(i[1])

pts_amount = len(bar_for_cave_pts)

for i in range(pts_amount - 1):
    ids = [i + 1, i, bar_for_cave_pts[i] + pts_amount]
    if bar_for_cave_pts[i] != bar_for_cave_pts[i + 1]:
        if i < 179:
            continue
        elif bar_for_cave_pts[i] != bar_for_cave_pts[i - 179]:
            continue
        ids = [i - 179, i, bar_for_cave_pts[i] + pts_amount]
    tr = vtk.vtkTriangle()
    for j in range(0, 3):
        tr.GetPointIds().SetId(j, ids[j])
    unstructuredGrid.InsertNextCell(tr.GetCellType(), tr.GetPointIds())


if pts_amount > 179:
    if bar_for_cave_pts[pts_amount - 1] == bar_for_cave_pts[pts_amount - 180]:
        ids = [pts_amount - 180, pts_amount - 1, bar_for_cave_pts[pts_amount - 1] + pts_amount]
        tr = vtk.vtkTriangle()
        for j in range(0, 3):
            tr.GetPointIds().SetId(j, ids[j])
        unstructuredGrid.InsertNextCell(tr.GetCellType(), tr.GetPointIds())

unstructuredGrid.SetPoints(points)
unstructuredGrid.GetPointData().AddArray(depth)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetInputDataObject(unstructuredGrid)
writer.SetFileName("Slices.vtu")
writer.Write()
