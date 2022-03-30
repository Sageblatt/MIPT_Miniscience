import numpy as np
import vtk
import open3d as o3d


cave_pts = np.load('cave_pts.npy')

pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(cave_pts)
pcd.estimate_normals()
distances = pcd.compute_nearest_neighbor_distance()
avg_dist = np.mean(distances)
r = 1.6
k = 6

cfs = [r]

for i in range(5):
    cfs.append(r * (1 - i/60))
for i in range(200):
    cfs.append(r * (1 + i**2/600))

bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(
    pcd,
    o3d.utility.DoubleVector(cfs))

o3d.visualization.draw_geometries([bpa_mesh, pcd], mesh_show_back_face=True)
# dec_mesh = bpa_mesh.simplify_quadric_decimation(100000)
#
# dec_mesh.remove_degenerate_triangles()
# dec_mesh.remove_duplicated_triangles()
# dec_mesh.remove_duplicated_vertices()
# dec_mesh.remove_non_manifold_edges()
# o3d.visualization.draw_geometries([dec_mesh])
