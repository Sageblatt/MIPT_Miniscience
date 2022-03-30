import numpy as np
import open3d as o3d


cave_pts = np.load('cave_pts.npy')

pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(cave_pts)
pcd.estimate_normals()

alpha = 30

mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)
mesh.compute_vertex_normals()
o3d.visualization.draw_geometries([mesh, pcd], mesh_show_back_face=True)

