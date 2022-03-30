import numpy as np
import open3d as o3d


cave_pts = np.load('cave_pts.npy')

pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(cave_pts)
pcd.estimate_normals()


mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(
    pcd, depth=14, scale=30, linear_fit=0, n_threads=6)[0]

bbox = pcd.get_axis_aligned_bounding_box()
mesh = mesh.crop(bbox)

o3d.visualization.draw_geometries([mesh, pcd], width=800, height=600, mesh_show_back_face=True)

