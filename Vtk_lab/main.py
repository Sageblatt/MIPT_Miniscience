import numpy as np
import gmsh
import vtk
import math
import os


class CalcMesh:
    def __init__(self, nodes_coords, tetrs_points):
        self.nodes = np.array([nodes_coords[0::3], nodes_coords[1::3], nodes_coords[2::3]])

        self.smth = (self.nodes[0, :] + self.nodes[1, :])*0

        self.velocity = np.zeros(shape=(3, int(len(nodes_coords) / 3)), dtype=np.double)

        self.tetrs = np.array([tetrs_points[0::4],tetrs_points[1::4],tetrs_points[2::4],tetrs_points[3::4]])
        self.tetrs -= 1

        tmp_hook_pts = []
        for i in range(len(self.nodes[0])):
            if (96.2 < self.nodes[0, i] < 105) and (78.8 < self.nodes[1, i] < 81) and (3.4 < self.nodes[2, i] < 9):
                tmp_hook_pts.append(i)

        self.hook_pts = np.array(tmp_hook_pts)
        self.init_vals = np.array([0.05, 1])
        self.comp_flag = False

    def move(self, tau):
        self.nodes += self.velocity * tau

    def set_velocity(self, time):
        if time < 6:
            self.velocity = np.zeros(self.velocity.shape)
            coord = 0
            new_vel = 0
            if time < 3:
                coord = 1
                new_vel = 4
            else:
                coord = 2
                new_vel = -1.3
            for i in self.hook_pts:
                self.velocity[coord, i] = new_vel
        elif (6 <= time < 10) and not self.comp_flag:
            self.velocity = np.zeros(self.velocity.shape)
            self.comp_flag = True
            sum_x = 0
            sum_z = 0
            z_min, z_max = self.nodes[2, 0], self.nodes[2, 0]
            for i in self.hook_pts:
                sum_x += self.nodes[0, i]
                sum_z += self.nodes[2, i]
                z_min = min(z_min, self.nodes[2, i])
                z_max = max(z_max, self.nodes[2, i])
            x_0 = sum_x / len(self.hook_pts)
            z_0 = 1.1*sum_z / len(self.hook_pts)+0.53
            k = (self.init_vals[1] - self.init_vals[0])/(z_max - z_min)
            b = self.init_vals[0] - k * z_min
            for i in self.hook_pts:
                vel_ratio = 0.07
                self.velocity[0, i] = (k * self.nodes[2, i] + b - 1) * (self.nodes[0, i] - x_0) / 5
                if self.nodes[2, i] < z_0:
                    self.velocity[2, i] -= 0.1 * (self.nodes[2, i] - z_0)
                if self.nodes[0, i] - x_0 > 0.1 and self.nodes[2, i] < z_0:
                    self.velocity[0, i] -= (3*k * self.nodes[2, i] + b) * vel_ratio
                elif self.nodes[0, i] - x_0 < -0.1 and self.nodes[2, i] < z_0:
                    self.velocity[0, i] += (3*k * self.nodes[2, i] + b) * vel_ratio

        # else:
        #     self.velocity = np.zeros(self.velocity.shape)

    def snapshot(self, snap_number, step):
        unstructuredGrid = vtk.vtkUnstructuredGrid()
        points = vtk.vtkPoints()

        smth = vtk.vtkDoubleArray()
        smth.SetName("Electric field")

        vel = vtk.vtkDoubleArray()
        vel.SetNumberOfComponents(3)
        vel.SetName("Velocity")

        self.move(step)
        self.set_velocity(step * snap_number)
        k = np.array([0.1, 0.5, 0.6])
        x_0, y_0, z_0 = self.nodes[0, self.hook_pts[0]], self.nodes[1, self.hook_pts[0]], self.nodes[
            2, self.hook_pts[0]]
        r = np.array([self.nodes[0, :] - x_0, self.nodes[1, :] - y_0, self.nodes[2, :] - y_0])
        mod_r = (np.sqrt(r[0, :] ** 2 + r[1, :] ** 2 + r[2, :] ** 2) + 1)
        self.smth = 1000 / mod_r * np.exp(-0.11*mod_r) * np.cos(1 * step * snap_number - 5*mod_r)**3

        for i in range(0, len(self.nodes[0])):
            x, y, z = self.nodes[0, i], self.nodes[1, i], self.nodes[2, i]
            points.InsertNextPoint(x, y, z)
            smth.InsertNextValue(self.smth[i])
            vel.InsertNextTuple((self.velocity[0, i], self.velocity[1, i], self.velocity[2, i]))

        unstructuredGrid.SetPoints(points)

        unstructuredGrid.GetPointData().AddArray(smth)
        unstructuredGrid.GetPointData().AddArray(vel)

        for i in range(0, len(self.tetrs[0])):
            tetr = vtk.vtkTetra()
            for j in range(0, 4):
                tetr.GetPointIds().SetId(j, self.tetrs[j, i])
            unstructuredGrid.InsertNextCell(tetr.GetCellType(), tetr.GetPointIds())

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputDataObject(unstructuredGrid)
        writer.SetFileName("building_process_-step-" + str(snap_number) + ".vtu")
        writer.Write()


gmsh.initialize()

try:
    path = os.path.dirname(os.path.abspath(__file__))
    gmsh.merge(os.path.join(path, 'CY.stl'))
except:
    print("Could not load STL mesh: bye!")
    gmsh.finalize()
    exit(-1)


angle = 10
forceParametrizablePatches = False
includeBoundary = True
curveAngle = 180
gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary, forceParametrizablePatches, curveAngle * math.pi / 180.)
gmsh.model.mesh.createGeometry()


s1 = gmsh.model.getEntities(2)
s2 = gmsh.model.geo.addSurfaceLoop([s1[i][1] for i in range(len(s1))])
gmsh.model.geo.addVolume([s2])

gmsh.model.geo.synchronize()

f = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(f, "F", "1")
gmsh.model.mesh.field.setAsBackgroundMesh(f)

gmsh.model.mesh.generate(3)

nodeTags, nodesCoord, parametricCoord = gmsh.model.mesh.getNodes()

GMSH_TETR_CODE = 4
tetrsNodesTags = None
elementTypes, elementTags, elementNodeTags = gmsh.model.mesh.getElements()
for i in range(0, len(elementTypes)):
    if elementTypes[i] != GMSH_TETR_CODE:
        continue
    tetrsNodesTags = elementNodeTags[i]

if tetrsNodesTags is None:
    print("Can not find tetra data. Exiting.")
    gmsh.finalize()
    exit(-2)

print("The model has %d nodes and %d tetrs" % (len(nodeTags), len(tetrsNodesTags) / 4))


for i in range(0, len(nodeTags)):
    assert (i == nodeTags[i] - 1)
assert(len(tetrsNodesTags) % 4 == 0)

mesh = CalcMesh(nodesCoord, tetrsNodesTags)

total_time = 10
steps = 100
for i in range(steps):
    mesh.snapshot(i, total_time/steps)

gmsh.finalize()
