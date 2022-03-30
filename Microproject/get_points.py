import numpy as np
from numpy import sqrt, pi, abs, sin, cos, tan
from pandas import read_csv


DBG = 0

if DBG:
    name = 'test'
else:
    name = 'main'

data = read_csv(name + ".csv").to_numpy()

data[:, 1] = data[:, 1] * pi / 180
data[:, 2] = data[:, 2] * pi / 180
data[:, 3] = data[:, 3] * pi / 180


def get_barrel_pts(l, alpha, beta, pts, n):
    if tan(alpha):
        a = abs(sin(alpha)) * np.array([cos(beta), -1 / tan(alpha), sin(beta)])
        x = pts[n] + l * a
    else:
        a = np.array([0, 1, 0])
        x = pts[n] - l * a
    return x, a


x0 = np.array([0, 0, 0], dtype=np.double)
barrel_pts = [x0]
[x1, a1] = get_barrel_pts(data[0, 0], data[0, 2], data[0, 3], barrel_pts, 0)
barrel_pts.append(x1)

guiding_vects = [a1]

for n in range(1, len(data)):
    l = data[n, 0] - data[n-1, 0]
    alpha = data[n, 2]
    beta = data[n, 3]
    (x, a) = get_barrel_pts(l, alpha, beta, barrel_pts, n)
    barrel_pts.append(x)
    guiding_vects.append(a)

barrel_pts = np.array(barrel_pts)
guiding_vects = np.array(guiding_vects)

cave_pts = []
bar_for_cave_pts = []


def rotation_matrtix(vect, angle):
    M = np.zeros((3, 3))
    x, y, z = vect
    co = cos(angle)
    si = sin(angle)
    M[0, 0] = co + (1 - co) * x**2
    M[0, 1] = (1 - co) * x * y - si * z
    M[0, 2] = (1 - co) * x * z + si * y
    M[1, 0] = (1 - co) * y * x + si * z
    M[1, 1] = co + (1 - co) * y**2
    M[1, 2] = (1 - co) * y * z - si * x
    M[2, 0] = (1 - co) * z * x - si * y
    M[2, 1] = (1 - co) * z * y + si * x
    M[2, 2] = co + (1 - co) * z**2
    return M


for i in range(0, len(data)):
    beta = data[i, 3]
    p0 = np.array([0, 0, 0], dtype=np.double)
    n = guiding_vects[i]
    n1, n2, n3 = n[0], n[1], n[2]
    if beta > pi:
        p0[2] = sqrt(1 / (1 + (n3 / n1) ** 2))
        p0[0] = -n3 * p0[2] / n1
    else:
        p0[2] = -sqrt(1 / (1 + (n3 / n1) ** 2))
        p0[0] = -n3 * p0[2] / n1
    for j in range(180):
        l = data[i, j + 4]
        if l > 0:
            ang = (2 * j - 90) * pi / 180
            M = rotation_matrtix(n, ang)
            p = M @ p0
            x0 = barrel_pts[i + 1]
            x = x0 + l * p
            cave_pts.append(x)
            bar_for_cave_pts.append(i + 1)

cave_pts = np.array(cave_pts)
bar_for_cave_pts = np.array(bar_for_cave_pts)

np.save('cave_pts', cave_pts)
np.save('bar_for_cave_pts', bar_for_cave_pts)
np.save('barrel_pts', barrel_pts)


# for debugging using mpl
if DBG == 1:
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(barrel_pts[:, 2], barrel_pts[:, 0], barrel_pts[:, 1])
    ax.scatter(cave_pts[:, 2], cave_pts[:, 0], cave_pts[:, 1])
    plt.show()
