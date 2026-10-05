import numpy as np
from pyqtgraph.opengl import MeshData, GLMeshItem
def create_cylinder(start, end, radius=0.1, color=(0.3, 0.3, 0.3, 1), sectors=32):
    """
    Creates a cylinder mesh between two 3D points.
    Maneja correctamente cilindros paralelos al eje Z.
    """
    from numpy.linalg import norm
    from numpy import cross, dot
    from math import acos

    start = np.array(start, dtype=float)
    end = np.array(end, dtype=float)
    
    height = norm(end - start)
    if height < 1e-6:
        return None
    
    direction = (end - start) / height

    # Create the base cylinder (along Z-axis)
    angle = np.linspace(0, 2 * np.pi, sectors, endpoint=False)
    x = radius * np.cos(angle)
    y = radius * np.sin(angle)
    z = np.zeros_like(x)

    vertices = []
    for i in range(sectors):
        vertices.append([x[i], y[i], 0])
        vertices.append([x[i], y[i], height])
    vertices = np.array(vertices)

    # Faces (triangles)
    faces = []
    for i in range(sectors):
        a = 2 * i
        b = (2 * i + 2) % (2 * sectors)
        c = a + 1
        d = b + 1
        faces.append([a, b, c])
        faces.append([c, b, d])
    faces = np.array(faces)

    # Rotate cylinder to match the direction
    z_axis = np.array([0.0, 0.0, 1.0])
    dot_product = dot(z_axis, direction)
    
    # Caso especial: paralelo al eje Z
    if abs(dot_product) > 0.9999:
        if dot_product < 0:  # apunta hacia -Z
            R = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]], dtype=float)
        else:  # apunta hacia +Z, no rotar
            R = np.eye(3)
    else:
        # Caso general
        axis = cross(z_axis, direction)
        axis = axis / norm(axis)
        angle_rad = acos(np.clip(dot_product, -1.0, 1.0))

        K = np.array([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0]
        ])
        R = np.eye(3) + np.sin(angle_rad) * K + (1 - np.cos(angle_rad)) * (K @ K)
    
    vertices = vertices @ R.T

    # Translate to start point
    vertices += start

    mesh = MeshData(vertexes=vertices, faces=faces)
    item = GLMeshItem(meshdata=mesh, smooth=True, color=color, drawFaces=True)
    return item
