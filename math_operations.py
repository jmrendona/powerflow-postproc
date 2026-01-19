import numpy as np

class rot_matrix:
    
    def rx(theta):
        theta = np.deg2rad(theta)
        c,s = np.cos(theta), np.sin(theta)
        return(np.array([[1,0,0],[0,c,-s],[0,s,c]],dtype=np.float64))

    def ry(theta):
        theta = np.deg2rad(theta)
        c,s = np.cos(theta), np.sin(theta)
        return(np.array([[c,0,s],[0,1,0],[-s,0,c]],dtype=np.float64))

    def rz(theta):
        theta = np.deg2rad(theta)
        c,s = np.cos(theta), np.sin(theta)
        return(np.array([[c,-s,0],[s,c,0],[0,0,1]],dtype=np.float64))


def rotate_plane(x, y, z, axis='y', angle=0):
    coords = np.vstack([x, y, z])

    if axis.lower() == 'x':
        R = rot_matrix.rx(angle)
    elif axis.lower() == 'y':
        R = rot_matrix.ry(angle)
    elif axis.lower() == 'z':
        R = rot_matrix.rz(angle)
    else:
        raise ValueError("axis must be 'x', 'y', or 'z'")

    rotated = R @ coords
    return rotated[0], rotated[1], rotated[2]
