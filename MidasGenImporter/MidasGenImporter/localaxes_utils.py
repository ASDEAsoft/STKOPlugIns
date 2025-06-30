'''
This file contains utility functions for handling local axes in the MidasGenImporter plugin for STKO.
'''

from PyMpc import *
from typing import List, Tuple

# This function compute the director vector for an element.
# for a 2-node line element it gives the tangent vector,
# for a 3-node and 4-node shell element it gives the normal vector.
def _director(points:List[Math.vec3]) -> Math.vec3:
    if len(points) == 2:
        dx = (points[1] - points[0]).normalized()
        return dx
    if len(points) == 3:
        dx = (points[1] - points[0]).normalized()
        dy = (points[2] - points[0]).normalized()
        dz = dx.cross(dy).normalized()
        return dz
    elif len(points) == 4:
        R = (points[1] + points[2]) / 2.0
        L = (points[0] + points[3]) / 2.0
        dx = (R - L).normalized()
        B = (points[1] + points[0]) / 2.0
        T = (points[2] + points[3]) / 2.0
        dy = (T - B).normalized()
        dz = dx.cross(dy).normalized()
        return dz
    else:
        raise Exception('Invalid number of points for shell local axes')

# returns the default local axes to be projected on the surface plane
def _surface_default() -> Math.mat3:
    x = Math.vec3(1.0,  0.0, 0.0)
    y = Math.vec3(0.0,  0.0, 1.0)
    z = Math.vec3(0.0, -1.0, 0.0)
    return Math.mat3(x, y, z)

# projects the input matrix on the surface defined by the normal vector dz
# if dz is zero, it returns the input matrix unchanged
def _project_on_surface(m : Math.mat3, dz : Math.vec3) -> Math.mat3:
    if dz.norm() < 1.0e-10:
        return m
    # start with the input Y axis (which by default is global Z)
    dy = m.col(1)
    if abs(dy.dot(dz)) > (1.0 - 1e-10):
        # if dy is aligned with normal vector, use trial dx (which by default is global X)
        dx = m.col(0)
        dy = dz.cross(dx).normalized()
        dx = dy.cross(dz).normalized()
        return Math.mat3(dx, dy, dz)
    else:
        # use trial dy
        dx = dy.cross(dz).normalized()
        dy = dz.cross(dx).normalized()
        return Math.mat3(dx, dy, dz)

class STKOLocalAxesConvention:
    @staticmethod
    def get_local_axes(points:List[Math.vec3]) -> Math.mat3:
        return _project_on_surface(_surface_default(), _director(points))

class MIDASLocalAxesConvention:
    @staticmethod
    def get_local_axes(points:List[Math.vec3]) -> Math.mat3:
        return _project_on_surface(_surface_default(), _director(points))