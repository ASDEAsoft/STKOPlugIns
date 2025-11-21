'''
This file contains utility functions for handling local axes in the MidasGenImporter plugin for STKO.
'''

from PyMpc import *
from typing import List, Tuple

# This class implements the local axes convention for MIDAS Gen
class MIDASLocalAxesConvention:

    @staticmethod
    def get_director(points:List[Math.vec3]) -> Math.mat3:
        """
        Computes a director vector or matrix based on the provided list of points.

        Args:
            points (List[Math.vec3]): A list of 2, 3, or 4 3D points (Math.vec3) defining a geometric entity.

        Returns:
            Math.mat3: The computed director vector or matrix, depending on the number of input points.

        Raises:
            Exception: If the number of points is not 2, 3, or 4.

        Notes:
            - For 2 points: Returns the normalized direction vector from the first to the second point.
            - For 3 points: Returns the normalized cross product of the vectors (point1 - point0) and (point2 - point0).
            - For 4 points: Computes midpoints to define local axes and returns the normalized cross product of the local axes.
        """
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

    @staticmethod
    def get_default(points:List[Math.vec3]) -> Math.mat3:
        """
        Returns the default local axes matrix for the given points.

        Args:
            points (List[Math.vec3]): A list of 2, 3, or 4 3D points (Math.vec3).

        Returns:
            Math.mat3: The default local axes matrix.
        """
        N = len(points)
        if N == 2:
            x = Math.vec3(-1.0,  0.0, 0.0)
            y = Math.vec3( 0.0, -1.0, 0.0)
            z = Math.vec3( 0.0,  0.0, 1.0)
            return Math.mat3(x, y, z)
        elif N == 3 or N == 4:
            #TODO: in midas it seems to be the same as in STKO
            x = Math.vec3(1.0,  0.0, 0.0)
            y = Math.vec3(0.0,  0.0, 1.0)
            z = Math.vec3(0.0, -1.0, 0.0)
            return Math.mat3(x, y, z)
        else:
            raise Exception('Invalid number of points for frame or shell local axes')

    @staticmethod
    def project_on_surface(m:Math.mat3, dz:Math.vec3) -> Math.mat3:
        """
        Projects the input matrix on the surface defined by the normal vector dz.

        Args:
            m (Math.mat3): The input matrix to be projected.
            dz (Math.vec3): The normal vector defining the surface.

        Returns:
            Math.mat3: The projected matrix.
        """
        if dz.norm() < 1.0e-4:
            return m
        # start with the input Y axis (which by default is global Z)
        dy = m.col(1)
        if abs(dy.dot(dz)) > (1.0 - 1e-4):
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
        
    @staticmethod
    def project_on_curve(m:Math.mat3, dx:Math.vec3) -> Math.mat3:
        """
        Projects the input matrix on the curve defined by the direction vector dx.

        Args:
            m (Math.mat3): The input matrix to be projected.
            dx (Math.vec3): The direction vector defining the curve.

        Returns:
            Math.mat3: The projected matrix.
        """
        if dx.norm() < 1.0e-4:
            return m
        # start with the local dz, and check if it is aligned with the direction vector (not for vertical frame)
        dz = m.col(2)
        if abs(dz.dot(dx)) > (1.0 - 1e-4):
            # this is a vertical frame
            # use trial dy
            dy = m.col(1)
            dz = dx.cross(dy).normalized()
            dy = dz.cross(dx).normalized()
            return Math.mat3(dx, dy, dz)
        else:
            # use trial dz
            dy = dz.cross(dx).normalized()
            dz = dx.cross(dy).normalized()
            return Math.mat3(dx, dy, dz)

    @staticmethod
    def project(m:Math.mat3, axis:Math.vec3, points:List[Math.vec3]) -> Math.mat3:
        """
        Projects the input matrix on the surface or curve defined by the axis and points.

        Args:
            m (Math.mat3): The input matrix to be projected.
            axis (Math.vec3): The axis vector defining the projection.
            points (List[Math.vec3]): A list of 2, 3, or 4 3D points (Math.vec3).

        Returns:
            Math.mat3: The projected matrix.
        """
        N = len(points)
        if N == 2:
            return MIDASLocalAxesConvention.project_on_curve(m, axis)
        elif N == 3 or N == 4:
            return MIDASLocalAxesConvention.project_on_surface(m, axis)
        else:
            raise Exception('Invalid number of points for frame or shell local axes projection')

    @staticmethod
    def get_local_axes(points:List[Math.vec3]) -> Math.mat3:
        return MIDASLocalAxesConvention.project(
            MIDASLocalAxesConvention.get_default(points),
            MIDASLocalAxesConvention.get_director(points),
            points)