from PyMpc import Math
import math
from MidasGenImporter.document import document, frame
from shapely.geometry import Point, Polygon, LineString, MultiLineString, GeometryCollection
from shapely.ops import linemerge, unary_union, polygonize
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection
from scipy.spatial import KDTree
import numpy as np
import itertools
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection

class floor_polygon:
    """
    A class to represent a polygon for floor loads.
    It is used to store the vertices of the polygon, assumed to be in the XY plane and
    ordered in a counter-clockwise manner.
    The vertices are stored as a list of tuples (ID, x, y).
    This class computes:
    - the total area of the polygon using the shoelace formula.
    - the total perimeter of the polygon.
    - the list of edges as a list of tuples of indices of the vertices.
    - for each edge, it computes the length and it's ratio to the total perimeter (this will be used to distribute the load along the edges).
    """

    class _edge:
        """
        A class to represent an edge of a polygon.
        It is used to store the start and end vertices of the edge,
        and compute the length of the edge.
        """
        def __init__(self, v1:int, v2:int, length:float, ratio:float):
            self.v1 = v1
            self.v2 = v2
            self.length = length
            self.ratio = ratio

    def __init__(self, vertices:List[Tuple[int, float, float]]):
        """
        Initializes the polygon with the given vertices.
        The vertices are assumed to be in the XY plane and ordered in a counter-clockwise manner.
        """
        self.vertices:List[Tuple[int, float, float]] = vertices
        self.edges:List[floor_polygon._edge] = []
        self.area = 0.0
        self.perimeter = 0.0
        self._compute_polygon_properties()

    def _compute_polygon_properties(self):
        """
        Computes the area, perimeter and edges of the polygon.
        The area is computed using the shoelace formula.
        The perimeter is computed as the sum of the lengths of the edges.
        The edges are stored as a list of tuples (v1, v2) where v1 and v2 are the indices of the vertices.
        For each edge, it computes the length and its ratio to the total perimeter.
        """
        n = len(self.vertices)
        if n < 3:
            raise Exception('A polygon must have at least 3 vertices')
        # compute area using the shoelace formula
        self.area = 0.0
        for i in range(n):
            id1, x1, y1 = self.vertices[i]
            id2, x2, y2 = self.vertices[(i + 1) % n]
            self.area += x1 * y2 - x2 * y1
        self.area = abs(self.area) / 2.0
        # compute total perimeter and edge structure
        self.perimeter = 0.0
        for i in range(n):
            v1 = i
            v2 = (i + 1) % n
            id1, x1, y1 = self.vertices[v1]
            id2, x2, y2 = self.vertices[v2]
            length = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            self.edges.append(floor_polygon._edge(v1, v2, length, 0.0))
            self.perimeter += length
        # compute the ratio of each edge length to the total perimeter
        for edge in self.edges:
            edge.ratio = edge.length / self.perimeter if self.perimeter > 0 else 0.0

    def __str__(self):
        # let's print the edges' ratios
        edges_str = ', '.join(f'({self.vertices[edge.v1][0]}, {self.vertices[edge.v2][0]}): {edge.ratio:.4f}' for edge in self.edges)
        return f'Polygon with {len(self.vertices)} vertices, Area: {self.area:.3g}, Perimeter: {self.perimeter:.3g}, Edges: [{edges_str}]'
    
    def __repr__(self):
        return self.__str__()

def floor_polygon_generator(
        corner_points:Dict[int, Math.vec3], 
        doc:document) -> List[floor_polygon]:
    
    # build the main poylgon
    poly = Polygon([(p.x, p.y) for p in corner_points.values()])
    tol = poly.length*1e-3

    # assume all at same z
    floor_z = corner_points.values().__iter__().__next__().z
    print(f'\nFloor polygon at z={floor_z} with area={poly.area} and tol={tol}')

    # build a lookup for vertices in doc which are inside the polygon
    def _is_inside(v : Math.vec3) -> bool:
        p = Point(v.x, v.y)
        return  abs(v.z - floor_z) < tol and (poly.contains(p) or poly.touches(p) or poly.exterior.distance(p) < tol)
    points_inside = {pid: (p.x, p.y) for pid, p in doc.vertices.items() if _is_inside(p)}
    print(f'Found {len(points_inside)} vertices inside the polygon')

    # Prepare cutting lines
    cutters = []
    for _, ele in doc.frames.items():
        try:
            x1, y1 = points_inside[ele.nodes[0]]
            x2, y2 = points_inside[ele.nodes[1]]
            cutters.append(LineString([(x1, y1), (x2, y2)]))
        except KeyError:
            continue
    print(f'Prepared {len(cutters)} cutting lines from frame elements')

    # combine cutters into 1 or more MultiLineString
    merged_cutters = linemerge(cutters)
    print(f'After union, have {len(cutters)} cutter geometries')

    cutter_faces = polygonize(unary_union(merged_cutters))
    print(f'After polygonize, have {len(cutter_faces)} cutter faces')

    pieces = [f.intersection(poly) for f in cutter_faces if f.intersects(poly)]
    print(f'After intersecting with main polygon, have {len(pieces)} pieces')

    # compare total area
    total_area = sum(p.area for p in pieces)
    print(f'Total area of pieces: {total_area}, original polygon area: {poly.area}, difference: {abs(total_area - poly.area)}')
    if abs(total_area - poly.area) > tol:

        # plot for debugging
        fig, ax = plt.subplots()
        mpl_poly = MplPolygon(list(poly.exterior.coords), closed=True, fill=None, edgecolor='black', linewidth=2)
        ax.add_patch(mpl_poly)
        # plot cutters
        for cutter in cutters:
            x, y = cutter.xy
            ax.plot(x, y, color='red', linewidth=1)
        # plot inside points
        x_inside = [p[0] for p in points_inside.values()]
        y_inside = [p[1] for p in points_inside.values()]
        ax.scatter(x_inside, y_inside, color='blue', s=5)
        ax.set_aspect('equal', 'box')
        plt.title('Floor Polygon Generation Debug Plot')
        plt.show()

        raise Exception('Area of pieces does not match original polygon area! Check tolerance or input geometry.')
    
    # create a tree for fast lookup of nearest vertex.
    # we need to map polygon vertices back to document vertex IDs
    # we use the tree for nearest neighbor search, because if 
    # shapely adds a new vertex at the intersection of two lines,
    # that vertex will not exist in the document, so we find the nearest one.
    # but this should never happen if the tolerance is set correctly.
    point_ids = list(itertools.chain(points_inside.keys(), corner_points.keys()))
    coords = np.array([doc.vertices[pid] for pid in point_ids])
    tree = KDTree(coords)

    # Now we need to convert pieces to floor_polygon instances
    floor_polygons = []
    for geom in pieces:
        exterior_coords = list(geom.exterior.coords)[:-1]  # skip the repeated last point
        vertices = []
        for x, y in exterior_coords:
            dist, idx = tree.query([x, y, floor_z])
            if dist > tol:
                raise Exception('Could not find matching vertex in document for polygon vertex')
            pid = point_ids[idx]
            vertices.append((pid, x, y))
        floor_polygons.append(floor_polygon(vertices))
    
    # done
    print(f'Generated {len(floor_polygons)} floor polygons:')
    return floor_polygons