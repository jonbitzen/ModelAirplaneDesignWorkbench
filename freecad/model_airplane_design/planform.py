import FreeCAD as App
import Part
from . import utilities
from typing import List, Tuple

class Planform():    
    """
    Helper class to locate chord length and rib XY center location given a list
    Edges that represent a wing planform

    Parameters
    ----------
    pf_edges: List[Part.Edge]
        A list of edges that represents the wing planform  The planform edges
        should meet the following requirements:
        - fully in the XY plane, Z=0
        - the wing root should be located at X=0 (on the Y axis)
        - the wing root points should be located such that they are symmetric
          with the X axis
        - the wing tip should be in -X
        - the wing leading edge is in -Y, the wing trailing edge is in +Y

    """
    def __init__(self, pf_edges: List[Part.Edge]) -> None:
        self.pf_edges = Part.__sortEdges__(pf_edges)

    def get_rib_chord_at(self, position: App.Vector) -> Tuple[float, App.Vector]:
        """
        Calculates the rib chord and XY center point at a given position within
        the planform.

        Parameters
        ----------
        position: App.Vector
            location at which we should search for intersections with the planform
            edges

        Return
        ------
        Tuple[float, App.Vector]
            float: the rib chord at position
            App.Vector: the rib center location in the XY plane at position

        """
        plane = Part.Plane(position, -utilities.x_axis)

        # the plane should intersect the wing planform in two places
        pf_intersections: List[Part.Vertex] = []

        for edge in self.pf_edges:

            trimmed_curve: Part.Curve = edge.Curve
            trimmed_curve = trimmed_curve.trim(*edge.ParameterRange)

            intersections = plane.intersect(trimmed_curve)
            # if there are no intersections, skip the rest of the loop
            if intersections is None:
                continue
            
            p: Part.Point
            for p in intersections[0]:
                pf_intersections.append(Part.Vertex(p))

            # if we have found the needed planform intersections, we're done
            if len(pf_intersections) == 2:
                p0: Part.Vertex = pf_intersections[0]
                p1: Part.Vertex = pf_intersections[1]
                ctr_pt = (p0.Point+p1.Point)/2
                dist_t = p0.distToShape(p1)
                dist = dist_t[0]
                return dist, ctr_pt