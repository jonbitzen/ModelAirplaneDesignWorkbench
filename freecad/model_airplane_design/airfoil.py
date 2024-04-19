import Draft
from enum import Enum
import FreeCAD as App
import numpy
import os.path
import Part
from pathlib import Path
import Sketcher
from typing import List
from . import utilities

class FileType(Enum):
    SELIG = 1
    LEDNICER = 2

class TrailingEdgeType(Enum):
    LINE = 1
    ROUNDED = 2

class AirfoilType(Enum):
    CLARKY = "~/Documents/airfoil-data/clarky.dat"
    GOE173 = "~/Documents/airfoil-data/goe173.dat"
    NACA2411 = "~/Documents/airfoil-data/naca2411.dat"
    NACA654221 = "~/Documents/airfoil-data/naca654221.dat"

class AirfoilData:
    """
    Container for a set of airfoil coordinates in Lednicer format, where
    coordinates originate from the upper airfoil surface trailing edge, and wind
    counter-clockwise through the leading edge and back again to the lower air-
    foil surface trailing edge
    The coordinates stored are the orignal coordinates from the file, and hence
    may not necessarily form a closed loop of edges
    """
    def __init__(self, 
            coords: List[App.Vector], 
            name: str,
            filename: str,
            type: FileType) -> None:
        self.coords = coords
        self.name = name
        self.filename = filename
        self.type =type
        self.master_sketch = self.__to_master_sketch()

    def to_sketch(self, chord: float) -> Sketcher.Sketch:

        scale_factor = chord / self.master_sketch.Shape.BoundBox.XLength
        scaled_af_shape = self.master_sketch.Shape.copy()
        scaled_af_shape.scale(scale_factor)
        return Draft.make_sketch(scaled_af_shape, autoconstraints=True, name=self.name+"-sketch")

    def __to_master_sketch(self, trailing_edge_type=TrailingEdgeType.LINE) -> Sketcher.Sketch:
        """
        Generate a fully-constrained sketch from the internally held Lednicer-type 
        coordinates
        """
        af_sk = Sketcher.Sketch()
        
        bsp = Part.BSplineCurve()
        bsp.interpolate(self.coords)
        constraints: List[Sketcher.Constraint] = []
        bsp_id = af_sk.addGeometry(bsp)
        constraints.append(Sketcher.Constraint("Block", bsp_id))

        # if the start and end point aren't coincident, then we need to close the
        # edge loop
        if self.coords[0] != self.coords[-1]:
            # get the start and end point of the line
            p1 = self.coords[-1]
            p2 = self.coords[0]
            match trailing_edge_type:
                case TrailingEdgeType.LINE:
                    line = Part.LineSegment(p1, p2)
                    line_id = af_sk.addGeometry(line)
                    constraints.append(Sketcher.Constraint("Coincident", bsp_id, 2, line_id, 1))
                    constraints.append(Sketcher.Constraint("Coincident", line_id, 2, bsp_id, 1))
                case TrailingEdgeType.ROUNDED:
                    mid_pt = (p1+p2)/2
                    radius = p1.distanceToPoint(p2)/2

                    # TODO: we use these arc start-stop angles in the hole generator
                    #       as well, maybe we should make a utilitiy for that
                    arc_top = Part.ArcOfCircle(
                        Part.Circle(
                            mid_pt,
                            utilities.z_axis,
                            radius
                        ),
                        numpy.deg2rad(0),
                        numpy.deg2rad(90)
                    )

                    id_arc_top = af_sk.addGeometry(arc_top)

                    constraints.append(Sketcher.Constraint("Tangent", bsp_id, 1, id_arc_top, 2))

                    arc_bottom = Part.ArcOfCircle(
                        Part.Circle(
                            mid_pt,
                            utilities.z_axis,
                            radius
                        ),
                        numpy.deg2rad(-90),
                        numpy.deg2rad(0)
                    )

                    id_arc_bottom = af_sk.addGeometry(arc_bottom)
                    
                    constraints.append(Sketcher.Constraint("Tangent", bsp_id, 2, id_arc_bottom, 1))
                    constraints.append(Sketcher.Constraint("Tangent", id_arc_top, 1, id_arc_bottom, 2))
                    af_sk.addConstraint(constraints)
                    af_sk.solve()
                    constraints.clear()

                    constraints.append(Sketcher.Constraint("Radius", id_arc_bottom, af_sk.Geometries[id_arc_bottom].Radius))
                    constraints.append(Sketcher.Constraint("Radius", id_arc_top, af_sk.Geometries[id_arc_top].Radius))
                    af_sk.addConstraint(constraints)
                    af_sk.solve()
                    constraints.clear()

        af_sk.addConstraint(constraints)

        return af_sk


def load(filename: str) -> AirfoilData:
    """
    Load airfoil coordinates from a .dat file.  Coordinates may be in either Selig
    or Lednicer format
    """
    airfoil_name = Path(filename).stem

    raw_coords = numpy.loadtxt(os.path.expanduser(filename), skiprows=1)

    # Selig contains an initial coordinate pair which represents the number of 
    # coordinates for the upper and lower airfoil surface, which will always 
    # be > 1.0.  Each set of coordinates is from the leading edge to the trailing 
    # edge.
    #
    # Lednicer contains coordinates going counter-clockwise from the upper surface
    # trailing edge, to the leading edge, and then back to the lower surface
    # trailing edge.  Hence it *never* has a coordinate > 1.0
    file_type = FileType.SELIG if raw_coords[0][0] > 1.0 else FileType.LEDNICER

    # rearrange the coordinates in Lednicer order - makes it easier to create
    # a bspline from the coordinates downstream
    match file_type:
        case FileType.SELIG:
            coords = __coords_from_selig(raw_coords)
        case FileType.LEDNICER:
            coords = __coords_from_lednicer(raw_coords)

    return AirfoilData(coords, airfoil_name, filename, file_type)


def __coords_from_selig(raw_coords: numpy.ndarray) -> List[App.Vector]:
    """
    Convert numpy coordinates in Selig format to Lednicer format in FreeCAD
    App.Vector coordinates
    """
    # the first pair of numbers is the number of coordinate rows for the upper
    # and lower airfoil surface
    num_rows = int(raw_coords[0][0])
    raw_coords = numpy.delete(raw_coords, 0, 0)

    # selig is leading edge to trailing edge; to put this in the desired order
    # we reverse the coordinate set for the upper surface of the airfoil
    raw_coords[0:num_rows] = raw_coords[0:num_rows][::-1]

    # the leading edge point will have a duplicate, remove it
    raw_coords = numpy.delete(raw_coords, num_rows, 0)

    coords: List[App.Vector] = [App.Vector(pt[0], pt[1], 0.0) for pt in raw_coords]
    return coords

def __coords_from_lednicer(raw_coords: numpy.ndarray) -> List[App.Vector]:
    """
    Convert numpy coordinates to a list of FreeCAD App.Vector coords
    """
    coords: List[App.Vector] = [App.Vector(pt[0], pt[1], 0.0) for pt in raw_coords]
    return coords