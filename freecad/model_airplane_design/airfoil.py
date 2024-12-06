import Draft
from contextlib import redirect_stdout
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


def get_intersections_with_curve(x_position: float, curve):
    intersections: List[App.Vector] = []
    curve: Part.BoundedCurve = curve
    cut_line = Part.Line(App.Vector(x_position, 0, 0), App.Vector(x_position, 1.0, 0.0))
    trimmed_curve = curve.trim(curve.FirstParameter, curve.LastParameter)
    normal_plane = Part.Plane(utilities.origin, utilities.z_axis)
    pts = trimmed_curve.intersect2d(cut_line, normal_plane)
    for pt in pts:
        intersections.append(App.Vector(pt[0], pt[1], 0))
    intersections.sort(key=lambda pt: pt.y, reverse=True)
    return intersections

# TODO: can we break this up in a way that we have a method to use for a single
#       curve, that way it can be used for trimming the rib
def get_intersections( 
        x_position: float, 
        profile: Part.Shape
    ) -> List[App.Vector]:
    """
        Locates the intersections between an infinite line drawn at the
        x_position, going up and down in the Y direction with the airfoil
        interior profile.

        Parameters
        ----------
        x_position: float
            location along the X axis where the line should be drawn in the
            Y direction, in order to cut the airfoil inner profile in two 
            places

        profile: Part.Shape
            An airfoil shape

        Return
        ------
        List[App.Vector]
            The locations where the infinite line intersects the airfoil shape.
            The upper surface intersection is at index 0, the lower surface
            intersection is at index 1

    """
    intersections: List[App.Vector] = []
    edge: Part.Edge
    normal_plane = Part.Plane(utilities.origin, utilities.z_axis)
    cut_line = Part.Line(App.Vector(x_position, 0, 0), App.Vector(x_position, 1.0, 0.0))
    for edge in profile.Edges:
        curve: Part.Curve = edge.Curve
        curve = curve.trim(*edge.ParameterRange)
        pts = curve.intersect2d(cut_line, normal_plane)
        for pt in pts:
            intersections.append(App.Vector(pt[0], pt[1], 0))
        if len(intersections) == 2:
            break
    
    # TODO: make a simple class that has members for the upper and lower intersection
    #       points explicitly
    intersections.sort(key=lambda pt: pt.y, reverse=True)
    return intersections

# TODO: I think one way for this to work is for the user to add some small set
#       of airfoils to the project.  That's probably a sane choice, since there
#       are over 1600, and it'll make subsequent choices easier since you probably
#       only work with a small subset anyway.  I think the trick may be, then,
#       to create the set of available airfoils using all-caps strings in the
#       FreeCAD enum set
class AirfoilType(Enum):
    CLARKY = 1
    GOE173 = 2
    NACA2411 = 3
    NACA654221 = 4
    AIRCO_DH1 = 5
    AIRCO_DH1_THICK_TE = 6
    AIRCO_DH1_CORRECTED = 7

    @classmethod
    def to_filename(self, airfoil_type: str) -> str:
        """
        Gets the filepath for a stringified version of an AirfoilType enum value

        Parameters
        ----------
        airfoil_type : str
            The AirfoilType enum value as a string (e.g., "CLARKY", GOE173, etc.)

        Return
        ------
        str
            A filepath for a dat file containing coordinates for the requested
            airfoil    
        
        """
        match airfoil_type:
            case AirfoilType.CLARKY.name:
                return "~/Documents/airfoil-data/clarky.dat"
            case AirfoilType.GOE173.name:
                return "~/Documents/airfoil-data/goe173.dat"
            case AirfoilType.NACA2411.name:
                return "~/Documents/airfoil-data/naca2411.dat"
            case AirfoilType.NACA654221.name:
                return "~/Documents/airfoil-data/naca654221.dat"
            case AirfoilType.AIRCO_DH1.name:
                return "~/Documents/airfoil-data/airco-dh1.dat" 
            case AirfoilType.AIRCO_DH1_THICK_TE.name:
                return "~/Documents/airfoil-data/airco-dh1-thick-TE.dat"
            case AirfoilType.AIRCO_DH1_CORRECTED.name:
                return "~/Documents/airfoil-data/airco-dh1-corrected.dat"
            case _:
                print("AirfoilType.to_filename - unknown airfoil type \"" + airfoil_type + "\"")
                return None

class AirfoilData:
    """
    Stores airfoil coordinate data loaded from a dat file, and provides factory-
    like methods to create useful FreeCAD representations of the airfoil data
    for solid modeling

    Parameters
    ----------
    coords : List[App.Vector]
        Original coordinate data loaded from file, in Lednicer format, counter-
        clockwise winding from the top surface trailing edge, through the leading-
        edge, and around back to the lower surface trailing edge.  Note that
        the coordinate from the original file may not necessarily form a closed
        loop

    name : str
        Human-friendly name of the airfoil data

    filename : str
        Orginal source filename for the dat file where the airfoil coordinates
        were obtained

    type : FileType
        Coordinate encoding in the original source data file.  May be one of
        {SELIG, LEDNICER}

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

    def to_shape(self, chord: float, body_centered: bool=False) -> Part.Shape:
        """
        Scales the master airfoil sketchto the desired chord size, returns a
        Shape

        Parameters
        ----------
        chord : float
            The rib chord is the distance from the airfoil leading edge to its
            trailing edge
        
        body_centered: bool
            The shape should be returned so that its body origin is located at
            the first moment of area

        Returns
        -------
        Part.Shape
            A version of the master airfoil shape scaled to have the user-
            directed chord length

        """
        af_sketch: Sketcher.Sketch
        bbox: App.BoundBox
        af_sketch, bbox = self.__to_master_sketch(chord)
        af_shape = af_sketch.Shape.copy()
        if body_centered:
            body_center = bbox.Center
            transform_mtx = App.Matrix()
            transform_mtx.move(-body_center)
            af_shape = af_shape.transformed(transform_mtx, True)

        return af_shape


    def to_sketch(self, chord: float, body_centered: bool=False) -> Sketcher.Sketch:
        """
        Scales the master airfoil sketch to the desired chord size, returns an
        editable SketchObject

        Parameters
        ----------
        chord: float
            The length of the airfoil from leading edge to trailing edge

        body_centered: bool
            The shape should be returned so that its body origin is located at
            the first moment of area

        Return
        ------
        Sketcher.Sketch
            A fully constrained sketch derived from the airfoil data file
            coordinates, scaled to the user-designated chord length
        """
        scaled_af_shape = self.to_shape(chord, body_centered)
        return Draft.make_sketch(scaled_af_shape, autoconstraints=True, name=self.name+"-sketch")

    def __to_master_sketch(self, chord: float, trailing_edge_type=TrailingEdgeType.LINE) -> Sketcher.Sketch:
        """
        Generate a fully-constrained sketch from the internally held Lednicer-type 
        coordinates.  The airfoil edges derived from the coordinates will be
        generated as a BSpline

        Parameters
        ----------
        trailing_edge_type : TrailingEdgeType
            An enum arg that indicates how to close the airfoil BSpline if the
            trailing edge coordinates for an open loop:
                LINE 
                    a line will be drawn between the open trailing edge coordinates. 
                ROUNDED
                    a pair of circular arcs will close the open trailing edge
                    coordinates.  Arcs will be constrained such that all joining
                    edges will be tangent to one another and to the airfoil
                    BSpline endoiints

        Return
        ------
        Sketcher.Sketch
            A fully constrained Sketch derived from the dat file coordinates.
            Note that the master sketch is not a DocumentObject, and hence is
            not a displayable entity

        """
        af_sk = Sketcher.Sketch()

        coords: List[App.Vector] = []
        bbox = App.BoundBox()
        bbox.XMin = 0.0; bbox.XMax = chord
        bbox.ZMin = 0.0; bbox.ZMax = 0.0
        bbox.YMin = 0.0; bbox.YMin = 0.0
        for coord in self.coords:
            scaled_coord = chord*coord
            coords.append(scaled_coord)
            bbox.YMin = min(bbox.YMin, scaled_coord.y)
            bbox.YMax = max(bbox.YMax, scaled_coord.y)

        bsp = Part.BSplineCurve()
        
        bsp.interpolate(coords)
        constraints: List[Sketcher.Constraint] = []
        bsp_id = af_sk.addGeometry(bsp)
        constraints.append(Sketcher.Constraint("Block", bsp_id))

        # if the start and end point aren't coincident, then we need to close the
        # edge loop
        if coords[0] != coords[-1]:
            # get the start and end point of the line
            p1 = coords[-1]
            p2 = coords[0]
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

        return af_sk, bbox

def load(filename: str) -> AirfoilData:
    """
    Load airfoil coordinates from a .dat file.  Coordinates may be in either Selig
    or Lednicer format

    Parameters
    ---------
    filename : str
        fully-qualified filepath of a dat file containing airfoil coordinates

    Return
    ------
    AirfoilData
        A factory object which wraps the data from the dat file.  May be used to
        create various entities useful for solid modeling.
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

    Parameters
    ----------
    raw_coords : numpy.ndarray
        Airfoil coordinates in Selig format

    Return
    ------
    List[App.Vector]
        Airfoil coordinates in Lednicer format, using FreeCAD types
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

    Parameters
    ----------
    raw_coords : numpy.ndarray
        Airfoil coordinates in Lednicer format

    Return
    ------
    List[App.Vector]
        Coordinates in Lednicer format, using FreeCAD types

    """
    coords: List[App.Vector] = [App.Vector(pt[0], pt[1], 0.0) for pt in raw_coords]
    return coords

def save(
        coord_sketch: Sketcher.Sketch, 
        filename: str,
        info: str) -> bool:
    """
    Saves airfoil coordinates from a sketch into a dat file, in Selig format

    Parameters
    ----------
    coord_sketch: Sketcher.SketchObject
        Sketch that contains airfoil coordinates in Selig arrangement, which entails
        the following on the part of the modeler:
        - the airfoil chord is normalized to have a length of 1 along the x axis
        - the airfoil leading edge is at the sketch origin point; it should be a
          single point where the leading edge is tangent to the y axis at the
          origin
        - the airfoil trailing edge is at coordinate (1,0) along the x axis; it 
          may be either a single point, or two points.  If the trailinge edge is
          two points, they should be centered on a vertical line passing through
          (1,0)
        - the remaining points should be x-coordinate sections that intersect
          the upper and lower surface exactly once at the same x coordinate
        - the airfoil upper surface y coordinates should be "up" (upper surface
          y coordinate is > the lower surface y coordinate)
    
    filename: str
        Fully qualified filepath to the location where the coordinates will be
        stored.  Should end with a ".dat" extension according to convention

    info: str
        Brief file header info used to describe the airfoil
    """
    af_coords: List[App.Vector] = []

    for idx in range(len(coord_sketch.Geometry)):
        if coord_sketch.getConstruction(idx):
            continue

        geom = coord_sketch.Geometry[idx]
        if geom.TypeId != "Part::GeomPoint":
            continue

        af_coords.append(App.Vector(geom.X, geom.Y, geom.Z))

    af_coords.sort(key=lambda pt: pt.x)

    upper_coords: List[App.Vector] = []
    lower_coords: List[App.Vector] = []

    # store the last coordinate in the list of it doesn't have a partner with
    # an identical x-coordinate
    extra_coord: App.Vector = None

    # if the first two coordinates dont have the same x-value, then we had a
    # single point at the leading edge
    if af_coords[0].x != af_coords[1].x:
        upper_coords.append(af_coords[0])
        lower_coords.append(af_coords[0])
        af_coords.pop(0)

    # if the last two coordinates dont have the same x-value, then we had a
    # single point at the trailing edge
    if af_coords[-1].x != af_coords[-2].x:
        extra_coord = af_coords[-1]
        af_coords.pop(-1)

    # now that we've handled the special cases, all the rest of the coordinates
    # should be upper surface/lower surface pairs, so walk them two at a time
    idx: int = 0
    while idx < len(af_coords):
        
        if af_coords[idx].x != af_coords[idx+1].x:
            print("Error: Upper and lower surface coordinate x values do not match")
            return False

        # order the coords in a list according to their y coordinate so the
        # upper coordinate is first in the list
        surf_coords: List[App.Vector] = []
        surf_coords.append(af_coords[idx])
        surf_coords.append(af_coords[idx+1])
        surf_coords.sort(key=lambda pt: pt.y, reverse=True)

        upper_coords.append(surf_coords[0])
        lower_coords.append(surf_coords[1])
        idx += 2
    
    if extra_coord is not None:
        upper_coords.append(extra_coord)
        lower_coords.append(extra_coord)
        

    def write_coords(coord_list: str):
        for idx in range(len(coord_list)):
            y_coord: str = "{0:.{1}f}".format(numpy.abs(coord_list[idx].y), 7)
            if coord_list[idx].y < 0:
                y_coord = '-' + y_coord[1:]
            print(" {0:.{1}f}".format(coord_list[idx].x, 7) + " " + y_coord)

    with open(filename, 'w') as out_file:
        with redirect_stdout(out_file):
            print(" " + info)
            print(' {:>8}.'.format(len(upper_coords)) + ' {:>8}.'.format(len(upper_coords)))
            print("")
            write_coords(upper_coords)
            print("")
            write_coords(lower_coords)

    return True