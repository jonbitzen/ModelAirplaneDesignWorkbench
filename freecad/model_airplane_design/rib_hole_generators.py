from . import utilities
from abc import ABC, abstractmethod
import Draft
import DraftGeomUtils
from enum import Enum
import FreeCAD as App
import math
import numpy
import Part
import Sketcher
from typing import List, Dict


class Interval():
    """
    Stores starting and ending interval values used to compute lightening hole
    placement
    """
    def __init__(self, start: float = 0.0, end: float = 0.0):
        self.start = start
        self.end = end

    def length(self) -> float:
        """
        Provides the length of the interval
        """
        return self.end-self.start

    def __str__(self):
        return "start=" + str(self.start) + ", end=" + str(self.end)
    
    def __lt__(self, other: 'Interval'):
        return self.start < other.start

class HoleExclusion():
    """
    Stores data that indicates where the hole generator may *not* generate a
    lightening hole in the wing rib.  Note that the hole exclusion must be in
    the same plane as the rib airfoil outline sketch

    Parameters
    ----------
    geometry: Part.Wire
        Wire outline of arbitrary excluded geometry, where a lightening hole may
        not be generated.  The geometry must be planar, and should be a closed
        wire as well
    
    standoff: float
        Measurement that will be applied to the geometry profile to form an
        exterior offset.  The exterior offset around the excluded geometry
        reflects structure that is necessary to support the structure in the rib
    """
    def __init__(
            self,
            geometry: Part.Wire,
            standoff: float
        ) -> None:
        
        if not DraftGeomUtils.is_planar(geometry):
            raise ValueError("Hole exclusion geometry is non-planar")

        if standoff < 0.0:
            raise ValueError("Hole exclusion standoff must not be negative")

        self.geometry = geometry
        self.standoff = standoff

    def get_excluded_region(self) -> App.BoundBox:
        region = self.geometry.makeOffset2D(self.standoff, join=2)
        return region.BoundBox
    
    def get_excluded_interval(self) -> Interval:
        bbox = self.get_excluded_region()
        return Interval(bbox.XMin, bbox.XMax)

class HoleBoundRegion():
    """
    A region where the hole generator may generate one or more holes in the rib
    structure

    Parameters
    ----------
    af_inner_contour: Part.Shape
        An interior contour offset from the rib's outer airfoil by some physical
        measurement

    interval: Interval
        A measurement along the X-axis which indicates that valid region where
        the hole generator may generate one or more holes
    """
    def __init__(
            self, 
            af_inner_contour: Part.Shape, 
            interval: Interval) -> None:
        if af_inner_contour is None:
            raise ValueError("Airfoil inner contour must not be None")
        
        if interval is None:
            raise ValueError("Hole bound interval must not be None")

        self.af_inner_contour = af_inner_contour
        self.interval = interval        

# TODO: rename to HoleIntervalGenerator, since it is a region for one or more
#       holes
class HoleBoundGenerator():
    """
        Generates one or more HoleBoundRegions inside a rib, such that the
        HoleBoundRegions do not intersect with a list of excluded regions where
        holes may not be drawn.  The hole generator may generate one or more
        holes within each HoleBoundRegion

        Parameters
        ----------

        airfoil_sk: Sketcher.Sketch
            an airfoil sketch object, which represents the outer physical boundary
            of the rib
        
        inset_width: float
            a measurement for the inner contour to be drawn from within the rib
            airfoil; this will form the outer boundary for all holes generated
            by the hole generator

        excluded_regions: List[HoleExclusion]
            a list of areas (including an optional keep-out margin) where holes
            may not be drawn by the hole generator
    """
    def __init__(
            self,
            airfoil_sk: Sketcher.Sketch,
            inset_width: float,
            excluded_regions: List[HoleExclusion]
        ) -> None:
        self.airfoil_sk = airfoil_sk
        self.inset_width = inset_width
        self.excluded_regions = excluded_regions

        # this is a place to keep and delete temporary doc objects for viewing
        self.show_items: List[App.DocumentObject] = []

        self.af_inner_profile = self.airfoil_sk.Shape.makeOffset2D(-self.inset_width, join=2)

        # create a list of intervals, and make sure they are sorted by starting
        # value so we can insert properly
        intf_intervals: List[Interval] = [excl.get_excluded_interval() for excl in self.excluded_regions]
        intf_intervals.sort()

        # first, calculate the nominal max hole width; this will be used to calculate
        # the number of holes in each of the space intervals between interferences
        inner_profile_bbox = self.af_inner_profile.BoundBox
    
        profile_interval = Interval(inner_profile_bbox.XMin, inner_profile_bbox.XMax)
        self.hole_intervals: List[HoleBoundRegion] = []

        current_interval = Interval(profile_interval.start, 0.0)
        for intf_interval in  intf_intervals:
            current_interval.end = intf_interval.start
            self.hole_intervals.append(HoleBoundRegion(self.af_inner_profile, current_interval))
            current_interval = Interval(intf_interval.end, 0.0)
        current_interval.end = profile_interval.end
        self.hole_intervals.append(HoleBoundRegion(self.af_inner_profile, current_interval))

    def get_hole_bounds(self) -> List[HoleBoundRegion]:
        """
        The list of valid regions where the hole generator may draw rib lightening
        holes
        """
        return self.hole_intervals
    
    def show(self) -> None:
        """
        This utility method draws reference geometry to display the rib inner
        profile, as well as the valid hole intervals in the X axis as boxes
        """
        af_profile = Part.show(self.af_inner_profile)
        af_profile.ViewObject.LineColor = utilities.RED(0.5)
        self.show_items.append(af_profile)

        af_prof_bbox: App.BoundBox = self.af_inner_profile.BoundBox

        y_coord: float = (af_prof_bbox.YMax + af_prof_bbox.YMin)/2
        box_height: float = af_prof_bbox.YLength
        for hbr in self.hole_intervals:
            width: float = hbr.interval.length()
            placement = utilities.xy_placement.copy()
            placement.Base.y = y_coord - box_height/2
            placement.Base.x = hbr.interval.start
            intv = Draft.make_rectangle(width, box_height, placement, face=False)
            intv.ViewObject.LineColor = utilities.GREEN(0.5)
            self.show_items.append(intv)

    def unshow(self) -> None:
        """
        This utility clears any previously drawn reference geometry
        """
        for doc_obj in self.show_items:
            App.ActiveDocument.removeObject(doc_obj.Name)

class RoundedTrapezoidEndcap():
    """
    Helper class which encapsulates control points for the end-caps of a rounded
    trapezoidal lightening hole.

    Parameters
    ----------
    pts: List[App.Vector]
        List of points indicating the geometry of the end-cap.  If one point is
        provided, then the end cap is a rounded.  If two points are provided, the
        end cap is a line
    """
    def __init__(self, pts: List[App.Vector]) -> None:
        self.pts = pts

    def is_line(self) -> bool:
        """
        Indicates whether the end cap represents a line
        """
        return len(self.pts) == 2
    
    def is_point(self) -> bool:
        """
        Indicates whether the end cap represents a rounded area
        """
        return len(self.pts) == 1

class EndcapSide(Enum):
    LEFT = 1
    RIGHT = 2

class ChamferPointLocation(Enum):
    UPPER = 1
    LOWER = 2
    CENTER = 3

class RoundedTrapezoidControlPoints():
    """
    Contains the control points used to generate a single rounded trapezoidal rib
    lightening hole

    Parameters
    ----------
    rib_inner_profile: Part.Shape
        the inner profile in the rib, which forms the outermost upper and lower
        bound for the rounded trapezoidal hole
    
    interval: Interval
        represents the X coordinate left and right bounds for the proposed hole

    chamfer_rad: float
        maximum radius of the rounded chamfers that join the left and right
        vertical hole lines to the upper and lower hole surfaces
    """
    def __init__(
            self,
            rib_inner_profile: Part.Shape,
            interval: Interval,
            chamfer_rad: float
        ) -> None:
        
        self.show_items: List[App.DocumentObject] = []

        # control points for the upper and lower surface of the airfoil inner profile
        self.upper_spline_pts: List[App.Vector] = []
        self.lower_spline_pts: List[App.Vector] = []

        self.left_chamfer_pts: Dict[ChamferPointLocation, App.Vector] = {}
        self.right_chamfer_pts: Dict[ChamferPointLocation, App.Vector] = {}

        self.left_endcap: RoundedTrapezoidEndcap = \
            self.__get_side_coords(interval.start, rib_inner_profile, chamfer_rad)
        self.right_endcap: RoundedTrapezoidEndcap = \
            self.__get_side_coords(interval.end, rib_inner_profile, chamfer_rad)

        # generate control points for the upper and lower splines, which will
        # hug the airfoil interior profile (airfoil_profile)
        # TODO: we need to make sure that right_bound-left_bound is greater than
        #       zero at a minimum, and probably also that it's greater than some
        #       practical minimum length as well, or the algo below will explode
        spl_left_bound = interval.start+chamfer_rad
        spl_right_bound = interval.end-chamfer_rad
        num_ctl_pts = 5
        for x_coord in numpy.linspace(spl_left_bound, spl_right_bound, num_ctl_pts):
            pts = self.__get_intersections(x_coord, rib_inner_profile)
            self.upper_spline_pts.append(pts[0])
            self.lower_spline_pts.append(pts[1])

        self.left_chamfer_pts = self.__set_chamfer_pts(self.left_endcap, EndcapSide.LEFT)
        self.right_chamfer_pts = self.__set_chamfer_pts(self.right_endcap, EndcapSide.RIGHT)

    def __set_chamfer_pts(
            self,
            endcap: RoundedTrapezoidEndcap, 
            side: EndcapSide
        ) -> Dict[ChamferPointLocation, App.Vector]:
        chamfer_pts: Dict[ChamferPointLocation, App.Vector] = {}
        
        # this helps us get spline points from either the left or right side
        # of the spline, to calculate the chamfer points
        side_idx: int
        match side:
            case EndcapSide.LEFT:
                side_idx = 0
            case EndcapSide.RIGHT:
                side_idx = -1

        # if we have a line-type endcap, then we have an upper and lower
        # chamfer point to calculate
        if endcap.is_line():
            chamfer_pts[ChamferPointLocation.UPPER] = \
                App.Vector(self.upper_spline_pts[side_idx].x, endcap.pts[0].y, 0)
            chamfer_pts[ChamferPointLocation.LOWER] = \
                App.Vector(self.lower_spline_pts[side_idx].x, endcap.pts[1].y, 0)
        # if the endcap was too small for a line, then find the midpoint
        # between the upper and lower spline endpoints
        else:
            upper_spline_pt = self.upper_spline_pts[side_idx]
            lower_spline_pt = self.lower_spline_pts[side_idx]
            chamfer_pts[ChamferPointLocation.CENTER] = (upper_spline_pt+lower_spline_pt)/2
        
        return chamfer_pts     

    def __get_side_coords(
            self, 
            x_position: float, 
            profile: Part.Shape,
            chamfer_rad: float
        ) -> RoundedTrapezoidEndcap:
        """
            Generates end cap coordinates for one side of the rounded trapezoidal
            hole.  The end cap side should include a line and two chamfers if 
            possible.  If the available space between the upper and lower surfaces
            of the airfoil inner profile is less than two chamfer radii, a single
            mid-point is generated instead, which is where two chamfer arcs will
            be joined without a line between them

            Parameters
            ----------
            x_position: float
                location along the X axis where the end cap coordinates should be
                generated

            profile: Part.Shape
                inner contour with respect to the rib's airfoil shape

            chamfer_rad: float
                the maximum radius allowed for the chamfer
        """

        corner_pts: List[App.Vector] = self.__get_intersections(x_position, profile)
        
        p1 = corner_pts[0]
        p2 = corner_pts[1]

        # ensure the distance between p1 and p2 is long enough to fit a non-zero
        # line segment between the chamfer reference points
        line_coord_pts: List[App.Vector] = []
        if p1.distanceToPoint(p2) > 2* chamfer_rad:
            b1 = App.Vector(p1)
            b2 = App.Vector(p2)
            b1.y = b1.y - chamfer_rad
            b2.y = b2.y + chamfer_rad
            line_coord_pts = [b1,b2]
        # if the distance between p1 and p2 is below the min threshold for the
        # chamfer offset, then just provide a single point midway between p1 and p2
        else:
            b1 = (p1 + p2)/2
            line_coord_pts = [b1]

        return RoundedTrapezoidEndcap(line_coord_pts)

    def __get_intersections(
            self, 
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
                airfoil inner profile
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
        
        intersections.sort(key=lambda pt: pt.y, reverse=True)
        return intersections
    
    # TODO: when the endcap is a center point rather than a line, you really want
    #       to keep the hole spacing uniform, you want to extend the spline out
    #       toward the edge, and then force the rounded part to go through the
    #       outer side point

    def show(self) -> None:
        self.show_items.extend(utilities.draw_points(self.upper_spline_pts))
        self.show_items.extend(utilities.draw_points(self.lower_spline_pts))
        self.show_items.extend(utilities.draw_points(self.left_endcap.pts))
        self.show_items.extend(utilities.draw_points(self.right_endcap.pts))
        for pt_loc in ChamferPointLocation:
            if pt_loc in self.left_chamfer_pts:
                self.show_items.extend(utilities.draw_points([self.left_chamfer_pts[pt_loc]]))

            if pt_loc in self.right_chamfer_pts:
                self.show_items.extend(utilities.draw_points([self.right_chamfer_pts[pt_loc]]))

    def unshow(self) -> None:
        for doc_obj in self.show_items:
            App.ActiveDocument.removeObject(doc_obj.Name)

class HoleGenerator(ABC):
    """
    Interface for classes that will draw lightening holes inside a wing rib
    """
    def generate_sketch(self, interval: HoleBoundRegion) -> Sketcher.Sketch:
        pass

class RoundedTrapezoidHoleGenerator(HoleGenerator):
    """
    Generates one or more rounded trapezoidal holes in a single HoleBoundingRegion

    Parameters
    ----------
    max_chamfer_rad: float
        maximum radius for chamfers at the corners between end cap lines and the
        upper or lower spline curves that form the hole

    max_hole_length: float
        maximum length that holes generated in the HoleBoundingRegion may be

    min_hole_spacing: float
        minumum material that separates each generated hole in the HoleBoundingRegion
    """
    def __init__(
            self,
            max_chamfer_rad: float,
            max_hole_length: float,
            min_hole_spacing: float
        ) -> None:
        super().__init__()
        
        self.max_chamfer_rad = max_chamfer_rad
        self.max_hole_length = max_hole_length
        self.min_hole_spacing = min_hole_spacing

    def generate_sketch(self, bdg_region: HoleBoundRegion) -> Sketcher.Sketch:
        """
        Generate a sketch representing one or more rounded trapezoidal lightening
        holes inside the given HoleboundingRegion

        Parameters
        ----------
        interval: HoleBoundRegion
            a region where the hole generator may generate one or more rib
            lightening holes
        """
                
        # create the sketch on the YZ plane
        sk: Sketcher.Sketch = App.ActiveDocument.addObject("Sketcher::SketchObject", "lightening-holes")
        sk.Placement = utilities.xy_placement    

        constraints: List[Sketcher.Constraint] = []

        def draw_line(endcap: RoundedTrapezoidEndcap) -> tuple[int, Part.LineSegment]:
            if endcap.is_line():
                line = Part.LineSegment(endcap.pts[0], endcap.pts[1])
                id = sk.addGeometry(line, False)
                constraints.append(Sketcher.Constraint("Block", id))
                return id, line
            return -1, None
        
        def draw_spline(pts: List[App.Vector]) -> tuple[int, Part.BSplineCurve]:
            bsp = Part.BSplineCurve()
            bsp.interpolate(pts)
            id = sk.addGeometry(bsp)
            constraints.append(Sketcher.Constraint("Block", id))
            return id, bsp

        # curve: Part.BoundedCurve
        def get_nearest_end_to_pt(target_pt: App.Vector, curve):
            pt_dist = [
                (1, curve.StartPoint.distanceToPoint(target_pt)),
                (2, curve.EndPoint.distanceToPoint(target_pt))
            ]
            pt_dist.sort(key=lambda x : x[1])
            return pt_dist[0][0]
        
        def draw_chamfers(
            chamfer_pts: Dict[ChamferPointLocation, App.Vector],
            max_chamfer_rad: float,
            line: Part.LineSegment, line_id: int,
            upper_bsp: Part.BSplineCurve, upper_bsp_id: int,
            lower_bsp: Part.BSplineCurve, lower_bsp_id: int,
            endcap_side: EndcapSide
        ) -> None:
            """
            Draw circular arcs that join a side boundary line with one of the bspline
            interior boundary sections
            """
                
            top_arc_start: float; top_arc_end: float
            bottom_arc_start: float; bottom_arc_end: float
            upper_spl_pt: App.Vector; lower_spl_pt: App.Vector

            match endcap_side:
                case EndcapSide.LEFT:
                    top_arc_start = numpy.deg2rad(90)
                    top_arc_end = numpy.deg2rad(180)
                    bottom_arc_start = numpy.deg2rad(-180)
                    bottom_arc_end = numpy.deg2rad(-90)
                    upper_spl_pt = upper_bsp.StartPoint
                    lower_spl_pt = lower_bsp.StartPoint
                    
                case EndcapSide.RIGHT:
                    top_arc_start = numpy.deg2rad(0)
                    top_arc_end = numpy.deg2rad(90)
                    bottom_arc_start = numpy.deg2rad(-90)
                    bottom_arc_end = numpy.deg2rad(0)
                    upper_spl_pt = upper_bsp.EndPoint
                    lower_spl_pt = lower_bsp.EndPoint
                    

            # if the space was too small to form a line with the target offset, then
            # we create two arcs and join them directly, without a line in the middle
            if line is None:
                # p1 = chamfer_pts[0]
                # p2 = chamfer_pts[1]
                p1 = upper_spl_pt
                p2 = lower_spl_pt

                mid_pt = (p1+p2)/2
                radius = p1.distanceToPoint(p2)/2  

                arc_top = Part.ArcOfCircle(
                    Part.Circle(
                        mid_pt,
                        utilities.z_axis,
                        radius
                    ),
                    top_arc_start,
                    top_arc_end
                )

                id_arc_top = sk.addGeometry(arc_top, False)
                bsp_pt_id = 1 if upper_bsp.StartPoint == p1 else 2
                arc_pt_id = get_nearest_end_to_pt(p1, arc_top)
                top_arc_other_id = 1 if arc_pt_id == 2  else 2

                # TODO: We really need to know why the order that we apply the 
                #       constraints matters so much; if you do the arc first, 
                #       everything is rejected as over-constrained 
                constraints.append(Sketcher.Constraint("Tangent", upper_bsp_id, bsp_pt_id, id_arc_top, arc_pt_id))

                arc_bottom = Part.ArcOfCircle(
                    Part.Circle(
                        mid_pt,
                        utilities.z_axis,
                        radius
                    ),
                    bottom_arc_start,
                    bottom_arc_end
                )

                id_arc_bottom = sk.addGeometry(arc_bottom, False)
                bsp_pt_id = 1 if lower_bsp.StartPoint == p2 else 2
                arc_pt_id = get_nearest_end_to_pt(p2, arc_bottom)
                bottom_arc_other_id = 1 if arc_pt_id == 2 else 2

                constraints.append(Sketcher.Constraint("Tangent", lower_bsp_id, bsp_pt_id, id_arc_bottom, arc_pt_id))
                constraints.append(Sketcher.Constraint("Tangent", id_arc_top, top_arc_other_id, id_arc_bottom, bottom_arc_other_id))
            # we have a line in the middle, so create arcs that join the line on each
            # end to the top and bottom spline sections
            else:
                arc_radius = max_chamfer_rad
                top_arc_center = App.Vector(upper_spl_pt.x, line.StartPoint.y, 0.0)

                arc_top = Part.ArcOfCircle(
                    Part.Circle(
                        top_arc_center,
                        utilities.z_axis,
                        arc_radius
                    ),
                    top_arc_start,
                    top_arc_end
                )

                id_arc_top = sk.addGeometry(arc_top, False)
                bsp_pt_id = 1 if upper_bsp.StartPoint == upper_spl_pt else 2
                arc_pt_id = get_nearest_end_to_pt(upper_spl_pt, arc_top)
                top_arc_other_id = 1 if arc_pt_id == 2  else 2
                
                constraints.append(Sketcher.Constraint("Tangent", upper_bsp_id, bsp_pt_id, id_arc_top, arc_pt_id))
                constraints.append(Sketcher.Constraint("Tangent", id_arc_top, top_arc_other_id, line_id, 1))

                btm_arc_ctr = App.Vector(lower_spl_pt.x, line.EndPoint.y, 0.0)

                arc_bottom = Part.ArcOfCircle(
                    Part.Circle(
                        btm_arc_ctr,
                        utilities.z_axis,
                        arc_radius
                    ),
                    bottom_arc_start,
                    bottom_arc_end
                )

                id_arc_bottom = sk.addGeometry(arc_bottom, False)
                bsp_pt_id = 1 if lower_bsp.StartPoint == lower_spl_pt else 2
                arc_pt_id = get_nearest_end_to_pt(lower_spl_pt, arc_bottom)
                btm_arc_other_id = 1 if arc_pt_id == 2  else 2

                constraints.append(Sketcher.Constraint("Tangent", lower_bsp_id, bsp_pt_id, id_arc_bottom, arc_pt_id))
                constraints.append(Sketcher.Constraint("Tangent", id_arc_bottom, btm_arc_other_id, line_id, 2))

        def get_hole_control_pts(
                bdg_region: HoleBoundRegion, 
                hole_spacing: float,
                max_chamfer_rad: float,
                max_hole_length: float
            ) -> List[RoundedTrapezoidControlPoints]:
            
            # create a number of hole segments from the bounding region interval
            # the segments do not (yet) include spacing in between holes
            num_segments: int = math.ceil(bdg_region.interval.length()/max_hole_length)
            if num_segments == 0:
                return []
            
            # create a set of control point groups where each group is separated
            # by at distance hole_spacing
            ctl_pts: List[RoundedTrapezoidControlPoints] = []
            
            total_wall_width = (num_segments-1)*hole_spacing
            total_hole_width = bdg_region.interval.length()-total_wall_width
            hole_width = total_hole_width/num_segments

            dist = bdg_region.interval.start
            # for idx in range(0,num_segments-1):
            for idx in range(0,num_segments):
                # get the endpoints for the hole
                start = dist
                dist += hole_width
                end = dist
                
                # get the starting point for the next hole, if there is one
                dist += hole_spacing

                # add a set of control points for the hole
                ctl_pts.append(
                    RoundedTrapezoidControlPoints(
                        bdg_region.af_inner_contour, 
                        Interval(start,end), 
                        max_chamfer_rad
                    )
                )
            return ctl_pts
        
        hole_ctl_pts_list = \
            get_hole_control_pts(
                bdg_region, 
                self.min_hole_spacing, 
                self.max_chamfer_rad, 
                self.max_hole_length
            )
        
        for ctl_pts in hole_ctl_pts_list:

            left_line_id, left_line = draw_line(ctl_pts.left_endcap)
            right_line_id, right_line = draw_line(ctl_pts.right_endcap)
            upper_bsp_id, upper_bsp = draw_spline(ctl_pts.upper_spline_pts)
            lower_bsp_id, lower_bsp = draw_spline(ctl_pts.lower_spline_pts)

            draw_chamfers(
                ctl_pts.left_chamfer_pts,
                self.max_chamfer_rad,
                left_line, 
                left_line_id, 
                upper_bsp, 
                upper_bsp_id, 
                lower_bsp, 
                lower_bsp_id, 
                EndcapSide.LEFT
            )
            
            draw_chamfers(
                ctl_pts.right_chamfer_pts, 
                self.max_chamfer_rad,
                right_line, 
                right_line_id, 
                upper_bsp, 
                upper_bsp_id, 
                lower_bsp, 
                lower_bsp_id, 
                EndcapSide.RIGHT
            )

        sk.addConstraint(constraints)

        return sk

class LighteningHoleBounds:
    """
    Stores a set of boundaries and control points to generate lightening holes
    for an airfoil
    """
    def __init__(
        self, 
        airfoil_shape: Part.Shape, 
        left_bound: Part.Line, 
        right_bound: Part.Line,
        chamfer_radius: float = 2.0
    ):
        
        # left and right boundary lines; these are sized to intersect the upper
        # and lower part of the internal profile, airfoil_shape
        self.left_bound = left_bound
        self.right_bound = right_bound
        self.airfoil_shape = airfoil_shape

        # the corners are the location where the left and right bound lines
        # intersect the airfoil_shape
        self.left_corner_pts: List[App.Vector] = self.get_intersections(self.left_bound)
        self.right_corner_pts: List[App.Vector] = self.get_intersections(self.right_bound)

        # left-hand and right-hand bounds of the upper and lower spline
        left_chamfer_ref: Part.Line = self.left_bound.copy()
        left_chamfer_ref.translate(App.Vector(0.0,chamfer_radius,0.0))
        self.spline_left_bnd_pts: List[App.Vector] = self.get_intersections(left_chamfer_ref)

        right_chamfer_ref: Part.Line = self.right_bound.copy()
        right_chamfer_ref.translate(App.Vector(0.0,-chamfer_radius,0.0))
        self.spline_right_bnd_pts: List[App.Vector] = self.get_intersections(right_chamfer_ref)

        # points to form the left and right line for the hole; these are usually
        # reduced by the chamfer size
        self.left_line_pts: List[App.Vector] = []
        self.right_line_pts: List[App.Vector] = []

        # set of control points for the upper and lower spline
        self.upper_spline_pts: List[App.Vector] = []
        self.lower_spline_pts: List[App.Vector] = []

        self.upper_spline_pts, self.lower_spline_pts = self.get_spline_control_pts()

        self.left_line_pts = self.get_hole_side_coords(self.left_corner_pts, chamfer_radius)
        self.right_line_pts = self.get_hole_side_coords(self.right_corner_pts, chamfer_radius)

    def get_spline_control_pts(self, num_ctl_pts: int = 5) -> tuple[List[App.Vector], List[App.Vector]]:
        """
        Gets spline control points to draw a lightening hole on the upper and
        lower parts of the inner profile
        """
        upper_spline_pts: List[App.Vector] = []
        lower_spline_pts: List[App.Vector] = []

        # so we need to get the left and right boundary Y coordinate to calculate
        # the extents of the bspline region
        y_L: float = self.spline_left_bnd_pts[0].y
        y_R: float =self.spline_right_bnd_pts[0].y
        dY: float = (y_R - y_L)/num_ctl_pts

        ref_line = Part.makeLine(
            App.Vector(0.0, y_L, self.left_bound.Vertexes[0].Point.z), 
            App.Vector(0.0, y_L, self.left_bound.Vertexes[1].Point.z)
            )

        for idx in range(0, num_ctl_pts+1):
            tmp_line = ref_line.copy()
            tmp_line.translate(App.Vector(0.0, dY*idx, 0.0))

            pts = self.get_intersections(tmp_line)
            pts.sort(key=lambda pt: pt.z, reverse=True)

            upper_spline_pts.append(pts[0])
            lower_spline_pts.append(pts[1])

        return upper_spline_pts, lower_spline_pts

    def get_hole_side_coords(self, pts: List[App.Vector], chamfer_radius: float) -> List[App.Vector]:
        """
        Get coordinates for the lines that form the left and right side of the
        lightening hole, leaving space for the chamfer radius on each end of the
        line.  If the line would be shorter than 2 chamfer radii, return the
        midpoint instead
        """
        line_coord_pts: List[App.Vector] = []
        # locate points on the left and right lines that intersect the contours
        # then back up and down using the offset dimension from each side
        p1 = pts[0]
        p2 = pts[1]

        # ensure the distance between p1 and p2 is long enough to fit a non-zero
        # line segment between the chamfer reference points
        if p1.distanceToPoint(p2) > 2* chamfer_radius:
            b1 = App.Vector(p1)
            b2 = App.Vector(p2)
            b1.z = b1.z - chamfer_radius
            b2.z = b2.z + chamfer_radius
            line_coord_pts = [b1,b2]
        # if the distance between p1 and p2 is below the min threshold for the
        # chamfer offset, then just provide a single point midway between p1 and p2
        else:
            b1 = (p1 + p2)/2
            line_coord_pts = [b1]

        line_coord_pts.sort(key=lambda pt: pt.z, reverse=True)

        return line_coord_pts


    def show(self) -> None:
        """
        Draws the hole boundary lines and all control points generated that are
        used to create the lightening hole geometry
        """
        def draw_points(pts: List[App.Vector], grp: App.DocumentObjectGroup, color: tuple = (0.0, 0.5, 0.0)) -> None:

            for pt in pts:
                pf = Draft.make_point(pt)
                pf.ViewObject.PointColor = color
                pf.ViewObject.PointSize = 10.0
                grp.addObject(pf)

        grp: App.DocumentObjectGroup = App.ActiveDocument.addObject("App::DocumentObjectGroup", "lhb-viz")

        left_f = Part.show(self.left_bound)
        left_f.ViewObject.LineColor = utilities.BLUE(0.5)
        grp.addObject(left_f)

        right_f = Part.show(self.right_bound)
        right_f.ViewObject.LineColor = utilities.BLUE(0.5)
        grp.addObject(right_f)

        draw_points(self.left_corner_pts, grp, utilities.RED(1.0))
        draw_points(self.right_corner_pts, grp, utilities.RED(1.0))

        draw_points(self.left_line_pts, grp, utilities.TEAL(1.0))
        draw_points(self.right_line_pts, grp, utilities.TEAL(1.0))

        draw_points(self.upper_spline_pts, grp, utilities.TEAL(1.0))
        draw_points(self.lower_spline_pts, grp, utilities.TEAL(1.0))

        App.ActiveDocument.recompute()
        
    def get_intersections(self, bound_line: Part.Line) -> List[App.Vector]:
        """
        Find the points on a line that intersect the top and bottom of the
        airfoil inner contour.  The points will be returned in sorted order so
        that the upper contour point is first, and lower contour point is last
        """
        intersections: List[App.Vector] = []

        edge: Part.Edge
        for edge in self.airfoil_shape.Edges:
            profile_curve: Part.Curve =  edge.Curve
            profile_curve = profile_curve.trim(*edge.ParameterRange)
            bl_curve: Part.Curve = bound_line.Edges[0].Curve.trim(*bound_line.Edges[0].ParameterRange)
            plane = Part.Plane(App.Vector(0,0,0), App.Vector(1,0,0))
            pts = profile_curve.intersect2d(bl_curve, plane)
            for pt in pts:
                intersections.append(App.Vector(0,-pt[1], pt[0]))

        # sort by highest z coordinate, so that intersections with the top of the
        # airfoil are first, and intersections with the bottom of the airfoil are
        # last
        intersections.sort(key=lambda pt: pt.z, reverse=True)

        return intersections                        

# TODO: This wants to be private, since its not a part of the API
def generate_hole_bounds(
        airfoil_sk: Sketcher.Sketch, 
        interferences: List[Sketcher.Sketch] = [],
        profile_offset: float = 2.0, 
        wall_thickness: float = 2.0,
        target_num_holes: int = 6
    ) -> List[LighteningHoleBounds]:
    """
    Generates a set of lightening hole bound objects, ensuring that such bounds
    do not intersect a set of interferences
    """
    # create the interior hole profile using an offset
    af_inner_profile = airfoil_sk.Shape.makeOffset2D(-profile_offset, join=2)

    orig_placements: List[App.Placement] = []
    for idx in range(0, len(interferences)):
        orig_placements.append(interferences[idx].Placement)
        interferences[idx].Placement = utilities.yz_placement

    # create a list of intervals, and make sure they are sorted by starting
    # value so we can insert properly
    intf_intervals: List[Interval] = [Interval(intf.Shape.BoundBox.YMin, intf.Shape.BoundBox.YMax) for intf in interferences]
    intf_intervals.sort(key=lambda intv: intv.start)

    # first, calculate the nominal max hole width; this will be used to calculate
    # the number of holes in each of the space intervals between interferences
    bbox = af_inner_profile.BoundBox
    max_hole_width: float = (bbox.YLength - wall_thickness*(target_num_holes+1))/target_num_holes

    profile_interval: Interval = Interval(bbox.YMin, bbox.YMax)
    hole_intervals: List[Interval] = []

    # second, calculate the valid set of intervals where we will generate holes
    current_interval: Interval = Interval(profile_interval.start, 0.0)
    for intf_interval in intf_intervals:
        current_interval.end = intf_interval.start
        hole_intervals.append(current_interval)
        current_interval = Interval(intf_interval.end, 0.0)
    current_interval.end = profile_interval.end
    hole_intervals.append(current_interval)

    # TODO: We often get extra wall spacing in areas where the line height is
    #       zero, but, basically, we ought to be able to compensate for this,
    #       since we can detect it when we generate the bound list
    def generate_bounds(interval: Interval) -> List[LighteningHoleBounds]:
        num_segments: int = math.ceil(interval.length()/max_hole_width)

        if num_segments == 0:
            return []

        segment_length: float = interval.length()/num_segments
  
        ref_lines: List[Part.Edge] = []
        # TODO: replace with linspace
        for y_coord in numpy.arange(interval.start, interval.end + segment_length, segment_length):
            top_pt = (0, y_coord, bbox.ZMax + 1)
            bottom_pt = (0, y_coord, bbox.ZMin - 1)
            
            ref_line = Part.makeLine(top_pt, bottom_pt)
            ref_lines.append(ref_line)
        # NOTE: LighteningHoleBounds actually stores control points
        bounds: List[LighteningHoleBounds] = []
        for ref_idx in range(0, num_segments):
            left = ref_lines[ref_idx].copy()

            left_offset = wall_thickness/2
            right_offset = -left_offset
            if ref_idx == 0:
                left_offset = wall_thickness
            if ref_idx == num_segments-1:
                right_offset = -wall_thickness
                pass
            
            left.translate(App.Vector(0.0, left_offset, 0.0))

            right = ref_lines[ref_idx+1].copy()
            right.translate(App.Vector(0.0, right_offset, 0.0))
        
            new_bound: LighteningHoleBounds = LighteningHoleBounds(af_inner_profile, left, right)
            bounds.append(new_bound)
        return bounds

    bnds: List[LighteningHoleBounds] = []
    for intv in hole_intervals:
        bnds.extend(generate_bounds(intv))

    for idx in range(0, len(interferences)):
        interferences[idx].Placement = orig_placements[idx]

    return bnds

# TODO: passing in entire sketches without actually performing any operations on
#       them is probably a bit heavy-weight, we should possibly lighten up these
#       args
def create_lightening_hole_sketch(
    airfoil_sk: Sketcher.Sketch,
    interferences: List[Sketcher.Sketch] = []
) -> Sketcher.Sketch:
    """
    Creates a sketch which contains geometry for a set of lightening holes for
    an airfoil.  A set of interferences may be provided, which represent areas
    withing the airfoil section where lightening holes should not be generated
    """

    af_orig_placement = airfoil_sk.Placement
    airfoil_sk.Placement = utilities.yz_placement

    bnds: List[LighteningHoleBounds] = \
            generate_hole_bounds(
                airfoil_sk,
                interferences
            )

    # create the sketch on the YZ plane
    sk: Sketcher.Sketch = App.ActiveDocument.addObject("Sketcher::SketchObject", "lightening-holes")
    sk.Placement = utilities.yz_placement    

    constraints: List[Sketcher.Constraint] = []

    # we need the placement of the sketch to convert the reference points in
    # global coordinate space into sketch coordinate space
    inv_placement = sk.Placement.inverse()

    def draw_line(pts: List[App.Vector]) -> tuple[int, Part.LineSegment]:
        if len(pts) > 1:
            p1 = inv_placement.multVec(pts[0])
            p2 = inv_placement.multVec(pts[1])
            line = Part.LineSegment(p1, p2)
            id = sk.addGeometry(line, False)
            constraints.append(Sketcher.Constraint("Block", id))
            return id, line
        return -1, None

    def draw_spline(pts: List[App.Vector]) -> tuple[int, Part.BSplineCurve]:
        xf_coords: List[App.Vector] = [inv_placement.multVec(pt) for pt in pts]
        bsp = Part.BSplineCurve()
        bsp.interpolate(xf_coords)
        
        id = sk.addGeometry(bsp)
        constraints.append(Sketcher.Constraint("Block", id))
        return id, bsp
    
    class ChamferSide(Enum):
        LEFT = 1
        RIGHT = 2

    # Note that curve must be an instance of Part.BoundedCurve; unfortunately
    # if you provide the abstract base class, the interpreter rejects it since
    # the abstract base probably has type hints, but isn't exposed directly to
    # the interpreter
    def get_nearest_end_to_pt(target_pt: App.Vector, curve):
        pt_dist = [
            (1, curve.StartPoint.distanceToPoint(target_pt)),
            (2, curve.EndPoint.distanceToPoint(target_pt))
        ]
        pt_dist.sort(key=lambda x : x[1])
        return pt_dist[0][0]

    def draw_chamfers(
            chamfer_pts: List[App.Vector], 
            line: Part.LineSegment, line_id: int,
            upper_bsp: Part.BSplineCurve, upper_bsp_id: int,
            lower_bsp: Part.BSplineCurve, lower_bsp_id: int,
            chamfer_side: ChamferSide) -> None:
        """
        Draw circular arcs that join a side boundary line with one of the bspline
        interior boundary sections
        """
        ch_pts: List[App.Vector] = [inv_placement.multVec(pt) for pt in chamfer_pts]
              

        top_arc_start: float; top_arc_end: float
        bottom_arc_start: float; bottom_arc_end: float

        # TODO: Really need to understand how arc start stop work here
        match chamfer_side:
            case ChamferSide.LEFT:
                top_arc_start = numpy.deg2rad(90)
                top_arc_end = numpy.deg2rad(180)
                bottom_arc_start = numpy.deg2rad(-180)
                bottom_arc_end = numpy.deg2rad(-90)
                
            case ChamferSide.RIGHT:
                top_arc_start = numpy.deg2rad(0)
                top_arc_end = numpy.deg2rad(90)
                bottom_arc_start = numpy.deg2rad(-90)
                bottom_arc_end = numpy.deg2rad(0)

        # if the space was too small to form a line with the target offset, then
        # we create two arcs and join them directly, without a line in the middle
        if line is None:
            p1 = ch_pts[0]
            p2 = ch_pts[1]
            mid_pt = (p1+p2)/2
            radius = p1.distanceToPoint(p2)/2  

            arc_top = Part.ArcOfCircle(
                Part.Circle(
                    mid_pt,
                    App.Vector(0,0,1),
                    radius
                ),
                top_arc_start,
                top_arc_end
            )

            id_arc_top = sk.addGeometry(arc_top, False)
            bsp_pt_id = 1 if upper_bsp.StartPoint == p1 else 2
            arc_pt_id = get_nearest_end_to_pt(p1, arc_top)
            top_arc_other_id = 1 if arc_pt_id == 2  else 2

            # TODO: We really need to know why the order that we apply the 
            #       constraints matters so much; if you do the arc first, 
            #       everything is rejected as over-constrained 
            constraints.append(Sketcher.Constraint("Tangent", upper_bsp_id, bsp_pt_id, id_arc_top, arc_pt_id))

            arc_bottom = Part.ArcOfCircle(
                Part.Circle(
                    mid_pt,
                    App.Vector(0,0,1),
                    radius
                ),
                bottom_arc_start,
                bottom_arc_end
            )

            id_arc_bottom = sk.addGeometry(arc_bottom, False)
            bsp_pt_id = 1 if lower_bsp.StartPoint == p2 else 2
            arc_pt_id = get_nearest_end_to_pt(p2, arc_bottom)
            bottom_arc_other_id = 1 if arc_pt_id == 2 else 2

            constraints.append(Sketcher.Constraint("Tangent", lower_bsp_id, bsp_pt_id, id_arc_bottom, arc_pt_id))
            constraints.append(Sketcher.Constraint("Tangent", id_arc_top, top_arc_other_id, id_arc_bottom, bottom_arc_other_id))
        # we have a line in the middle, so create arcs that join the line on each
        # end to the top and bottom spline sections
        else:
            arc_radius = numpy.abs(ch_pts[0].x - line.StartPoint.x) 
            top_arc_center = App.Vector(ch_pts[0].x, line.StartPoint.y, 0.0)

            arc_top = Part.ArcOfCircle(
                Part.Circle(
                    top_arc_center,
                    App.Vector(0,0,1),
                    arc_radius
                ),
                top_arc_start,
                top_arc_end
            )

            id_arc_top = sk.addGeometry(arc_top, False)
            bsp_pt_id = 1 if upper_bsp.StartPoint == ch_pts[0] else 2
            arc_pt_id = get_nearest_end_to_pt(ch_pts[0], arc_top)
            top_arc_other_id = 1 if arc_pt_id == 2  else 2
            
            constraints.append(Sketcher.Constraint("Tangent", upper_bsp_id, bsp_pt_id, id_arc_top, arc_pt_id))
            constraints.append(Sketcher.Constraint("Tangent", id_arc_top, top_arc_other_id, line_id, 1))

            btm_arc_ctr = App.Vector(ch_pts[1].x, line.EndPoint.y, 0.0)

            arc_bottom = Part.ArcOfCircle(
                Part.Circle(
                    btm_arc_ctr,
                    App.Vector(0,0,1),
                    arc_radius
                ),
                bottom_arc_start,
                bottom_arc_end
            )

            id_arc_bottom = sk.addGeometry(arc_bottom, False)
            bsp_pt_id = 1 if lower_bsp.StartPoint == ch_pts[1] else 2
            arc_pt_id = get_nearest_end_to_pt(ch_pts[1], arc_bottom)
            btm_arc_other_id = 1 if arc_pt_id == 2  else 2

            constraints.append(Sketcher.Constraint("Tangent", lower_bsp_id, bsp_pt_id, id_arc_bottom, arc_pt_id))
            constraints.append(Sketcher.Constraint("Tangent", id_arc_bottom, btm_arc_other_id, line_id, 2))

    for bnd in bnds:
        left_line_id, left_line = draw_line(bnd.left_line_pts)
        right_line_id, right_line = draw_line(bnd.right_line_pts)
        upper_bsp_id, upper_bsp = draw_spline(bnd.upper_spline_pts)
        lower_bsp_id, lower_bsp = draw_spline(bnd.lower_spline_pts)

        draw_chamfers(bnd.spline_left_bnd_pts, left_line, left_line_id, upper_bsp, upper_bsp_id, lower_bsp, lower_bsp_id, ChamferSide.LEFT)
        draw_chamfers(bnd.spline_right_bnd_pts, right_line, right_line_id, upper_bsp, upper_bsp_id, lower_bsp, lower_bsp_id, ChamferSide.RIGHT)

    sk.addConstraint(constraints)

    airfoil_sk.Placement = af_orig_placement
    sk.Placement = airfoil_sk.Placement

    return sk



    
            