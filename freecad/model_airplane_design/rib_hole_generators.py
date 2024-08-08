from . import utilities
from . import airfoil
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

class HoleGeneratorType(Enum):
    NONE = 1
    ROUNDED_TRAPEZOID = 2

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
            bbox: App.BoundBox,
            standoff: float
        ) -> None:
        
        if standoff < 0.0:
            raise ValueError("Hole exclusion standoff must not be negative")

        self.bbox = bbox
        self.standoff = standoff
    
    def get_excluded_interval(self) -> Interval:
        return Interval(self.bbox.XMin-self.standoff, self.bbox.XMax+self.standoff)

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
        
        # TODO: removing join=2 made this work, but why?
        self.af_inner_profile = self.airfoil_sk.Shape.makeOffset2D(-self.inset_width)
        self.af_inner_profile.tessellate(0.01)

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

            # if the interference is past the end of the internal profile, then
            # skip it
            if intf_interval.start >= self.af_inner_profile.BoundBox.XMax:
                continue

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
            pts = airfoil.get_intersections(x_coord, rib_inner_profile)
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

        corner_pts: List[App.Vector] = airfoil.get_intersections(x_position, profile)
        
        # print("------------------")
        # print("range: " + str(profile.BoundBox.XMin) + " - " + str(profile.BoundBox.XMax))
        # print("x_position=" + str(x_position))

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
    
        # create the sketch on the XY plane
        sk: Sketcher.Sketch = App.ActiveDocument.addObject("Sketcher::SketchObject", "lightening-holes")
        sk.Placement = utilities.xy_placement    

        constraints: List[Sketcher.Constraint] = []

        def draw_line(endcap: RoundedTrapezoidEndcap) -> tuple[int, Part.LineSegment]:
            if endcap.is_line():
                line = Part.LineSegment(endcap.pts[0], endcap.pts[1])
                id = sk.addGeometry(line, False)
                constraints.append(Sketcher.Constraint("Vertical", id))
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
            # we have a line in the middle, so create arcs that join the line on
            # each end to the top and bottom spline sections
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
        

        with utilities.TempDocObjectHelper() as tmp_obj_helper:
            
            tmp_obj_helper.addObject(sk)
        
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

            # if we got here without an exception, then it is safe to remove
            # the object from the helper so we can return it
            tmp_obj_helper.removeObject(sk)

        return sk
    
class HoleGeneratorFactory():
    @staticmethod
    def create_generator(generator_type: str, chord_param: float) -> HoleGenerator:
        match generator_type:
            case HoleGeneratorType.NONE.name:
                return None
            case HoleGeneratorType.ROUNDED_TRAPEZOID.name:
                scale_factor = chord_param/175.0
                chamfer_rad: float = scale_factor*2.0
                max_hole_width: float = scale_factor*(100.0/6)
                min_hole_spacing: float = scale_factor*2.0
                return RoundedTrapezoidHoleGenerator(chamfer_rad, max_hole_width, min_hole_spacing)
            case _:
                print("HoleGeneratorFactory.create_generator: unknown generator")
                return None