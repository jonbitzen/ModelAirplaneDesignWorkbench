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
from typing import List


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
    def __init__(
            self, 
            af_inner_contour: Part.Shape, 
            interval: Interval) -> None:
        if af_inner_contour is None:
            raise ValueError("Airfoil inner contour must not be None")
        
        if interval is None:
            raise ValueError("Hole bound interval must not be None")

        self.af_inner_counter = af_inner_contour
        self.interval = interval        

class HoleBoundGenerator():
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
        return self.hole_intervals
    
    def show(self) -> None:
        
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
        for doc_obj in self.show_items:
            App.ActiveDocument.removeObject(doc_obj.Name)

class HoleGenerator(ABC):
    def generate_sketch(self, interval: HoleBoundRegion) -> Sketcher.Sketch:
        pass

class RoundedTrapezoidSidePoints():
    def __init__(self, pts: List[App.Vector]) -> None:
        self.pts = pts

    def is_line(self) -> bool:
        return len(self.pts) == 2
    
    def is_point(self) -> bool:
        return len(self.pts) == 1

class RoundedTrapezoidControlPoints():
    def __init__(
            self,
            airfoil_profile: Part.Shape,
            interval: Interval,
            chamfer_rad: float
        ) -> None:
        
        self.show_items: List[App.DocumentObject] = []

        # control points for the upper and lower surface of the airfoil inner profile
        self.upper_spline_pts: List[App.Vector] = []
        self.lower_spline_pts: List[App.Vector] = []

        self.left_side_pts: RoundedTrapezoidSidePoints = \
            self.__get_side_coords(interval.start, airfoil_profile, chamfer_rad)
        self.right_side_pts: RoundedTrapezoidSidePoints = \
            self.__get_side_coords(interval.end, airfoil_profile, chamfer_rad)


        # generate control points for the upper and lower splines, which will
        # hug the airfoil interior profile (airfoil_profile)
        # TODO: we need to make sure that right_bound-left_bound is greater than
        #       zero at a minimum, and probably also that it's greater than some
        #       practical minimum length as well, or the algo below will explode
        spl_left_bound = interval.start+chamfer_rad
        spl_right_bound = interval.end-chamfer_rad
        num_ctl_pts = 5
        for x_coord in numpy.linspace(spl_left_bound, spl_right_bound, num_ctl_pts):
            pts = self.__get_intersections(x_coord, airfoil_profile)
            self.upper_spline_pts.append(pts[0])
            self.lower_spline_pts.append(pts[1])

    def __get_side_coords(
            self, 
            x_position: float, 
            profile: Part.Shape,
            chamfer_rad: float
        ) -> RoundedTrapezoidSidePoints:

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

        return RoundedTrapezoidSidePoints(line_coord_pts)

    def __get_intersections(self, x_position: float, profile: Part.Shape) -> List[App.Vector]:
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
    
    def show(self) -> None:
        # pts: List[Draft.Point] = utilities.draw_points(self.upper_spline_pts)
        self.show_items.extend(utilities.draw_points(self.upper_spline_pts))
        self.show_items.extend(utilities.draw_points(self.lower_spline_pts))
        self.show_items.extend(utilities.draw_points(self.left_side_pts.pts))
        self.show_items.extend(utilities.draw_points(self.right_side_pts.pts))

    def unshow(self) -> None:
        for doc_obj in self.show_items:
            App.ActiveDocument.removeObject(doc_obj.Name)

class RoundedTrapezoidHoleGenerator(HoleGenerator):
    def __init__(
            self,
            max_chamfer_rad: float,
            max_hole_length: float,
            min_hole_spacing: float
        ) -> None:
        super().__init__()
        
        # self.max_chamfer_rad = max_chamfer_rad
        # self.max_hole_length = max_hole_length
        # self.min_hole_spacing = min_hole_spacing


        



    def generate_sketch(self, interval: HoleBoundRegion) -> Sketcher.Sketch:
        return super().generate_sketch(interval)

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



    
            