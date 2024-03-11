import FreeCAD as App
import Part
import Draft
import Sketcher
import math
import numpy
from enum import Enum
from typing import List

# Helper methods, used throughout
def create_group(name: str, doc: App.Document = App.ActiveDocument) -> App.DocumentObjectGroup:
    return doc.addObject("App::DocumentObjectGroup", name)

def TEAL(level: float) -> tuple:
    level = numpy.clip(level, 0.1, 1.0)
    return (0.0, level, level)

def BLACK() -> tuple:
    return (0.0, 0.0, 0.0)

def GREY(grey: float = 0.5):
    return (grey, grey, grey)

def RED(red: float = 1.0) -> tuple:
    return (red, 0.0, 0.0)

def GREEN(green: float = 1.0) -> tuple:
    return (0.0, green, 0.0)

def BLUE(blue: float = 1.0) -> tuple:
    return (0.0, 0.0, blue)

class LighteningHoleBounds:
    """
    Stores a set of boundaries and control points to generate lightening holes
    for an airfoil
    """
    def __init__(self, airfoil_shape: Part.Shape, left_bound: Part.Line, right_bound: Part.Line):
        self.left_bound = left_bound
        self.right_bound = right_bound
        self.airfoil_shape = airfoil_shape

        # TODO: improve and document these names before I forget what they all
        #       mean; probably need to make a general explanation of how LHB works
        # boundary points
        self.left_bnd_pts: List[App.Vector] = []
        self.right_bnd_pts: List[App.Vector] = []

        self.left_chamfer_pts: List[App.Vector] = []
        self.right_chamfer_pts: List[App.Vector] = []

        # points on the bound line
        self.left_l_chamfer_pts: List[App.Vector] = []
        self.right_l_chamfer_pts: List[App.Vector] = []

        # set of control points for the upper and lower spline
        self.upper_spline_pts: List[App.Vector] = []
        self.lower_spline_pts: List[App.Vector] = []

        self.get_chamfer_corners()
        self.upper_spline_pts, self.lower_spline_pts = self.get_spline_control_pts()

        self.offset = 2

        self.left_l_chamfer_pts = self.get_chamfer_line_pts(self.left_bnd_pts, self.offset)
        self.right_l_chamfer_pts = self.get_chamfer_line_pts(self.right_bnd_pts, self.offset)

    def get_spline_control_pts(self, num_ctl_pts: int = 5) -> tuple[List[App.Vector], List[App.Vector]]:
        """
        Gets spline control points to draw a lightening hole on the upper and
        lower parts of the inner profile
        """
        upper_spline_pts: List[App.Vector] = []
        lower_spline_pts: List[App.Vector] = []

        # so we need to get the left and right boundary Y coordinate to calculate
        # the extents of the bspline region
        y_L: float = self.left_chamfer_pts[0].y
        y_R: float =self.right_chamfer_pts[0].y
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

    # This really finds the hole edge bounds, subtracting the chamfer radius
    # from each set of bounds.  So, we find the upper and lower points of the
    # left and right line-bound, and the right and left bounds for the upper
    # spline and lower spline
    # TODO: It might be more intuitive to arrange the spline bounds as left-right
    #       pairs for the upper and lower splines, rather than the current
    #       unintuitive arrangement.  In future we might want to gather info that
    #       way in the first place to keep the algo simple to understand
    # TODO: we should refactor corner_offset to chamfer_radius
    
    def get_chamfer_corners(self, corner_offset: float = 2) -> List[App.Vector]:
        """
        
        """
        # Get reference points on the left and right boundary line
        self.left_bnd_pts = self.get_intersections(self.left_bound)
        self.right_bnd_pts = self.get_intersections(self.right_bound)


        # TODO: We should move this off someplace else to whatever point we
        #       start creating the hole's edge geometry, since we're essentially
        #       creating the bspline bounds here, but not the correspond left
        #       and right edge bounds
        # create lines that are offset from the left and right, locate where they
        # intersect with the top and bottom curve, these will join with the top
        # bspline for the interior contour
        left_chamfer_ref: Part.Line = self.left_bound.copy()
        left_chamfer_ref.translate(App.Vector(0.0,corner_offset,0.0))
        self.left_chamfer_pts = self.get_intersections(left_chamfer_ref)

        right_chamfer_ref: Part.Line = self.right_bound.copy()
        right_chamfer_ref.translate(App.Vector(0.0,-corner_offset,0.0))
        self.right_chamfer_pts = self.get_intersections(right_chamfer_ref)

    def get_chamfer_line_pts(self, pts: List[App.Vector], offset_mm: float) -> List[App.Vector]:
        """
        Find points on the left and right bounding lines that are backed from
        the intersection with the top and bottom part of the interior profile a
        fixed distance.  If the line so formed would be of zero length, return a
        single point midway between the intersection points
        """
        chm_pts: List[App.Vector] = []
        # locate points on the left and right lines that intersect the contours
        # then back up and down using the offset dimension from each side
        p1 = pts[0]
        p2 = pts[1]

        # ensure the distance between p1 and p2 is long enough to fit a non-zero
        # line segment between the chamfer reference points
        if p1.distanceToPoint(p2) > 2* offset_mm:
            b1 = App.Vector(p1)
            b2 = App.Vector(p2)
            b1.z = b1.z - offset_mm
            b2.z = b2.z + offset_mm
            chm_pts = [b1,b2]
        # if the distance between p1 and p2 is below the min threshold for the
        # chamfer offset, then just provide a single point midway between p1 and p2
        else:
            b1 = (p1 + p2)/2
            chm_pts = [b1]

        chm_pts.sort(key=lambda pt: pt.z, reverse=True)

        return chm_pts


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

        grp: App.DocumentObjectGroup = create_group("lhb-viz")

        left_f = Part.show(self.left_bound)
        left_f.ViewObject.LineColor = BLUE(0.5)
        grp.addObject(left_f)

        right_f = Part.show(self.right_bound)
        right_f.ViewObject.LineColor = BLUE(0.5)
        grp.addObject(right_f)

        # TODO: maybe rename these left boundary corners and right boundary
        #       corners
        draw_points(self.left_bnd_pts, grp, RED(1.0))
        draw_points(self.right_bnd_pts, grp, RED(1.0))

        # TODO: these are the right and left bounds for the spline profile
        #       parts.  I have a feeling we can build this into the part where
        #       calculate the splines, and eliminate this member entirely
        # draw_points(self.left_chamfer_pts, grp)
        # draw_points(self.right_chamfer_pts, grp)

        # TODO: these are the endpoints for the left and right hole line; we
        #       should just rename them to that
        draw_points(self.left_l_chamfer_pts, grp, TEAL(1.0))
        draw_points(self.right_l_chamfer_pts, grp, TEAL(1.0))

        draw_points(self.upper_spline_pts, grp, TEAL(1.0))
        draw_points(self.lower_spline_pts, grp, TEAL(1.0))

        App.ActiveDocument.recompute()
        
    # find the points on a line that intersect the top and bottom of the airfoil
    # inner contour
    def get_intersections(self, bound_line: Part.Line) -> List[App.Vector]:
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

class Interval:
    """
    Stores starting and ending interval values used to compute lightening hole
    placement
    """
    start: float
    end: float
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
        segment_length: float = interval.length()/num_segments
  
        ref_lines: List[Part.Edge] = []
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
    bnds: List[LighteningHoleBounds] = \
            generate_hole_bounds(
                airfoil_sk,
                interferences
            )

    # create the sketch on the YZ plane
    sk: Sketcher.Sketch = App.ActiveDocument.addObject("Sketcher::SketchObject", "lightening-holes")
    sk.Placement = App.Placement(
        App.Vector(0.000000, 0.000000, 0.000000), 
        App.Rotation(0.500000, 0.500000, 0.500000, 0.500000)
    )

    # we need the placement of the sketch to convert the reference points in
    # global coordinate space into sketch coordinate space
    inv_placement = sk.Placement.inverse()

    def draw_line(pts: List[App.Vector]) -> tuple[int, Part.LineSegment]:
        if len(pts) > 1:
            p1 = inv_placement.multVec(pts[0])
            p2 = inv_placement.multVec(pts[1])
            line = Part.LineSegment(p1, p2)
            id = sk.addGeometry(line, False)
            sk.addConstraint(Sketcher.Constraint("Block", id))
            return id, line
        return -1, None

    def draw_spline(pts: List[App.Vector]) -> tuple[int, Part.BSplineCurve]:
        xf_coords: List[App.Vector] = [inv_placement.multVec(pt) for pt in pts]
        bsp = Part.BSplineCurve()
        bsp.interpolate(xf_coords)
        
        id = sk.addGeometry(bsp)
        sk.addConstraint(Sketcher.Constraint("Block", id))
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
            sk.addConstraint(Sketcher.Constraint("Tangent", upper_bsp_id, bsp_pt_id, id_arc_top, arc_pt_id))

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

            sk.addConstraint(Sketcher.Constraint("Tangent", lower_bsp_id, bsp_pt_id, id_arc_bottom, arc_pt_id))
            sk.addConstraint(Sketcher.Constraint("Tangent", id_arc_top, top_arc_other_id, id_arc_bottom, bottom_arc_other_id))
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
            
            sk.addConstraint(Sketcher.Constraint("Tangent", upper_bsp_id, bsp_pt_id, id_arc_top, arc_pt_id))
            sk.addConstraint(Sketcher.Constraint("Tangent", id_arc_top, top_arc_other_id, line_id, 1))

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

            sk.addConstraint(Sketcher.Constraint("Tangent", lower_bsp_id, bsp_pt_id, id_arc_bottom, arc_pt_id))
            sk.addConstraint(Sketcher.Constraint("Tangent", id_arc_bottom, btm_arc_other_id, line_id, 2))

    for bnd in bnds:
        left_line_id, left_line = draw_line(bnd.left_l_chamfer_pts)
        right_line_id, right_line = draw_line(bnd.right_l_chamfer_pts)
        upper_bsp_id, upper_bsp = draw_spline(bnd.upper_spline_pts)
        lower_bsp_id, lower_bsp = draw_spline(bnd.lower_spline_pts)

        draw_chamfers(bnd.left_chamfer_pts, left_line, left_line_id, upper_bsp, upper_bsp_id, lower_bsp, lower_bsp_id, ChamferSide.LEFT)
        draw_chamfers(bnd.right_chamfer_pts, right_line, right_line_id, upper_bsp, upper_bsp_id, lower_bsp, lower_bsp_id, ChamferSide.RIGHT)

    return sk



    
            