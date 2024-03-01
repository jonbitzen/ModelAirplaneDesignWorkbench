import FreeCAD as App
import Part
import Draft
import Sketcher
import numpy
from typing import List

class LighteningHoleBounds:
    def __init__(self, airfoil_shape: Part.Shape, left: Part.Line, right: Part.Line):
        self.left = left
        self.right = right
        self.airfoil_shape = airfoil_shape

        # boundary points
        self.left_bnd_pts: List[App.Vector] = []
        self.right_bnd_pts: List[App.Vector] = []

        self.left_chamfer_pts: List[App.Vector] = []
        self.right_chamfer_pts: List[App.Vector] = []

        self.left_l_chamfer_pts: List[App.Vector] = []
        self.right_l_chamfer_pts: List[App.Vector] = []

        self.upper_spline_pts: List[App.Vector] = []
        self.lower_spline_pts: List[App.Vector] = []


    def get_spline_control_pts(self, num_ctl_pts: int = 5) -> tuple[List[App.Vector], List[App.Vector]]:
        
        upper_spline_pts: List[App.Vector] = []
        lower_spline_pts: List[App.Vector] = []

        # so we need to get the left and right boundary Y coordinate to calculate
        # the extents of the bspline region
        y_L: float = self.left_chamfer_pts[0].y
        y_R: float =self.right_chamfer_pts[0].y
        dY: float = (y_R - y_L)/num_ctl_pts

        ref_line = Part.makeLine(
            App.Vector(0.0, y_L, self.left.Vertexes[0].Point.z), 
            App.Vector(0.0, y_L, self.left.Vertexes[1].Point.z)
            )

        for idx in range(0, num_ctl_pts+1):
            tmp_line = ref_line.copy()
            tmp_line.translate(App.Vector(0.0, dY*idx, 0.0))

            pts = self.get_intersections(tmp_line)
            pts.sort(key=lambda pt: pt.z, reverse=True)

            self.upper_spline_pts.append(pts[0])
            self.lower_spline_pts.append(pts[1])

        return (self.upper_spline_pts, self.lower_spline_pts)

    # locate the points that will be used to form the chamfer coordinate
    def get_chamfer_corners(self, corner_offset_mm: float = 2) -> List[App.Vector]:

        # Get reference points on the left and right boundary line
        self.left_bnd_pts = self.get_intersections(self.left)
        self.right_bnd_pts = self.get_intersections(self.right)

        # create lines that are offset from the left and right, locate where they
        # intersect with the top and bottom curve, these will join with the top
        # bspline for the interior contour
        left_chamfer_ref: Part.Line = self.left.copy()
        left_chamfer_ref.translate(App.Vector(0.0,corner_offset_mm,0.0))
        self.left_chamfer_pts = self.get_intersections(left_chamfer_ref)

        right_chamfer_ref: Part.Line = self.right.copy()
        right_chamfer_ref.translate(App.Vector(0.0,-corner_offset_mm,0.0))
        self.right_chamfer_pts = self.get_intersections(right_chamfer_ref)

        self.left_l_chamfer_pts = self.get_chamfer_line_pts(self.left_bnd_pts, corner_offset_mm)
        self.right_l_chamfer_pts = self.get_chamfer_line_pts(self.right_bnd_pts, corner_offset_mm)

    def get_chamfer_line_pts(self, pts: List[App.Vector], offset_mm: float) -> List[App.Vector]:

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


    def draw_points(self, pts: List[App.Vector], grp: App.DocumentObjectGroup, color: tuple = (0.0, 0.5, 0.0)) -> None:
        for pt in pts:
            pf = Draft.make_point(pt)
            pf.ViewObject.PointColor = color
            pf.ViewObject.PointSize = 10.0
            grp.addObject(pf)

    def show(self) -> None:

        grp: App.DocumentObjectGroup = create_group("lhb-viz")

        left_f = Part.show(self.left)
        left_f.ViewObject.LineColor = BLUE(0.5)
        grp.addObject(left_f)

        right_f = Part.show(self.right)
        right_f.ViewObject.LineColor = BLUE(0.5)
        grp.addObject(right_f)

        self.draw_points(self.left_bnd_pts, grp)
        self.draw_points(self.right_bnd_pts, grp)

        self.draw_points(self.left_chamfer_pts, grp)
        self.draw_points(self.right_chamfer_pts, grp)

        self.draw_points(self.left_l_chamfer_pts, grp, TEAL(1.0))
        self.draw_points(self.right_l_chamfer_pts, grp, TEAL(1.0))

        self.draw_points(self.upper_spline_pts, grp, TEAL(1.0))
        self.draw_points(self.lower_spline_pts, grp, TEAL(1.0))

        App.ActiveDocument.recompute()
        

    # locate the points on the top and bottom of the inner contour that will be
    # the knot locations to form the top and bottom bspline
    def get_inner_contour_control_points(self):
        pass

    def generate_lightening_holes(self):
        pass

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

def create_airfoil_inner_profile(
        airfoil_sk : Sketcher.Sketch, 
        offset_mm: float = 2.0, 
        structure_thk_mm: float = 2.0
):
    print("creating inner profile for airfoil using airfoil sketch " + airfoil_sk.Label)
    offset = airfoil_sk.Shape.makeOffset2D(-offset_mm, join=2)
    f = Part.show(offset, "airfoil-offset")
    f.ViewObject.LineColor = RED(0.3)
    bbox = offset.BoundBox
    num_segments = 5
    segment_length = bbox.YLength/num_segments
    
    # create a set of reference lines the demarcate the rib without thinking explictly of structure
    ref_lines: List[Part.Edge] = []
    for y_coord in numpy.arange(bbox.YMin, bbox.YMax + segment_length, segment_length):
        top_pt = (0, y_coord, bbox.ZMax + 1)
        bottom_pt = (0, y_coord, bbox.ZMin - 1)
        
        ref_line = Part.makeLine(top_pt, bottom_pt)
        ref_lines.append(ref_line)
        
    # create a set of lightening hole bound areas
    bounds: List[LighteningHoleBounds] = []
    for ref_idx in range(0, num_segments):
        left = ref_lines[ref_idx].copy()
        if ref_idx != 0:
            left.translate(App.Vector(0.0, structure_thk_mm/2, 0.0))

        right = ref_lines[ref_idx+1].copy()
        right.translate(App.Vector(0.0, -structure_thk_mm/2, 0.0))
    
        new_bound: LighteningHoleBounds = LighteningHoleBounds(offset, left, right)
        bounds.append(new_bound)

    return bounds

def create_lightening_hole(int_profile: Part.Shape, left_bound: Part.Line, right_bound: Part.Line):
    print("create_lightening_hole")


