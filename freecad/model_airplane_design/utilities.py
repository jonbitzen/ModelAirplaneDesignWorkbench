import FreeCAD as App
import Part
import Draft
import Sketcher
import numpy
from typing import List

class LighteningHoleBounds:
    def __init__(self, shape: Part.Shape, left: Part.Line, right: Part.Line):
        self.left = left
        self.right = right
        self.shape = shape

    def get_right_intersections(self):
        self.get_intersections(self.right)

    def get_left_intersections(self):
        self.get_intersections(self.left)

    def get_intersections(self, bound_line: Part.Line) -> List[Draft.Point]:
        intersections: List[Draft.Point]

        l_f = Part.show(bound_line)
        l_f.ViewObject.LineColor = BLUE(0.3)

        edge: Part.Edge
        for edge in self.shape.Edges:
            profile_curve: Part.Curve =  edge.Curve
            profile_curve = profile_curve.trim(*edge.ParameterRange)
            bl_curve: Part.Curve = bound_line.Edges[0].Curve.trim(*bound_line.Edges[0].ParameterRange)
            plane = Part.Plane(App.Vector(0,0,0), App.Vector(1,0,0))
            pts = profile_curve.intersect2d(bl_curve, plane)
            for pt in pts:
                Draft.make_point(0, -pt[1], pt[0], color=GREEN(0.7), point_size=10.0)

        App.ActiveDocument.recompute()

    def show_bounds(self) -> None:

        left_f = Part.show(self.left)
        left_f.ViewObject.LineColor = BLUE(0.5)

        right_f = Part.show(self.right)
        right_f.ViewObject.LineColor = BLUE(0.5)

def create_group(name: str, doc: App.Document = App.ActiveDocument) -> App.DocumentObjectGroup:
    return doc.addObject("App::DocumentObjectGroup", name)

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
    ref_grp = create_group("ref_lines")
    # ref_lines = [Part.Edge]
    ref_lines: List[Part.Edge] = []
    for y_coord in numpy.arange(bbox.YMin, bbox.YMax + segment_length, segment_length):
        top_pt = (0, y_coord, bbox.ZMax + 1)
        bottom_pt = (0, y_coord, bbox.ZMin - 1)
        
        ref_line = Part.makeLine(top_pt, bottom_pt)
        ref_lines.append(ref_line)

        # ref_feature = Part.show(ref_line, "ref_line")
        # ref_feature.ViewObject.LineColor = GREY(0.5)
        # ref_grp.addObject(ref_feature)
        
    # create a set of lightening hole bound areas
    bounds: List[LighteningHoleBounds] = []
    for ref_idx in range(0, num_segments):
        left = ref_lines[ref_idx].copy()
        if ref_idx != 0:
            left.translate(App.Vector(0.0, structure_thk_mm/2, 0.0))
        
        # sl = Part.show(left, "str_line_L")
        # sl.ViewObject.LineColor = BLUE(0.3)
        # ref_grp.addObject(sl)

        right = ref_lines[ref_idx+1].copy()
        right.translate(App.Vector(0.0, -structure_thk_mm/2, 0.0))
        
        # sr = Part.show(right, "str_line_R")
        # sr.ViewObject.LineColor = BLUE(0.3)
        # ref_grp.addObject(sr)

        new_bound: LighteningHoleBounds = LighteningHoleBounds(offset, left, right)
        bounds.append(new_bound)

        bounds[0].get_left_intersections()
        bounds[0].get_right_intersections()

    return bounds[0]

    # # run over each the bounding areas, and locate control points to generate the
    # # hole geometry
    # for bounds_set in bounds:
    #     # iterate over every edge in the shape, see where it intersects the curent bounds
    #     print("---------------------------")
    #     for edge_idx in range(0, len(bounds_set.shape.Edges)):
    #         curve: Part.Curve = bounds_set.shape.Edges[edge_idx].Curve

    #         i_left = curve.intersect(left.Edges[0].Curve)

    #         pt: Draft.Point
    #         for pt in i_left:
    #             Draft.make_point(pt.X, pt.Y, pt.Z, point_size=5.0, color=GREEN(0.7))
            
    #         i_right = curve.intersect(right.Edges[0].Curve)
            
           


    # for rl_idx in range(1, len(ref_lines)):

    #     left = ref_lines[rl_idx].copy()
    #     if rl_idx != 1:
    #         left.translate(App.Vector(0.0, -structure_thk_mm/2, 0.0))
    #     sl = Part.show(left, "str_line")
    #     sl.ViewObject.LineColor = BLUE(0.3)

    #     right = ref_lines[rl_idx].copy()
    #     right.translate(App.Vector(0.0, structure_thk_mm/2, 0.0))
    #     sr = Part.show(right, "str_line")
    #     sr.ViewObject.LineColor = BLUE(0.3)

    #     create_lightening_hole(offset, left, right)

def create_lightening_hole(int_profile: Part.Shape, left_bound: Part.Line, right_bound: Part.Line):
    print("create_lightening_hole")


