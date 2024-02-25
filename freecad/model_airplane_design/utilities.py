import FreeCAD as App
import Part
import Sketcher
import numpy

class LighteningHoleBounds:
    def __init__(self, shape: Part.Shape, left: Part.Line, right: Part.Line):
        self.left = left
        self.right = right
        self.shape = shape

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

def create_airfoil_inner_profile(airfoil_sk : Sketcher.Sketch, offset_mm: float = 2.0, structure_thk_mm: float = 2.0):
    print("creating inner profile for airfoil using airfoil sketch " + airfoil_sk.Label)
    offset = airfoil_sk.Shape.makeOffset2D(-offset_mm, join=2)
    f = Part.show(offset, "airfoil-offset")
    f.ViewObject.LineColor = RED(1.0)
    bbox = offset.BoundBox
    num_segments = 5
    segment_length = bbox.YLength/num_segments
    
    ref_grp = create_group("ref_lines")
    ref_lines = [Part.Edge]
    for y_coord in numpy.arange(bbox.YMin, bbox.YMax + segment_length, segment_length):
        top_pt = (0, y_coord, bbox.ZMax + 1)
        bottom_pt = (0, y_coord, bbox.ZMin - 1)
        
        l = Part.makeLine(top_pt, bottom_pt)
        f = Part.show(l, "ref_line")
        f.ViewObject.LineColor = GREY(0.5)
        ref_grp.addObject(f)
        ref_lines.append(l)

    for rl_idx in range(2, len(ref_lines)-1):

        left = ref_lines[rl_idx].copy()
        left.translate(App.Vector(0.0, -structure_thk_mm/2, 0.0))
        sl = Part.show(left, "str_line")
        sl.ViewObject.LineColor = BLUE(0.3)

        right = ref_lines[rl_idx].copy()
        right.translate(App.Vector(0.0, structure_thk_mm/2, 0.0))
        sr = Part.show(right, "str_line")
        sr.ViewObject.LineColor = BLUE(0.3)
        
        create_lightening_hole(offset, left, right)

def create_lightening_hole(int_profile: Part.Shape, left_bound: Part.Line, right_bound: Part.Line):
    print("create_lightening_hole")


