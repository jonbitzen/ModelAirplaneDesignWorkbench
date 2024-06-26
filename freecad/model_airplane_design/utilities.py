import Draft
import FreeCAD as App
import Part
import numpy
from typing import List

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

origin = App.Vector(0,0,0)
x_axis = App.Vector(1,0,0)
y_axis = App.Vector(0,1,0)
z_axis = App.Vector(0,0,1)

xy_placement = \
    App.Placement(
        App.Vector(0,0,0),
        App.Vector(0,0,0),
        0
    )

yz_placement = \
    App.Placement(
        App.Vector(0.000000, 0.000000, 0.000000), 
        App.Rotation(0.500000, 0.500000, 0.500000, 0.500000)
    )

xz_placement = \
    App.Placement(
        App.Vector(0.000000, 0.000000, 0.000000), 
        App.Rotation(0.707107, 0.000000, 0.000000, 0.707107)
    )

# TODO: I think we need a utility that automates change of placement and geometry
#       that I seem to be doing everywhere on an ad-hoc basis

def makePointV(point: App.Vector, color: tuple = BLACK(), point_size: float = 5.0) -> Part.Feature:
    p_t = Part.Point(point)
    p_t = Part.show(p_t.toShape())
    p_t.ViewObject.PointColor = color
    p_t.ViewObject.PointSize = point_size
    return p_t

def makePointP(point: Part.Point, color: tuple = BLACK(), point_size: float = 5.0) -> Part.Feature:
    p_t = Part.show(point.toShape())
    p_t.ViewObject.PointColor = color
    p_t.ViewObject.PointSize = point_size
    return p_t

def draw_points(pts: List[App.Vector], point_size=10.0, color: tuple = GREEN(0.5)) -> None:
    drawn_pts: List[Draft.Point] = []
    for pt in pts:
        pf = Draft.make_point(pt)
        pf.ViewObject.PointColor = color
        pf.ViewObject.PointSize = 10.0
        drawn_pts.append(pf)
    return drawn_pts
    
