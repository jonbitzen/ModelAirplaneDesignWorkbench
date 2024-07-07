import Draft
import FreeCAD as App
from freecad.model_airplane_design import ASSETPATH
import os
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

epsilon: float = 1.0e-6

def save_feature_asset(doc_obj: Part.Feature) -> None:
    """
    Save a document object into the plugin resources/assets folder.  This is a
    feature intended for developer use, to author default assets that may be
    loaded during plugin operation

    Parameters
    ----------
    doc_obj: Part.Feature
        The document object to store as an asset in the plugin resources/assets
        folder.  Will be saved with a filename of the form "${doc_obj.Name}.asset"
    """
    filename = doc_obj.Name + ".asset"
    filepath = os.path.join(ASSETPATH, filename)
    data_file = open(filepath, "wb")
    data = doc_obj.dumpContent()
    data_file.write(data)
    data_file.close() 

def load_feature_asset(asset_name: str, obj_type: str, obj_name: str = None) -> Part.Feature:
    """
    Load a document object from the plugin resources/assets folder

    Parameters
    ----------
    asset_name: str
        Name of the asset as it would appear in the FreeCAD document where it was
        originally authored.  For example,  if a sketch object could be accessed
        as "App.ActiveDocument.my_sketch", the asset name will be "my_sketch"

    obj_type: str
        DocumentObject class (e.g., Part::Feature, Sketcher::SketchObject, etc)
        that would normally be applied when the asset was originally authored

    obj_name: str
        Rename the DocumentObject "Name" and "Label" field to use the user-provided
        name; if not specified the "asset_name" variable value will be used
        instead 
    """
    if obj_name is None:
        obj_name = asset_name
    filename = asset_name + ".asset"
    filepath = os.path.join(ASSETPATH, filename)
    data_file = open(filepath, "rb")
    data = data_file.read()
    data_file.close()
    obj: Part.Feature = App.ActiveDocument.addObject(obj_type, obj_name)
    obj.restoreContent(data)
    obj.Label = obj_name
    return obj

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

def draw_points(pts: List[App.Vector], point_size=10.0, color: tuple = GREEN(0.5)) -> List[Draft.Point]:
    drawn_pts: List[Draft.Point] = []
    for pt in pts:
        pf = draw_point(pt, point_size=point_size, color=color)
        drawn_pts.append(pf)
    return drawn_pts

def draw_point(pt: App.Vector, point_size=10.0, color: tuple = GREEN(0.5)) -> Draft.Point:
    pf = Draft.make_point(pt)
    pf.ViewObject.PointColor = color
    pf.ViewObject.PointSize = point_size
    return pf

def draw_line(start: App.Vector, dir: App.Vector, line_len=1.0) -> Part.Feature:
    norm_dir = dir.normalize()
    pl = Part.makeLine(start, start+line_len*norm_dir)
    pl = Part.show(pl)
    return pl
    
