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

class TempDocObjectHelper():
    """
    This class holds temporary objects in a with-as context block, to ensure that
    the held objects are removed even if some of the geometry generation fails
    """
    def __init__(self) -> None:
        self.doc = App.ActiveDocument
        self.tmp_objects : List[App.DocumentObject] = []

    def addObject(self, tmp_obj: App.DocumentObject, do_delete: bool = True) -> App.DocumentObject:
        """
        Adds a document object to the internal cache to be released when the
        with-as block exits

        Parameters
        ----------
        tmp_obj: App.DocumentObject
            the document object to be cleaned up when the with-as block ends

        Return
        ------
        A reference to the original document object that will be managed

        """
        if do_delete:
            self.tmp_objects.append(tmp_obj)
        return tmp_obj
    
    def removeObject(self, tmp_obj: App.DocumentObject) -> None:
        """
        Removes a document object from the internal cache to be released when
        the with-as block exits.

        Parameters
        ----------
        tmp_obj: App.DocumentObject
            the document object that will be released from the cache; the object
            will not be cleaned up when the with-as block ends
        """
        for idx in range(len(self.tmp_objects)):
            if self.tmp_objects[idx].Name == tmp_obj.Name:
                self.tmp_objects.pop(idx)
                break

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        for tmp in self.tmp_objects:
            self.doc.removeObject(tmp.Name)
        self.tmp_objects = []

    def __del__(self):
        for tmp in self.tmp_objects:
            self.doc.removeObject(tmp.Name)
        self.tmp_objects = []
            

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
    
