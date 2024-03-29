# from color import *
import FreeCAD as App
import math
import numpy
import Part
import Sketcher
from typing import Tuple
from typing import List

def create(
    obj_name: str, 
    root_airfoil: Sketcher.Sketch,
    planform: Sketcher.Sketch,
    path: Sketcher.Sketch,
    num_sections: int = 5
):    

    obj: App.DocumentObject = App.ActiveDocument.addObject(
        "Part::FeaturePython",
        obj_name
    )

    Wing(obj, root_airfoil, planform, path, num_sections),
    WingViewProvider(obj.ViewObject),

    App.ActiveDocument.recompute()
    return obj

class RibPose():
    position: App.Vector
    direction: App.Vector
    def __init__(self, position: App.Vector, direction: App.Vector) -> None:
        self.position = position
        self.direction = direction

class PathInterval():
    start: float
    end: float
    offset: float
    edge: Part.Edge

    def __init__(self, offset: float, edge: Part.Edge) -> None:
        self.edge = edge
        self.offset = offset
        # TODO: is there a universe where FirstParameter isn't zero, if so
        #       there is a discontinuity to cope with
        self.start = edge.FirstParameter + offset
        self.end = edge.LastParameter + offset

    def contains(self, location: float):
        return location >= self.start and location <= self.end

    def get_pose_at(self, location: float):
        if not self.contains(location):
            return None
        loc: float = location - self.offset
        pos: App.Vector = self.edge.valueAt(loc)
        tan: App.Vector = self.edge.tangentAt(loc)
        return RibPose(pos, tan)

class PathHelper():
    path: List[PathInterval] = []
    length: float = 0
    def __init__(self, path: List[Part.Edge]) -> None:

        origin: App.Vector = App.Vector(0,0,0)
        path.sort(key=lambda edge: edge.valueAt(edge.FirstParameter).distanceToPoint(origin))

        self.length = 0
        for edge in path:
            path_segment = \
                PathInterval(
                    self.length,
                    edge
                )

            self.path.append(path_segment)
            self.length += edge.Length

    def get_rib_poses(self, num_poses: int):
        pose_list: List[RibPose] = []
        for edge_dist in numpy.linspace(0.0, self.length, num_poses):
            pose_generated: bool = False
            for edge_int in self.path:
                
                rib_pose: RibPose = edge_int.get_pose_at(edge_dist)
                if rib_pose is not None:
                    pose_list.append(rib_pose)
                    pose_generated = True
                    continue

            if not pose_generated:
                print("PathHelper.get_rib_poses: failed to generate pose at " +  str(edge_dist))
        return pose_list                    
            
            
            
class Wing():

    path_helper: PathHelper

    def __init__(
            self, 
            obj: App.DocumentObject,
            root_airfoil: Sketcher.Sketch,
            planform: Sketcher.Sketch,
            path: Sketcher.Sketch,
            num_sections: int = 5
        ) -> None:

        self.attach(obj)

        obj.addProperty(
            "App::PropertyLink", 
            "root_airfoil", 
            "Wing", 
            "The airfoil to use for the wing root"
        ).root_airfoil = root_airfoil
        obj.addObject(obj.root_airfoil)

        obj.addProperty(
            "App::PropertyLink",
            "planform",
            "Wing",
            "A sketch of the wing planform"
        ).planform = planform
        obj.addObject(obj.planform)

        obj.addProperty(
            "App::PropertyLink",
            "path",
            "Wing",
            "A sketch giving the wing extent through the planform"
        ).path = path
        obj.addObject(obj.path)

        obj.addProperty(
            "App::PropertyQuantity",
            "num_sections",
            "Wing",
            "The number of wing sections to generate along the path"
        ).num_sections = num_sections

        self.path_helper = PathHelper(path.Shape.Edges)

        # TODO: this wants to be encasulated somewhere
        rib_poses: List[RibPose] = self.path_helper.get_rib_poses(num_sections)
        for pose in rib_poses:
            next_af: Sketcher.SketchObject = App.ActiveDocument.copyObject(root_airfoil)

            # store off the original placement - we need it to get back into the original plane later
            orig_placement = next_af.Placement.Base
            orig_axis = next_af.Placement.Rotation.Axis
            orig_angle = math.degrees(next_af.Placement.Rotation.Angle)

            # clear the placement, so we can rotate the rib section properly
            next_af.Placement.Base = App.Vector(0,0,0)
            next_af.Placement.Rotation.Axis = App.Vector(0,0,0)
            next_af.Placement.Rotation.Angle = 0

            # get the angle needed to rotate a rib section so its normal aligns with the path tangent
            root_norm = App.Vector(-1,0,0)
            dir: App.Vector = pose.direction.normalize()
            rot_axis: App.Vector = root_norm.cross(dir).normalize()
            angle: float = math.degrees(root_norm.getAngle(dir))

            # rotate the rib section into the path tangent
            next_af.Placement.rotate(
                next_af.Shape.BoundBox.Center,
                rot_axis,
                angle
            )

            # rotate once more to get the rib section back into its original frame
            next_af.Placement.rotate(
                App.Vector(0,0,0),
                orig_axis,
                orig_angle
            )

            # locate the rib in the appropriate place along the path
            next_af.Placement.Base = pose.position + orig_placement
            obj.addObject(next_af)
        
        # Add this last, or chaos ensues
        obj.Proxy = self

    def onChanged(self, obj: App.DocumentObject, property: str) -> None:
        property_list: List[str] = ["root_airfoil", "planform", "path", "num_sections"]

        for prop in property_list:
            if not hasattr(obj, prop):
                return

        if property in property_list:
            self.execute(obj)

    def attach(self, obj: App.DocumentObject) -> None:
        obj.addExtension("App::OriginGroupExtensionPython")
        obj.Origin = App.ActiveDocument.addObject("App::Origin", "Origin")

    def execute(self, obj: App.DocumentObject) -> None:
        print("Wing.execute")

class WingViewProvider():
    def __init__(self, vobj: App.Gui.ViewProviderDocumentObject) -> None:
        vobj.Proxy = self
        self.Object = vobj.Object
        self.attach(vobj)

    def getIcon(self) -> str:
        return None
    
    def attach(self, vobj: App.Gui.ViewProviderDocumentObject) -> None:
        vobj.addExtension("Gui::ViewProviderOriginGroupExtensionPython")
        vobj.Proxy = self
        self.Object = vobj.Object
        self.ViewObject = vobj
        
    def onDelete(self, vobj: App.Gui.ViewProviderDocumentObject, subelements: Tuple[str]) -> bool:
        return True
    
    def onChanged(self, vobj: App.Gui.ViewProviderDocumentObject, prop: str) -> None:
        pass
    
