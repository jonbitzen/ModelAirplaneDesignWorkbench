# from color import *
import Draft
import FreeCAD as App
import math
import numpy
import Part
import Sketcher
from typing import Tuple
from typing import List
from . import utilities

def create(
    obj_name: str, 
    root_airfoil: Sketcher.Sketch,
    planform: Sketcher.Sketch,
    path: Sketcher.Sketch,
    num_sections: int = 5
):    

    # TODO: Probably create this as an App::FeaturePython, since the wing is
    #       a container, and isn't expected to have a Shape
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

# TODO: I think eventually the Rib is going to need to encapsulate its airfoil
#       as well as a structure generator.  I think it'll also need to help out
#       when we go to compute an intersection with a spar
#
#       This will also likely have to be a Part::FeaturePython thing, so that it 
#       shows up in the part container with everything it owns.  Yeah, I guess 
#       whenever it calls execute, it'll walk its own sketches, and then it'll
#       
class Rib():

    airfoil:    Sketcher.Sketch
    structure:  Sketcher.Sketch

    # so we could give this the root airfoil and a pose, and I think it could
    # encapsulate all the calculations related to 
    def __init__(self, airfoil: Sketcher.Sketch, structure: Sketcher.Sketch = None) -> None:
        self.airfoil = airfoil
        self.structure = structure
            
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

        xy_placement = App.Placement(
            App.Vector(0,0,0),
            App.Vector(0,0,0),
            0
        )

        # TODO: this wants to be encapsulated somewhere
        rib_list: List[Rib] = []
        rib_poses: List[RibPose] = self.path_helper.get_rib_poses(num_sections)
        for pose in rib_poses:
            
            # get the angle needed to rotate a rib section so its normal aligns
            # with the path tangent
            root_norm = App.Vector(-1,0,0)
            wp_tangent: App.Vector = pose.direction.normalize()
            rot_axis: App.Vector = root_norm.cross(wp_tangent)
            if rot_axis.Length != 0:
                rot_axis = rot_axis.normalize()
            angle: float = math.degrees(root_norm.getAngle(wp_tangent))

            # get the scale factor and position offset to apply to the rib section
            # that keeps it inside the planform bounds at this position
            pf_dist, min_y = \
                self.__get_planform_dist(
                    planform.Shape.Edges,
                    pose.position, 
                    root_norm
                )

            base_dist = root_airfoil.Shape.BoundBox.YLength

            scale_factor = pf_dist/base_dist

            ra_orig_placement = root_airfoil.Placement
            root_airfoil.Placement = xy_placement
            af_shape: Part.Shape = root_airfoil.Shape.copy()
            af_shape.scale(scale_factor)
            next_af: Sketcher.SketchObject = Draft.make_sketch(af_shape, True)
            
            root_airfoil.Placement = ra_orig_placement

            next_af.Placement = ra_orig_placement

            # store off the original placement - we need it to get back into the
            # original plane later
            orig_placement = next_af.Placement.Base
            orig_axis = next_af.Placement.Rotation.Axis
            orig_angle = math.degrees(next_af.Placement.Rotation.Angle)

            # clear the placement, so we can rotate the rib section properly
            next_af.Placement = xy_placement

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

            # adjust the placement Y coordinate to keep the airfoil inside the
            # wing planform
            next_af.Placement.Base.y += min_y - root_airfoil.Shape.BoundBox.YMin * scale_factor
            next_af.recompute()
            rib_list.append(Rib(next_af))

            obj.addObject(next_af)

        for rib in rib_list:
            lh_sk = utilities.create_lightening_hole_sketch(rib.airfoil, [])
            rib.structure = lh_sk
            obj.addObject(lh_sk)

        # Add this last, or chaos ensues
        obj.Proxy = self

    # TODO:  I have a feeling the Planform is going to want to be its own object
    #        with data and operations, eventually
    def __get_planform_dist(
            self, 
            pf_edges: List[Part.Edge], 
            position: App.Vector, 
            direction: App.Vector
        ) -> Tuple[float, float]:
        pf_edges = Part.__sortEdges__(pf_edges)
        plane = Part.Plane(position, direction)

        # the plane should intersect the wing planform in two places
        pf_intersections: List[Part.Vertex] = []

        for edge in pf_edges:

            trimmed_curve: Part.Curve = edge.Curve
            trimmed_curve = trimmed_curve.trim(*edge.ParameterRange)

            intersections = plane.intersect(trimmed_curve)
            # if there are no intersections, skip the rest of the loop
            if intersections is None:
                continue
            
            p: Part.Point
            for p in intersections[0]:
                pf_intersections.append(Part.Vertex(p))

            # if we have found the needed planform intersections, we're done
            if len(pf_intersections) == 2:
                p0: Part.Vertex = pf_intersections[0]
                p1: Part.Vertex = pf_intersections[1]
                dist_t = p0.distToShape(p1)
                dist = dist_t[0]
                return dist, min(p0.Y, p1.Y)
        
        print("Wing.__get_planform_dist - no intersections to planform found at position = " + str(position))
        return None

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
    
