from . import rib_hole_generators
from . import utilities
import Draft
from enum import Enum
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

def make_spar_cuboid_profile(width: float, height: float, name: str) -> Sketcher.Sketch:
    spar_sk: Sketcher.Sketch = App.activeDocument().addObject('Sketcher::SketchObject', name)
    spar_sk.Placement = utilities.xy_placement
    spar_sk.MapMode = "Deactivated"

    geoList = []
    geoList.append(Part.LineSegment(App.Vector(-1,-1,0),App.Vector(-1,1,0)))
    geoList.append(Part.LineSegment(App.Vector(-1,1,0),App.Vector(1,1,0)))
    geoList.append(Part.LineSegment(App.Vector(1,1,0),App.Vector(1,-1,0)))
    geoList.append(Part.LineSegment(App.Vector(1,-1,0),App.Vector(-1,-1,0)))
    geoList.append(Part.Point(App.Vector(0.000000,0.000000,0)))
    spar_sk.addGeometry(geoList, False)

    conList = []
    conList.append(Sketcher.Constraint('Coincident',0,2,1,1))
    conList.append(Sketcher.Constraint('Coincident',1,2,2,1))
    conList.append(Sketcher.Constraint('Coincident',2,2,3,1))
    conList.append(Sketcher.Constraint('Coincident',3,2,0,1))
    conList.append(Sketcher.Constraint('Horizontal',1))
    conList.append(Sketcher.Constraint('Horizontal',3))
    conList.append(Sketcher.Constraint('Vertical',0))
    conList.append(Sketcher.Constraint('Vertical',2))
    conList.append(Sketcher.Constraint('Symmetric',1,2,0,1,4,1))
    spar_sk.addConstraint(conList)
    del geoList, conList

    spar_sk.addConstraint(Sketcher.Constraint('Coincident',4,1,-1,1))

    spar_sk.addConstraint(Sketcher.Constraint('Distance',1,width))
    spar_sk.setDatum(10,App.Units.Quantity(width))
    spar_sk.renameConstraint(10, u'width')

    spar_sk.addConstraint(Sketcher.Constraint('Distance',2,height))
    spar_sk.setDatum(11,App.Units.Quantity(height))
    spar_sk.renameConstraint(11, u'height')

    return spar_sk
    
class RibPose():
    def __init__(self, position: App.Vector, direction: App.Vector) -> None:
        self.position = position
        self.direction = direction

class PathInterval():
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
    def __init__(self, path: List[Part.Edge]) -> None:

        self.path: List[PathInterval] = []
        self.length: float = 0

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
#       whenever it calls execute, it'll walk its own sketches
#       
class Rib():
    def __init__(
            self, 
            root_airfoil: Sketcher.Sketch,
            root_norm: App.Vector, 
            pose: RibPose,
            pf_dist: float, 
            min_y: float
        ) -> None:

        self.airfoil: Sketcher.Sketch = None
        self.structure: Sketcher.Sketch = None
        self.interferences: List[Sketcher.Sketch] = []

        wp_tangent: App.Vector = pose.direction.normalize()
        rot_axis: App.Vector = root_norm.cross(wp_tangent)
        if rot_axis.Length != 0:
            rot_axis = rot_axis.normalize()
        angle: float = math.degrees(root_norm.getAngle(wp_tangent))

        base_dist = root_airfoil.Shape.BoundBox.YLength

        scale_factor = pf_dist/base_dist

        ra_orig_placement = root_airfoil.Placement
        root_airfoil.Placement = utilities.xy_placement
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
        next_af.Placement = utilities.xy_placement

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

        self.airfoil = next_af

class Planform():    
    def __init__(self, pf_edges: List[Part.Edge]) -> None:
        self.pf_edges = Part.__sortEdges__(pf_edges)

    def get_rib_length_at(self, position: App.Vector, direction: App.Vector) -> Tuple[float, float]:
        '''
        use position and direction to construct a plane, then calculate the
        intersection between the plane and the planform to get rib length, as 
        well as the leading edge y coordinate
        '''
        plane = Part.Plane(position, direction)

        # the plane should intersect the wing planform in two places
        pf_intersections: List[Part.Vertex] = []

        for edge in self.pf_edges:

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

class WingEdge(Enum):
    LEADING = 1
    TRAILING = 2

class Wing():
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

        planform_helper = Planform(planform.Shape.Edges)

        # TODO: this wants to be encapsulated somewhere
        rib_list: List[Rib] = []
        rib_poses: List[RibPose] = self.path_helper.get_rib_poses(num_sections)
        for pose in rib_poses:
            
            # get the angle needed to rotate a rib section so its normal aligns
            # with the path tangent
            root_norm = App.Vector(-1,0,0)
            pf_dist, min_y = planform_helper.get_rib_length_at(pose.position, root_norm)

            next_rib = \
                Rib(root_airfoil,
                    root_norm,
                    pose,
                    pf_dist,
                    min_y
                )

            next_af = next_rib.airfoil

            rib_list.append(next_rib)

            obj.addObject(next_af)

        r0 = rib_list[0]
        r1 = rib_list[num_sections-2]
        r2 = rib_list[num_sections-1]

        fs_loft = self.__make_spar("front-spar", r0, 0.2, r1, 0.2, 4.0, 8.0, WingEdge.LEADING)
        self.__make_spar_penetrations(rib_list, fs_loft)

        rs_loft = self.__make_spar("rear-spar", r0, 0.33, r2, 0.5, 3.0, 3.0, WingEdge.TRAILING)
        self.__make_spar_penetrations(rib_list, rs_loft)
        
        # generate the lightening hole sketches
        for rib in rib_list:
            lh_sk = rib_hole_generators.create_lightening_hole_sketch(rib.airfoil, rib.interferences)
            rib.structure = lh_sk
            obj.addObject(lh_sk)

        # Add this last, or chaos ensues
        obj.Proxy = self

    def __make_spar_penetrations(self, rib_list: List[Rib], fs_loft: Part.Shape) -> List[Sketcher.Sketch]:
        for rib in rib_list:
            rib_face: Part.Feature = Part.makeFace([rib.airfoil.Shape.copy()], "Part::FaceMakerSimple")
            rib_face = Part.show(rib_face, "rib_face")
            rib_face.Shape = rib_face.Shape.transformGeometry(rib.airfoil.Placement.Matrix.inverse())
            rib_face.Placement = rib.airfoil.Placement
            
            sec: Part.Feature = App.ActiveDocument.addObject("Part::Section", "spar-section")
            sec.Base = fs_loft
            sec.Tool = rib_face
            sec.recompute()
            sec.Shape = sec.Shape.transformGeometry(rib_face.Placement.Matrix.inverse())
            spar_hole: Sketcher.Sketch = Draft.make_sketch([sec.Shape], autoconstraints=True, name="spar_hole")
            if spar_hole is not None:
                spar_hole.Placement = rib_face.Placement
                spar_hole.recompute()
                rib.interferences.append(spar_hole)
                
            App.ActiveDocument.removeObject(sec.Name)

    def __make_spar(self, 
            name: str,
            first_rib: Rib, first_frac: float, 
            last_rib: Rib, last_frac: float,
            width: float, height: float,
            edge: WingEdge) -> Part.Shape:
        
        p0_r = self.__get_spar_point(first_rib.airfoil, first_frac, edge)
        p0_r = Part.show(p0_r)
        p0_r.ViewObject.PointColor = utilities.BLUE(1.0)
        p0_r.ViewObject.PointSize = 10.0

        p0_t = self.__get_spar_point(last_rib.airfoil, last_frac, edge)
        p0_t = Part.show(p0_t)
        p0_t.ViewObject.PointColor = utilities.BLUE(1.0)
        p0_t.ViewObject.PointSize = 10.0

        fs_root_sk = make_spar_cuboid_profile(width, height, name+"-root")
        fs_root_sk.Placement.rotate(
            App.Vector(0,0,0),
            p0_r.Placement.Rotation.Axis,
            math.degrees(p0_r.Placement.Rotation.Angle)
        )
        fs_root_sk.Placement.Base += p0_r.Placement.Base
        fs_root_sk.recompute()

        fs_tip_sk = make_spar_cuboid_profile(width, height, name+"-tip")
        fs_tip_sk.Placement.rotate(
            App.Vector(0,0,0),
            p0_t.Placement.Rotation.Axis,
            math.degrees(p0_t.Placement.Rotation.Angle)
        )
        fs_tip_sk.Placement.Base += p0_t.Placement.Base
        fs_tip_sk.recompute()

        fs_loft: Part.Feature = App.ActiveDocument.addObject("Part::Loft", name)
        fs_loft.Sections = [fs_root_sk, fs_tip_sk]
        fs_loft.Solid = True
        fs_loft.Ruled = False
        fs_loft.Closed = False

        return fs_loft

    def __get_spar_point(
            self, 
            airfoil: Sketcher.Sketch, 
            edge_fraction: float,
            wing_edge: WingEdge
        ) -> Part.Vertex:
        orig_pose: Part.Placement = airfoil.Placement
        airfoil.Placement = utilities.xy_placement

        af_bbox = airfoil.Shape.BoundBox

        x_fs: float
        if wing_edge == WingEdge.LEADING:
            x_fs = af_bbox.XMin + af_bbox.XLength*edge_fraction
        else:
            x_fs = af_bbox.XMax - af_bbox.XLength*edge_fraction

        fs_line = Part.makeLine(App.Vector(x_fs, af_bbox.YMax+10, 0), App.Vector(x_fs, af_bbox.YMin-10, 0))
        e_l: Part.Curve = fs_line.Edges[0].Curve

        e_list = airfoil.Shape.Edges

        pts: List[Tuple[float,float]]
        for e in e_list:
            pts = e_l.intersect2d(e.Curve, Part.Plane(App.Vector(0,0,0), App.Vector(0,0,1)))
            if len(pts) > 0:
                break

        pts.sort(key=lambda y: y[1], reverse=True)
        
        p_top = pts[0]
        p_bot = pts[1]

        mid_pt = App.Vector(x_fs, (p_top[1] + p_bot[1])/2, 0)
        p_m = Part.Point(App.Vector(0,0,0)).toShape()
        p_m.Placement.Base = mid_pt

        airfoil.Placement = orig_pose

        ob = p_m.Placement.Base

        p_m.Placement.rotate(
            -ob, 
            orig_pose.Rotation.Axis, 
            math.degrees(orig_pose.Rotation.Angle)
        )
        
        p_m.Placement.Base += orig_pose.Base

        return p_m

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
    
