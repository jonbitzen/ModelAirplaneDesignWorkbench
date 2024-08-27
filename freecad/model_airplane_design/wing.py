from . import utilities
from . import elevation_path
from . import planform
from . import rib
import FreeCAD as App
import PartDesign
import math
from typing import List, Tuple

def create(obj_name: str) -> App.DocumentObject:

    obj: App.DocumentObject = App.ActiveDocument.addObject(
        "Part::FeaturePython",
        obj_name
    )

    Wing(obj)
    WingViewProvider(obj.ViewObject)

    App.ActiveDocument.recompute()

    return obj

class Wing():
    def __init__(self, obj: App.DocumentObject) -> None:
        
        # attach the origin group extension, so that this wing is created as a
        # container for other doc objects
        self.attach(obj)

        obj.addProperty(
            "App::PropertyLink",
            "elevation_path",
            "Wing",
            "A sketch giving the wing extent through the planform"
        ).elevation_path = \
            utilities.load_feature_asset(
                "elevation_path_default", 
                "Sketcher::SketchObject",
                obj_name="elevation_path"
            )
        obj.addObject(obj.elevation_path)
   
        obj.addProperty(
            "App::PropertyLink",
            "planform",
            "Wing",
            "A sketch of the wing planform"
        ).planform = \
            utilities.load_feature_asset(
                "planform_default",
                "Sketcher::SketchObject",
                obj_name="planform"
            )
        obj.addObject(obj.planform)

        obj.addProperty(
            "App::PropertyLinkList",
            "rib_list",
            "Wing",
            "The number of ribs to generate"
        ).rib_list = []
        

        obj.addProperty(
            "App::PropertyInteger",
            "num_ribs",
            "Wing",
            "The number of ribs to generate"
        ).num_ribs = 6

        obj.addProperty(
            "App::PropertyAngle",
            "root_cant_angle",
            "Wing",
            "Cant angle of the root rib, range -15/+15 degrees"
        ).root_cant_angle = 0

        # Add this last, or chaos ensues
        obj.Proxy = self

        self.rebuild_wing(obj)

    def onChanged(self, obj: App.DocumentObject, property: str) -> None:
        do_exec = False
        match property:
            case "num_ribs":
                if obj.num_ribs < 2:
                    print("Wing.onChanged: number of ribs may not be less than 2")
                    obj.num_ribs = 2
                do_exec = True

            case "root_cant_angle":
                if obj.root_cant_angle > 15.0 or obj.root_cant_angle < -15.0:
                    print("Wing.onChanged: root_cant_angle must be between -15.0/+15.0 degrees")
                    obj.root_cant_angle = 0.0
                do_exec = True

            case _:
                pass

    def attach(self, obj: App.DocumentObject) -> None:
        obj.addExtension("App::OriginGroupExtensionPython")
        obj.Origin = App.ActiveDocument.addObject("App::Origin", "Origin")

    def rebuild_wing(self, obj: App.DocumentObject) -> None:
        rib_list: List[PartDesign.Body] = obj.rib_list
        for r in rib_list:
            r.removeObjectsFromDocument()
            obj.removeObject(r)
            App.ActiveDocument.removeObject(r.Name)
        obj.rib_list = []
        rib_list = []

        path_helper = elevation_path.PathHelper(obj.elevation_path.Shape.Edges)
        planform_helper = planform.Planform(obj.planform.Shape.Edges)
        rib_poses = path_helper.get_poses(obj.num_ribs)
        
        idx: int = 0 
        for pose in rib_poses:
            chord, ctr_pt = planform_helper.get_rib_chord_at(pose.position)
            ctr_pt.z = pose.position.z

            rib_name = "rib"+str(idx)
            rib_base_name = rib_name+"_base"
            rib_body: PartDesign.Body = rib.create(rib_name)
            r = rib_body.getObject(rib_base_name)
            r.chord = chord

            path_tan = pose.direction.normalize()
            rib_norm = -utilities.x_axis
            rot_axis = rib_norm.cross(path_tan)
            angle = math.degrees(rib_norm.getAngle(path_tan))
            
            # TODO: Why does this work when we use the x_axis, but not when we
            #       use rot_axis?  I have a bad feeling using x_axis may not be
            #       general
            p = utilities.yz_placement.copy()
            p.rotate(
                App.Vector(0,0,0),
                utilities.x_axis,
                angle
            )

            p.Base = ctr_pt
            rib_body.Placement = p
            rib_body.Placement.Base = ctr_pt
            rib_body.recompute(True)
            rib_list.append(rib_body)
            obj.addObject(rib_body)
            idx += 1
        obj.rib_list = rib_list


    def execute(self, obj: App.DocumentObject) -> None:
        pass
        
class WingViewProvider():
    def __init__(self, vobj: App.Gui.ViewProviderDocumentObject) -> None:
        vobj.Proxy = self
        self.Object = vobj.Object
        self.attach(vobj)

    def doubleClicked(self, vobj: App.Gui.ViewProviderDocumentObject) -> None:        
        w_obj: Wing = vobj.Object.Proxy
        w_obj.rebuild_wing(vobj.Object)

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