from . import utilities
import FreeCAD as App
import Part
import Sketcher
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

        if do_exec is True:    
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