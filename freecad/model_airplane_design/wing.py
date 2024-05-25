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

        # Add this last, or chaos ensues
        obj.Proxy = self

    def onChanged(self, obj: App.DocumentObject, property: str) -> None:
        print("Wing.onChanged: " + property)
        
    def attach(self, obj: App.DocumentObject) -> None:
        obj.addExtension("App::OriginGroupExtensionPython")
        obj.Origin = App.ActiveDocument.addObject("App::Origin", "Origin")

    def execute(self, obj: App.DocumentObject) -> None:
        pass

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