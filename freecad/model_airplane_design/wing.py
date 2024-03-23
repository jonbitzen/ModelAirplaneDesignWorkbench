import Draft
import FreeCAD as App
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

        num_sections: int = int(obj.num_sections)
        
        path_shape: Part.Shape = obj.path.Shape

        edge_list: List[Part.Edge] = path_shape.Edges

        e: Part.Edge = edge_list[0]

        path_len = e.Length

        for u in numpy.linspace(0.0, path_len, num_sections):
            pt: App.Vector = e.valueAt(u)

            next_af: Sketcher.SketchObject = App.ActiveDocument.copyObject(root_airfoil)
            next_af.Placement.Base = root_airfoil.Shape.Placement.Base + pt - path_shape.Placement.Base
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
    
