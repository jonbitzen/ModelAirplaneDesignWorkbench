import FreeCAD as App
from . import airfoil
from typing import Tuple
from typing import List
import PartDesign
import Sketcher

def create(
    obj_name: str,
) -> App.DocumentObject:
    obj: App.DocumentObject = App.ActiveDocument.addObject(
        "Part::FeaturePython",
        obj_name
    )

    Rib(obj)
    RibViewProvider(obj.ViewObject)

    App.ActiveDocument.recompute()

    return obj

class Rib():
    def __init__(
            self,
            obj: App.DocumentObject
        ) -> None:

        self.attach(obj)

        obj.addProperty(
            "App::PropertyFloat",
            "chord",
            "Rib",
            "Rib length from leading edge to trailing edge"
        ).chord = 100.0

        obj.addProperty(
            "App::PropertyFloat",
            "thickness",
            "Rib",
            "Rib thickness"
        ).thickness = 2.0

        obj.addProperty(
            "App::PropertyEnumeration",
            "airfoil",
            "Rib",
            "Airfoil type"
        ).airfoil = [e.name for e in airfoil.AirfoilType]

        obj.Proxy = self

    def onChanged(self, obj: App.DocumentObject, property: str) -> None:
        """
        Called when the FeaturePython object properties change.  Note that this
        includes any properties of the object, in addition to the App::Property*
        objects that may be added explictly by the user
        """  

        # TODO: why is it that when we change a number like the chord we dont
        #       need to explicitly call execute, but when we change the airfoil
        #       we do?     
        if property == "airfoil":
            self.execute(obj)

    def attach(self, obj: App.DocumentObject) -> None:
        pass

    def execute(self, obj: App.DocumentObject) -> None:
        self.airfoil_data = airfoil.load(airfoil.AirfoilType.to_filename(obj.getPropertyByName("airfoil")))
        tmp_placement: App.Placement = obj.Placement

        tmp_af: Sketcher.SketchObject = self.airfoil_data.to_sketch(obj.getPropertyByName("chord"))
        tmp_af.recompute()
        tmp_body: PartDesign.Body = App.ActiveDocument.addObject("PartDesign::Body", "tmp_body")
        tmp_pad = tmp_body.newObject("PartDesign::Pad", "tmp_pad")
        tmp_pad.Profile = tmp_af
        tmp_pad.Length = obj.getPropertyByName("thickness")
        tmp_pad.recompute()
        tmp_body.recompute()

        obj.Shape = tmp_body.Shape.copy()
        obj.Placement = tmp_placement

        App.ActiveDocument.removeObject(tmp_body.Name)
        App.ActiveDocument.removeObject(tmp_af.Name)
        App.ActiveDocument.removeObject(tmp_pad.Name)

class RibViewProvider():
    def __init__(self, vobj: App.Gui.ViewProviderDocumentObject) -> None:
        vobj.Proxy = self
        self.Object = vobj.Object
        self.attach(vobj)

    def getIcon(self) -> str:
        return None
    
    def attach(self, vobj: App.Gui.ViewProviderDocumentObject) -> None:
        vobj.Proxy = self
        self.Object = vobj.Object
        self.ViewObject = vobj

    def onDelete(self, vobj: App.Gui.ViewProviderDocumentObject, subelements: Tuple[str]) -> bool:
        return True
    
    def claimChildren(self) -> List[App.DocumentObject]:
        """
        Claim children of this ViewObject.  Claimed children will be displayed
        as children underneath this ViewObject.  Each child object should be an
        App.DocumentObject
        """
        return []

    def onChanged(self, vobj: App.Gui.ViewProviderDocumentObject, prop: str) -> None:
        pass
