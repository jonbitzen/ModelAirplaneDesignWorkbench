import FreeCAD as App
from . import airfoil
from enum import Enum
from typing import Tuple
from typing import Optional

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
        ).chord = 1.0

        obj.addProperty(
            "App::PropertyPythonObject",
            "airfoil_data",
            "Rib",
            "Airfoil definition data"
        ).airfoil_data = None

        # self._airfoil_data = None
        obj.Proxy = self

    # @property
    # def airfoil_data(self) -> Optional[airfoil.AirfoilData]:
    #     return self._airfoil_data
    
    # @airfoil_data.setter
    # def airfoil_data(self, value: Optional[airfoil.AirfoilData]) -> None:
    #     self._airfoil_data = value

    def onChanged(self, obj: App.DocumentObject, property: str) -> None:

        match property:
            case "Label":
                pass
            case "Proxy":
                pass
            case _:
                print(property + " = " + str(obj.getPropertyByName(property)))

    def attach(self, obj: App.DocumentObject) -> None:
        pass

    def execute(self, obj: App.DocumentObject) -> None:
        print("Rib.execute")



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
    
    def onChanged(self, vobj: App.Gui.ViewProviderDocumentObject, prop: str) -> None:
        pass
