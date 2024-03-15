import FreeCAD as App
import Part
import Sketcher
from typing import List

# def create(obj_name: str, root_airfoil: Sketcher.Sketch):
def create(obj_name: str):

    fp_obj = App.ActiveDocument.addObject("Part::FeaturePython", obj_name)
    # Wing(fp_obj, root_airfoil)
    Wing(fp_obj)
    WingViewProvider(fp_obj.ViewObject)

    App.ActiveDocument.recompute()


    return fp_obj


class Wing():
    # def __init__(self, fp_obj, root_airfoil: Sketcher.Sketch) -> None:
    def __init__(self, fp_obj):
        self.Type = 'Wing'
        fp_obj.Proxy = self

        fp_obj.addExtension("Part::AttachExtensionPython")

        fp_obj.addProperty(
            "App::PropertyLink", 
            "root_airfoil", 
            "Wing", 
            "The airfoil to use for the wing root"
        ).root_airfoil = None

    def onChanged(self, fp_obj, property: str) -> None:
        print("-- onChanged")
        print("     property = " + property)

    def execute(self, fp_obj) -> None:
        print("execute")

class WingViewProvider():
    def __init__(self, fp_view_obj: App.Gui.ViewProvider) -> None:
        fp_view_obj.Proxy = self
        self.Object = fp_view_obj.Object

    def getIcon(self):
        return None
    
    def attach(self, fp_obj):
        self.Object = fp_obj
        self.onChanged(fp_obj, "root_airfoil")

    def claimChildren(self):

        return [self.Object.root_airfoil]
    
    def onDelete(self, feature, subelements):
        print("--- onDelete ---")
        print(type(feature))
        print(type(subelements))
        return True
    
    def onChanged(self, fp_obj, prop):
        
        pass

    def __get_state__(self):
        return None

    def __set_state__(self, state):
        return None
    
