from . import utilities
import FreeCAD as App
import Part
import Sketcher
from typing import Dict, List, Tuple

def create_cut_tool(obj_name: str) -> App.DocumentObject:

    obj: App.DocumentObject = \
        App.ActiveDocument.addObject(
            "Part::FeaturePython",
            obj_name
        )
    
    AileronWedgeCutter(obj)
    AileronWedgeCutterViewProvider(obj.ViewObject)

    App.ActiveDocument.recompute()

    return obj


class AileronWedgeCutter():
    def __init__(self, obj: App.DocumentObject) -> None:
        self.attach(obj)

        obj.addProperty(
            "App::PropertyFloat",
            "wing_edge_height",
            "Aileron Cut Tool",
            "Wing thickness at the aileron-wing joint"
        ).wing_edge_height = 15.0

        obj.addProperty(
            "App::PropertyFloat",
            "top_cut_width",
            "Aileron Cut Tool",
            "Width of the aileron cut at the upper airfoil surface"
        ).top_cut_width = 5.0

        obj.addProperty(
            "App::PropertyFloat",
            "cut_angle",
            "Aileron Cut Tool",
            "Cut angle into the aileron member"
        ).cut_angle = 60.0

        obj.addProperty(
            "App::PropertyFloat",
            "tool_length",
            "Aileron Cut Tool",
            "Length of the aileron cut tool"
        ).tool_length = 40.0

        self.Object = obj
        obj.Proxy = self
        
    def onChanged(self, obj: App.DocumentObject, property: str) -> None:
        pass

    def attach(self, obj: App.DocumentObject) -> None:
        pass

    # TODO: this belongs in utilities
    def __get_constraint_index_map(self, sketch: Sketcher.Sketch) -> Dict[str, int]:
        constraint_map: Dict[str,int] = {}
        for idx in range(len(sketch.Constraints)):
            constraint = sketch.Constraints[idx]
            if constraint.Name != '':
                constraint_map[constraint.Name] = idx
        return constraint_map

    def execute(self, obj: App.DocumentObject) -> None:
        with utilities.TempDocObjectHelper() as tmp_helper:
             # TODO: need to do some range checking here
            profile: Sketcher.SketchObject = tmp_helper.addObject(utilities.load_feature_asset("aileron_cut_profile", "Sketcher::SketchObject"), do_delete=True)
            profile_constraint_map = self.__get_constraint_index_map(profile)
            profile.setDatum(profile_constraint_map["wing_edge_height"], obj.wing_edge_height)
            profile.setDatum(profile_constraint_map["top_cut_width"], obj.top_cut_width)
            profile.setDatum(profile_constraint_map["cut_angle"], App.Units.Quantity(str(obj.cut_angle) + ' deg'))
            profile.recompute()

            path: Sketcher.SketchObject = tmp_helper.addObject(utilities.load_feature_asset("aileron_cut_path", "Sketcher::SketchObject"),do_delete=True)
            path_constraint_map = self.__get_constraint_index_map(path)
            path.setDatum(path_constraint_map["tool_length"], obj.tool_length)
            path.recompute()

            sweep: Part.Feature = tmp_helper.addObject(App.ActiveDocument.addObject("Part::Sweep", "tmp_sweep"), do_delete=True)
            sweep.Sections = [profile]
            sweep.Spine = path
            sweep.Solid = True
            sweep.recompute()
            obj.Shape = sweep.Shape.copy()

class AileronWedgeCutterViewProvider():

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