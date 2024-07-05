from . import airfoil
from . import rib_hole_generators as rhg
from . import utilities
from BOPTools import BOPFeatures
import FreeCAD as App
import Part
import Sketcher
from typing import Tuple
from typing import List

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
    """
        Proxy class that handles wing rib generation logic on behalf of a 
        Part::FeaturePython instance
    """
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

        obj.addProperty(
            "App::PropertyEnumeration",
            "hole_type",
            "Rib",
            "Lightening hole type"
        ).hole_type = [e.name for e in rhg.HoleGeneratorType]

        obj.addProperty(
            "App::PropertyLinkList",
            "interferences",
            "Rib",
            "Other solids that intersect the rib"
        ).interferences = []

        self.Object = obj
        obj.Proxy = self

    def onChanged(self, obj: App.DocumentObject, property: str) -> None:
        """
        Called when the FeaturePython object properties change.  Note that this
        includes any properties of the object, in addition to the App::Property*
        objects that may be added explictly by the user
        """  

        match property:
            case "airfoil" | "hole_type" | "interferences":
                self.execute(obj)
            case _:
                pass

    def attach(self, obj: App.DocumentObject) -> None:
        pass


    def __make_body_centered_rib_feature(self, rib_sketch: Sketcher.Sketch, thickness: float) -> Part.Feature:
        rib_extr: Part.Feature = App.ActiveDocument.addObject("Part::Extrusion", "rib_extr")
        rib_extr.Base = rib_sketch
        rib_extr.Dir = App.Vector(0, 0, thickness)
        rib_extr.Solid = True
        rib_extr.recompute()

        body_center = rib_sketch.Shape.BoundBox.Center
        transform_mtx = App.Matrix()
        transform_mtx.move(-body_center)
        
        rib_ftr: Part.Feature = App.ActiveDocument.addObject("Part::Feature", "rib_ftr")
        rib_ftr.Shape = rib_extr.Shape.transformed(transform_mtx, copy=True)

        App.ActiveDocument.removeObject(rib_extr.Name)

        return rib_ftr

    def execute(self, obj: App.DocumentObject) -> None:
        self.airfoil_data = airfoil.load(airfoil.AirfoilType.to_filename(obj.getPropertyByName("airfoil")))

        # register temporary doc objects here, so they can be easily disposed
        # when we're done with them
        tmp_to_delete: List[Part.Feature] = []

        orig_placement: App.Placement = obj.Placement

        rib_sketch: Sketcher.SketchObject = self.airfoil_data.to_sketch(obj.getPropertyByName("chord"))
        rib_sketch.Shape.tessellate(0.01)
        rib_sketch.recompute()
        tmp_to_delete.append(rib_sketch)

        # create an extruded rib form, we're going to use it to find common bool
        # shapes with the intersecting objects
        rib_extr = self.__make_body_centered_rib_feature(rib_sketch, obj.getPropertyByName("thickness"))
        tmp_to_delete.append(rib_extr)
        rib_extr.recompute()

        # create a set of hole exclusion regions by finding the common bool
        # between each of the intersecting objects and the rib and applying a
        # standoff to the bbox of the resulting shape
        boolean_tool = BOPFeatures.BOPFeatures(App.ActiveDocument)
        hole_exclusions: List[rhg.HoleExclusion] = []
        scale_factor: float = obj.chord / 175.0
        rib_pen_standoff: float = scale_factor * 2.0
        for interference in obj.interferences:
            rib_intersection: Part.Feature = boolean_tool.make_common([rib_extr.Name, interference.Name])
            rib_intersection.recompute()
            hole_exclusions.append(rhg.HoleExclusion(rib_intersection.Shape.BoundBox, rib_pen_standoff))
            tmp_to_delete.append(rib_intersection)

        inner_profile_standoff: float = scale_factor * 2.0
        hbg =rhg.HoleBoundGenerator(rib_sketch, inner_profile_standoff, hole_exclusions)

        hbr_list = hbg.get_hole_bounds()
 
        # Calculate lightening holes for each hole bound region, and add the
        # geometry to the rib sketch 
        hbg = rhg.HoleGeneratorFactory.create_generator(obj.hole_type, obj.chord)
        if hbg is not None:
            for hbr in hbr_list:
                lh_sk = hbg.generate_sketch(hbr)
                lh_sk.recompute()
                tmp_to_delete.append(lh_sk)
                rib_sketch.addGeometry(lh_sk.Geometry)
                
        rib_sketch.recompute()        

        # create a solid for the final rib shape, starting with the sketch that
        # has the lightening holes in it
        rib_final_extr = self.__make_body_centered_rib_feature(rib_sketch, obj.getPropertyByName("thickness"))
        rib_final_extr.Placement = orig_placement
        rib_final_extr.recompute()
        tmp_to_delete.append(rib_final_extr)

        # make cuts for each of the intereferences in the rib solid
        for interference in obj.interferences:
            rib_final_extr: Part.Feature = boolean_tool.make_cut([rib_final_extr.Name, interference.Name])
            rib_final_extr.recompute()
            tmp_to_delete.append(rib_final_extr)

        obj.Shape = rib_final_extr.Shape.transformed(orig_placement.Matrix.inverse(), True)

        for tmp_ft in tmp_to_delete:
            App.ActiveDocument.removeObject(tmp_ft.Name)

class RibViewProvider():
    """
        Proxy class that handles wing rib view operations on behalf of a 
        Part::FeaturePython class' ViewProvider member
    """
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
