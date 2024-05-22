from . import airfoil
from . import rib_hole_generators as rhg
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
            "intersections",
            "Rib",
            "Other solids that intersect the rib"
        ).intersections = []

        self.Object = obj
        obj.Proxy = self

    def onChanged(self, obj: App.DocumentObject, property: str) -> None:
        """
        Called when the FeaturePython object properties change.  Note that this
        includes any properties of the object, in addition to the App::Property*
        objects that may be added explictly by the user
        """  

        match property:
            case "airfoil" | "hole_type" | "intersections":
                self.execute(obj)
            case _:
                pass

    def attach(self, obj: App.DocumentObject) -> None:
        pass

    def execute(self, obj: App.DocumentObject) -> None:
        self.airfoil_data = airfoil.load(airfoil.AirfoilType.to_filename(obj.getPropertyByName("airfoil")))

        # register temporary doc objects here, so they can be easily disposed
        # when we're done with them
        tmp_to_delete: List[Part.Feature] = []

        tmp_placement: App.Placement = obj.Placement

        tmp_af: Sketcher.SketchObject = self.airfoil_data.to_sketch(obj.getPropertyByName("chord"))
        tmp_af.recompute()
        tmp_to_delete.append(tmp_af)

        tmp_ext = App.ActiveDocument.addObject("Part::Extrusion", "tmp_ext")
        tmp_ext.Base = tmp_af
        tmp_ext.Dir = App.Vector(0, 0, obj.getPropertyByName("thickness"))
        tmp_ext.Solid = True
        tmp_to_delete.append(tmp_ext)
        tmp_ext.recompute()

        bp = BOPFeatures.BOPFeatures(App.ActiveDocument)
        exclusions: List[rhg.HoleExclusion] = []
        scale_factor: float = obj.chord / 175.0
        rib_pen_standoff: float = scale_factor * 2.0
        for intersection in obj.intersections:
            c: Part.Feature = bp.make_common([tmp_ext.Name, intersection.Name])
            c.recompute()
            exclusions.append(rhg.HoleExclusion(c.Shape.BoundBox, rib_pen_standoff))
            tmp_to_delete.append(c)

        inner_profile_standoff: float = scale_factor * 2.0
        hbg =rhg.HoleBoundGenerator(tmp_af, inner_profile_standoff, exclusions)

        hbr_list = hbg.get_hole_bounds()

        # TODO: eventually what we're going to want is the enum field which allows
        #       you to select from one of several hole generators, and then a
        #       python object field which can raise a custom form for each
        #       generators parameters.  
        hbg = rhg.HoleGeneratorFactory.create_generator(obj.hole_type, obj.chord)
        if hbg is not None:
            for hbr in hbr_list:
                lh_sk = hbg.generate_sketch(hbr)
                lh_sk.recompute()
                tmp_to_delete.append(lh_sk)
                tmp_af.addGeometry(lh_sk.Geometry)
                
        tmp_af.recompute()        

        ext = App.ActiveDocument.addObject("Part::Extrusion", "tmp_ext")
        ext.Base = tmp_af
        ext.Dir = App.Vector(0, 0, obj.getPropertyByName("thickness"))
        ext.Solid = True
        ext.recompute()
        tmp_to_delete.append(ext)
        
        for intersection in obj.intersections:
            ext: Part.Feature = bp.make_cut([ext.Name, intersection.Name])
            ext.recompute()
            tmp_to_delete.append(ext)

        obj.Shape = ext.Shape.copy()
        obj.Placement = tmp_placement

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
