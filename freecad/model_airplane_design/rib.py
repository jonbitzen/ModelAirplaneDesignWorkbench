from . import airfoil
from . import rib_hole_generators as rhg
import Draft
import FreeCAD as App
import Part
import PartDesign
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

        self.intersections: List[Part.Shape] = []

        self.Object = obj
        obj.Proxy = self

    def add_cut(self, other_solid: Part.Feature) -> None:
        """
        Adds a cut into the rib solid, and registers a section derived from the
        cut as an exclusion region so that other rib structures may not be
        generated there

        Parameters
        ----------
        other_solid: Part.Feature
            The solid feature that will intersect with the wing rib and generate
            a cut and an exclusion region

        """

        # TODO: store the name of the feature used to do the cut, and make a
        #       mapping to the originating feature and its cut, so that we can
        #       - remove cut features
        #       - check whether the cut features have changed in execute, so we
        #         can skip the hole generation if they havent
        #       Also I think this is a natural place to add a parameter for the
        #       standoff to use for the cut feature, since you may need different
        #       amounts of structure around the cut, depending on what it is

        # TODO: better yet, we might even like to store a reference to the
        #       interfering object, that way when we move things around a little
        #       we don't need to explicitly go and re-create all the interference
        #       sections, we can just recompute them by iterating over the list
        #       of things again.  Plus it might be easier to keep a list of
        #       object's we're currently checking, and evict them from the UI
        #       later

        # these start out in the xy plane
        af_sketch = self.airfoil_data.to_sketch(self.Object.chord)
        af_sketch.Placement = self.Object.Placement
        af_sketch.recompute()
        rib_face: Part.Feature = Part.makeFace([af_sketch.Shape], "Part::FaceMakerSimple")
        rib_face = Part.show(rib_face)
        
        sec: Part.Feature = App.ActiveDocument.addObject("Part::Section", "sec")
        sec.Base = other_solid
        sec.Tool = rib_face
        sec.recompute()

        # TODO: Since these exclusion/intersections MUST be in the XY plane, we
        #       probably ought to force them to be initialized in a way such that
        #       is *always* true, or that we at least guide ourselves to make it
        #       so
        # store intersections in the XY plane - when we calculate the lightening
        # holes, we will do everything on the "base" airfoil template scaled to
        # the correct chord, and then transform the hole sketch back into the
        # Rib's Placement
        new_intersection = sec.Shape.transformGeometry(self.Object.Placement.Matrix.inverse())
        self.intersections.append(new_intersection)

        App.ActiveDocument.removeObject(af_sketch.Name)
        App.ActiveDocument.removeObject(sec.Name)
        App.ActiveDocument.removeObject(rib_face.Name)

        self.execute(self.Object)


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

        # register temporary doc objects here, so they can be easily disposed
        # when we're done with them
        tmp_to_delete: List[Part.Feature] = []

        tmp_placement: App.Placement = obj.Placement

        tmp_af: Sketcher.SketchObject = self.airfoil_data.to_sketch(obj.getPropertyByName("chord"))
        tmp_af.recompute()
        tmp_to_delete.append(tmp_af)

        tmp_body: PartDesign.Body = App.ActiveDocument.addObject("PartDesign::Body", "tmp_body")
        tmp_to_delete.append(tmp_body)

        tmp_pad = tmp_body.newObject("PartDesign::Pad", "tmp_pad")
        tmp_pad.Profile = tmp_af
        tmp_pad.Length = obj.getPropertyByName("thickness")
        tmp_pad.recompute()
        tmp_to_delete.append(tmp_pad)

        tmp_pkt = None
        tmp_sketch = None
        if len(self.intersections) > 0:
            tmp_sketch = Draft.make_sketch(self.intersections, autoconstraints=True, name="tmp_sk")
            tmp_sketch.recompute()
            tmp_to_delete.append(tmp_sketch)

            tmp_pkt = tmp_body.newObject("PartDesign::Pocket", "tmp_pkt")
            tmp_pkt.Profile = tmp_sketch
            tmp_pkt.Length = 2*obj.getPropertyByName("thickness")
            tmp_pkt.Midplane = True
            tmp_pkt.recompute()
            tmp_to_delete.append(tmp_pkt)

        exclusions: List[rhg.HoleExclusion] = []
        for iw in self.intersections:
            exclusions.append(rhg.HoleExclusion(Part.Wire(iw.Edges), 2.0))

        hbg =rhg.HoleBoundGenerator(tmp_af, 2.0, exclusions)

        hbr_list = hbg.get_hole_bounds()

        max_hole_length = obj.getPropertyByName("chord") / 6

        hbg = rhg.RoundedTrapezoidHoleGenerator(2.0, max_hole_length, 2.0)
        for hbr in hbr_list:
            lh_sk = hbg.generate_sketch(hbr)
            lh_sk.recompute()
            tmp_to_delete.append(lh_sk)

            lh_pkt = tmp_body.newObject("PartDesign::Pocket", "lh_pkt")
            lh_pkt.Profile = lh_sk
            lh_pkt.Length = 2*obj.getPropertyByName("thickness")
            lh_pkt.Midplane = True
            lh_pkt.recompute()
            tmp_to_delete.append(lh_pkt)

        tmp_body.recompute()

        obj.Shape = tmp_body.Shape.copy()
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
