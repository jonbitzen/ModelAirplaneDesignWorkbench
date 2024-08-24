from . import airfoil
from . import rib_hole_generators as rhg
from . import utilities
from BOPTools import BOPFeatures
import FreeCAD as App
import Part, PartDesign
import Sketcher
from typing import List, Dict, Tuple

def create(
    obj_name: str
) -> App.DocumentObject:
    
    body: PartDesign.Body = App.ActiveDocument.addObject("PartDesign::Body", obj_name)
    rib_feature = App.ActiveDocument.addObject("PartDesign::FeaturePython", obj_name + "_base")

    Rib(rib_feature)
    RibViewProvider(rib_feature.ViewObject)

    body.addObject(rib_feature)

    return body

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
            "trim_le",
            "Rib",
            "Trim the rib some distance from the leading edge"
        ).trim_le = 0.0

        obj.addProperty(
            "App::PropertyFloat",
            "trim_te",
            "Rib",
            "Trim the rib some distance from the trailing edge"
        ).trim_te = 0.0

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

        obj.addProperty(
            "App::PropertyLink",
            "le_trim_line",
            "Rib",
            "Line reference used to trim the rib from the leading edge back toward the trailing edge"
        ).le_trim_line = None

        obj.addProperty(
            "App::PropertyLink",
            "te_trim_line",
            "Rib",
            "Line reference used to trim the rib from the trailing edge toward the leading edge"
        ).te_trim_line = None

        # TODO: These don't get copied when the object is copied
        self.Object = obj
        obj.Proxy = self

    def onChanged(self, obj: App.DocumentObject, property: str) -> None:
        """
        Called when the FeaturePython object properties change.  Note that this
        includes any properties of the object, in addition to the App::Property*
        objects that may be added explictly by the user
        """  

        match property:
            case "te_trim_line":
                if obj.te_trim_line is None:
                    obj.trim_te = 0.0
                self.execute(obj)
            case "le_trim_line":
                if obj.le_trim_line is None:
                    obj.trim_le = 0.0
                self.execute(obj)
            case "airfoil" | "hole_type" | "interferences":
                    self.execute(obj)
            case _:
                pass

    def attach(self, obj: App.DocumentObject) -> None:
        pass


    def __make_body_centered_rib_feature(self, rib_sketch: Sketcher.Sketch, thickness: float) -> Part.Feature:
        """
        Creates a Part::Feature solid object from the given rib airfoil sketch,
        and the target thickness property
        """
        with utilities.TempDocObjectHelper() as tmp_obj_helper:
            rib_extr: Part.Feature = tmp_obj_helper.addObject(App.ActiveDocument.addObject("Part::Extrusion", "rib_extr"), do_delete=True)
            rib_extr.Base = rib_sketch
            rib_extr.Dir = App.Vector(0, 0, thickness)
            rib_extr.Solid = True
            rib_extr.recompute()

            body_center = rib_sketch.Shape.BoundBox.Center
            transform_mtx = App.Matrix()
            transform_mtx.move(-body_center)
            
            rib_ftr: Part.Feature = tmp_obj_helper.addObject(App.ActiveDocument.addObject("Part::Feature", "rib_ftr"), do_delete=True)
            rib_ftr.Shape = rib_extr.Shape.transformed(transform_mtx, copy=True)
            tmp_obj_helper.removeObject(rib_ftr)
            
            return rib_ftr

    def __get_te_distance(self, edges: List[Part.Edge], v_te: float, matrix: App.Matrix) -> float:

        if not edges:
            return None

        v_0 = utilities.origin

        v_te_rib = matrix.multVec(v_te)
        v_0_rib = matrix.multVec(v_0)

        plane_normal_xf = matrix.multVec(utilities.z_axis)
        plane = Part.Plane(v_te_rib, plane_normal_xf)

        p: App.Vector = None

        # search the list of edges for the first point intersection with the rib plane
        for edge in edges:
            curve: Part.Curve = edge.Curve
            trimmed_curve = curve.trim(curve.FirstParameter, curve.LastParameter)

            intersections = plane.intersect(trimmed_curve)

            for inter in intersections:
                if not inter:
                    continue

                i = inter[0]
                if type(i) is Part.Point:
                    p = App.Vector(i.X, i.Y, i.Z)
                break

            if p is not None:
                break

        # if there was no intersection between the rib plane and the the reference
        # line, then return an empty value
        if p is None:
            return None
        
        # if we get here, then we've got all we need to get the projected length
        # from the reference line down on the rib
        A: App.Vector = p - v_te_rib
        B: App.Vector = v_0_rib - v_te_rib
        
        projection = A.dot(B)/B.Length

        return projection

    def execute(self, obj: App.DocumentObject) -> None:

        with utilities.TempDocObjectHelper() as tmp_obj_helper:

            intf_vis: Dict[str, bool] = {}
            for intf in obj.interferences:
                intf_vis[intf.Name] = intf.ViewObject.Visibility

            bdy: PartDesign.Body = obj.getParent()

            orig_placement: App.Placement = bdy.Placement

            self.airfoil_data = airfoil.load(airfoil.AirfoilType.to_filename(obj.airfoil))

            # create rib sketchin in the body local coord system (LCS)
            rib_sketch: Sketcher.SketchObject = tmp_obj_helper.addObject(self.airfoil_data.to_sketch(obj.chord), do_delete=True)
            rib_sketch.Shape.tessellate(0.01)
            rib_sketch.recompute()

            # create body-centered solid in body LCS
            rib_extr = tmp_obj_helper.addObject(self.__make_body_centered_rib_feature(rib_sketch, obj.thickness), do_delete=True)

            # TODO: make this a method somehow
            # calculate trim distances to reference linesin the body world 
            # coordinate system (WCS)
            v_te = App.Vector(rib_extr.Shape.BoundBox.XMax, 0, 0)
            trim_le = 0.0
            trim_te = 0.0
            if obj.te_trim_line is not None:
                ref_edges = obj.te_trim_line.Shape.Edges
                te_distance = self.__get_te_distance(ref_edges, v_te, orig_placement.Matrix)
                if te_distance < obj.chord and te_distance > 0:
                    obj.trim_te = te_distance
                    trim_te = te_distance
                else:
                    obj.trim_te = 0.0
            else:
                trim_te = obj.trim_te

            if obj.le_trim_line is not None:
                ref_edges = obj.le_trim_line.Shape.Edges
                te_distance = self.__get_te_distance(ref_edges, v_te, orig_placement.Matrix)
                if te_distance < obj.chord and te_distance > 0:
                    obj.trim_le = obj.chord - te_distance
                    trim_le = obj.chord - te_distance
                else:
                    obj.trim_le = 0.0
            else:
                trim_le = obj.trim_le

            # TODO: we only need to calculate interferences if the hole generator is
            #       valid, so we could skip all this in many cases

            # translate the rib solid from body LCS to body WCS so that we can
            # find boolean intersection with the interferences
            rib_extr.Placement = orig_placement
            rib_extr.recompute()

            # create a set of hole exclusion regions by finding the common bool
            # between each of the intersecting objects and the rib and applying a
            # standoff to the bbox of the resulting shape
            boolean_tool = BOPFeatures.BOPFeatures(App.ActiveDocument)
            hole_exclusions: List[rhg.HoleExclusion] = []
            scale_factor: float = obj.chord / 175.0
            rib_pen_standoff: float = scale_factor * 2.0
            for interference in obj.interferences:

                # find the intersection between the rib solid and the intersection
                # in the body WCS
                rib_intersection: Part.Feature = tmp_obj_helper.addObject(boolean_tool.make_common([rib_extr.Name, interference.Name]), do_delete=True)
                rib_intersection.recompute()

                if rib_intersection.Shape.isNull():
                    continue

                # NOTE: for some reason isNull() sometimes returns False, even when
                #       there is no possible intersection, but the shape volume is
                #       at least zero so hopefully we can rely on that
                if rib_intersection.Shape.Volume < utilities.epsilon:
                    continue

                # the boolean intersection's *shape* coords are in world coords
                # and its placement is at the origin; we apply the body inverse
                # world transform to the boolean *shape* coordinates to get it
                # back into the body LCS, so we can calculate interference
                # intervals along the x axis for each; however after c
                rib_i_ft: Part.Feature = tmp_obj_helper.addObject(App.ActiveDocument.addObject("Part::Feature", "rib_i_ft"), do_delete=True)
                rib_i_ft.Shape = rib_intersection.Shape.transformed(rib_extr.Placement.Matrix.inverse(), copy=True)
                rib_i_ft.Placement.Matrix.move(rib_sketch.Shape.BoundBox.Center)
                rib_i_ft.recompute()

                hole_exclusions.append(rhg.HoleExclusion(rib_i_ft.Shape.BoundBox, rib_pen_standoff))            

            # create trim boxes in the rib's body-centered coordinate space to
            # trim the rib some distance with respect to the leading edge and/or
            # trailing edge
            trim_le = trim_le if trim_le >= 0.0 else 0.0
            trim_te = trim_te if trim_te >= 0.0 else 0.0
            
            do_trim = False if (trim_te+trim_le) > 0.8 * obj.chord else True

            le_trim_bbox: App.BoundBox = None
            te_trim_bbox: App.BoundBox = None
            if do_trim:
                if trim_le > 0.0:
                    le_trim_bbox = App.BoundBox()
                    le_trim_bbox.XMin = rib_sketch.Shape.BoundBox.XMin
                    le_trim_bbox.XMax = le_trim_bbox.XMin + trim_le
                    le_trim_bbox.YMax = rib_sketch.Shape.BoundBox.YMax + utilities.epsilon
                    le_trim_bbox.YMin = rib_sketch.Shape.BoundBox.YMin - utilities.epsilon
                    hole_exclusions.append(rhg.HoleExclusion(le_trim_bbox, rib_pen_standoff))

                if trim_te > 0.0:
                    te_trim_bbox = App.BoundBox()
                    te_trim_bbox.XMin = rib_sketch.Shape.BoundBox.XMax - trim_te
                    te_trim_bbox.XMax = rib_sketch.Shape.BoundBox.XMax
                    te_trim_bbox.YMax = rib_sketch.Shape.BoundBox.YMax + utilities.epsilon
                    te_trim_bbox.YMin = rib_sketch.Shape.BoundBox.YMin - utilities.epsilon
                    hole_exclusions.append(rhg.HoleExclusion(te_trim_bbox, rib_pen_standoff))

            inner_profile_standoff: float = scale_factor * 2.0
            hbg =rhg.HoleBoundGenerator(rib_sketch, inner_profile_standoff, hole_exclusions)

            hbr_list = hbg.get_hole_bounds()
    
            # Calculate lightening holes for each hole bound region, and add the
            # geometry to the rib sketch 
            hbg = rhg.HoleGeneratorFactory.create_generator(obj.hole_type, obj.chord)
            if hbg is not None:
                for hbr in hbr_list:
                    lh_sk = tmp_obj_helper.addObject(hbg.generate_sketch(hbr), do_delete=True)
                    lh_sk.recompute()
                    if lh_sk.Shape.isNull() or not lh_sk.Shape.isClosed():
                        print("lightening hole sketch geometry is either null, or is not closed, skipping")
                        continue
                    rib_sketch.addGeometry(lh_sk.Geometry)
                    
            rib_sketch.recompute()        

            # create a solid for the final rib shape, starting with the sketch that
            # has the lightening holes in it; this final shape is in the body LCS
            rib_final_extr = tmp_obj_helper.addObject(self.__make_body_centered_rib_feature(rib_sketch, obj.thickness), do_delete=True)
            
            # make boolean trim cuts here before we change the rib transform
            if le_trim_bbox is not None:
                le_trim_tool = tmp_obj_helper.addObject(App.ActiveDocument.addObject("Part::Box", "le_trim_tool"), do_delete=True)
                # NOTE: we add margins to the tool size, and compensate with the
                #       placement to ensure that the tool entirely covers the
                #       airfoil extrusion, and doesn't leave slivers behind
                le_trim_tool.Length = le_trim_bbox.XLength + 1
                le_trim_tool.Width = le_trim_bbox.YLength + 2
                le_trim_tool.Height = obj.thickness * 2
                
                le_trim_tool.Placement.move(
                    App.Vector(
                        rib_final_extr.Shape.BoundBox.XMin - 1, 
                        rib_final_extr.Shape.BoundBox.YMin - 1, 
                        -0.5*obj.thickness
                    )
                )
                le_trim_tool.recompute()
                rib_final_extr: Part.Feature = tmp_obj_helper.addObject(boolean_tool.make_cut([rib_final_extr.Name, le_trim_tool.Name]), do_delete=True)
                rib_final_extr.recompute()

            if te_trim_bbox is not None:
                te_trim_tool = tmp_obj_helper.addObject(App.ActiveDocument.addObject("Part::Box", "te_trim_tool"), do_delete=True)
                te_trim_tool.Length = te_trim_bbox.XLength + 1
                te_trim_tool.Width = te_trim_bbox.YLength + 2
                te_trim_tool.Height = obj.thickness * 2
                
                te_trim_tool.Placement.move(
                    App.Vector(
                        rib_final_extr.Shape.BoundBox.XMax - trim_te, 
                        rib_final_extr.Shape.BoundBox.YMin - 1, 
                        -0.5*obj.thickness
                    )
                )
                te_trim_tool.recompute()
                rib_final_extr: Part.Feature = tmp_obj_helper.addObject(boolean_tool.make_cut([rib_final_extr.Name, te_trim_tool.Name]), do_delete=True)
                rib_final_extr.recompute()


            obj.Shape = rib_final_extr.Shape.copy()

            # move the rib into its final transform space
            # TODO: we might not need to do this, since everything is in the
            #       parent's coordinate space now, which ends up being what we
            #       always wanted 
            rib_final_extr.Shape = rib_final_extr.Shape.transformed(orig_placement.Matrix,True)

            # make cuts for each of the interferences in the rib solid; note that
            # after these cuts are made, the new final rib solid's *shape* is in
            # world coordinates, and its placement is at the origin
            for interference in obj.interferences:
                rib_final_extr: Part.Feature = tmp_obj_helper.addObject(boolean_tool.make_cut([rib_final_extr.Name, interference.Name]), do_delete=True)
                rib_final_extr.recompute()

            obj.Shape = rib_final_extr.Shape.transformed(orig_placement.Matrix.inverse(), True)

            for intf in obj.interferences:
                intf.ViewObject.Visibility = intf_vis[intf.Name]


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
