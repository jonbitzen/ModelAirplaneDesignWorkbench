from . import airfoil
from . import rib_hole_generators as rhg
from . import utilities
from BOPTools import BOPFeatures
import FreeCAD as App
import Part, PartDesign
import Sketcher
from typing import List, Tuple

def create(
    obj_name: str
) -> App.DocumentObject:
    
    body: PartDesign.Body = App.ActiveDocument.addObject("PartDesign::Body", obj_name)
    rib_feature = App.ActiveDocument.addObject("PartDesign::FeaturePython", obj_name + "_base")

    Rib(rib_feature)
    RibViewProvider(rib_feature.ViewObject)

    body.addObject(rib_feature)

    return body

def make_copy_in_global_coords(feature: Part.Feature):
        tmp_base: Part.Feature = App.ActiveDocument.addObject("Part::Feature", "tmp_base")
        s = feature.Shape.copy()
        s.Placement = utilities.xy_placement
        tmp_base.Shape = s
        tmp_base.Placement = feature.getGlobalPlacement()
        return tmp_base

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

    class TrimDist():
        def __init__(self, le_dist: float, te_dist: float) -> None:
            self.le_dist = le_dist
            self.te_dist = te_dist

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

    def __get_rib_ref_line_trim_distance(
            self,
            body: Part.Feature, 
            rib_bbox: App.BoundBox, 
            le_trim_ref: Part.Feature, 
            te_trim_ref: Part.Feature) -> Tuple[float]:
                
        chord = rib_bbox.XLength
        le_trim: float = self.__get_te_trim_distance(body, rib_bbox, le_trim_ref)
        if le_trim is not None:
            le_trim = chord - le_trim
        te_trim: float = self.__get_te_trim_distance(body, rib_bbox, te_trim_ref)

        return Rib.TrimDist(le_trim, te_trim)

    def __get_te_trim_distance(
            self,
            body: Part.Feature, 
            rib_bbox: App.BoundBox, 
            trim_ref: Part.Feature) -> float:
        
        if trim_ref is None:
            return None
        
        edges: List[Part.Edge] = trim_ref.Shape.Edges
        if not edges:
            return None

        pt_te = App.Vector(rib_bbox.XMax, 0, 0)
        pt_le = App.Vector(rib_bbox.XMin, 0, 0)

        plane_normal = body.getGlobalPlacement().Matrix.multVec(utilities.z_axis)
        plane = Part.Plane(pt_te, plane_normal)

        pt_on_ref: App.Vector = None
        for edge in edges:
            curve: Part.Curve = edge.Curve
            trimmed_curve = curve.trim(curve.FirstParameter, curve.LastParameter)
            intersections = plane.intersect(trimmed_curve)
            for intersection in intersections:
                if not intersection:
                    continue

                i = intersection[0]
                if type(i) is Part.Point:
                    pt_on_ref = App.Vector(i.X, i.Y, i.Z)
                    break

            if pt_on_ref is not None:
                break

        if pt_on_ref is None:
            return None
        
        A: App.Vector = pt_on_ref - body.getGlobalPlacement().Matrix.multVec(pt_te)
        B: App.Vector = body.getGlobalPlacement().Matrix.multVec(pt_le-pt_te)

        te_dist = A.dot(B)/B.Length
        rib_chord = rib_bbox.XLength

        # only return a value if the projected distance is positive, and if it
        # is within the rib chord length; anything else indicates that the ref
        # line was either in front the LE, or behind the TE
        if te_dist < rib_chord and te_dist > 0:
            return te_dist
        
        return None


    def execute(self, obj: App.DocumentObject) -> None:

        with utilities.TempDocObjectHelper() as tmp_obj_helper:

            # somehow we can get here before the rib feature is attached to the
            # body, which causes this to be None
            bdy: PartDesign.Body = obj.getParent()
            if bdy is None:
                return

            self.airfoil_data = airfoil.load(airfoil.AirfoilType.to_filename(obj.airfoil))

            # create rib sketch in the body local coord system (LCS)
            rib_sketch: Sketcher.SketchObject = tmp_obj_helper.addObject(self.airfoil_data.to_sketch(obj.chord), do_delete=True)
            rib_sketch.Shape.tessellate(0.01)
            rib_sketch.recompute()

            # this is initially at the origin, it will be moved to the global
            # coordinate system later
            rib_intf_global = \
                tmp_obj_helper.addObject(
                    self.__make_body_centered_rib_feature(rib_sketch, obj.thickness), 
                    do_delete=True
                )

            # get the body-centered rib bbox while it's still at the origin 
            rib_bbox_bc = rib_intf_global.Shape.BoundBox

            # translate the rib solid from body LCS to body WCS so that we can
            # find boolean intersection with the interferences
            rib_intf_global.Placement = bdy.getGlobalPlacement()
            rib_intf_global.recompute()

            trim_dist: Rib.TrimDist = \
                self.__get_rib_ref_line_trim_distance(
                    bdy,
                    rib_bbox_bc,
                    obj.le_trim_line,
                    obj.te_trim_line
                )
            
            trim_te = obj.trim_te
            if obj.te_trim_line is not None:
                if trim_dist.te_dist is not None:
                    obj.trim_te = trim_dist.te_dist
                    trim_te = trim_dist.te_dist
                else:
                    obj.trim_te = 0.0

            trim_le = obj.trim_le
            if obj.le_trim_line is not None:
                if trim_dist.le_dist is not None:
                    obj.trim_le = trim_dist.le_dist
                    trim_le = trim_dist.le_dist
                else:
                    obj.trim_le = 0.0

            # TODO: we only need to calculate interferences if the hole generator is
            #       valid, so we could skip all this in many cases
            # create a set of hole exclusion regions by finding the common bool
            # between each of the intersecting objects and the rib and applying a
            # standoff to the bbox of the resulting shape
            boolean_tool = BOPFeatures.BOPFeatures(App.ActiveDocument)
            hole_exclusions: List[rhg.HoleExclusion] = []
            scale_factor: float = obj.chord / 175.0
            rib_pen_standoff: float = scale_factor * 2.0
            for interference in obj.interferences:

                intf_global = \
                    tmp_obj_helper.addObject(
                        make_copy_in_global_coords(interference),
                        do_delete=True
                    )

                # find the intersection between the rib solid and the intersection
                # in the body WCS
                rib_intersection: Part.Feature = \
                    tmp_obj_helper.addObject(
                        boolean_tool.make_common([rib_intf_global.Name, intf_global.Name]), 
                        do_delete=True
                    )
                rib_intersection.recompute()

                # NOTE: for some reason isNull() sometimes returns False, even when
                #       there is no possible intersection, but the shape volume is
                #       at least zero so hopefully we can rely on that
                if rib_intersection.Shape.isNull() or rib_intersection.Shape.Volume < utilities.epsilon:
                    print("No intersection between " + bdy.Name + " and " + interference.Name)
                    continue

                # the boolean intersection's *shape* coords are in world coords
                # and its placement is at the origin; we apply the body inverse
                # world transform to the boolean *shape* coordinates to get it
                # back into the body LCS, so we can calculate interference
                # intervals along the x axis for each; however after c
                rib_i_ft: Part.Feature = \
                    tmp_obj_helper.addObject(
                        App.ActiveDocument.addObject("Part::Feature", "rib_i_ft"), 
                        do_delete=True
                    )
                rib_i_ft.Shape = \
                    rib_intersection.Shape.transformed(
                        bdy.getGlobalPlacement().Matrix.inverse(), 
                        copy=True
                    )
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

            # move the working shape from the body LCS into body's WCS, so we can 
            # calculate make holes for the interferences; do it to the shape rather
            # than the transform, so we only have one thing to undo if there was
            # no interference calcuation
            rib_final_extr.Shape = \
                rib_final_extr.Shape.transformed(
                    bdy.getGlobalPlacement().Matrix, 
                    copy=True
                )
            # make cuts for each of the interferences in the rib solid; note that
            # after these cuts are made, the new final rib solid's *shape* is in
            # world coordinates, and its placement is at the origin
            for interference in obj.interferences:

                intf_global = \
                    tmp_obj_helper.addObject(
                        make_copy_in_global_coords(interference), 
                        do_delete=True
                    )

                rib_final_extr: Part.Feature = \
                    tmp_obj_helper.addObject(
                        boolean_tool.make_cut([rib_final_extr.Name, intf_global.Name]), 
                        do_delete=True
                    )
                rib_final_extr.recompute()

            obj.Shape = \
                rib_final_extr.Shape.transformed(
                    bdy.getGlobalPlacement().Matrix.inverse(), 
                    copy=True
                )


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
